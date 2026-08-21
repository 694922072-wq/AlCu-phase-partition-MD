#!/usr/bin/env python
"""Build M3/M4 PPP Al-Cu main models with Atomsk and validate type mapping."""

from __future__ import annotations

import csv
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\leng\AlCu_MAIN_PPP_Zhou")
OLD_STRUCTURE = Path(r"D:\leng\AlCu_SAFE_SC\structure")
STRUCTURE = ROOT / "structure"
DOCS = ROOT / "docs"
VALIDATION = ROOT / "validation"
TMP = STRUCTURE / "_build_tmp"

TARGET_L = 145.608
REMOVE_DOUBLES_CUTOFF = 1.00

AL_MASS = 26.9815385
CU_MASS = 63.546

PHASES = {
    "Al": {"a_xy": 4.050, "c_z": 4.050, "nx": 36, "ny": 36},
    "Cu": {"a_xy": 3.615, "c_z": 3.615, "nx": 40, "ny": 40},
    "Al2Cu": {"a_xy": 6.067, "c_z": 4.877, "nx": 24, "ny": 24},
    "Al4Cu9": {"a_xy": 8.61773054, "c_z": 8.61773054, "nx": 17, "ny": 17},
}


def run(cmd: list[str]) -> None:
    print(" ".join(str(part) for part in cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)


def ensure_dirs() -> None:
    for directory in [
        ROOT / "inputs",
        STRUCTURE,
        ROOT / "potentials",
        ROOT / "scripts",
        ROOT / "logs",
        ROOT / "dumps",
        ROOT / "outputs",
        ROOT / "figures",
        DOCS,
        VALIDATION,
        TMP,
    ]:
        directory.mkdir(parents=True, exist_ok=True)
    if TMP.exists():
        for item in TMP.glob("*"):
            if item.is_file():
                item.unlink()


def check_sources() -> bool:
    missing = []
    for filename in ["Al2Cu.cif", "Al4Cu9_full.cif"]:
        source = OLD_STRUCTURE / filename
        if not source.exists():
            missing.append(str(source))
        else:
            shutil.copy2(source, STRUCTURE / filename)

    if missing:
        text = [
            "# Missing Structure Files",
            "",
            "The main PPP model build was stopped because required structure files were not found.",
            "",
            "Missing files:",
        ]
        text.extend(f"- `{item}`" for item in missing)
        text.extend(
            [
                "",
                "Please provide valid CIF/data files for Al2Cu and Al4Cu9. The script will not fabricate IMC unit cells.",
            ]
        )
        (DOCS / "missing_structure_files.md").write_text("\n".join(text), encoding="utf-8")
        return False
    return True


def write_type_mapping() -> None:
    (STRUCTURE / "type_mapping_main.txt").write_text(
        "type 1 = Al\n"
        "type 2 = Cu\n"
        "pair_coeff order = Al Cu\n",
        encoding="utf-8",
    )


def strain_for(phase: str) -> float:
    info = PHASES[phase]
    return TARGET_L / (info["nx"] * info["a_xy"]) - 1.0


def write_mismatch_tables() -> None:
    rows = []
    warnings = []
    for phase, info in PHASES.items():
        raw = info["nx"] * info["a_xy"]
        strain = TARGET_L / raw - 1.0
        rows.append(
            {
                "phase": phase,
                "nx": info["nx"],
                "ny": info["ny"],
                "raw_Lx_A": raw,
                "target_Lx_A": TARGET_L,
                "in_plane_strain": strain,
                "in_plane_strain_percent": 100.0 * strain,
            }
        )
        if abs(strain) > 0.02:
            warnings.append(f"- {phase}: in-plane strain is {100.0 * strain:.4f}%, larger than 2%.")

    for name in ["lattice_mismatch_M3.csv", "lattice_mismatch_M4.csv"]:
        with (VALIDATION / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    warning_path = VALIDATION / "warnings.md"
    if warnings:
        warning_path.write_text("# Build Warnings\n\n" + "\n".join(warnings) + "\n", encoding="utf-8")
    elif not warning_path.exists():
        warning_path.write_text("# Build Warnings\n\nNo in-plane strain exceeds 2%.\n", encoding="utf-8")


def write_thickness_table() -> None:
    al2cu_c = PHASES["Al2Cu"]["c_z"]
    al4cu9_c = PHASES["Al4Cu9"]["c_z"]
    m4_single_side_target = 4 * al4cu9_c + 6 * al2cu_c
    m3_repeat = round(m4_single_side_target / al2cu_c)
    rows = [
        {
            "model": "M3_SYM",
            "slab": "single_side_Al2Cu_thick",
            "target_thickness_A": m4_single_side_target,
            "actual_thickness_A": m3_repeat * al2cu_c,
            "repeat_z": m3_repeat,
            "notes": "Chosen to match one side of M4_SYM Al4Cu9 + Al2Cu IMC thickness.",
        },
        {
            "model": "M4_SYM",
            "slab": "single_side_Al4Cu9",
            "target_thickness_A": "",
            "actual_thickness_A": 4 * al4cu9_c,
            "repeat_z": 4,
            "notes": "Part of one-side composite IMC.",
        },
        {
            "model": "M4_SYM",
            "slab": "single_side_Al2Cu",
            "target_thickness_A": "",
            "actual_thickness_A": 6 * al2cu_c,
            "repeat_z": 6,
            "notes": "Part of one-side composite IMC.",
        },
        {
            "model": "M4_SYM",
            "slab": "single_side_Al4Cu9_plus_Al2Cu",
            "target_thickness_A": "60-65",
            "actual_thickness_A": m4_single_side_target,
            "repeat_z": "4 + 6",
            "notes": "Reference IMC thickness used to choose M3_SYM Al2Cu(thick).",
        },
    ]
    with (VALIDATION / "thickness_control.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def create_slab(phase: str, nz: int, tag: str) -> Path:
    raw = TMP / f"{tag}_raw.cfg"
    matched = TMP / f"{tag}_match.cfg"
    if phase == "Al":
        run(
            [
                "atomsk",
                "--create",
                "fcc",
                "4.05",
                "Al",
                "orient",
                "[100]",
                "[010]",
                "[001]",
                "-duplicate",
                "36",
                "36",
                str(nz),
                str(raw),
            ]
        )
    elif phase == "Cu":
        run(
            [
                "atomsk",
                "--create",
                "fcc",
                "3.615",
                "Cu",
                "orient",
                "[100]",
                "[010]",
                "[001]",
                "-duplicate",
                "40",
                "40",
                str(nz),
                str(raw),
            ]
        )
    elif phase == "Al2Cu":
        run(
            [
                "atomsk",
                str(STRUCTURE / "Al2Cu.cif"),
                "-orthogonal-cell",
                "-duplicate",
                "24",
                "24",
                str(nz),
                str(raw),
            ]
        )
    elif phase == "Al4Cu9":
        run(
            [
                "atomsk",
                str(STRUCTURE / "Al4Cu9_full.cif"),
                "-orthogonal-cell",
                "-duplicate",
                "17",
                "17",
                str(nz),
                str(raw),
            ]
        )
    else:
        raise ValueError(f"Unknown phase: {phase}")

    strain = strain_for(phase)
    if abs(strain) > 1e-10:
        run(
            [
                "atomsk",
                str(raw),
                "-deform",
                "x",
                f"{strain:.12g}",
                "0",
                "-deform",
                "y",
                f"{strain:.12g}",
                "0",
                "-wrap",
                str(matched),
            ]
        )
    else:
        run(["atomsk", str(raw), "-wrap", str(matched)])
    return matched


def parse_lammps_data(path: Path) -> tuple[dict[str, float], dict[int, str], list[tuple[int, int, float, float, float]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    box: dict[str, float] = {}
    type_elements: dict[int, str] = {}
    atoms: list[tuple[int, int, float, float, float]] = []

    for line in lines:
        parts = line.split()
        if len(parts) >= 4 and parts[-2:] == ["xlo", "xhi"]:
            box["xlo"], box["xhi"] = float(parts[0]), float(parts[1])
        elif len(parts) >= 4 and parts[-2:] == ["ylo", "yhi"]:
            box["ylo"], box["yhi"] = float(parts[0]), float(parts[1])
        elif len(parts) >= 4 and parts[-2:] == ["zlo", "zhi"]:
            box["zlo"], box["zhi"] = float(parts[0]), float(parts[1])

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped.startswith("Masses"):
            i += 1
            while i < len(lines):
                mass_line = lines[i].strip()
                if not mass_line:
                    i += 1
                    continue
                if re.match(r"^[A-Za-z]", mass_line):
                    break
                fields = mass_line.split("#", 1)
                tokens = fields[0].split()
                if len(tokens) >= 2 and tokens[0].isdigit():
                    old_type = int(tokens[0])
                    mass = float(tokens[1])
                    comment = fields[1].strip().lower() if len(fields) > 1 else ""
                    if "al" in comment and "cu" not in comment:
                        element = "Al"
                    elif "cu" in comment:
                        element = "Cu"
                    elif abs(mass - AL_MASS) < abs(mass - CU_MASS):
                        element = "Al"
                    else:
                        element = "Cu"
                    type_elements[old_type] = element
                i += 1
        elif stripped.startswith("Atoms"):
            i += 1
            while i < len(lines):
                atom_line = lines[i].strip()
                if not atom_line:
                    i += 1
                    continue
                if re.match(r"^[A-Za-z]", atom_line):
                    break
                tokens = atom_line.split()
                if len(tokens) >= 5 and tokens[0].lstrip("-").isdigit():
                    atoms.append((int(tokens[0]), int(tokens[1]), float(tokens[2]), float(tokens[3]), float(tokens[4])))
                i += 1
        else:
            i += 1

    if not type_elements:
        raise ValueError(f"Could not read Masses/type mapping from {path}")
    if not atoms:
        raise ValueError(f"Could not read Atoms section from {path}")
    return box, type_elements, atoms


def write_mapped_data(source: Path, target: Path, title: str) -> dict[str, int]:
    box, type_elements, atoms = parse_lammps_data(source)
    element_to_new_type = {"Al": 1, "Cu": 2}
    type_counts = {1: 0, 2: 0}
    with target.open("w", encoding="utf-8") as handle:
        handle.write(f"# {title}\n\n")
        handle.write(f"{len(atoms):12d}  atoms\n")
        handle.write("           2  atom types\n\n")
        handle.write(f"{box['xlo']:20.12f} {box['xhi']:20.12f}  xlo xhi\n")
        handle.write(f"{box['ylo']:20.12f} {box['yhi']:20.12f}  ylo yhi\n")
        handle.write(f"{box['zlo']:20.12f} {box['zhi']:20.12f}  zlo zhi\n\n")
        handle.write("Masses\n\n")
        handle.write(f"1 {AL_MASS:.7f} # Al\n")
        handle.write(f"2 {CU_MASS:.6f} # Cu\n\n")
        handle.write("Atoms # atomic\n\n")
        for atom_id, old_type, x, y, z in atoms:
            element = type_elements[old_type]
            new_type = element_to_new_type[element]
            type_counts[new_type] += 1
            handle.write(f"{atom_id:10d} {new_type:3d} {x:20.12f} {y:20.12f} {z:20.12f}\n")
    return {"Al": type_counts[1], "Cu": type_counts[2], "total": len(atoms)}


def write_data_from_atoms(
    box: dict[str, float],
    atoms: list[tuple[int, int, float, float, float]],
    target: Path,
    title: str,
) -> dict[str, int]:
    type_counts = {1: 0, 2: 0}
    with target.open("w", encoding="utf-8") as handle:
        handle.write(f"# {title}\n\n")
        handle.write(f"{len(atoms):12d}  atoms\n")
        handle.write("           2  atom types\n\n")
        handle.write(f"{box['xlo']:20.12f} {box['xhi']:20.12f}  xlo xhi\n")
        handle.write(f"{box['ylo']:20.12f} {box['yhi']:20.12f}  ylo yhi\n")
        handle.write(f"{box['zlo']:20.12f} {box['zhi']:20.12f}  zlo zhi\n\n")
        handle.write("Masses\n\n")
        handle.write(f"1 {AL_MASS:.7f} # Al\n")
        handle.write(f"2 {CU_MASS:.6f} # Cu\n\n")
        handle.write("Atoms # atomic\n\n")
        for new_id, (_, atom_type, x, y, z) in enumerate(atoms, start=1):
            type_counts[atom_type] += 1
            handle.write(f"{new_id:10d} {atom_type:3d} {x:20.12f} {y:20.12f} {z:20.12f}\n")
    return {"Al": type_counts[1], "Cu": type_counts[2], "total": len(atoms)}


def find_close_pairs(
    box: dict[str, float],
    atoms: list[tuple[int, int, float, float, float]],
    cutoff: float,
) -> list[tuple[int, int, float]]:
    xlo, xhi = box["xlo"], box["xhi"]
    ylo, yhi = box["ylo"], box["yhi"]
    zlo, zhi = box["zlo"], box["zhi"]
    lx, ly, lz = xhi - xlo, yhi - ylo, zhi - zlo
    cell = max(cutoff * 1.2, 1.0)
    nx = max(1, int(lx / cell))
    ny = max(1, int(ly / cell))
    nz = max(1, int(lz / cell))
    bins: dict[tuple[int, int, int], list[int]] = defaultdict(list)

    for idx, (_, _, x, y, z) in enumerate(atoms):
        ix = int((x - xlo) / lx * nx) % nx
        iy = int((y - ylo) / ly * ny) % ny
        iz = int((z - zlo) / lz * nz) % nz
        bins[(ix, iy, iz)].append(idx)

    cutoff2 = cutoff * cutoff
    pairs: list[tuple[int, int, float]] = []
    for (ix, iy, iz), local_indices in bins.items():
        for idx_a in local_indices:
            id_a, _, xa, ya, za = atoms[idx_a]
            for dx_cell in (-1, 0, 1):
                for dy_cell in (-1, 0, 1):
                    for dz_cell in (-1, 0, 1):
                        key = ((ix + dx_cell) % nx, (iy + dy_cell) % ny, (iz + dz_cell) % nz)
                        for idx_b in bins.get(key, []):
                            if idx_b <= idx_a:
                                continue
                            id_b, _, xb, yb, zb = atoms[idx_b]
                            dx = xa - xb
                            dy = ya - yb
                            dz = za - zb
                            dx -= round(dx / lx) * lx
                            dy -= round(dy / ly) * ly
                            dz -= round(dz / lz) * lz
                            r2 = dx * dx + dy * dy + dz * dz
                            if r2 < cutoff2:
                                pairs.append((id_a, id_b, math.sqrt(r2)))
    return pairs


def remove_obvious_overlaps(raw_data: Path, ready_data: Path, cutoff: float) -> tuple[dict[str, int], int, int]:
    box, _, atoms = parse_lammps_data(raw_data)
    initial = len(atoms)
    total_removed = 0

    for _ in range(8):
        close_pairs = find_close_pairs(box, atoms, cutoff)
        if not close_pairs:
            break
        remove_ids: set[int] = set()
        for id_a, id_b, _ in close_pairs:
            if id_a not in remove_ids and id_b not in remove_ids:
                remove_ids.add(max(id_a, id_b))
        if not remove_ids:
            break
        atoms = [atom for atom in atoms if atom[0] not in remove_ids]
        total_removed += len(remove_ids)

    counts = write_data_from_atoms(
        box,
        atoms,
        ready_data,
        f"{ready_data.stem}; type 1 = Al, type 2 = Cu; Python overlap cutoff {cutoff} A",
    )
    return counts, initial, total_removed


def min_distance_check(data_path: Path, csv_path: Path) -> None:
    box, _, atoms = parse_lammps_data(data_path)
    xlo, xhi = box["xlo"], box["xhi"]
    ylo, yhi = box["ylo"], box["yhi"]
    zlo, zhi = box["zlo"], box["zhi"]
    lx, ly, lz = xhi - xlo, yhi - ylo, zhi - zlo
    cell = 3.0
    nx = max(1, int(lx / cell))
    ny = max(1, int(ly / cell))
    nz = max(1, int(lz / cell))
    bins: dict[tuple[int, int, int], list[tuple[int, float, float, float]]] = defaultdict(list)

    for atom_id, _, x, y, z in atoms:
        ix = int((x - xlo) / lx * nx) % nx
        iy = int((y - ylo) / ly * ny) % ny
        iz = int((z - zlo) / lz * nz) % nz
        bins[(ix, iy, iz)].append((atom_id, x, y, z))

    min_r2 = float("inf")
    min_pair = ("", "")
    below_1 = 0
    below_15 = 0
    for (ix, iy, iz), local_atoms in bins.items():
        for atom_a in local_atoms:
            id_a, xa, ya, za = atom_a
            for dx_cell in (-1, 0, 1):
                for dy_cell in (-1, 0, 1):
                    for dz_cell in (-1, 0, 1):
                        neighbor_key = ((ix + dx_cell) % nx, (iy + dy_cell) % ny, (iz + dz_cell) % nz)
                        for id_b, xb, yb, zb in bins.get(neighbor_key, []):
                            if id_b <= id_a:
                                continue
                            dx = xa - xb
                            dy = ya - yb
                            dz = za - zb
                            dx -= round(dx / lx) * lx
                            dy -= round(dy / ly) * ly
                            dz -= round(dz / lz) * lz
                            r2 = dx * dx + dy * dy + dz * dz
                            if r2 < 1.0:
                                below_1 += 1
                            if r2 < 2.25:
                                below_15 += 1
                            if r2 < min_r2:
                                min_r2 = r2
                                min_pair = (id_a, id_b)

    rows = [
        {
            "data_file": str(data_path),
            "min_distance_A": math.sqrt(min_r2),
            "atom_i": min_pair[0],
            "atom_j": min_pair[1],
            "pairs_below_1p0_A": below_1,
            "pairs_below_1p5_A": below_15,
            "cell_size_A": cell,
            "pbc": "p p p",
        }
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_model(model: str, slabs: list[Path], raw_out: Path, ready_out: Path) -> None:
    raw_atomsk = TMP / f"{model}_raw_atomsk.lmp"
    run(["atomsk", "--merge", "Z", str(len(slabs)), *[str(path) for path in slabs], str(raw_atomsk)])
    raw_counts = write_mapped_data(raw_atomsk, raw_out, f"{model} PPP raw data; type 1 = Al, type 2 = Cu")
    ready_counts, initial_atoms, removed = remove_obvious_overlaps(raw_out, ready_out, REMOVE_DOUBLES_CUTOFF)
    removed_pct = 100.0 * removed / raw_counts["total"] if raw_counts["total"] else 0.0

    (VALIDATION / f"{model}_atom_count.txt").write_text(
        f"raw_atoms = {raw_counts['total']}\n"
        f"ready_atoms = {ready_counts['total']}\n"
        f"removed_atoms = {removed}\n"
        f"removed_percent = {removed_pct:.6f}\n"
        f"ready_Al_atoms = {ready_counts['Al']}\n"
        f"ready_Cu_atoms = {ready_counts['Cu']}\n",
        encoding="utf-8",
    )
    (VALIDATION / f"{model}_build_summary.txt").write_text(
        f"{model} build summary\n"
        f"raw data: {raw_out}\n"
        f"ready data: {ready_out}\n"
        f"type 1 = Al\n"
        f"type 2 = Cu\n"
        f"Python obvious-overlap cutoff: {REMOVE_DOUBLES_CUTOFF} A\n"
        f"removed atoms: {removed} ({removed_pct:.6f}%)\n"
        f"Al atoms: {ready_counts['Al']}\n"
        f"Cu atoms: {ready_counts['Cu']}\n",
        encoding="utf-8",
    )
    min_distance_check(ready_out, VALIDATION / f"{model}_min_distance_check.csv")

    if removed_pct > 0.5:
        with (VALIDATION / "warnings.md").open("a", encoding="utf-8") as handle:
            handle.write(f"\n- {model}: remove-doubles deleted {removed_pct:.4f}% atoms (>0.5%). Inspect raw and ready files.\n")


def main() -> int:
    ensure_dirs()
    if not check_sources():
        return 2
    write_type_mapping()
    write_mismatch_tables()
    write_thickness_table()

    cu_half = create_slab("Cu", 12, "Cu_half_z12")
    cu_full = create_slab("Cu", 24, "Cu_full_z24")
    al = create_slab("Al", 24, "Al_z24")
    al2cu_thick = create_slab("Al2Cu", 13, "Al2Cu_thick_z13")
    al2cu_m4 = create_slab("Al2Cu", 6, "Al2Cu_M4_z6")
    al4cu9_m4 = create_slab("Al4Cu9", 4, "Al4Cu9_M4_z4")

    build_model(
        "M3_SYM",
        [cu_half, al2cu_thick, al, al2cu_thick, cu_half],
        STRUCTURE / "M3_SYM_Cu_Al2CuThick_Al_PPP_raw.data",
        STRUCTURE / "M3_SYM_Cu_Al2CuThick_Al_PPP_ready.data",
    )
    build_model(
        "M4_SYM",
        [cu_half, al4cu9_m4, al2cu_m4, al, al2cu_m4, al4cu9_m4, cu_half],
        STRUCTURE / "M4_SYM_Cu_Al4Cu9_Al2Cu_Al_PPP_raw.data",
        STRUCTURE / "M4_SYM_Cu_Al4Cu9_Al2Cu_Al_PPP_ready.data",
    )
    build_model(
        "M4_LIT",
        [cu_full, al4cu9_m4, al2cu_m4, al],
        STRUCTURE / "M4_LIT_Cu_Al4Cu9_Al2Cu_Al_PPP_raw.data",
        STRUCTURE / "M4_LIT_Cu_Al4Cu9_Al2Cu_Al_PPP_ready.data",
    )

    print("Build completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
