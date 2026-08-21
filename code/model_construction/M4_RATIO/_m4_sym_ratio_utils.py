from __future__ import annotations

import csv
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\leng\AlCu_MAIN_PPP_Zhou")
MODEL = "M4_SYM_RATIO"
MODEL_DIR = ROOT / "models" / MODEL
STRUCTURE_DIR = MODEL_DIR / "structure"
SCRIPT_DIR = MODEL_DIR / "scripts"
METADATA_DIR = MODEL_DIR / "metadata"
DEBUG_DIR = MODEL_DIR / "debug"
PREVIEW_DIR = MODEL_DIR / "preview"
TMP_DIR = STRUCTURE_DIR / "_build_tmp"

AL_MASS = 26.9815385
CU_MASS = 63.546
TARGET_L = 145.608
REMOVE_DOUBLES_CUTOFF = 1.00

PHASES = {
    "Al": {"a_xy": 4.050, "c_z": 4.050, "nx": 36, "ny": 36},
    "Cu": {"a_xy": 3.615, "c_z": 3.615, "nx": 40, "ny": 40},
    "Al2Cu": {"a_xy": 6.067, "c_z": 4.877, "nx": 24, "ny": 24},
    "Al4Cu9": {"a_xy": 8.61773054, "c_z": 8.61773054, "nx": 17, "ny": 17},
}

REF_REPEATS = {
    "Cu_half": 12,
    "Al4Cu9": 4,
    "Al2Cu": 6,
    "Al": 24,
}

RATIO_REPEATS = {
    "Cu_half": 12,
    "Al4Cu9": 2,
    "Al2Cu": 10,
    "Al": 24,
}


def ensure_model_dirs() -> None:
    for directory in [MODEL_DIR, STRUCTURE_DIR, SCRIPT_DIR, METADATA_DIR, DEBUG_DIR, PREVIEW_DIR, TMP_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def strain_for(phase: str) -> float:
    info = PHASES[phase]
    return TARGET_L / (info["nx"] * info["a_xy"]) - 1.0


def phase_thickness(phase: str, repeat_z: int) -> float:
    return PHASES[phase]["c_z"] * repeat_z


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def parse_lammps_data(path: Path):
    lines = read_text(path).splitlines()
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
                    atom_type = int(tokens[0])
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
                    type_elements[atom_type] = element
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

    required = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi"}
    if set(box) != required:
        missing = sorted(required - set(box))
        raise ValueError(f"Missing box fields in {path}: {missing}")
    if not type_elements:
        raise ValueError(f"Could not read type mapping from {path}")
    if not atoms:
        raise ValueError(f"Could not read atoms from {path}")
    return box, type_elements, atoms


def box_lengths(box: dict[str, float]) -> dict[str, float]:
    return {
        "lx": box["xhi"] - box["xlo"],
        "ly": box["yhi"] - box["ylo"],
        "lz": box["zhi"] - box["zlo"],
    }


def count_elements(type_elements: dict[int, str], atoms) -> dict[str, int]:
    counts = {"Al": 0, "Cu": 0}
    for _, atom_type, *_ in atoms:
        element = type_elements.get(atom_type)
        if element in counts:
            counts[element] += 1
    counts["total"] = len(atoms)
    return counts


def write_mapped_data(source: Path, target: Path, title: str) -> dict[str, int]:
    box, type_elements, atoms = parse_lammps_data(source)
    element_to_type = {"Al": 1, "Cu": 2}
    counts = {1: 0, 2: 0}
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
            new_type = element_to_type[element]
            counts[new_type] += 1
            handle.write(f"{atom_id:10d} {new_type:3d} {x:20.12f} {y:20.12f} {z:20.12f}\n")
    return {"Al": counts[1], "Cu": counts[2], "total": len(atoms)}


def write_data_from_atoms(box: dict[str, float], atoms, target: Path, title: str) -> dict[str, int]:
    counts = {1: 0, 2: 0}
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
            counts[atom_type] += 1
            handle.write(f"{new_id:10d} {atom_type:3d} {x:20.12f} {y:20.12f} {z:20.12f}\n")
    return {"Al": counts[1], "Cu": counts[2], "total": len(atoms)}


def find_close_pairs(box: dict[str, float], atoms, cutoff: float):
    lengths = box_lengths(box)
    lx, ly, lz = lengths["lx"], lengths["ly"], lengths["lz"]
    xlo, ylo, zlo = box["xlo"], box["ylo"], box["zlo"]
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
    pairs = []
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


def remove_obvious_overlaps(raw_data: Path, ready_data: Path, cutoff: float = REMOVE_DOUBLES_CUTOFF):
    box, _, atoms = parse_lammps_data(raw_data)
    initial = len(atoms)
    total_removed = 0
    iterations = 0
    for _ in range(8):
        close_pairs = find_close_pairs(box, atoms, cutoff)
        if not close_pairs:
            break
        iterations += 1
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
        f"{ready_data.stem}; type 1 = Al, type 2 = Cu; overlap cutoff {cutoff} A",
    )
    return counts, initial, total_removed, iterations


def min_distance_summary(data_path: Path, cutoff_small: float = 1.0, cutoff_warn: float = 1.5) -> dict[str, object]:
    box, _, atoms = parse_lammps_data(data_path)
    lengths = box_lengths(box)
    lx, ly, lz = lengths["lx"], lengths["ly"], lengths["lz"]
    xlo, ylo, zlo = box["xlo"], box["ylo"], box["zlo"]
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
    below_small = 0
    below_warn = 0
    for (ix, iy, iz), local_atoms in bins.items():
        for id_a, xa, ya, za in local_atoms:
            for dx_cell in (-1, 0, 1):
                for dy_cell in (-1, 0, 1):
                    for dz_cell in (-1, 0, 1):
                        key = ((ix + dx_cell) % nx, (iy + dy_cell) % ny, (iz + dz_cell) % nz)
                        for id_b, xb, yb, zb in bins.get(key, []):
                            if id_b <= id_a:
                                continue
                            dx = xa - xb
                            dy = ya - yb
                            dz = za - zb
                            dx -= round(dx / lx) * lx
                            dy -= round(dy / ly) * ly
                            dz -= round(dz / lz) * lz
                            r2 = dx * dx + dy * dy + dz * dz
                            if r2 < cutoff_small * cutoff_small:
                                below_small += 1
                            if r2 < cutoff_warn * cutoff_warn:
                                below_warn += 1
                            if r2 < min_r2:
                                min_r2 = r2
                                min_pair = (id_a, id_b)
    return {
        "data_file": str(data_path),
        "min_distance_A": math.sqrt(min_r2),
        "atom_i": min_pair[0],
        "atom_j": min_pair[1],
        "pairs_below_1p0_A": below_small,
        "pairs_below_1p5_A": below_warn,
        "cell_size_A": cell,
        "pbc": "p p p",
    }


def z_profile(data_path: Path, csv_path: Path, png_path: Path | None = None, bins: int = 300):
    box, type_elements, atoms = parse_lammps_data(data_path)
    zlo, zhi = box["zlo"], box["zhi"]
    dz = (zhi - zlo) / bins
    rows = []
    for bin_index in range(bins):
        rows.append(
            {
                "bin": bin_index,
                "z_min": zlo + bin_index * dz,
                "z_max": zlo + (bin_index + 1) * dz,
                "z_center": zlo + (bin_index + 0.5) * dz,
                "Al_atoms": 0,
                "Cu_atoms": 0,
                "total_atoms": 0,
                "Al_fraction": "",
                "Cu_fraction": "",
            }
        )
    for _, atom_type, _, _, z in atoms:
        idx = int((z - zlo) / dz)
        if idx == bins:
            idx -= 1
        idx = max(0, min(bins - 1, idx))
        element = type_elements.get(atom_type)
        if element == "Al":
            rows[idx]["Al_atoms"] += 1
        elif element == "Cu":
            rows[idx]["Cu_atoms"] += 1
        rows[idx]["total_atoms"] += 1
    for row in rows:
        total = row["total_atoms"]
        if total:
            row["Al_fraction"] = row["Al_atoms"] / total
            row["Cu_fraction"] = row["Cu_atoms"] / total
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    if png_path is not None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        xs = [row["z_center"] for row in rows if row["total_atoms"]]
        al = [row["Al_fraction"] for row in rows if row["total_atoms"]]
        cu = [row["Cu_fraction"] for row in rows if row["total_atoms"]]
        plt.figure(figsize=(9, 4.8))
        plt.plot(xs, al, label="Al fraction", color="#2c7fb8", linewidth=1.8)
        plt.plot(xs, cu, label="Cu fraction", color="#d95f0e", linewidth=1.8)
        plt.xlabel("z (A)")
        plt.ylabel("Atomic fraction")
        plt.ylim(-0.05, 1.05)
        plt.grid(True, alpha=0.25)
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(png_path, dpi=220)
        plt.close()
    return rows


def run_command(cmd: list[str], cwd: Path = ROOT, log_path: Path | None = None) -> subprocess.CompletedProcess:
    completed = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write("$ " + " ".join(cmd) + "\n")
            handle.write("--- stdout ---\n")
            handle.write(completed.stdout or "")
            handle.write("\n--- stderr ---\n")
            handle.write(completed.stderr or "")
            handle.write(f"\n--- returncode: {completed.returncode} ---\n\n")
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(cmd)}")
    return completed


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def expected_ratio_layer_sequence() -> list[tuple[str, str, int]]:
    return [
        ("Cu_half_left", "Cu", RATIO_REPEATS["Cu_half"]),
        ("Al4Cu9_left_thin", "Al4Cu9", RATIO_REPEATS["Al4Cu9"]),
        ("Al2Cu_left_adjusted", "Al2Cu", RATIO_REPEATS["Al2Cu"]),
        ("Al_center", "Al", RATIO_REPEATS["Al"]),
        ("Al2Cu_right_adjusted", "Al2Cu", RATIO_REPEATS["Al2Cu"]),
        ("Al4Cu9_right_thin", "Al4Cu9", RATIO_REPEATS["Al4Cu9"]),
        ("Cu_half_right", "Cu", RATIO_REPEATS["Cu_half"]),
    ]

