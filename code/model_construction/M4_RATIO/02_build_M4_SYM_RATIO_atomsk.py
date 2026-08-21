from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _m4_sym_ratio_utils import (  # noqa: E402
    AL_MASS,
    CU_MASS,
    DEBUG_DIR,
    METADATA_DIR,
    MODEL,
    PREVIEW_DIR,
    RATIO_REPEATS,
    REMOVE_DOUBLES_CUTOFF,
    ROOT,
    STRUCTURE_DIR,
    TARGET_L,
    TMP_DIR,
    box_lengths,
    count_elements,
    ensure_model_dirs,
    expected_ratio_layer_sequence,
    min_distance_summary,
    parse_lammps_data,
    phase_thickness,
    remove_obvious_overlaps,
    run_command,
    strain_for,
    write_csv,
    write_mapped_data,
    z_profile,
)


def atomsk_log() -> Path:
    return DEBUG_DIR / "atomsk_build_stdout_stderr.log"


def create_slab(phase: str, nz: int, tag: str, warnings: list[str]) -> Path:
    existing = ROOT / "structure" / "_build_tmp" / f"{tag}_match.cfg"
    matched = TMP_DIR / f"{tag}_match.cfg"
    raw = TMP_DIR / f"{tag}_raw.cfg"
    if existing.exists():
        shutil.copy2(existing, matched)
        return matched

    if phase == "Al":
        cmd = [
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
    elif phase == "Cu":
        cmd = [
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
    elif phase == "Al2Cu":
        source = ROOT / "structure" / "Al2Cu.cif"
        cmd = ["atomsk", str(source), "-orthogonal-cell", "-duplicate", "24", "24", str(nz), str(raw)]
    elif phase == "Al4Cu9":
        source = ROOT / "structure" / "Al4Cu9_full.cif"
        cmd = ["atomsk", str(source), "-orthogonal-cell", "-duplicate", "17", "17", str(nz), str(raw)]
    else:
        raise ValueError(f"Unknown phase: {phase}")

    run_command(cmd, ROOT, atomsk_log())
    strain = strain_for(phase)
    if abs(strain) > 1.0e-10:
        run_command(
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
            ],
            ROOT,
            atomsk_log(),
        )
    else:
        run_command(["atomsk", str(raw), "-wrap", str(matched)], ROOT, atomsk_log())

    if phase == "Al4Cu9" and "Al4Cu9_full.cif" not in " ".join(cmd):
        warnings.append("Al4Cu9_full.cif was not used for the Al4Cu9 slab.")
    return matched


def write_input_files(root_ready_name: str) -> None:
    inputs = ROOT / "inputs"
    inputs.mkdir(exist_ok=True)
    pair_order = "Al Cu"
    potential = "potentials/CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr_Zhou04.eam.alloy"

    min_text = f"""clear

units           metal
dimension       3
boundary        p p p
atom_style      atomic

read_data       structure/{root_ready_name}

mass            1 {AL_MASS}
mass            2 {CU_MASS}

pair_style      eam/alloy
pair_coeff      * * {potential} {pair_order}

neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes one 3000 page 30000

compute         peatom all pe/atom
compute         csym all centro/atom fcc

thermo          1000
thermo_style    custom step temp pe etotal press pxx pyy pzz lx ly lz vol atoms

fix             BR all box/relax aniso 0.0 vmax 0.001

min_style       cg
minimize        1.0e-12 1.0e-12 30000 300000

unfix           BR

write_data      outputs/M4_SYM_RATIO_min.data
write_restart   outputs/M4_SYM_RATIO_min.restart
"""
    (inputs / "in.M4_SYM_RATIO.01_min.lmp").write_text(min_text, encoding="utf-8")

    # Inherit the latest clean M4_SYM 02c style: NPT xy-coupled + short NVT polish, 50000 + 20000 steps.
    eq_text = f"""clear

units           metal
dimension       3
boundary        p p p
atom_style      atomic

read_data       outputs/M4_SYM_RATIO_min.data

mass            1 {AL_MASS}
mass            2 {CU_MASS}

pair_style      eam/alloy
pair_coeff      * * {potential} {pair_order}

neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes one 3000 page 30000

timestep        0.001
reset_timestep  0

compute         peatom all pe/atom
compute         csym all centro/atom fcc

thermo          1000
thermo_style    custom step temp pe etotal press pxx pyy pzz lx ly lz vol

dump            deq all custom 1000 dumps/M4_SYM_RATIO_eq_300K_xy_z.lammpstrj id type x y z c_peatom c_csym
dump_modify     deq sort id flush yes

velocity        all create 300.0 4928459 mom yes rot yes dist gaussian

fix             EQ1 all npt temp 300.0 300.0 0.1 x 0.0 0.0 1.0 y 0.0 0.0 1.0 z 0.0 0.0 1.0 couple xy
run             50000
unfix           EQ1

fix             EQ2 all nvt temp 300.0 300.0 0.1
run             20000
unfix           EQ2

undump          deq

write_data      outputs/M4_SYM_RATIO_eq_300K_xy_z.data
write_restart   outputs/M4_SYM_RATIO_eq_300K_xy_z.restart
"""
    (inputs / "in.M4_SYM_RATIO.02c_eq_300K_xy_z.lmp").write_text(eq_text, encoding="utf-8")

    tensile_text = f"""clear

units           metal
dimension       3
boundary        p p p
atom_style      atomic

read_data       outputs/M4_SYM_RATIO_eq_300K_xy_z.data

mass            1 {AL_MASS}
mass            2 {CU_MASS}

pair_style      eam/alloy
pair_coeff      * * {potential} {pair_order}

neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes one 3000 page 30000

timestep        0.001
reset_timestep  0

compute         peatom all pe/atom
compute         csym all centro/atom fcc

variable        L0 equal $(lx)
variable        erate equal 1.0e-3
variable        target_strain equal 0.25
variable        total_steps equal 250000

variable        strain equal (lx-v_L0)/v_L0
variable        strain_percent equal 100.0*v_strain
variable        stress_xx_GPa equal -pxx/10000.0

variable        vstep equal step
variable        vtemp equal temp
variable        vpe equal pe
variable        vetotal equal etotal
variable        vpxx equal pxx
variable        vpyy equal pyy
variable        vpzz equal pzz
variable        vlx equal lx
variable        vly equal ly
variable        vlz equal lz

thermo          1000
thermo_style    custom step temp pe etotal press pxx pyy pzz lx ly lz vol v_strain v_strain_percent v_stress_xx_GPa

fix             NVT all nvt temp 300.0 300.0 0.1
fix             DEF all deform 1 x erate ${{erate}} units box remap x

fix             SS all print 100 "${{vstep}} ${{strain}} ${{strain_percent}} ${{stress_xx_GPa}} ${{vtemp}} ${{vpe}} ${{vetotal}} ${{vpxx}} ${{vpyy}} ${{vpzz}} ${{vlx}} ${{vly}} ${{vlz}}" file outputs/M4_SYM_RATIO_stress_strain_x25_xy_z_clean.dat screen no title "step strain strain_percent stress_xx_GPa temp pe etotal pxx pyy pzz lx ly lz"

dump            D1 all custom 1000 dumps/M4_SYM_RATIO_tensile_x25_xy_z_clean_light.lammpstrj id type x y z c_peatom c_csym
dump_modify     D1 sort id flush yes

run             ${{total_steps}}

unfix           SS
unfix           DEF
unfix           NVT

undump          D1

write_data      outputs/M4_SYM_RATIO_tensile_x25_xy_z_clean_final.data
write_restart   outputs/M4_SYM_RATIO_tensile_x25_xy_z_clean_final.restart
"""
    (inputs / "in.M4_SYM_RATIO.03d_tensile_x_25_xy_z_clean.lmp").write_text(tensile_text, encoding="utf-8")


def write_server_command_files() -> None:
    commands = [
        "mpiexec -np 32 lmp -in inputs/in.M4_SYM_RATIO.01_min.lmp -log logs/log.M4_SYM_RATIO.01_min.lammps",
        "mpiexec -np 32 lmp -in inputs/in.M4_SYM_RATIO.02c_eq_300K_xy_z.lmp -log logs/log.M4_SYM_RATIO.02c_eq_300K_xy_z.lammps",
        "mpiexec -np 32 lmp -in inputs/in.M4_SYM_RATIO.03d_tensile_x_25_xy_z_clean.lmp -log logs/log.M4_SYM_RATIO.03d_tensile_x_25_xy_z_clean.lammps",
    ]
    (ROOT / "models" / MODEL / "scripts" / "run_M4_SYM_RATIO_on_server_commands.txt").write_text(
        "\n".join(commands) + "\n", encoding="utf-8"
    )
    bat = ["@echo off", "setlocal", "D:", f"cd {ROOT}"] + commands + ["endlocal", ""]
    (ROOT / "models" / MODEL / "scripts" / "run_M4_SYM_RATIO_on_server_32core.bat").write_text(
        "\n".join(bat), encoding="utf-8"
    )


def layer_bounds_from_design(data_path: Path) -> list[dict[str, object]]:
    box, type_elements, atoms = parse_lammps_data(data_path)
    z = box["zlo"]
    rows = []
    for region_name, phase, repeat_z in expected_ratio_layer_sequence():
        thickness = phase_thickness(phase, repeat_z)
        z_min = z
        z_max = z + thickness
        is_last_region = region_name == "Cu_half_right"
        region_atoms = [
            atom
            for atom in atoms
            if (z_min <= atom[4] < z_max) or (is_last_region and z_min <= atom[4] <= z_max)
        ]
        counts = count_elements(type_elements, region_atoms)
        total = counts["total"]
        al_fraction = counts["Al"] / total if total else ""
        cu_fraction = counts["Cu"] / total if total else ""
        if phase == "Cu":
            region_type = "Cu"
        elif phase == "Al":
            region_type = "Al"
        else:
            region_type = phase
        rows.append(
            {
                "region_name": region_name,
                "region_type": region_type,
                "z_min": z_min,
                "z_max": z_max,
                "thickness_A": thickness,
                "Al_fraction": al_fraction,
                "Cu_fraction": cu_fraction,
                "notes": "Bounds from Atomsk merge order and discrete design thickness.",
            }
        )
        z = z_max
    return rows


def write_final_box_atoms(data_path: Path, root_copy: Path, copy_note: str) -> dict[str, object]:
    box, type_elements, atoms = parse_lammps_data(data_path)
    lengths = box_lengths(box)
    counts = count_elements(type_elements, atoms)
    row = {
        "model": MODEL,
        "data_file": str(data_path),
        "atoms": counts["total"],
        "xlo": box["xlo"],
        "xhi": box["xhi"],
        "lx": lengths["lx"],
        "ylo": box["ylo"],
        "yhi": box["yhi"],
        "ly": lengths["ly"],
        "zlo": box["zlo"],
        "zhi": box["zhi"],
        "lz": lengths["lz"],
        "type1_element": "Al",
        "type2_element": "Cu",
        "Al_atoms": counts["Al"],
        "Cu_atoms": counts["Cu"],
        "Al_fraction": counts["Al"] / counts["total"],
        "Cu_fraction": counts["Cu"] / counts["total"],
        "notes": f"type 1 = Al, type 2 = Cu; root copy: {root_copy.relative_to(ROOT)}; {copy_note}",
    }
    write_csv(METADATA_DIR / "M4_SYM_RATIO_final_box_atoms.csv", [row])
    return row


def main() -> int:
    ensure_model_dirs()
    if atomsk_log().exists():
        atomsk_log().unlink()

    warnings: list[str] = []
    required = [
        ROOT / "potentials" / "CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr_Zhou04.eam.alloy",
        ROOT / "structure" / "Al4Cu9_full.cif",
        ROOT / "structure" / "Al2Cu.cif",
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        (DEBUG_DIR / "build_warnings.md").write_text(
            "# Build Warnings\n\nBuild stopped because required files are missing:\n"
            + "\n".join(f"- `{path}`" for path in missing)
            + "\n",
            encoding="utf-8",
        )
        return 2

    for old in TMP_DIR.glob("*"):
        if old.is_file():
            old.unlink()

    slabs = [
        create_slab("Cu", RATIO_REPEATS["Cu_half"], "Cu_half_z12", warnings),
        create_slab("Al4Cu9", RATIO_REPEATS["Al4Cu9"], "Al4Cu9_ratio_z2", warnings),
        create_slab("Al2Cu", RATIO_REPEATS["Al2Cu"], "Al2Cu_ratio_z10", warnings),
        create_slab("Al", RATIO_REPEATS["Al"], "Al_z24", warnings),
        create_slab("Al2Cu", RATIO_REPEATS["Al2Cu"], "Al2Cu_ratio_z10_right", warnings),
        create_slab("Al4Cu9", RATIO_REPEATS["Al4Cu9"], "Al4Cu9_ratio_z2_right", warnings),
        create_slab("Cu", RATIO_REPEATS["Cu_half"], "Cu_half_z12_right", warnings),
    ]

    raw_atomsk = TMP_DIR / "M4_SYM_RATIO_raw_atomsk.lmp"
    run_command(["atomsk", "--merge", "Z", str(len(slabs)), *[str(path) for path in slabs], str(raw_atomsk)], ROOT, atomsk_log())

    raw_data = STRUCTURE_DIR / "M4_SYM_RATIO_raw.data"
    ready_data = STRUCTURE_DIR / "M4_SYM_RATIO_ready.data"
    raw_counts = write_mapped_data(raw_atomsk, raw_data, "M4_SYM_RATIO raw; type 1 = Al, type 2 = Cu")
    ready_counts, initial_atoms, removed, iterations = remove_obvious_overlaps(raw_data, ready_data, REMOVE_DOUBLES_CUTOFF)
    removed_pct = 100.0 * removed / initial_atoms if initial_atoms else 0.0

    root_ready = ROOT / "structure" / "M4_SYM_RATIO_ready.data"
    copy_note = "standard requested filename"
    if root_ready.exists():
        root_ready = ROOT / "structure" / "M4_SYM_RATIO_ready_v2.data"
        copy_note = "M4_SYM_RATIO_ready.data already existed; copied to v2 name"
    shutil.copy2(ready_data, root_ready)
    (METADATA_DIR / "M4_SYM_RATIO_root_structure_copy.txt").write_text(
        f"model_ready_data={ready_data}\nroot_ready_data={root_ready}\n{copy_note}\n",
        encoding="utf-8",
    )

    final_row = write_final_box_atoms(ready_data, root_ready, copy_note)
    write_csv(METADATA_DIR / "M4_SYM_RATIO_layer_bounds_auto.csv", layer_bounds_from_design(ready_data))
    z_profile(
        ready_data,
        METADATA_DIR / "M4_SYM_RATIO_z_profile.csv",
        PREVIEW_DIR / "M4_SYM_RATIO_z_profile.png",
        bins=300,
    )
    min_summary = min_distance_summary(ready_data)
    write_csv(METADATA_DIR / "M4_SYM_RATIO_min_distance_check.csv", [min_summary])

    ref_box, _, _ = parse_lammps_data(ROOT / "structure" / "M4_SYM_Cu_Al4Cu9_Al2Cu_Al_PPP_ready.data")
    ref_lengths = box_lengths(ref_box)
    lx_diff = final_row["lx"] - ref_lengths["lx"]
    ly_diff = final_row["ly"] - ref_lengths["ly"]
    if abs(lx_diff) > 0.05 or abs(ly_diff) > 0.05:
        warnings.append(f"x/y dimensions differ from M4_SYM by lx={lx_diff:.6f} A, ly={ly_diff:.6f} A.")
    if not (ROOT / "structure" / "Al4Cu9_full.cif").exists():
        warnings.append("Al4Cu9_full.cif was not found after build.")
    log_text = atomsk_log().read_text(encoding="utf-8", errors="replace").lower()
    if "partial occupancy" in log_text:
        warnings.append("Atomsk partial occupancy warning detected in build log.")
    if min_summary["pairs_below_1p0_A"]:
        warnings.append(f"Remaining pairs below 1.0 A: {min_summary['pairs_below_1p0_A']}.")
    if removed_pct > 0.5:
        warnings.append(f"Overlap removal deleted {removed_pct:.4f}% atoms, above the 0.5% watch threshold.")
    if not warnings:
        warnings.append("No blocking build warnings. Note: pairs below 1.5 A may occur at abrupt merged interfaces; inspect if needed.")

    build_report = [
        "# Build Warnings",
        "",
        f"- Raw atoms: {raw_counts['total']}",
        f"- Ready atoms: {ready_counts['total']}",
        f"- Removed atoms by {REMOVE_DOUBLES_CUTOFF:.2f} A overlap cleanup: {removed} ({removed_pct:.6f}%), iterations={iterations}",
        f"- Minimum distance after cleanup: {float(min_summary['min_distance_A']):.6f} A",
        f"- Pairs below 1.0 A after cleanup: {min_summary['pairs_below_1p0_A']}",
        f"- Pairs below 1.5 A after cleanup: {min_summary['pairs_below_1p5_A']}",
        f"- Root ready data copy: `{root_ready.relative_to(ROOT)}`",
        "",
        "## Notes / warnings",
    ]
    build_report.extend(f"- {item}" for item in warnings)
    (DEBUG_DIR / "build_warnings.md").write_text("\n".join(build_report) + "\n", encoding="utf-8")

    write_input_files(root_ready.name)
    write_server_command_files()

    print(raw_data)
    print(ready_data)
    print(root_ready)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
