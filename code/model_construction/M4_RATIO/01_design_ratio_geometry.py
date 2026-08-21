from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _m4_sym_ratio_utils import (  # noqa: E402
    DEBUG_DIR,
    METADATA_DIR,
    MODEL,
    PHASES,
    RATIO_REPEATS,
    REF_REPEATS,
    ROOT,
    ensure_model_dirs,
    phase_thickness,
    write_csv,
)


def main() -> int:
    ensure_model_dirs()

    t_ref_al4 = phase_thickness("Al4Cu9", REF_REPEATS["Al4Cu9"])
    t_ref_al2 = phase_thickness("Al2Cu", REF_REPEATS["Al2Cu"])
    t_ref_imc = t_ref_al4 + t_ref_al2
    target_al4 = 0.5 * t_ref_al4

    best = None
    for al4_repeat in range(1, REF_REPEATS["Al4Cu9"] + 1):
        t_al4 = phase_thickness("Al4Cu9", al4_repeat)
        if t_al4 >= t_ref_al4:
            continue
        target_al2 = t_ref_imc - t_al4
        for al2_repeat in range(1, 30):
            t_al2 = phase_thickness("Al2Cu", al2_repeat)
            t_imc = t_al4 + t_al2
            score = (abs(t_al4 - target_al4), abs(t_imc - t_ref_imc), abs(t_al2 - target_al2))
            candidate = (score, al4_repeat, al2_repeat, t_al4, t_al2, t_imc)
            if best is None or candidate[0] < best[0]:
                best = candidate

    if best is None:
        raise RuntimeError("No valid ratio geometry candidate found.")

    _, al4_repeat, al2_repeat, t_al4, t_al2, t_imc = best
    RATIO_REPEATS["Al4Cu9"] = al4_repeat
    RATIO_REPEATS["Al2Cu"] = al2_repeat

    rows = [
        {
            "model": MODEL,
            "region_name": "Cu_half_left",
            "target_thickness_A": phase_thickness("Cu", REF_REPEATS["Cu_half"]),
            "actual_thickness_A": phase_thickness("Cu", RATIO_REPEATS["Cu_half"]),
            "source_phase": "Cu",
            "source_file": "atomsk --create fcc 3.615 Cu",
            "z_repeat_or_cut": RATIO_REPEATS["Cu_half"],
            "thickness_error_A": 0.0,
            "notes": "Same as M4_SYM Cu_half.",
        },
        {
            "model": MODEL,
            "region_name": "Al4Cu9_left_thin",
            "target_thickness_A": target_al4,
            "actual_thickness_A": t_al4,
            "source_phase": "Al4Cu9",
            "source_file": "structure/Al4Cu9_full.cif",
            "z_repeat_or_cut": al4_repeat,
            "thickness_error_A": t_al4 - target_al4,
            "notes": "Al4Cu9 repeat reduced from M4_SYM z4 to z2.",
        },
        {
            "model": MODEL,
            "region_name": "Al2Cu_left_adjusted",
            "target_thickness_A": t_ref_imc - t_al4,
            "actual_thickness_A": t_al2,
            "source_phase": "Al2Cu",
            "source_file": "structure/Al2Cu.cif",
            "z_repeat_or_cut": al2_repeat,
            "thickness_error_A": t_al2 - (t_ref_imc - t_al4),
            "notes": "Adjusted to keep single-side IMC total close to M4_SYM.",
        },
        {
            "model": MODEL,
            "region_name": "Al_center",
            "target_thickness_A": phase_thickness("Al", REF_REPEATS["Al"]),
            "actual_thickness_A": phase_thickness("Al", RATIO_REPEATS["Al"]),
            "source_phase": "Al",
            "source_file": "atomsk --create fcc 4.05 Al",
            "z_repeat_or_cut": RATIO_REPEATS["Al"],
            "thickness_error_A": 0.0,
            "notes": "Same as M4_SYM Al layer.",
        },
        {
            "model": MODEL,
            "region_name": "Al2Cu_right_adjusted",
            "target_thickness_A": t_ref_imc - t_al4,
            "actual_thickness_A": t_al2,
            "source_phase": "Al2Cu",
            "source_file": "structure/Al2Cu.cif",
            "z_repeat_or_cut": al2_repeat,
            "thickness_error_A": t_al2 - (t_ref_imc - t_al4),
            "notes": "Symmetric counterpart of left Al2Cu_adjusted.",
        },
        {
            "model": MODEL,
            "region_name": "Al4Cu9_right_thin",
            "target_thickness_A": target_al4,
            "actual_thickness_A": t_al4,
            "source_phase": "Al4Cu9",
            "source_file": "structure/Al4Cu9_full.cif",
            "z_repeat_or_cut": al4_repeat,
            "thickness_error_A": t_al4 - target_al4,
            "notes": "Symmetric counterpart of left Al4Cu9_thin.",
        },
        {
            "model": MODEL,
            "region_name": "Cu_half_right",
            "target_thickness_A": phase_thickness("Cu", REF_REPEATS["Cu_half"]),
            "actual_thickness_A": phase_thickness("Cu", RATIO_REPEATS["Cu_half"]),
            "source_phase": "Cu",
            "source_file": "atomsk --create fcc 3.615 Cu",
            "z_repeat_or_cut": RATIO_REPEATS["Cu_half"],
            "thickness_error_A": 0.0,
            "notes": "Same as M4_SYM Cu_half.",
        },
    ]
    design_csv = METADATA_DIR / "M4_SYM_RATIO_geometry_design.csv"
    write_csv(design_csv, rows)

    imc_error = t_imc - t_ref_imc
    imc_error_pct = 100.0 * imc_error / t_ref_imc
    al4_fraction = t_al4 / t_imc
    warning = ""
    if abs(imc_error) > 3.0 and abs(imc_error_pct) > 5.0:
        warning = "\nWARNING: single-side IMC thickness error exceeds both 3 A and 5%.\n"

    notes = [
        "# M4_SYM_RATIO Design Notes",
        "",
        "## M4_SYM reference",
        "",
        f"- Al4Cu9 thickness: {t_ref_al4:.8f} A (z repeat {REF_REPEATS['Al4Cu9']})",
        f"- Al2Cu thickness: {t_ref_al2:.8f} A (z repeat {REF_REPEATS['Al2Cu']})",
        f"- Single-side IMC total: {t_ref_imc:.8f} A",
        "",
        "## M4_SYM_RATIO",
        "",
        f"- Al4Cu9 thickness: {t_al4:.8f} A (z repeat {al4_repeat})",
        f"- Al2Cu thickness: {t_al2:.8f} A (z repeat {al2_repeat})",
        f"- Single-side IMC total: {t_imc:.8f} A",
        f"- Single-side IMC error vs M4_SYM: {imc_error:+.8f} A ({imc_error_pct:+.4f}%)",
        f"- Al4Cu9 share of single-side IMC: {100.0 * al4_fraction:.4f}%",
        "",
        "## Difference from M4_SYM",
        "",
        "- Al4Cu9 is thinned from z4 to z2, exactly half of the M4_SYM Al4Cu9 thickness.",
        "- Al2Cu is increased from z6 to z10 to preserve the single-side IMC total.",
        "- Cu_half, Al center, x/y target size, type mapping, and clean LAMMPS settings are inherited from M4_SYM.",
        warning,
    ]
    (DEBUG_DIR / "M4_SYM_RATIO_design_notes.md").write_text("\n".join(notes).strip() + "\n", encoding="utf-8")
    print(design_csv)
    print(DEBUG_DIR / "M4_SYM_RATIO_design_notes.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

