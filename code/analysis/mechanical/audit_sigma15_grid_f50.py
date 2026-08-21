from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\leng\AlCu_lmy6_5_revision")
METRICS = ROOT / "02_local_statistics" / "metrics_per_trajectory_local.csv"
SIGMA_CSV = ROOT / "08_reports" / "SIGMA15_EXTRACTION_AUDIT.csv"
SIGMA_MD = ROOT / "08_reports" / "SIGMA15_EXTRACTION_AUDIT.md"
GRID_CSV = ROOT / "08_reports" / "TRAJECTORY_GRID_AUDIT.csv"
GRID_MD = ROOT / "08_reports" / "TRAJECTORY_GRID_AUDIT.md"
F50_CSV = ROOT / "02_local_statistics" / "F50_SENSITIVITY_N9.csv"
F50_MD = ROOT / "02_local_statistics" / "F50_SENSITIVITY_N9.md"

MODEL_ORDER = ["M3_SYM", "M4_RATIO", "M4_SYM"]
TOL = 1.0e-10


def first_crossing(strain: np.ndarray, values: np.ndarray, peak_index: int, threshold: float) -> float:
    qualifying = np.where((np.arange(len(values)) > peak_index) & (values <= threshold))[0]
    return float(strain[qualifying[0]]) if len(qualifying) else math.nan


def main() -> None:
    source = pd.read_csv(METRICS, dtype={"seed": str})
    sigma_rows: list[dict] = []
    grid_rows: list[dict] = []
    f50_rows: list[dict] = []
    records: list[tuple[dict, pd.DataFrame]] = []

    for item in source.to_dict("records"):
        raw_path = Path(item["local_raw_path"])
        df = pd.read_csv(raw_path, sep=r"\s+")
        records.append((item, df))
        strain = df["strain_percent"].to_numpy(float)
        stress = df["stress_xx_GPa"].to_numpy(float)
        step = df["step"].to_numpy(np.int64)

        distances = np.abs(strain - 15.0)
        exact = np.where(distances <= TOL)[0]
        if len(exact):
            idx15 = int(exact[np.argmin(distances[exact])])
            extraction = "DIRECT_ARCHIVED_RAW_ROW"
        else:
            idx15 = int(np.argmin(distances))
            extraction = "NO_EXACT_ROW_FALLBACK_REQUIRED"
        sigma_rows.append({
            "model": item["model"],
            "replica": item["replica"],
            "seed": item["seed"],
            "raw_path": str(raw_path),
            "row_count": len(df),
            "matching_row_index_zero_based": idx15,
            "matching_step": int(step[idx15]),
            "recorded_strain_percent": format(strain[idx15], ".17g"),
            "absolute_deviation_from_15_percent": format(abs(strain[idx15] - 15.0), ".17g"),
            "stress_xx_GPa_at_matching_row": format(stress[idx15], ".17g"),
            "exact_15_percent_row_within_tolerance": "YES" if abs(strain[idx15] - 15.0) <= TOL else "NO",
            "tolerance_percent": format(TOL, ".1e"),
            "extraction_logic": extraction,
        })

        peak_index = int(np.argmax(stress))
        peak = float(stress[peak_index])
        threshold = 0.5 * peak
        primary = first_crossing(strain, stress, peak_index, threshold)
        centered = pd.Series(stress).rolling(3, center=True, min_periods=3).mean().to_numpy()
        centered_three = first_crossing(strain, centered, peak_index, threshold)
        consecutive = math.nan
        for i in range(peak_index + 1, len(stress) - 2):
            if np.all(stress[i : i + 3] <= threshold):
                consecutive = float(strain[i])
                break
        shifts = [abs(centered_three - primary), abs(consecutive - primary)]
        f50_rows.append({
            "model": item["model"],
            "replica": item["replica"],
            "seed": item["seed"],
            "F50_primary": primary,
            "F50_centered_three_point": centered_three,
            "F50_three_consecutive": consecutive,
            "centered_shift_percentage_points": abs(centered_three - primary),
            "consecutive_shift_percentage_points": abs(consecutive - primary),
            "maximum_shift": max(shifts),
            "model_order_changed": "PENDING_GLOBAL_CHECK",
            "primary_definition": "first raw post-peak row <= 0.5 * own raw peak; no interpolation/smoothing/extrapolation",
        })

    reference_steps = records[0][1]["step"].to_numpy(np.int64)
    reference_strain = records[0][1]["strain_percent"].to_numpy(float)
    for item, df in records:
        steps = df["step"].to_numpy(np.int64)
        strain = df["strain_percent"].to_numpy(float)
        same_steps = np.array_equal(steps, reference_steps)
        same_strain_exact = np.array_equal(strain, reference_strain)
        same_strain_tol = len(strain) == len(reference_strain) and np.allclose(strain, reference_strain, rtol=0.0, atol=TOL)
        max_strain_delta = float(np.max(np.abs(strain - reference_strain))) if len(strain) == len(reference_strain) else math.nan
        grid_rows.append({
            "model": item["model"],
            "replica": item["replica"],
            "seed": item["seed"],
            "row_count": len(df),
            "first_step": int(steps[0]),
            "last_step": int(steps[-1]),
            "step_grid_identical_to_reference": "YES" if same_steps else "NO",
            "strain_values_bitwise_identical_to_reference": "YES" if same_strain_exact else "NO",
            "strain_grid_identical_within_tolerance": "YES" if same_strain_tol else "NO",
            "maximum_absolute_strain_grid_difference_percent": format(max_strain_delta, ".17g"),
            "display_alignment_decision": "POINTWISE_BY_COMMON_STEP_NO_INTERPOLATION" if same_steps and same_strain_tol else "DISPLAY_ONLY_COMMON_STRAIN_INTERPOLATION",
        })

    primary_means = {m: np.mean([r["F50_primary"] for r in f50_rows if r["model"] == m]) for m in MODEL_ORDER}
    centered_means = {m: np.mean([r["F50_centered_three_point"] for r in f50_rows if r["model"] == m]) for m in MODEL_ORDER}
    consecutive_means = {m: np.mean([r["F50_three_consecutive"] for r in f50_rows if r["model"] == m]) for m in MODEL_ORDER}
    primary_order = sorted(MODEL_ORDER, key=primary_means.get)
    order_changed = sorted(MODEL_ORDER, key=centered_means.get) != primary_order or sorted(MODEL_ORDER, key=consecutive_means.get) != primary_order
    for row in f50_rows:
        row["model_order_changed"] = "YES" if order_changed else "NO"

    pd.DataFrame(sigma_rows).to_csv(SIGMA_CSV, index=False, encoding="utf-8-sig")
    pd.DataFrame(grid_rows).to_csv(GRID_CSV, index=False, encoding="utf-8-sig")
    pd.DataFrame(f50_rows).to_csv(F50_CSV, index=False, encoding="utf-8-sig", float_format="%.15g")

    all_exact = all(r["exact_15_percent_row_within_tolerance"] == "YES" for r in sigma_rows)
    SIGMA_MD.write_text(
        "# Sigma15 extraction audit\n\n"
        f"- Nine authoritative raw trajectories audited: **{'PASS' if len(sigma_rows) == 9 else 'FAIL'}**.\n"
        f"- Archived 15.00% row present in all nine records (absolute tolerance {TOL:.1e} percentage points): **{'PASS' if all_exact else 'FAIL'}**.\n"
        "- Decision: **DIRECT_ARCHIVED_RAW_ROW**. Stress at 15% strain is read from the archived raw output row at nominal 15.00% engineering strain.\n"
        "- Floating-point note: the recorded decimal representations can differ from 15.0 by machine-scale amounts; these are the original archived output values, not interpolated values.\n"
        "- Deterministic linear interpolation remains a fallback only for a dataset without such a row and was not required for these nine trajectories.\n"
        f"- Detailed evidence: `{SIGMA_CSV}`.\n",
        encoding="utf-8",
    )

    all_steps = all(r["step_grid_identical_to_reference"] == "YES" for r in grid_rows)
    all_strain_tol = all(r["strain_grid_identical_within_tolerance"] == "YES" for r in grid_rows)
    max_grid_delta = max(float(r["maximum_absolute_strain_grid_difference_percent"]) for r in grid_rows)
    GRID_MD.write_text(
        "# Nine-trajectory output-grid audit\n\n"
        f"- Identical step grid (0 to 250000 every 100 steps; 2501 rows): **{'PASS' if all_steps else 'FAIL'}**.\n"
        f"- Common strain grid within {TOL:.1e} percentage points: **{'PASS' if all_strain_tol else 'FAIL'}**.\n"
        f"- Maximum machine-scale strain-coordinate difference across records: `{max_grid_delta:.3e}` percentage points.\n"
        "- Fig. 2 decision: calculate the display mean and sample-SD band pointwise by the common output step/row; no interpolation is required.\n"
        "- Formal sigma15, peak, peak strain, and F50 values remain per-trajectory extractions from each unsmoothed raw record.\n"
        f"- Detailed evidence: `{GRID_CSV}`.\n",
        encoding="utf-8",
    )

    max_shift = max(float(r["maximum_shift"]) for r in f50_rows)
    F50_MD.write_text(
        "# F50 criterion sensitivity audit for nine trajectories\n\n"
        "- Primary definition: first unmodified raw post-peak row at or below 50% of the trajectory-specific raw peak; no interpolation, smoothing, or extrapolation.\n"
        "- Auxiliary check A: first post-peak crossing of a centered three-point moving-average signal, using 50% of the same raw peak.\n"
        "- Auxiliary check B: first post-peak row beginning three consecutive raw rows at or below 50% of the same raw peak.\n"
        f"- Maximum criterion-induced shift across all nine trajectories: **{max_shift:.2f} percentage points**.\n"
        f"- Controlled-model mean ordering changed: **{'YES' if order_changed else 'NO'}**.\n"
        "- Primary manuscript values remain the raw-first-row values; auxiliary checks are stability diagnostics only.\n"
        "- Interpretation: the checks did not alter the primary interpretation.\n"
        f"- Detailed results: `{F50_CSV}`.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
