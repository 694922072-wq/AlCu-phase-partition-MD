#!/usr/bin/env python3
"""Extract formal per-trajectory and n=3 mechanical metrics from raw data.

F50_DEFINITION_LOCK = PAPER_RAW_FIRST_ROW. F50 is the first unmodified raw
post-peak row at or below 50% of that trajectory's own raw peak. No smoothing,
interpolation, or extrapolation is performed for F50.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from f50_contract import (
    F50_DEFINITION_LOCK,
    FORMAL_F50_FIELD,
    MODEL_ORDER,
    csv_value,
    summarize_formal_f50,
    write_rows,
)

REQUIRED = ("step", "strain_percent", "stress_xx_GPa")


def read_curve(path: Path) -> list[tuple[int, float, float]]:
    """Read raw rows in file order without sorting, filtering, or smoothing."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = handle.readline().strip().split()
        missing = [name for name in REQUIRED if name not in header]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}; header={header}")
        indices = {name: header.index(name) for name in REQUIRED}
        rows: list[tuple[int, float, float]] = []
        for line_number, line in enumerate(handle, start=2):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            try:
                row = (
                    int(float(parts[indices["step"]])),
                    float(parts[indices["strain_percent"]]),
                    float(parts[indices["stress_xx_GPa"]]),
                )
            except (ValueError, IndexError) as exc:
                raise ValueError(f"{path}:{line_number}: malformed numeric row") from exc
            if not all(math.isfinite(value) for value in row):
                raise ValueError(f"{path}:{line_number}: non-finite raw value")
            rows.append(row)
    if len(rows) < 2:
        raise ValueError(f"{path}: insufficient raw rows")
    return rows


def value_at_strain(rows: list[tuple[int, float, float]], target: float) -> float:
    """Linear interpolation for sigma15 only; this function is never used for F50."""
    for index in range(1, len(rows)):
        x1, x2 = rows[index - 1][1], rows[index][1]
        if (x1 <= target <= x2) or (x2 <= target <= x1):
            y1, y2 = rows[index - 1][2], rows[index][2]
            return y2 if x2 == x1 else y1 + (target - x1) * (y2 - y1) / (x2 - x1)
    raise ValueError(f"target strain {target} is not bracketed")


def extract(rows: list[tuple[int, float, float]]) -> dict[str, object]:
    steps = [row[0] for row in rows]
    strains = [row[1] for row in rows]
    stresses = [row[2] for row in rows]
    if steps[0] != 0 or steps[-1] != 250000:
        raise ValueError(f"incomplete step range: {steps[0]}..{steps[-1]}")
    if abs(strains[-1] - 25.0) > 1e-6:
        raise ValueError(f"final strain is {strains[-1]}, expected 25%")
    if any(steps[index] <= steps[index - 1] for index in range(1, len(steps))):
        raise ValueError("raw steps are not strictly increasing")
    if any(strains[index] < strains[index - 1] for index in range(1, len(strains))):
        raise ValueError("raw strain is not monotonically nondecreasing")

    peak_index = max(range(len(rows)), key=lambda index: stresses[index])
    peak_step, peak_strain, peak_stress = rows[peak_index]
    threshold = 0.5 * peak_stress
    f50_row: tuple[int, float, float] | None = None
    for row in rows[peak_index + 1 :]:
        if row[2] <= threshold:
            f50_row = row
            break

    return {
        "n_rows": len(rows),
        "first_step": steps[0],
        "last_step": steps[-1],
        "final_strain_percent": strains[-1],
        "sigma15_GPa": value_at_strain(rows, 15.0),
        "peak_stress_GPa": peak_stress,
        "peak_strain_percent": peak_strain,
        "peak_step": peak_step,
        "F50_threshold_GPa": threshold,
        FORMAL_F50_FIELD: None if f50_row is None else f50_row[1],
        "F50_raw_first_row_step": None if f50_row is None else f50_row[0],
        "F50_raw_first_row_stress_GPa": None if f50_row is None else f50_row[2],
        "F50_status": "NA_no_crossing_by_25pct" if f50_row is None else "FOUND_raw_first_post_peak_row",
        "F50_definition_lock": F50_DEFINITION_LOCK,
    }


def trajectory_path(project_root: Path, manifest_row: dict[str, str]) -> Path:
    if "data_path" in manifest_row and manifest_row["data_path"].strip():
        return project_root / manifest_row["data_path"]
    required = ("job_relpath", "outprefix")
    missing = [name for name in required if not manifest_row.get(name, "").strip()]
    if missing:
        raise ValueError(f"manifest row lacks {missing}: {manifest_row}")
    return project_root / manifest_row["job_relpath"] / "outputs" / f"{manifest_row['outprefix']}_stress_strain.dat"


def metric_summary(per_trajectory: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in per_trajectory:
        grouped[str(row["model"])].append(row)

    f50_rows = [
        {"model": row["model"], "replica": row["replica"], FORMAL_F50_FIELD: row[FORMAL_F50_FIELD]}
        for row in per_trajectory
    ]
    f50_by_model = {row["model"]: row for row in summarize_formal_f50(f50_rows)}
    summary: list[dict[str, object]] = []
    for model in MODEL_ORDER:
        group = grouped[model]
        if len(group) != 3:
            raise ValueError(f"{model}: expected exactly three trajectories, got {len(group)}")
        item: dict[str, object] = {"model": model, "n_trajectories": 3}
        for metric in ("sigma15_GPa", "peak_stress_GPa", "peak_strain_percent"):
            values = [float(row[metric]) for row in group]
            item[f"{metric}_mean"] = statistics.mean(values)
            item[f"{metric}_sample_SD"] = statistics.stdev(values)
            item[f"{metric}_min"] = min(values)
            item[f"{metric}_max"] = max(values)
        item.update({key: value for key, value in f50_by_model[model].items() if key != "model"})
        item["F50_definition_lock"] = F50_DEFINITION_LOCK
        summary.append(item)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    project_default = Path(__file__).resolve().parents[1]
    parser.add_argument("--project-root", type=Path, default=project_default)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=project_default / "05_active_all_new_9_v2" / "seed_manifest_all_new_9_v2.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=project_default / "06_results_n3")
    parser.add_argument("--mode", default="all_new_9_v2")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    manifest = args.manifest.resolve()
    output_dir = args.output_dir.resolve()
    per_trajectory: list[dict[str, object]] = []
    missing: list[Path] = []
    with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        for manifest_row in csv.DictReader(handle):
            path = trajectory_path(project_root, manifest_row)
            if not path.exists():
                missing.append(path)
                continue
            metrics = extract(read_curve(path))
            identity = {
                "model": manifest_row["model"],
                "replica": manifest_row["replica"],
                "seed": manifest_row["seed"],
                "raw_data_path": str(path),
            }
            per_trajectory.append({**identity, **metrics})

    if missing:
        print("Missing formal trajectory files; no summary was generated:", file=sys.stderr)
        for path in missing:
            print(f"  - {path}", file=sys.stderr)
        return 2
    if len(per_trajectory) != 9:
        raise ValueError(f"expected 9 formal trajectories, got {len(per_trajectory)}")

    per_path = output_dir / f"trajectory_metrics_{args.mode}.csv"
    summary_path = output_dir / f"model_summary_{args.mode}.csv"
    write_rows(per_path, per_trajectory)
    write_rows(summary_path, metric_summary(per_trajectory))
    print(per_path)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

