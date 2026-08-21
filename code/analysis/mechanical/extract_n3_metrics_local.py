#!/usr/bin/env python3
"""Independent local re-extraction of the nine authoritative Al-Cu trajectories.

Only paths listed in AUTHORITATIVE_NINE_TRAJECTORIES.csv are read.  No glob or
recursive stress-file discovery is used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import statistics
import sys
from pathlib import Path


EXPECTED = [
    ("M3_SYM", "R1", "3847291"),
    ("M3_SYM", "R2", "28471"),
    ("M3_SYM", "R3", "39581"),
    ("M4_RATIO", "R1", "4928459"),
    ("M4_RATIO", "R2", "41777"),
    ("M4_RATIO", "R3", "52891"),
    ("M4_SYM", "R1", "9520510"),
    ("M4_SYM", "R2", "63997"),
    ("M4_SYM", "R3", "74131"),
]

METRICS = (
    "sigma15_GPa",
    "peak_stress_GPa",
    "peak_strain_percent",
    "F50_raw_first_row_percent",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def map_archived_path(transfer_root: Path, archived_path: str) -> Path:
    normalized = archived_path.replace("/", "\\")
    prefix = "D:\\AlCu_n3_128core_20260801\\"
    if not normalized.lower().startswith(prefix.lower()):
        raise ValueError(f"Unexpected archived path prefix: {archived_path}")
    relative = normalized[len(prefix) :]
    return transfer_root.joinpath(*relative.split("\\"))


def read_raw(path: Path) -> tuple[list[str], list[list[float]]]:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip().split()
        rows = [[float(value) for value in line.split()] for line in handle if line.strip()]
    required = ["step", "strain_percent", "stress_xx_GPa"]
    missing = [name for name in required if name not in header]
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")
    if not rows:
        raise ValueError(f"No raw rows in {path}")
    if any(len(row) != len(header) for row in rows):
        raise ValueError(f"Column-count mismatch in {path}")
    return header, rows


def linear_at_target(xs: list[float], ys: list[float], target: float) -> float:
    for index, value in enumerate(xs):
        if value == target:
            return ys[index]
        if value > target:
            if index == 0:
                raise ValueError(f"Target {target} is below the raw range")
            x0, x1 = xs[index - 1], value
            y0, y1 = ys[index - 1], ys[index]
            if not (x0 < target < x1):
                raise ValueError(f"Target {target} is not bracketed monotonically")
            return y0 + (target - x0) * (y1 - y0) / (x1 - x0)
    raise ValueError(f"Target {target} is above the raw range")


def extract_metrics(path: Path) -> dict[str, object]:
    header, rows = read_raw(path)
    col = {name: index for index, name in enumerate(header)}
    strain = [row[col["strain_percent"]] for row in rows]
    stress = [row[col["stress_xx_GPa"]] for row in rows]
    if any(b <= a for a, b in zip(strain, strain[1:])):
        raise ValueError(f"strain_percent is not strictly increasing in {path}")

    sigma15 = linear_at_target(strain, stress, 15.0)
    peak_index = max(range(len(stress)), key=stress.__getitem__)
    peak_stress = stress[peak_index]
    peak_strain = strain[peak_index]
    threshold = 0.5 * peak_stress
    f50_index = next(
        (index for index in range(peak_index + 1, len(stress)) if stress[index] <= threshold),
        None,
    )
    if f50_index is None:
        raise ValueError(f"No raw-first-row F50 event after the peak in {path}")

    return {
        "row_count": len(rows),
        "sigma15_GPa": sigma15,
        "peak_stress_GPa": peak_stress,
        "peak_strain_percent": peak_strain,
        "F50_raw_first_row_percent": strain[f50_index],
        "final_strain_percent": strain[-1],
        "peak_row_index_zero_based": peak_index,
        "F50_row_index_zero_based": f50_index,
    }


def summarize(local_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for model in ("M3_SYM", "M4_RATIO", "M4_SYM"):
        selected = [row for row in local_rows if row["model"] == model]
        if len(selected) != 3:
            raise ValueError(f"Expected n=3 for {model}, found {len(selected)}")
        result: dict[str, object] = {"model": model, "n": 3, "ddof": 1}
        for metric in METRICS:
            values = [float(row[metric]) for row in selected]
            mean = statistics.mean(values)
            sd = statistics.stdev(values)
            result[f"{metric}_mean"] = mean
            result[f"{metric}_sample_SD"] = sd
            result[f"{metric}_min"] = min(values)
            result[f"{metric}_max"] = max(values)
            result[f"{metric}_CV"] = sd / mean
        output.append(result)
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--transfer-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    transfer_root = args.transfer_root.resolve()
    output_dir = args.output_dir.resolve()
    acceptance = transfer_root / "11_final_acceptance"
    authoritative_path = acceptance / "reports" / "AUTHORITATIVE_NINE_TRAJECTORIES.csv"
    result_manifest_path = acceptance / "reports" / "RESULT_FILE_MANIFEST_SHA256.csv"
    compute_per_path = acceptance / "reference_metrics" / "metrics_per_trajectory_compute_pc.csv"
    compute_summary_path = acceptance / "reference_metrics" / "n3_summary_compute_pc.csv"

    authoritative = read_csv(authoritative_path)
    actual_keys = [(row["model"], row["replica"], row["seed"]) for row in authoritative]
    if actual_keys != EXPECTED:
        raise ValueError(f"Authoritative manifest rows do not match the locked nine: {actual_keys}")
    if any(row["selection"] != "AUTHORITATIVE_FINAL" for row in authoritative):
        raise ValueError("One or more manifest rows are not AUTHORITATIVE_FINAL")

    result_manifest = read_csv(result_manifest_path)
    raw_manifest = {
        (row["model"], row["replica"], row["seed"]): row
        for row in result_manifest
        if row["file_role"] == "raw_stress_strain"
    }
    compute_per = {
        (row["model"], row["replica"], row["seed"]): row for row in read_csv(compute_per_path)
    }
    if set(raw_manifest) != set(EXPECTED) or set(compute_per) != set(EXPECTED):
        raise ValueError("Raw manifest or compute-PC table does not contain exactly the locked nine rows")

    local_rows: list[dict[str, object]] = []
    compare_rows: list[dict[str, object]] = []
    all_pass = True

    for row in authoritative:
        key = (row["model"], row["replica"], row["seed"])
        raw_path = map_archived_path(transfer_root, row["raw_path"])
        if not raw_path.is_file():
            raise FileNotFoundError(raw_path)
        digest = sha256(raw_path)
        manifest_row = raw_manifest[key]
        compute_row = compute_per[key]
        hash_pass = digest.lower() == manifest_row["sha256"].lower() == compute_row["raw_sha256"].lower()
        size_pass = raw_path.stat().st_size == int(manifest_row["length"])
        eligible_pass = manifest_row["formal_statistics_eligible"].lower() == "true"

        metrics = extract_metrics(raw_path)
        completion_pass = (
            int(metrics["row_count"]) == int(compute_row["row_count"])
            and float(metrics["final_strain_percent"]) >= 24.99
            and compute_row["completion_status"] == "PASS"
        )
        status = hash_pass and size_pass and eligible_pass and completion_pass
        all_pass = all_pass and status
        local = {
            "model": row["model"],
            "replica": row["replica"],
            "seed": row["seed"],
            "authoritative_raw_path_archived": row["raw_path"],
            "local_raw_path": str(raw_path),
            "raw_sha256": digest,
            **metrics,
            "raw_hash_match": hash_pass,
            "manifest_size_match": size_pass,
            "formal_statistics_eligible": eligible_pass,
            "completion_status": "PASS" if completion_pass else "FAIL",
        }
        local_rows.append(local)

        for metric in METRICS:
            local_value = float(metrics[metric])
            compute_value = float(compute_row[metric])
            delta = abs(local_value - compute_value)
            passed = delta <= 1.0e-9
            all_pass = all_pass and passed
            compare_rows.append(
                {
                    "scope": "per_trajectory",
                    "model": row["model"],
                    "replica": row["replica"],
                    "seed": row["seed"],
                    "metric": metric,
                    "compute_pc_value": compute_value,
                    "local_value": local_value,
                    "absolute_difference": delta,
                    "tolerance": 1.0e-9,
                    "status": "PASS" if passed else "FAIL",
                }
            )

    summary_rows = summarize(local_rows)
    compute_summary = {row["model"]: row for row in read_csv(compute_summary_path)}
    for summary in summary_rows:
        model = str(summary["model"])
        reference = compute_summary[model]
        for metric in METRICS:
            for suffix in ("mean", "sample_SD", "min", "max", "CV"):
                column = f"{metric}_{suffix}"
                local_value = float(summary[column])
                compute_value = float(reference[column])
                delta = abs(local_value - compute_value)
                passed = delta <= 1.0e-6
                all_pass = all_pass and passed
                compare_rows.append(
                    {
                        "scope": "n3_summary",
                        "model": model,
                        "replica": "",
                        "seed": "",
                        "metric": column,
                        "compute_pc_value": compute_value,
                        "local_value": local_value,
                        "absolute_difference": delta,
                        "tolerance": 1.0e-6,
                        "status": "PASS" if passed else "FAIL",
                    }
                )

    local_fields = [
        "model", "replica", "seed", "authoritative_raw_path_archived", "local_raw_path",
        "raw_sha256", "row_count", *METRICS, "final_strain_percent",
        "peak_row_index_zero_based", "F50_row_index_zero_based", "raw_hash_match",
        "manifest_size_match", "formal_statistics_eligible", "completion_status",
    ]
    summary_fields = ["model", "n", "ddof"] + [
        f"{metric}_{suffix}"
        for metric in METRICS
        for suffix in ("mean", "sample_SD", "min", "max", "CV")
    ]
    compare_fields = [
        "scope", "model", "replica", "seed", "metric", "compute_pc_value",
        "local_value", "absolute_difference", "tolerance", "status",
    ]
    write_csv(output_dir / "metrics_per_trajectory_local.csv", local_rows, local_fields)
    write_csv(output_dir / "n3_summary_local.csv", summary_rows, summary_fields)
    write_csv(output_dir / "compute_vs_local.csv", compare_rows, compare_fields)

    print(f"LOCAL_EXTRACTION_STATUS={'PASS' if all_pass else 'FAIL'}")
    print(f"AUTHORITATIVE_ROWS={len(authoritative)}")
    print(f"COMPARE_ROWS={len(compare_rows)}")
    print(f"OUTPUT_DIR={output_dir}")
    for summary in summary_rows:
        print(
            f"{summary['model']} sigma15={summary['sigma15_GPa_mean']!r} "
            f"sd={summary['sigma15_GPa_sample_SD']!r} "
            f"F50={summary['F50_raw_first_row_percent_mean']!r} "
            f"sd={summary['F50_raw_first_row_percent_sample_SD']!r}"
        )
    return 0 if all_pass else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"LOCAL_EXTRACTION_ERROR={type(exc).__name__}: {exc}", file=sys.stderr)
        raise
