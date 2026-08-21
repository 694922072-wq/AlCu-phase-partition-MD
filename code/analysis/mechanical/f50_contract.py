#!/usr/bin/env python3
"""Shared F50 data contract for formal n=3 summaries and figure updates."""
from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

F50_DEFINITION_LOCK = "PAPER_RAW_FIRST_ROW"
FORMAL_F50_FIELD = "F50_raw_first_row_percent"
FORBIDDEN_FORMAL_F50_FIELD = "F50_interpolated_percent"
MODEL_ORDER = ("M3_SYM", "M4_RATIO", "M4_SYM")


def parse_optional_float(value: object) -> float | None:
    text = "" if value is None else str(value).strip()
    if not text or text.upper() == "NA":
        return None
    number = float(text)
    if not math.isfinite(number):
        raise ValueError(f"non-finite F50 value: {value!r}")
    return number


def read_formal_f50_rows(path: Path) -> list[dict[str, object]]:
    """Read only the locked formal F50 field from a trajectory metric table."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or ())
        required = {"model", "replica", FORMAL_F50_FIELD}
        missing = required - fields
        if missing:
            raise ValueError(f"{path}: missing required fields {sorted(missing)}")
        rows = []
        for source in reader:
            rows.append(
                {
                    "model": source["model"],
                    "replica": source["replica"],
                    FORMAL_F50_FIELD: parse_optional_float(source[FORMAL_F50_FIELD]),
                }
            )
    return rows


def summarize_formal_f50(
    rows: list[dict[str, object]],
    models: tuple[str, ...] = MODEL_ORDER,
) -> list[dict[str, object]]:
    """Compute mean and sample SD only when all three raw-row values exist."""
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["model"])].append(row)

    result: list[dict[str, object]] = []
    for model in models:
        group = grouped[model]
        if len(group) != 3:
            raise ValueError(f"{model}: expected exactly three trajectories, got {len(group)}")
        replicas = [str(row["replica"]) for row in group]
        if len(set(replicas)) != 3:
            raise ValueError(f"{model}: replica labels are not unique: {replicas}")
        values = [row[FORMAL_F50_FIELD] for row in group]
        valid = [float(value) for value in values if value is not None]
        complete = len(valid) == 3
        result.append(
            {
                "model": model,
                "n_trajectories": 3,
                "n_valid_F50": len(valid),
                f"{FORMAL_F50_FIELD}_mean": statistics.mean(valid) if complete else None,
                f"{FORMAL_F50_FIELD}_sample_SD": statistics.stdev(valid) if complete else None,
                f"{FORMAL_F50_FIELD}_min": min(valid) if complete else None,
                f"{FORMAL_F50_FIELD}_max": max(valid) if complete else None,
            }
        )
    return result


def csv_value(value: object) -> object:
    return "NA" if value is None else value


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty table")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows({key: csv_value(value) for key, value in row.items()} for row in rows)

