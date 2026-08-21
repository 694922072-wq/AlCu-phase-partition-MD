from __future__ import annotations

import csv
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
import ovito
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, PythonScriptModifier


ROOT = Path(r".")
BASE_SCRIPT = Path(r"[external archival path omitted]")
ANALYSIS = Path(r"[external archival path omitted]")
OUT = ROOT / "figures" / "Fig6_stage2A_schemeA_revision"
CSV_OUT = OUT / "Fig6_stage2A_dense_metric_schemeA.csv"
MD_OUT = OUT / "Fig6_stage2A_dense_metric_schemeA.md"
JSON_OUT = OUT / "Fig6_stage2A_dense_metric_schemeA.json"

SAMPLE_STEP_FRAME = 10  # 1.0% engineering strain in the existing clean tensile trajectories.


def load_base_module():
    spec = importlib.util.spec_from_file_location("fig6_stage1b_renderer", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def attach_metric_properties(frame, data, phase_by_id: np.ndarray, base):
    ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
    structures = np.asarray(data.particles["Structure Type"], dtype=np.int32)
    phases = phase_by_id[ids].astype(np.int32, copy=False)
    matrix_mask = np.isin(phases, (base.PHASE_CODE["Cu"], base.PHASE_CODE["Al"]))
    selected_hcp = matrix_mask & (structures == 2)
    data.particles_.create_property("Phase Code", data=phases)
    data.particles_.create_property("Selected Matrix HCP-like", data=selected_hcp.astype(np.int8))


def frame_to_actual_strain_percent(frame_index: int) -> float:
    return frame_index / 10.0


def sample_frames_for_model(f50_frame: int) -> list[int]:
    frames = list(range(0, f50_frame + 1, SAMPLE_STEP_FRAME))
    if f50_frame not in frames:
        frames.append(f50_frame)
    return sorted(set(frames))


def main() -> int:
    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    base = load_base_module()
    analysis = json.loads(ANALYSIS.read_text(encoding="utf-8"))
    model_analysis = {item["model"]: item for item in analysis["models"]}
    rows: list[dict] = []

    for model, config in base.MODELS.items():
        f50_frame = int(model_analysis[model]["keyframes"]["F50"]["frame_index"])
        frames = sample_frames_for_model(f50_frame)
        print(f"[{model}] sampling {len(frames)} frames through F50 frame {f50_frame}", flush=True)

        phase_by_id = base.build_phase_map(config["construction"], config["layers"])
        phase_codes = phase_by_id[1:]
        matrix_atom_count = int(
            np.count_nonzero(
                np.isin(phase_codes, (base.PHASE_CODE["Cu"], base.PHASE_CODE["Al"]))
            )
        )

        pipeline = import_file(str(config["trajectory"]), multiple_frames=True)
        ptm = PolyhedralTemplateMatchingModifier(
            rmsd_cutoff=base.RMSD_CUTOFF,
            color_by_type=False,
        )
        for structure in ptm.structures:
            structure.enabled = structure.id in (
                PolyhedralTemplateMatchingModifier.Type.FCC,
                PolyhedralTemplateMatchingModifier.Type.HCP,
            )
        pipeline.modifiers.append(ptm)
        pipeline.modifiers.append(
            PythonScriptModifier(
                function=lambda frame, data, phase_by_id=phase_by_id, base=base: attach_metric_properties(
                    frame, data, phase_by_id, base
                )
            )
        )

        for frame_index in frames:
            data = pipeline.compute(frame_index)
            selected = np.asarray(data.particles["Selected Matrix HCP-like"], dtype=np.int8)
            selected_count = int(np.count_nonzero(selected))
            rows.append(
                {
                    "model": model,
                    "frame_index": frame_index,
                    "timestep": int(frame_index) * 1000,
                    "actual_strain_percent": frame_to_actual_strain_percent(frame_index),
                    "matrix_atom_count": matrix_atom_count,
                    "matrix_hcp_like_count": selected_count,
                    "matrix_hcp_like_fraction": selected_count / matrix_atom_count,
                    "is_key_state": frame_index
                    in {
                        int(model_analysis[model]["keyframes"]["initial"]["frame_index"]),
                        int(model_analysis[model]["keyframes"]["15pct"]["frame_index"]),
                        int(model_analysis[model]["keyframes"]["F50"]["frame_index"]),
                    },
                    "is_model_specific_F50": frame_index == f50_frame,
                }
            )
            print(f"  frame {frame_index:03d}: {100.0 * selected_count / matrix_atom_count:.3f}%", flush=True)

    fieldnames = [
        "model",
        "frame_index",
        "timestep",
        "actual_strain_percent",
        "matrix_atom_count",
        "matrix_hcp_like_count",
        "matrix_hcp_like_fraction",
        "is_key_state",
        "is_model_specific_F50",
    ]
    with CSV_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    JSON_OUT.write_text(
        json.dumps(
            {
                "software": f"OVITO Pro {ovito.version_string}",
                "sample_step_frame": SAMPLE_STEP_FRAME,
                "sample_step_strain_percent": 1.0,
                "rows": rows,
                "elapsed_seconds": time.time() - started,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    MD_OUT.write_text(
        "# Fig. 6 Stage 2A dense metric data\n\n"
        "Dense-sampled selected-matrix PTM-HCP-like fraction was computed from the same clean tensile trajectories used for the Scheme A snapshots. "
        "Sampling used a 1.0% engineering-strain step from 0% to each model-specific F50 region, with the exact F50 frame retained when it does not fall on the 1.0% grid.\n\n"
        f"- Output CSV: `{CSV_OUT.name}`\n"
        f"- Total sampled states: {len(rows)}\n"
        "- PTM settings: same RMSD cutoff and enabled FCC/HCP templates as the snapshot rendering.\n"
        "- Interpretation: selected-matrix PTM-HCP-like fraction is used only as a structural-activity indicator.\n",
        encoding="utf-8",
    )

    print(f"Wrote {CSV_OUT}")
    print(f"Finished in {time.time() - started:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
