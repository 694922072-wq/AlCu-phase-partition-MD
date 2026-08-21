from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import numpy as np
import ovito
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, PythonScriptModifier
from ovito.vis import OpenGLRenderer


ROOT = Path(r".")
BASE_SCRIPT = Path(r"[external archival path omitted]")
ANALYSIS = Path(r"[external archival path omitted]")
OUT = ROOT / "figures" / "Fig6_stage2A_schemeA_revision"
RAW = OUT / "_ovito_renders_schemeA_50A"
METADATA = OUT / "Fig6_stage2A_schemeA_render_metadata.json"
DISPLAY_DEPTH_A = 50.0


SCHEME_A_HEX = {
    "Cu_context": "#4F6DFF",
    "Al_context": "#8AA8FF",
    "Al4Cu9_context": "#2DB7C9",
    "Al2Cu_context": "#74D66A",
    "PTM_HCP_like": "#E64B35",
    "matrix_other": "#C7D4E8",
    "guide": "#555555",
    "background": "#FFFFFF",
}


def hex_to_rgb01(value: str) -> tuple[float, float, float]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


SCHEME_A_RGB = {key: hex_to_rgb01(value) for key, value in SCHEME_A_HEX.items()}


def load_base_module():
    spec = importlib.util.spec_from_file_location("fig6_stage1b_renderer", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def apply_scheme_a_context_and_overlay(frame, data, phase_by_id: np.ndarray, base):
    ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
    structures = np.asarray(data.particles["Structure Type"], dtype=np.int32)
    phases = phase_by_id[ids].astype(np.int32, copy=False)

    colors = np.empty((len(ids), 3), dtype=np.float32)
    radii = np.full(len(ids), 0.70, dtype=np.float32)
    transparency = np.full(len(ids), 0.0, dtype=np.float32)

    colors[phases == base.PHASE_CODE["Cu"]] = SCHEME_A_RGB["Cu_context"]
    colors[phases == base.PHASE_CODE["Al"]] = SCHEME_A_RGB["Al_context"]
    colors[phases == base.PHASE_CODE["Al2Cu"]] = SCHEME_A_RGB["Al2Cu_context"]
    colors[phases == base.PHASE_CODE["Al4Cu9"]] = SCHEME_A_RGB["Al4Cu9_context"]

    matrix_mask = np.isin(phases, (base.PHASE_CODE["Cu"], base.PHASE_CODE["Al"]))
    matrix_other = matrix_mask & ~np.isin(structures, (1, 2))
    selected_hcp = matrix_mask & (structures == 2)

    colors[matrix_other] = SCHEME_A_RGB["matrix_other"]
    transparency[matrix_other] = 0.0
    radii[matrix_other] = 0.68

    colors[selected_hcp] = SCHEME_A_RGB["PTM_HCP_like"]
    transparency[selected_hcp] = 0.0
    radii[selected_hcp] = 1.03

    data.particles_.create_property("Phase Code", data=phases)
    data.particles_.create_property("Selected Matrix HCP-like", data=selected_hcp.astype(np.int8))
    data.particles_.create_property("Color", data=colors)
    data.particles_.create_property("Radius", data=radii)
    data.particles_.create_property("Transparency", data=transparency)


def make_flat_renderer():
    renderer = OpenGLRenderer()
    renderer.antialiasing_level = 4
    renderer.order_independent_transparency = True
    return renderer


def main() -> int:
    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)

    base = load_base_module()
    base.DISPLAY_DEPTH_A = DISPLAY_DEPTH_A
    base.PALETTE = SCHEME_A_RGB
    base.make_renderer = make_flat_renderer

    analysis = json.loads(ANALYSIS.read_text(encoding="utf-8"))
    model_analysis = {item["model"]: item for item in analysis["models"]}

    metadata = {
        "software": f"OVITO Pro {ovito.version_string}",
        "stage": "Fig. 6 Stage 2A Scheme A revision",
        "strategy": "50 A central-Y-slab phase-context background plus selected Cu/Al-matrix PTM-HCP-like overlay",
        "palette_hex": SCHEME_A_HEX,
        "display_depth_A": DISPLAY_DEPTH_A,
        "rendering_adjustments": {
            "renderer": "OpenGLRenderer flat-style",
            "context_radius_A": 0.70,
            "context_transparency": 0.0,
            "other_unclassified_radius_A": 0.68,
            "other_unclassified_transparency": 0.0,
            "hcp_like_radius_A": 1.03,
            "hcp_like_transparency": 0.0,
        },
        "camera": {
            "projection": "orthographic",
            "direction": [-0.14, 0.990, 0.0],
            "up": [0.0, 0.0, 1.0],
            "fov": base.ORTHO_FOV,
        },
        "models": [],
    }

    for model, config in base.MODELS.items():
        print(f"[{model}] loading trajectory", flush=True)
        phase_by_id = base.build_phase_map(config["construction"], config["layers"])
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
                function=lambda frame, data, phase_by_id=phase_by_id, base=base: apply_scheme_a_context_and_overlay(
                    frame, data, phase_by_id, base
                )
            )
        )

        model_record = {
            "model": model,
            "trajectory": str(config["trajectory"]),
            "construction": str(config["construction"]),
            "layers": config["layers"],
            "frames": [],
        }

        for state in ("initial", "15pct", "F50"):
            frame_map = model_analysis[model]["keyframes"][state]
            frame_index = int(frame_map["frame_index"])
            output = RAW / f"Fig6_stage2A_schemeA_{model}_{state}_50A.png"
            print(
                f"[{model}] {state}: frame={frame_index}, "
                f"strain={frame_map['actual_strain_percent']:.4f}%",
                flush=True,
            )
            render_info = base.render_frame(pipeline, frame_index, output)
            model_record["frames"].append(
                {
                    "state": state,
                    "raw_path": str(output),
                    **frame_map,
                    **render_info,
                }
            )
        metadata["models"].append(model_record)

    metadata["elapsed_seconds"] = time.time() - started
    METADATA.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"Finished in {metadata['elapsed_seconds']:.1f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
