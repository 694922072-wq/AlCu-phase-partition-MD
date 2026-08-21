from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import ovito
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, ExpressionSelectionModifier, PythonScriptModifier
from ovito.vis import OSPRayRenderer, Viewport


PROJECT = Path(r"D:\leng\AlCu_MAIN_PPP_Zhou")
OUT = PROJECT / (
    r"cms_manuscript_finalization\09_text_references_figures_integration"
    r"\stage3_figure_rebuild\01_fig1_fig2_rebuild_revision2\03_Fig1_R2"
)
RAW = OUT / "ovito_renders"

SLAB_TESTS_A = (20.0, 25.0, 30.0)
SELECTED_SLAB_A = 20.0
TEST_SIZE = (900, 1440)
FINAL_SIZE = (2000, 3000)
PARTICLE_RADIUS_A = 1.12
ORTHO_FOV_A = 204.0

# Element identity remains explicit: all Al atoms are warm coral and all Cu
# atoms are blue. The small phase-conditioned shifts reproduce the approved
# source renderer and do not introduce artificial phase blocks.
COLORS = {
    "Al_matrix": (0.94, 0.49, 0.52),
    "Cu_matrix": (0.36, 0.52, 0.91),
    "Al_in_Al2Cu": (0.98, 0.57, 0.59),
    "Cu_in_Al2Cu": (0.31, 0.61, 0.94),
    "Al_in_Al4Cu9": (0.91, 0.40, 0.46),
    "Cu_in_Al4Cu9": (0.47, 0.46, 0.90),
}

MODELS = {
    "M3_SYM": {
        "equilibrated": PROJECT / r"outputs\M3_SYM_eq_300K_xy_z.data",
        "construction": PROJECT / r"structure\M3_SYM_Cu_Al2CuThick_Al_PPP_ready.data",
        "layers": [
            ("Cu", 43.380), ("Al2Cu", 63.401), ("Al", 97.200),
            ("Al2Cu", 63.401), ("Cu", 43.380),
        ],
    },
    "M4_RATIO": {
        "equilibrated": PROJECT / r"outputs\M4_SYM_RATIO_eq_300K_xy_z.data",
        "construction": PROJECT / r"structure\M4_SYM_RATIO_ready.data",
        "layers": [
            ("Cu", 43.380), ("Al4Cu9", 17.235461), ("Al2Cu", 48.770),
            ("Al", 97.200), ("Al2Cu", 48.770), ("Al4Cu9", 17.235461),
            ("Cu", 43.380),
        ],
    },
    "M4_SYM": {
        "equilibrated": PROJECT / r"outputs\M4_SYM_eq_300K_xy_z.data",
        "construction": PROJECT / r"structure\M4_SYM_Cu_Al4Cu9_Al2Cu_Al_PPP_ready.data",
        "layers": [
            ("Cu", 43.380), ("Al4Cu9", 34.470922), ("Al2Cu", 29.262),
            ("Al", 97.200), ("Al2Cu", 29.262), ("Al4Cu9", 34.470922),
            ("Cu", 43.380),
        ],
    },
    "M4_LIT": {
        "equilibrated": PROJECT / r"outputs\M4_LIT_eq_300K.data",
        "construction": PROJECT / r"structure\M4_LIT_Cu_Al4Cu9_Al2Cu_Al_PPP_ready.data",
        "layers": [
            ("Cu", 86.760), ("Al4Cu9", 34.470922),
            ("Al2Cu", 29.262), ("Al", 97.200),
        ],
    },
}

PHASE_CODE = {"Cu": 1, "Al": 2, "Al2Cu": 3, "Al4Cu9": 4}


def read_atom_count(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2 and fields[1] == "atoms":
                return int(fields[0])
    raise RuntimeError(f"Could not read atom count from {path}")


def build_phase_map(path: Path, layers: list[tuple[str, float]]) -> np.ndarray:
    atom_count = read_atom_count(path)
    phase_by_id = np.zeros(atom_count + 1, dtype=np.int8)
    boundaries = np.cumsum([thickness for _, thickness in layers], dtype=float)
    codes = np.asarray([PHASE_CODE[phase] for phase, _ in layers], dtype=np.int8)
    in_atoms = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not in_atoms:
                if stripped.startswith("Atoms"):
                    in_atoms = True
                continue
            if not stripped:
                continue
            if stripped[0].isalpha():
                break
            fields = stripped.split()
            if len(fields) < 5:
                continue
            atom_id = int(fields[0])
            z = float(fields[4])
            index = min(int(np.searchsorted(boundaries, z, side="right")), len(codes) - 1)
            phase_by_id[atom_id] = codes[index]
    if int(np.count_nonzero(phase_by_id[1:])) != atom_count:
        raise RuntimeError(f"Incomplete phase assignment for {path.name}")
    return phase_by_id


def color_particles(frame, data, phase_by_id: np.ndarray) -> None:
    ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
    elements = np.asarray(data.particles["Particle Type"], dtype=np.int32)
    phases = phase_by_id[ids].astype(np.int32, copy=False)
    if np.any(phases == 0):
        raise RuntimeError("Construction/equilibrated atom identifiers do not match")
    colors = np.empty((len(ids), 3), dtype=np.float32)
    colors[phases == 1] = COLORS["Cu_matrix"]
    colors[phases == 2] = COLORS["Al_matrix"]
    mask = phases == 3
    colors[mask & (elements == 1)] = COLORS["Al_in_Al2Cu"]
    colors[mask & (elements == 2)] = COLORS["Cu_in_Al2Cu"]
    mask = phases == 4
    colors[mask & (elements == 1)] = COLORS["Al_in_Al4Cu9"]
    colors[mask & (elements == 2)] = COLORS["Cu_in_Al4Cu9"]
    data.particles_.create_property("Phase", data=phases)
    data.particles_.create_property("Color", data=colors)


def disable_display_cell(frame, data) -> None:
    data.cell_.pbc = (False, False, False)
    data.cell_.vis.enabled = False
    data.cell_.vis.render_cell = False


def make_renderer(final: bool) -> OSPRayRenderer:
    renderer = OSPRayRenderer()
    renderer.samples_per_pixel = 32 if final else 10
    renderer.refinement_iterations = 2 if final else 1
    renderer.denoising_enabled = True
    renderer.ambient_light_enabled = True
    renderer.ambient_brightness = 0.96
    renderer.direct_light_enabled = True
    renderer.direct_light_intensity = 0.30
    renderer.direct_light_angular_diameter = 30.0
    renderer.material_specular_brightness = 0.015
    renderer.material_shininess = 8.0
    renderer.outlines_enabled = False
    return renderer


def render_one(model: str, config: dict, slab_a: float, final: bool) -> dict:
    phase_by_id = build_phase_map(config["construction"], config["layers"])
    pipeline = import_file(str(config["equilibrated"]))
    pipeline.source.data.particles.vis.radius = PARTICLE_RADIUS_A
    pipeline.modifiers.append(
        PythonScriptModifier(function=lambda frame, data: color_particles(frame, data, phase_by_id))
    )

    full_data = pipeline.compute()
    cell = np.asarray(full_data.cell)
    xlo, ylo, zlo = cell[0, 3], cell[1, 3], cell[2, 3]
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xcenter, ycenter, zcenter = xlo + lx / 2, ylo + ly / 2, zlo + lz / 2
    ymin, ymax = ycenter - slab_a / 2, ycenter + slab_a / 2
    pipeline.modifiers.append(
        ExpressionSelectionModifier(expression=f"Position.Y < {ymin:.12f} || Position.Y > {ymax:.12f}")
    )
    pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.modifiers.append(PythonScriptModifier(function=disable_display_cell))
    display_data = pipeline.compute()

    pipeline.add_to_scene()
    viewport = Viewport(type=Viewport.Type.Ortho)
    viewport.camera_dir = (0.0, 1.0, 0.0)
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.camera_pos = (xcenter, ycenter - 600.0, zcenter)
    viewport.fov = ORTHO_FOV_A
    size = FINAL_SIZE if final else TEST_SIZE
    if final:
        out_file = RAW / "final" / f"{model}_verticalZ_slab{int(slab_a):02d}_highres.png"
    else:
        out_file = RAW / "slab_tests" / f"{model}_verticalZ_slab{int(slab_a):02d}.png"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    viewport.render_image(
        filename=str(out_file),
        size=size,
        renderer=make_renderer(final),
        background=(1.0, 1.0, 1.0),
        alpha=False,
    )
    pipeline.remove_from_scene()

    phase_counts = np.bincount(phase_by_id[1:], minlength=5)
    return {
        "model": model,
        "equilibrated_file": str(config["equilibrated"]),
        "construction_identity_file": str(config["construction"]),
        "full_atom_count": int(full_data.particles.count),
        "display_atom_count": int(display_data.particles.count),
        "phase_identity_counts": {
            "Cu": int(phase_counts[1]), "Al": int(phase_counts[2]),
            "Al2Cu": int(phase_counts[3]), "Al4Cu9": int(phase_counts[4]),
        },
        "box_A": {"x": float(lx), "y": float(ly), "z": float(lz)},
        "central_slab_A": slab_a,
        "camera": {
            "projection": "orthographic",
            "direction": [0.0, 1.0, 0.0],
            "up": [0.0, 0.0, 1.0],
            "orientation": "Z vertical, X horizontal, Y projection depth; no side face",
        },
        "particle_radius_A": PARTICLE_RADIUS_A,
        "render_size_px": list(size),
        "output": str(out_file),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Display-only OVITO rerender for Fig.1 Revision 2")
    parser.add_argument("--mode", choices=("tests", "final", "all"), default="all")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    runs: list[dict] = []
    if args.mode in ("tests", "all"):
        for slab in SLAB_TESTS_A:
            for model, config in MODELS.items():
                item = render_one(model, config, slab, final=False)
                runs.append(item)
                print(f"TEST {model} slab={slab:g} A atoms={item['display_atom_count']}")
    if args.mode in ("final", "all"):
        for model, config in MODELS.items():
            item = render_one(model, config, SELECTED_SLAB_A, final=True)
            runs.append(item)
            print(f"FINAL {model} slab={SELECTED_SLAB_A:g} A atoms={item['display_atom_count']}")

    metadata = {
        "ovito_version": ovito.version_string,
        "render_engine": "OSPRayRenderer",
        "scientific_operation": "none; display-only central-slab rerender",
        "projection": "orthographic",
        "slab_tests_A": list(SLAB_TESTS_A),
        "selected_slab_A": SELECTED_SLAB_A,
        "colors": COLORS,
        "runs": runs,
    }
    (RAW / f"render_metadata_{args.mode}.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
