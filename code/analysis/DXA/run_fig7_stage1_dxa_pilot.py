from __future__ import annotations

import csv
import gc
import json
import math
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ovito
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.modifiers import DislocationAnalysisModifier, PythonScriptModifier
from ovito.vis import OSPRayRenderer, Viewport


ROOT = Path(r".")
OUT = ROOT / "figures" / "Fig7_stage1_DXA_pilot"
RAW = OUT / "raw_snapshots"
SOURCE = OUT / "Fig7_stage1_source_data"
FIG6_FRAME_MAP = Path(r"[external archival path omitted]")

RENDER_SIZE = (1700, 2050)
CONTACT_THUMB = (560, 676)
ORTHO_FOV_PAD = 1.08
DXA_LATTICE = "FCC"

PHASE_CODE = {"Cu": 1, "Al": 2, "Al2Cu": 3, "Al4Cu9": 4}
MATRIX_PHASES = (PHASE_CODE["Cu"], PHASE_CODE["Al"])

PHASE_COLORS = {
    "Cu": (0.48, 0.55, 0.72),
    "Al": (0.68, 0.78, 0.88),
    "Al2Cu": (0.72, 0.84, 0.70),
    "Al4Cu9": (0.56, 0.78, 0.82),
}

MODELS = {
    "M3_SYM": {
        "construction": ROOT / "structure" / "M3_SYM_Cu_Al2CuThick_Al_PPP_ready.data",
        "trajectory": ROOT / "dumps" / "M3_SYM_tensile_x25_xy_z_clean_light.lammpstrj",
        "stress_strain": ROOT / "outputs" / "M3_SYM_stress_strain_x25_xy_z_clean.dat",
        "f50": 19.20,
        "layers": [
            ("Cu", 43.380),
            ("Al2Cu", 63.401),
            ("Al", 97.200),
            ("Al2Cu", 63.401),
            ("Cu", 43.380),
        ],
    },
    "M4_RATIO": {
        "construction": ROOT / "structure" / "M4_SYM_RATIO_ready.data",
        "trajectory": ROOT / "dumps" / "M4_SYM_RATIO_tensile_x25_xy_z_clean_light.lammpstrj",
        "stress_strain": ROOT / "outputs" / "M4_SYM_RATIO_stress_strain_x25_xy_z_clean.dat",
        "f50": 17.98,
        "layers": [
            ("Cu", 43.380),
            ("Al4Cu9", 17.235461),
            ("Al2Cu", 48.770),
            ("Al", 97.200),
            ("Al2Cu", 48.770),
            ("Al4Cu9", 17.235461),
            ("Cu", 43.380),
        ],
    },
    "M4_SYM": {
        "construction": ROOT / "structure" / "M4_SYM_Cu_Al4Cu9_Al2Cu_Al_PPP_ready.data",
        "trajectory": ROOT / "dumps" / "M4_SYM_tensile_x25_xy_z_clean_light.lammpstrj",
        "stress_strain": ROOT / "outputs" / "M4_SYM_stress_strain_x25_xy_z_clean.dat",
        "f50": 17.74,
        "layers": [
            ("Cu", 43.380),
            ("Al4Cu9", 34.470922),
            ("Al2Cu", 29.262),
            ("Al", 97.200),
            ("Al2Cu", 29.262),
            ("Al4Cu9", 34.470922),
            ("Cu", 43.380),
        ],
    },
}

STATE_ORDER = ["initial", "15pct", "F50"]
STATE_LABEL = {"initial": "Initial", "15pct": "15%", "F50": "F50"}
MODEL_LABEL = {"M3_SYM": "M3_SYM", "M4_RATIO": "M4_RATIO", "M4_SYM": "M4_SYM"}

LENGTH_KEYS = [
    "DislocationAnalysis.length.1/6<112>",
    "DislocationAnalysis.length.1/2<110>",
    "DislocationAnalysis.length.1/6<110>",
    "DislocationAnalysis.length.1/3<100>",
    "DislocationAnalysis.length.1/3<111>",
    "DislocationAnalysis.length.other",
]


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        r"[external archival path omitted]" if bold else r"[external archival path omitted]",
        r"[external archival path omitted]" if bold else r"[external archival path omitted]",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size)
    return ImageFont.load_default()


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

    assigned = int(np.count_nonzero(phase_by_id[1:]))
    if assigned != atom_count:
        raise RuntimeError(f"Incomplete phase assignment for {path}: {assigned}/{atom_count}")
    return phase_by_id


def build_frame_maps() -> dict[str, dict[str, dict]]:
    if FIG6_FRAME_MAP.exists():
        data = json.loads(FIG6_FRAME_MAP.read_text(encoding="utf-8"))
        return {
            model["model"]: {state: dict(model["keyframes"][state]) for state in STATE_ORDER}
            for model in data["models"]
            if model["model"] in MODELS
        }
    raise RuntimeError(f"Missing frame map: {FIG6_FRAME_MAP}")


def matrix_thickness(layers: list[tuple[str, float]]) -> float:
    return float(sum(thickness for phase, thickness in layers if phase in ("Cu", "Al")))


def prep_particles(data, phase_by_id: np.ndarray):
    ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
    phases = phase_by_id[ids].astype(np.int8, copy=False)
    if np.any(phases == 0):
        raise RuntimeError("Construction and trajectory atom identifiers do not match.")

    colors = np.empty((len(ids), 3), dtype=np.float32)
    colors[phases == PHASE_CODE["Cu"]] = PHASE_COLORS["Cu"]
    colors[phases == PHASE_CODE["Al"]] = PHASE_COLORS["Al"]
    colors[phases == PHASE_CODE["Al2Cu"]] = PHASE_COLORS["Al2Cu"]
    colors[phases == PHASE_CODE["Al4Cu9"]] = PHASE_COLORS["Al4Cu9"]

    matrix_mask = np.isin(phases, MATRIX_PHASES)
    radii = np.full(len(ids), 0.32, dtype=np.float32)
    radii[matrix_mask] = 0.26
    transparency = np.full(len(ids), 0.90, dtype=np.float32)
    transparency[matrix_mask] = 0.82

    data.particles_.create_property("Phase Code", data=phases)
    data.particles_.create_property("Selection", data=matrix_mask.astype(np.int8))
    data.particles_.create_property("Color", data=colors)
    data.particles_.create_property("Radius", data=radii)
    data.particles_.create_property("Transparency", data=transparency)


def make_renderer() -> OSPRayRenderer:
    renderer = OSPRayRenderer()
    renderer.samples_per_pixel = 12
    renderer.refinement_iterations = 1
    renderer.denoising_enabled = True
    renderer.ambient_light_enabled = True
    renderer.ambient_brightness = 1.0
    renderer.direct_light_enabled = True
    renderer.direct_light_intensity = 0.14
    renderer.direct_light_angular_diameter = 45.0
    renderer.material_specular_brightness = 0.0
    renderer.material_shininess = 1.0
    renderer.outlines_enabled = False
    return renderer


def configure_scene_objects(data) -> int:
    segment_count = 0
    for obj in data.objects:
        name = type(obj).__name__
        if name == "DislocationNetwork":
            segment_count = len(obj.segments)
            obj.vis.enabled = True
            obj.vis.line_width = 3.0
            obj.vis.show_burgers_vectors = False
            obj.vis.show_line_directions = False
            obj.vis.coloring_mode = obj.vis.ColoringMode.ByDislocationType
        elif name == "SurfaceMesh" and getattr(obj, "identifier", "") == "dxa-defect-mesh":
            obj.vis.enabled = False
        elif name == "SimulationCell":
            obj.vis.enabled = True
            obj.vis.render_cell = True
            obj.vis.line_width = 0.18
            obj.vis.rendering_color = (0.60, 0.63, 0.70)
    return segment_count


def render_snapshot(pipeline, frame_index: int, output_path: Path) -> dict:
    data = pipeline.compute(frame_index)
    segment_count = configure_scene_objects(data)

    cell = np.asarray(data.cell)
    xlo, ylo, zlo = cell[0, 3], cell[1, 3], cell[2, 3]
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xcenter, ycenter, zcenter = xlo + lx / 2, ylo + ly / 2, zlo + lz / 2

    pipeline.add_to_scene()
    viewport = Viewport(type=Viewport.Type.Ortho)
    viewport.camera_dir = (-0.12, 0.992, 0.0)
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.camera_pos = (xcenter + 78.0, ycenter - 610.0, zcenter)
    viewport.fov = max(lx, lz) * ORTHO_FOV_PAD
    viewport.render_image(
        filename=str(output_path),
        size=RENDER_SIZE,
        renderer=make_renderer(),
        background=(1.0, 1.0, 1.0),
        alpha=False,
        frame=frame_index,
    )
    pipeline.remove_from_scene()
    return {
        "box_x_A": float(lx),
        "box_y_A": float(ly),
        "box_z_A": float(lz),
        "line_segment_count": int(segment_count),
    }


def extract_metrics(model: str, state: str, frame_map: dict, data, config: dict, render_info: dict) -> dict:
    attrs = data.attributes
    cell = np.asarray(data.cell)
    lx, ly, lz = float(cell[0, 0]), float(cell[1, 1]), float(cell[2, 2])
    selected_matrix_volume = lx * ly * matrix_thickness(config["layers"])
    ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
    phases = np.asarray(data.particles["Phase Code"], dtype=np.int8)
    selected_matrix_atoms = int(np.count_nonzero(np.isin(phases, MATRIX_PHASES)))
    total_length = float(attrs.get("DislocationAnalysis.total_line_length", math.nan))
    density_nm2 = total_length / selected_matrix_volume * 100.0 if selected_matrix_volume > 0 else math.nan
    normalized_A_per_100k_atoms = (
        total_length / selected_matrix_atoms * 100000.0 if selected_matrix_atoms > 0 else math.nan
    )
    row = {
        "model": model,
        "state": state,
        "target_strain_percent": float(frame_map["target_strain_percent"]),
        "actual_strain_percent": float(frame_map["actual_strain_percent"]),
        "frame_index": int(frame_map["frame_index"]),
        "timestep": int(frame_map["timestep"]),
        "selected_matrix_atoms": selected_matrix_atoms,
        "selected_matrix_volume_A3": selected_matrix_volume,
        "total_DXA_line_length_A": total_length,
        "DXA_line_density_nm^-2": density_nm2,
        "normalized_line_length_A_per_100k_matrix_atoms": normalized_A_per_100k_atoms,
        "line_segment_count": int(render_info["line_segment_count"]),
        "DXA_count_FCC": int(attrs.get("DislocationAnalysis.counts.FCC", 0)),
        "DXA_count_HCP": int(attrs.get("DislocationAnalysis.counts.HCP", 0)),
        "DXA_count_OTHER": int(attrs.get("DislocationAnalysis.counts.OTHER", 0)),
        "box_x_A": float(render_info["box_x_A"]),
        "box_y_A": float(render_info["box_y_A"]),
        "box_z_A": float(render_info["box_z_A"]),
    }
    for key in LENGTH_KEYS:
        label = key.replace("DislocationAnalysis.length.", "length_")
        row[f"{label}_A"] = float(attrs.get(key, 0.0))
    return row


def write_csv(rows: list[dict], path: Path):
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_metric_plot(rows: list[dict], path: Path):
    fig, ax = plt.subplots(figsize=(6.0, 4.55), dpi=320)
    colors = {"M3_SYM": "#356AA0", "M4_RATIO": "#D9822B", "M4_SYM": "#8D4AAE"}
    markers = {"initial": "o", "15pct": "s", "F50": "^"}
    for model in MODELS:
        subset = [r for r in rows if r["model"] == model]
        subset.sort(key=lambda r: r["actual_strain_percent"])
        ax.plot(
            [r["actual_strain_percent"] for r in subset],
            [r["normalized_line_length_A_per_100k_matrix_atoms"] for r in subset],
            color=colors[model],
            lw=1.8,
            label=model,
        )
        for row in subset:
            ax.scatter(
                row["actual_strain_percent"],
                row["normalized_line_length_A_per_100k_matrix_atoms"],
                color=colors[model],
                marker=markers[row["state"]],
                s=34,
                edgecolor="white",
                linewidth=0.7,
                zorder=3,
            )
    ax.set_xlabel("Engineering strain (%)")
    ax.set_ylabel("DXA line length / 100k matrix atoms (A)")
    ax.set_title("Stage 1 selected-matrix DXA metric preview", fontsize=10)
    ax.grid(True, color="#D8DDE6", lw=0.6)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    fig.savefig(path.with_suffix(".tif"), dpi=600)
    plt.close(fig)


def make_contact_sheet(rows: list[dict], plot_path: Path, output_path: Path):
    margin = 44
    gutter = 24
    title_h = 92
    footer_h = 138
    cols = 3
    rows_count = 3
    panel_w, panel_h = CONTACT_THUMB
    plot_w, plot_h = 780, 590
    sheet_w = margin * 2 + cols * panel_w + (cols - 1) * gutter + 28 + plot_w
    sheet_h = margin + title_h + rows_count * panel_h + (rows_count - 1) * gutter + footer_h
    canvas = Image.new("RGB", (sheet_w, sheet_h), "white")
    draw = ImageDraw.Draw(canvas)

    draw.text((margin, 28), "Fig. 7 Stage 1 DXA pilot feasibility", font=font(34, True), fill="#111827")
    draw.text(
        (margin, 68),
        "DXA restricted to selected FCC Cu/Al matrix regions; pilot snapshots only.",
        font=font(20),
        fill="#4B5563",
    )

    top = margin + title_h
    by_key = {(r["model"], r["state"]): r for r in rows}
    for row_idx, model in enumerate(MODELS):
        y = top + row_idx * (panel_h + gutter)
        draw.text((margin - 8, y + 4), MODEL_LABEL[model], font=font(21, True), fill="#111827")
        for col_idx, state in enumerate(STATE_ORDER):
            x = margin + col_idx * (panel_w + gutter)
            record = by_key[(model, state)]
            img = Image.open(record["snapshot_path"]).convert("RGB")
            img.thumbnail((panel_w, panel_h - 58), Image.Resampling.LANCZOS)
            px = x + (panel_w - img.width) // 2
            py = y + 40
            draw.rounded_rectangle(
                (x, y + 34, x + panel_w, y + panel_h),
                radius=8,
                outline="#D7DEE8",
                width=1,
                fill="#FFFFFF",
            )
            canvas.paste(img, (px, py))
            label = f"{STATE_LABEL[state]} | eps={record['actual_strain_percent']:.2f}%"
            draw.text((x + 12, y + 10), label, font=font(18, True), fill="#111827")
            draw.text(
                (x + 12, y + panel_h - 31),
                f"L={record['total_DXA_line_length_A']:.1f} A, segments={record['line_segment_count']}",
                font=font(15),
                fill="#4B5563",
            )

    plot_x = margin + cols * panel_w + (cols - 1) * gutter + 28
    plot_y = top + 42
    plot = Image.open(plot_path).convert("RGB")
    plot.thumbnail((plot_w, plot_h), Image.Resampling.LANCZOS)
    draw.rounded_rectangle(
        (plot_x, plot_y - 34, plot_x + plot_w, plot_y + plot_h + 120),
        radius=8,
        outline="#D7DEE8",
        width=1,
        fill="#FFFFFF",
    )
    draw.text((plot_x + 12, plot_y - 26), "Metric preview", font=font(22, True), fill="#111827")
    canvas.paste(plot, (plot_x + (plot_w - plot.width) // 2, plot_y))
    note_y = plot_y + plot_h + 20
    notes = [
        "Metric shown: normalized DXA line length.",
        "Density estimate is also exported in CSV.",
        "Stage 1 is a feasibility gate, not final layout.",
        "No crack/void/fracture claim is made.",
    ]
    for i, note in enumerate(notes):
        draw.text((plot_x + 16, note_y + i * 27), f"- {note}", font=font(17), fill="#374151")

    footer_y = sheet_h - footer_h + 18
    draw.line((margin, footer_y - 18, sheet_w - margin, footer_y - 18), fill="#E5E7EB", width=1)
    draw.text((margin, footer_y), "Line coloring uses OVITO DXA dislocation types. Context atoms are semi-transparent phase identities.", font=font(17), fill="#374151")
    draw.text((margin, footer_y + 28), "Interpretation boundary: line-defect activity in selected FCC matrix/interface-neighboring matrix only.", font=font(17, True), fill="#1F2937")
    canvas.save(output_path)
    canvas.save(output_path.with_suffix(".tif"), dpi=(600, 600), compression="tiff_lzw")


def write_markdown(rows: list[dict], metadata: dict):
    metric_md = OUT / "Fig7_stage1_DXA_metric_preview.md"
    lines = [
        "# Fig. 7 Stage 1 DXA metric preview",
        "",
        "Stage 1 tests whether DXA can be run reproducibly in selected FCC Cu/Al matrix regions.",
        "It does not establish fracture, void nucleation, or IMC-internal dislocation mechanisms.",
        "",
        "| Model | State | Strain (%) | Total line length (A) | Segments | Normalized length (A / 100k matrix atoms) | Density estimate (nm^-2) | Dominant length class |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        class_values = {
            key.replace("length_", "").replace("_A", ""): row[key]
            for key in row
            if key.startswith("length_") and key.endswith("_A")
        }
        dominant = max(class_values.items(), key=lambda item: item[1])[0]
        lines.append(
            f"| {row['model']} | {STATE_LABEL[row['state']]} | {row['actual_strain_percent']:.2f} | "
            f"{row['total_DXA_line_length_A']:.2f} | {row['line_segment_count']} | "
            f"{row['normalized_line_length_A_per_100k_matrix_atoms']:.2f} | "
            f"{row['DXA_line_density_nm^-2']:.5f} | {dominant} |"
        )
    lines.extend(
        [
            "",
            "Metric definitions:",
            "",
            "- `total_DXA_line_length_A`: OVITO `DislocationAnalysis.total_line_length`, restricted with `only_selected=True`.",
            "- `normalized_line_length_A_per_100k_matrix_atoms`: total DXA line length divided by selected Cu/Al matrix atom count and scaled to 100,000 atoms.",
            "- `DXA_line_density_nm^-2`: total DXA line length divided by an approximate selected-matrix layer volume from phase-layer thickness and current x-y cell size; this is a pilot estimate.",
            "",
            "Selection:",
            "",
            "- Atom IDs are mapped back to the construction-source phase layers.",
            "- Selected matrix phases are Cu and Al only.",
            "- Al2Cu and Al4Cu9 atoms are excluded from the DXA selection.",
            "",
        ]
    )
    metric_md.write_text("\n".join(lines), encoding="utf-8")

    generation_log = OUT / "Fig7_stage1_generation_log.md"
    generation_log.write_text(
        "\n".join(
            [
                "# Fig. 7 Stage 1 generation log",
                "",
                f"- Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
                f"- OVITO version: {metadata['software']['ovito_version']}",
                f"- Python executable: {metadata['software']['python_executable']}",
                "- Analysis: OVITO DislocationAnalysisModifier, input crystal structure FCC, only_selected=True.",
                "- Models: M3_SYM, M4_RATIO, M4_SYM.",
                "- States: initial, 15%, model-specific F50.",
                "- Frame mapping source: Fig. 6 Stage 1 color-corrected metadata.",
                "- Output type: Stage 1 feasibility pilot only.",
                "- No manuscript, SI, final Fig. 6, source trajectory, or model file was modified.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    quality = OUT / "Fig7_stage1_quality_check.md"
    quality.write_text(
        "\n".join(
            [
                "# Fig. 7 Stage 1 quality check",
                "",
                "| Check | Status | Note |",
                "|---|---|---|",
                "| Three core model trajectories readable | PASS | M3_SYM, M4_RATIO, M4_SYM were processed. |",
                "| Nine key frames located | PASS | initial, 15%, F50 for all three models. |",
                "| Selected Cu/Al matrix definition | PASS | Construction-source atom-ID phase map; Cu/Al included, IMCs excluded. |",
                "| DXA only on selected matrix | PASS | `DislocationAnalysisModifier.only_selected=True`. |",
                "| Total line length exported | PASS | `total_DXA_line_length_A` in CSV. |",
                "| Segment count exported | PASS | `line_segment_count` in CSV. |",
                "| Normalized metric exported | PASS | Normalized line length and pilot density estimate both exported. |",
                "| Pilot snapshots generated | PASS | 9 individual PNGs and one contact sheet. |",
                "| Claim-boundary language | PASS | No crack/void/fracture interpretation is used. |",
                "| Human review still required | YES | Stage 1 is not the final Fig. 7 layout. |",
                "",
            ]
        ),
        encoding="utf-8",
    )

    boundary = OUT / "claim_boundary_check_stage1.md"
    boundary.write_text(
        "\n".join(
            [
                "# Fig. 7 Stage 1 claim-boundary check",
                "",
                "PASS.",
                "",
                "Allowed interpretation:",
                "",
                "> DXA line activity is used as a line-defect-level structural indicator in selected FCC matrix regions.",
                "",
                "Not claimed:",
                "",
                "- DXA proves fracture.",
                "- DXA confirms void nucleation.",
                "- DXA proves damage inside Al2Cu or Al4Cu9.",
                "- DXA confirms the failure mechanism.",
                "",
                "All metrics and snapshots are restricted to selected Cu/Al matrix atoms through `only_selected=True`.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    unresolved = OUT / "unresolved_items_stage1.md"
    unresolved.write_text(
        "\n".join(
            [
                "# Fig. 7 Stage 1 unresolved items",
                "",
                "- Stage 1 snapshots are pilot-quality and should not be submitted as final artwork.",
                "- Visual density/class coloring may need simplification in Stage 2 if the line network is visually busy.",
                "- `DXA_line_density_nm^-2` uses an approximate selected-matrix volume; normalized line length per selected matrix atom is the safer metric if strict volume definition is challenged.",
                "- M4_LIT is not included in the main pilot chain; it may be added later only as SI/reference if requested.",
                "",
                "No blocking item prevents ChatGPT/human review of Stage 1.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    manifest = OUT / "Fig7_stage1_source_data_manifest.md"
    manifest.write_text(
        "\n".join(
            [
                "# Fig. 7 Stage 1 source-data manifest",
                "",
                "- `Fig7_stage1_DXA_metric_preview.csv`: primary Stage 1 DXA metrics.",
                "- `Fig7_stage1_DXA_metric_preview.md`: readable metric summary and method note.",
                "- `Fig7_stage1_pilot_contact_sheet.png`: pilot snapshot contact sheet.",
        "- `Fig7_stage1_DXA_metric_plot.png/.tif`: metric preview plot.",
        "- `Fig7_stage1_pilot_contact_sheet.tif`: high-resolution contact sheet.",
                "- `raw_snapshots/*.png`: nine individual pilot snapshots.",
                "- `Fig7_stage1_metadata.json`: full run metadata.",
                "- `run_fig7_stage1_dxa_pilot.py`: reproducible generation script.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def copy_source_data():
    for file_name in [
        "Fig7_stage1_DXA_metric_preview.csv",
        "Fig7_stage1_DXA_metric_preview.md",
        "Fig7_stage1_metadata.json",
        "claim_boundary_check_stage1.md",
        "Fig7_stage1_quality_check.md",
        "unresolved_items_stage1.md",
    ]:
        src = OUT / file_name
        if src.exists():
            (SOURCE / file_name).write_bytes(src.read_bytes())


def main() -> int:
    start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)

    frame_maps = build_frame_maps()
    rows: list[dict] = []
    metadata = {
        "software": {
            "ovito_version": ovito.version_string,
            "python_executable": str(Path(__import__("sys").executable)),
        },
        "method": {
            "modifier": "DislocationAnalysisModifier",
            "input_crystal_structure": DXA_LATTICE,
            "only_selected": True,
            "selected_phases": ["Cu", "Al"],
            "excluded_phases": ["Al2Cu", "Al4Cu9"],
            "stage": "Stage 1 pilot feasibility",
        },
        "models": [],
    }

    for model, config in MODELS.items():
        model_start = time.time()
        print(f"[{model}] loading", flush=True)
        phase_by_id = build_phase_map(config["construction"], config["layers"])
        pipeline = import_file(str(config["trajectory"]), multiple_frames=True)
        frame_count = int(pipeline.source.num_frames)

        def prep(frame, data, phase_by_id=phase_by_id):
            prep_particles(data, phase_by_id)

        pipeline.modifiers.append(PythonScriptModifier(function=prep))
        dxa = DislocationAnalysisModifier()
        dxa.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
        dxa.only_selected = True
        dxa.color_by_type = True
        dxa.line_coarsening_enabled = True
        dxa.line_smoothing_enabled = True
        dxa.line_smoothing_level = 1
        pipeline.modifiers.append(dxa)

        model_record = {
            "model": model,
            "construction": str(config["construction"]),
            "trajectory": str(config["trajectory"]),
            "stress_strain": str(config["stress_strain"]),
            "frame_count": frame_count,
            "f50_target_percent": config["f50"],
            "frames": [],
        }

        for state in STATE_ORDER:
            frame_map = frame_maps[model][state]
            frame_index = int(frame_map["frame_index"])
            snapshot_path = RAW / f"Fig7_stage1_DXA_{model}_{state}.png"
            print(
                f"[{model}] {state}: frame={frame_index}, strain={frame_map['actual_strain_percent']:.4f}%",
                flush=True,
            )
            render_info = render_snapshot(pipeline, frame_index, snapshot_path)
            data = pipeline.compute(frame_index)
            configure_scene_objects(data)
            row = extract_metrics(model, state, frame_map, data, config, render_info)
            row["snapshot_path"] = str(snapshot_path)
            rows.append(row)
            model_record["frames"].append(dict(row))
            print(
                f"[{model}] {state}: length={row['total_DXA_line_length_A']:.2f} A, "
                f"segments={row['line_segment_count']}",
                flush=True,
            )

        model_record["elapsed_seconds"] = time.time() - model_start
        metadata["models"].append(model_record)
        pipeline = None
        gc.collect()

    csv_path = OUT / "Fig7_stage1_DXA_metric_preview.csv"
    write_csv(rows, csv_path)
    plot_path = OUT / "Fig7_stage1_DXA_metric_plot.png"
    make_metric_plot(rows, plot_path)
    contact_path = OUT / "Fig7_stage1_pilot_contact_sheet.png"
    make_contact_sheet(rows, plot_path, contact_path)

    metadata["elapsed_seconds"] = time.time() - start
    metadata["outputs"] = {
        "metric_csv": str(csv_path),
        "metric_plot": str(plot_path),
        "contact_sheet": str(contact_path),
        "raw_snapshots": [row["snapshot_path"] for row in rows],
    }
    (OUT / "Fig7_stage1_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_markdown(rows, metadata)
    copy_source_data()

    print(f"Finished Fig. 7 Stage 1 in {metadata['elapsed_seconds']:.1f} s", flush=True)
    print(contact_path, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
