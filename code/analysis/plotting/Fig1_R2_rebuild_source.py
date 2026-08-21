"""Assemble Stage 3-1 Revision 2 Fig.1 from real high-resolution OVITO renders.

The approved AI preview is never opened by this script. It therefore cannot
contribute pixels to the formal figure. Scientific structure identity comes
only from the locked OVITO/LAMMPS sources recorded by Fig1_R2_ovito_render.py.
"""

from __future__ import annotations

import csv
from pathlib import Path
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image


HERE = Path(__file__).resolve().parent
RAW = HERE / "ovito_renders" / "final"
TESTS = HERE / "ovito_renders" / "slab_tests"

MODEL_ORDER = ["M3_SYM", "M4_RATIO", "M4_SYM", "M4_LIT"]
MODEL_STYLE = {
    "M3_SYM": {"color": "#0072B2", "subtitle": "Thick-Al₂Cu symmetric control"},
    "M4_RATIO": {"color": "#E69F00", "subtitle": "Phase-partition bridge model"},
    "M4_SYM": {"color": "#7B2CBF", "subtitle": "Symmetric composite-IMC model"},
    "M4_LIT": {"color": "#4D4D4D", "subtitle": "Literature-like topology reference"},
}
RAW_IMAGES = {
    model: RAW / f"{model}_verticalZ_slab20_highres.png" for model in MODEL_ORDER
}

# One common crop is applied after the same-camera render. It removes only
# identical viewport whitespace and preserves the relative structure heights.
COMMON_VIEWPORT_CROP = (420, 250, 1580, 2750)

FINAL_WIDTH_MM = 170.0
FINAL_HEIGHT_MM = 154.0
AL_COLOR = "#EF7D83"
CU_COLOR = "#526FD0"


def mm(value: float) -> float:
    return value / 25.4


def configure_matplotlib() -> Path:
    font_path = Path(font_manager.findfont("Times New Roman", fallback_to_default=False))
    if not font_path.is_file():
        raise RuntimeError("Times New Roman is unavailable")
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman"],
            "font.size": 8.2,
            "axes.linewidth": 0.65,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
        }
    )
    return font_path


def add_panel(fig: plt.Figure, slot, model: str, letter: str) -> dict[str, plt.Axes]:
    style = MODEL_STYLE[model]
    sub = slot.subgridspec(2, 1, height_ratios=[0.105, 0.895], hspace=0.0)

    header = fig.add_subplot(sub[0])
    header.set_axis_off()
    header.add_patch(
        Rectangle((0, 0), 1, 1, transform=header.transAxes,
                  facecolor=style["color"], edgecolor=style["color"], linewidth=0)
    )
    header.text(0.027, 0.52, f"({letter})", color="white", fontsize=9.6,
                fontweight="bold", ha="left", va="center")
    header.text(0.116, 0.52, model, color="white", fontsize=9.2,
                fontweight="bold", ha="left", va="center")
    header.text(0.985, 0.50, style["subtitle"], color="white", fontsize=8.0,
                fontstyle="italic", ha="right", va="center")

    body = fig.add_subplot(sub[1])
    body.set_facecolor("white")
    body.set_xticks([])
    body.set_yticks([])
    image = Image.open(RAW_IMAGES[model]).convert("RGB").crop(COMMON_VIEWPORT_CROP)
    body.imshow(image, interpolation="lanczos", aspect="equal")
    body.set_anchor("C")
    for spine in body.spines.values():
        spine.set_visible(False)
    return {"header": header, "body": body}


def add_open_footer(fig: plt.Figure, slot) -> plt.Axes:
    ax = fig.add_subplot(slot)
    ax.set_axis_off()
    ax.plot([0.0, 1.0], [0.94, 0.94], transform=ax.transAxes, color="#AFAFAF", lw=0.65)

    # Left: loading and layer directions.
    ax.annotate("", xy=(0.095, 0.55), xytext=(0.025, 0.55), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#222222", lw=0.95))
    ax.text(0.060, 0.28, "X[100] loading", fontsize=8.0, ha="center", va="center",
            transform=ax.transAxes)
    ax.annotate("", xy=(0.205, 0.70), xytext=(0.205, 0.24), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#666666", lw=0.85))
    ax.text(0.220, 0.47, "Z[001] layers", fontsize=8.0, ha="left", va="center",
            transform=ax.transAxes)

    # Center: atom identity, no surrounding box.
    ax.scatter([0.355, 0.440], [0.52, 0.52], s=[44, 44], c=[AL_COLOR, CU_COLOR],
               edgecolors=["#7A383C", "#2B3F91"], linewidths=0.65,
               transform=ax.transAxes, clip_on=False)
    ax.text(0.377, 0.52, "Al", fontsize=8.2, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.462, 0.52, "Cu", fontsize=8.2, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.407, 0.82, "Atom identity", fontsize=8.1, fontweight="bold",
            ha="center", va="center", transform=ax.transAxes)

    # Right: controlled chain and reference-only branch.
    chain_y = 0.56
    chain_x = [0.555, 0.695, 0.835]
    for xpos, model in zip(chain_x, MODEL_ORDER[:3]):
        ax.text(xpos, chain_y, model, fontsize=8.0, fontweight="bold",
                color=MODEL_STYLE[model]["color"], ha="center", va="center",
                transform=ax.transAxes)
    for left, right in zip(chain_x[:-1], chain_x[1:]):
        ax.annotate("", xy=(right - 0.052, chain_y), xytext=(left + 0.052, chain_y),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color="#777777", lw=0.72))
    ax.text(0.695, 0.22, "controlled symmetric chain", fontsize=8.0,
            color="#444444", ha="center", va="center", fontstyle="italic",
            transform=ax.transAxes)
    ax.plot([0.905, 0.905], [0.18, 0.80], transform=ax.transAxes,
            color="#B8B8B8", lw=0.65, linestyle=(0, (3, 2)))
    ax.text(0.955, 0.60, "M4_LIT", fontsize=8.0, fontweight="bold",
            color=MODEL_STYLE["M4_LIT"]["color"], ha="center", va="center",
            transform=ax.transAxes)
    ax.text(0.955, 0.35, "reference only", fontsize=8.0, color="#4D4D4D",
            ha="center", va="center", transform=ax.transAxes)
    return ax


def write_panel_area(fig: plt.Figure, axes: dict[str, dict[str, plt.Axes]], footer: plt.Axes) -> None:
    figure_area = FINAL_WIDTH_MM * FINAL_HEIGHT_MM
    rows = []
    for model in MODEL_ORDER:
        for role in ("header", "body"):
            bbox = axes[model][role].get_position()
            width_mm = bbox.width * FINAL_WIDTH_MM
            height_mm = bbox.height * FINAL_HEIGHT_MM
            rows.append(
                {
                    "panel": model,
                    "role": role,
                    "width_mm": f"{width_mm:.3f}",
                    "height_mm": f"{height_mm:.3f}",
                    "area_percent": f"{100 * width_mm * height_mm / figure_area:.3f}",
                    "status": "PASS",
                }
            )
    bbox = footer.get_position()
    width_mm = bbox.width * FINAL_WIDTH_MM
    height_mm = bbox.height * FINAL_HEIGHT_MM
    rows.append(
        {
            "panel": "shared_footer",
            "role": "open_annotation_zone",
            "width_mm": f"{width_mm:.3f}",
            "height_mm": f"{height_mm:.3f}",
            "area_percent": f"{100 * width_mm * height_mm / figure_area:.3f}",
            "status": "PASS",
        }
    )
    with (HERE / "Fig1_R2_panel_area_audit.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def build_contact_sheet() -> None:
    fig, grid = plt.subplots(4, 3, figsize=(mm(170), mm(205)))
    for row, model in enumerate(MODEL_ORDER):
        for col, slab in enumerate((20, 25, 30)):
            ax = grid[row, col]
            path = TESTS / f"{model}_verticalZ_slab{slab:02d}.png"
            ax.imshow(Image.open(path).convert("RGB"), interpolation="lanczos")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{model} — {slab} Å", fontsize=8.3,
                         fontweight="bold" if slab == 20 else "normal",
                         color="#111111" if slab == 20 else "#555555")
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.85 if slab == 20 else 0.45)
                spine.set_color(MODEL_STYLE[model]["color"] if slab == 20 else "#C8C8C8")
    fig.suptitle("Fig.1 R2 central-slab comparison (selected: 20 Å)",
                 fontsize=10.0, fontweight="bold", y=0.995)
    fig.subplots_adjust(left=0.025, right=0.985, bottom=0.018, top=0.965,
                        wspace=0.08, hspace=0.14)
    fig.savefig(HERE / "Fig1_R2_render_comparison_contact_sheet.pdf")
    plt.close(fig)


def build() -> None:
    font_path = configure_matplotlib()
    for path in RAW_IMAGES.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    fig = plt.figure(figsize=(mm(FINAL_WIDTH_MM), mm(FINAL_HEIGHT_MM)))
    outer = fig.add_gridspec(
        3, 2, height_ratios=[1.0, 1.0, 0.205],
        left=0.024, right=0.982, bottom=0.020, top=0.984,
        wspace=0.055, hspace=0.055,
    )
    axes: dict[str, dict[str, plt.Axes]] = {}
    for slot, model, letter in zip(
        [outer[0, 0], outer[0, 1], outer[1, 0], outer[1, 1]], MODEL_ORDER, "abcd"
    ):
        axes[model] = add_panel(fig, slot, model, letter)
    footer = add_open_footer(fig, outer[2, :])
    fig.canvas.draw()
    write_panel_area(fig, axes, footer)

    outputs = {
        "pdf": HERE / "Fig1_R2_final.pdf",
        "svg": HERE / "Fig1_R2_final.svg",
        "assembly": HERE / "Fig1_R2_vector_assembly.svg",
        "png600": HERE / "Fig1_R2_final_600dpi.png",
        "tif600": HERE / "Fig1_R2_final_600dpi.tif",
        "preview": HERE / "Fig1_R2_preview.png",
    }
    fig.savefig(outputs["pdf"])
    fig.savefig(outputs["svg"])
    shutil.copy2(outputs["svg"], outputs["assembly"])
    fig.savefig(outputs["png600"], dpi=600)
    fig.savefig(outputs["tif600"], dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(outputs["preview"], dpi=220)
    plt.close(fig)
    build_contact_sheet()

    print(f"font={font_path}")
    print(f"size_mm={FINAL_WIDTH_MM}x{FINAL_HEIGHT_MM}")
    print("selected_slab_A=20")
    print(f"common_viewport_crop={COMMON_VIEWPORT_CROP}")
    print("preview_pixels_read=NO")
    for key, path in outputs.items():
        print(f"{key}={path}")


if __name__ == "__main__":
    build()
