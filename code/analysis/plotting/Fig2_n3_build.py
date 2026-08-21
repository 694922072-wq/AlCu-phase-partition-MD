"""Build manuscript Fig. 2 from the authoritative nine-trajectory manifest.

Individual trajectories and all event metrics use unsmoothed raw records.
The nine controlled records share one output-step grid. Display means and
sample-SD bands are therefore computed pointwise by common row without
interpolation. Only the display mean before 15.5% strain receives the legacy
centered 31-row moving average. M4_LIT remains an n=1 dashed reference.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REVISION_ROOT = HERE.parent
METRICS_PATH = REVISION_ROOT / "02_local_statistics" / "metrics_per_trajectory_local.csv"
SUMMARY_PATH = REVISION_ROOT / "02_local_statistics" / "n3_summary_local.csv"
M4_LIT_PATH = Path(r"D:\leng\AlCu_MAIN_PPP_Zhou\outputs\M4_LIT_stress_strain_x25_clean.dat")

MODEL_ORDER = ["M3_SYM", "M4_RATIO", "M4_SYM", "M4_LIT"]
CONTROLLED = MODEL_ORDER[:3]
MODEL_STYLE = {
    "M3_SYM": {"color": "#0072B2", "marker": "o", "linestyle": "-"},
    "M4_RATIO": {"color": "#E69F00", "marker": "s", "linestyle": "-"},
    "M4_SYM": {"color": "#7B2CBF", "marker": "^", "linestyle": "-"},
    "M4_LIT": {"color": "#4D4D4D", "marker": "D", "linestyle": "--"},
}
M4_LIT_METRICS = {
    "sigma15_GPa": 7.64304336869481,
    "peak_stress_GPa": 8.62305963183287,
    "peak_strain_percent": 18.24,
    "F50_raw_first_row_percent": 18.66,
}

FINAL_WIDTH_MM = 170.0
FINAL_HEIGHT_MM = 102.0
SMOOTH_WINDOW = 31
SMOOTH_END_PERCENT = 15.5


def mm(value: float) -> float:
    return value / 25.4


def configure_matplotlib() -> Path:
    font_path = Path(font_manager.findfont("Times New Roman", fallback_to_default=False))
    if not font_path.is_file():
        raise RuntimeError("Times New Roman is unavailable")
    mpl.rcParams.update({
        "font.family": "Times New Roman", "font.serif": ["Times New Roman"],
        "font.size": 8.2, "axes.labelsize": 9.1, "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0, "legend.fontsize": 7.6, "axes.linewidth": 0.78,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "mathtext.fontset": "custom", "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic", "mathtext.bf": "Times New Roman:bold",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "savefig.transparent": False,
    })
    return font_path


def centered_moving_average(values: np.ndarray, window: int) -> np.ndarray:
    half = window // 2
    padded = np.pad(values, (half, half), mode="reflect")
    return np.convolve(padded, np.full(window, 1.0 / window), mode="valid")


def panel_label(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.145, 1.005, f"({letter})", transform=ax.transAxes, ha="left",
            va="bottom", fontsize=9.5, fontweight="bold", clip_on=False)


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=2.7, width=0.72, pad=2.0)


def load_inputs():
    metrics = pd.read_csv(METRICS_PATH)
    summary = pd.read_csv(SUMMARY_PATH).set_index("model")
    if len(metrics) != 9 or set(metrics["model"]) != set(CONTROLLED):
        raise RuntimeError("Authoritative n=3 metric table mismatch")
    if not (metrics["completion_status"] == "PASS").all():
        raise RuntimeError("A trajectory failed completion audit")
    trajectories: dict[str, list[pd.DataFrame]] = {m: [] for m in CONTROLLED}
    for row in metrics.itertuples(index=False):
        path = Path(row.local_raw_path)
        frame = pd.read_csv(path, sep=r"\s+")
        if len(frame) != 2501 or abs(float(frame["strain_percent"].iloc[-1]) - 25.0) > 1e-8:
            raise RuntimeError(f"Incomplete raw trajectory: {path}")
        trajectories[row.model].append(frame)
    for model in CONTROLLED:
        reference_steps = trajectories[model][0]["step"].to_numpy(np.int64)
        for frame in trajectories[model][1:]:
            if not np.array_equal(frame["step"].to_numpy(np.int64), reference_steps):
                raise RuntimeError(f"Non-common step grid detected for {model}")
    lit = pd.read_csv(M4_LIT_PATH, sep=r"\s+")
    return metrics, summary, trajectories, lit


def metric_panel(ax: plt.Axes, letter: str, title: str, xlabel: str,
                 metric_col: str, summary: pd.DataFrame, metrics: pd.DataFrame,
                 xlim: tuple[float, float], decimals: int, suffix: str = "",
                 peak_labels: bool = False) -> None:
    y_positions = np.arange(len(MODEL_ORDER))[::-1]
    span = xlim[1] - xlim[0]
    for ypos, model in zip(y_positions, MODEL_ORDER):
        style = MODEL_STYLE[model]
        if model in CONTROLLED:
            vals = metrics.loc[metrics["model"].eq(model), metric_col].to_numpy(float)
            mean = float(summary.loc[model, f"{metric_col}_mean"])
            sd = float(summary.loc[model, f"{metric_col}_sample_SD"])
            jitter = np.array([-0.12, 0.0, 0.12])
            ax.scatter(vals, ypos + jitter, s=13, marker=style["marker"],
                       facecolor="white", edgecolor=style["color"], linewidth=0.72,
                       alpha=0.92, zorder=4)
            ax.errorbar(mean, ypos, xerr=sd, fmt=style["marker"], markersize=5.5,
                        color=style["color"], markerfacecolor=style["color"],
                        markeredgecolor="white", markeredgewidth=0.45,
                        elinewidth=1.15, capsize=2.4, capthick=1.0, zorder=6)
            if peak_labels:
                pm = float(summary.loc[model, "peak_strain_percent_mean"])
                ps = float(summary.loc[model, "peak_strain_percent_sample_SD"])
                label = f"{mean:.{decimals}f}±{sd:.{decimals}f} ({pm:.2f}±{ps:.2f}%)"
            else:
                label = f"{mean:.{decimals}f}±{sd:.{decimals}f}{suffix}"
        else:
            value = M4_LIT_METRICS[metric_col]
            ax.scatter(value, ypos, s=27, marker=style["marker"], facecolor="white",
                       edgecolor=style["color"], linewidth=0.85, zorder=5)
            if peak_labels:
                label = f"{value:.{decimals}f} ({M4_LIT_METRICS['peak_strain_percent']:.2f}%)"
            else:
                label = f"{value:.{decimals}f}{suffix}"
            mean = value
        ax.hlines(ypos, xlim[0], mean, color="#D7D7D7", linewidth=0.72, zorder=1)
        ax.text(xlim[1] + 0.015 * span, ypos, label, fontsize=6.9, color="#3D3D3D",
                ha="left", va="center", clip_on=False)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.7, 3.7)
    ax.set_yticks(y_positions, MODEL_ORDER)
    for tick, model in zip(ax.get_yticklabels(), MODEL_ORDER):
        tick.set_color(MODEL_STYLE[model]["color"])
        tick.set_fontweight("bold")
    ax.set_title(title, loc="left", fontsize=8.7, fontweight="bold", pad=4.0)
    ax.set_xlabel(xlabel, labelpad=1.8)
    ax.grid(axis="x", color="#ECECEC", linewidth=0.55, zorder=0)
    clean_axes(ax)
    panel_label(ax, letter)


def write_figure_manifest(metrics: pd.DataFrame, summary: pd.DataFrame) -> None:
    rows = []
    for row in metrics.itertuples(index=False):
        rows.append({"model": row.model, "replica": row.replica, "seed": row.seed,
                     "raw_path": row.local_raw_path, "raw_sha256": row.raw_sha256,
                     "sigma15_GPa": row.sigma15_GPa, "peak_stress_GPa": row.peak_stress_GPa,
                     "peak_strain_percent": row.peak_strain_percent,
                     "F50_percent": row.F50_raw_first_row_percent})
    with (HERE / "Fig2_n3_source_manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    summary.reset_index().to_csv(HERE / "Fig2_n3_summary_values.csv", index=False, encoding="utf-8-sig")


def build() -> None:
    font_path = configure_matplotlib()
    metrics, summary, trajectories, lit = load_inputs()
    write_figure_manifest(metrics, summary)

    fig = plt.figure(figsize=(mm(FINAL_WIDTH_MM), mm(FINAL_HEIGHT_MM)))
    grid = fig.add_gridspec(3, 2, width_ratios=[2.18, 1.0], height_ratios=[1, 1, 1],
                           left=0.077, right=0.790, bottom=0.108, top=0.950,
                           wspace=0.30, hspace=0.52)
    main = fig.add_subplot(grid[:, 0])
    for model in CONTROLLED:
        style = MODEL_STYLE[model]
        pointwise_stress = []
        pointwise_strain = []
        for frame in trajectories[model]:
            x = frame["strain_percent"].to_numpy(float)
            y = frame["stress_xx_GPa"].to_numpy(float)
            main.plot(x, y, color=style["color"], linewidth=0.55, alpha=0.28, zorder=2)
            pointwise_strain.append(x)
            pointwise_stress.append(y)
        display_strain = np.mean(np.vstack(pointwise_strain), axis=0)
        stacked = np.vstack(pointwise_stress)
        mean_raw = np.mean(stacked, axis=0)
        sd_raw = np.std(stacked, axis=0, ddof=1)
        main.fill_between(display_strain, mean_raw - sd_raw, mean_raw + sd_raw,
                          color=style["color"], alpha=0.10, linewidth=0.0, zorder=1)
        mean_smooth = centered_moving_average(mean_raw, SMOOTH_WINDOW)
        pre = display_strain <= SMOOTH_END_PERCENT
        post = display_strain >= SMOOTH_END_PERCENT
        main.plot(display_strain[pre], mean_smooth[pre], color=style["color"],
                  linewidth=1.80, alpha=0.98, zorder=5)
        main.plot(display_strain[post], mean_raw[post], color=style["color"],
                  linewidth=1.45, alpha=0.90, zorder=4)
    lx = lit["strain_percent"].to_numpy(float)
    ly = lit["stress_xx_GPa"].to_numpy(float)
    main.plot(lx, ly, color=MODEL_STYLE["M4_LIT"]["color"], linestyle="--",
              linewidth=1.25, alpha=0.88, zorder=3)
    main.axvline(15.0, color="#A5A5A5", linestyle=(0, (3, 2)), linewidth=0.78, zorder=0)
    main.text(15.18, 0.22, "15%", fontsize=8.0, color="#666666", rotation=90, va="bottom")
    main.set_xlim(0, 25); main.set_ylim(-0.5, 9.4)
    main.set_xlabel("Engineering strain (%)")
    main.set_ylabel(r"Tensile stress, $\sigma_{xx}$ (GPa)")
    main.grid(color="#EBEBEB", linewidth=0.52, zorder=0)
    clean_axes(main); panel_label(main, "a")
    handles = [Line2D([0], [0], color=MODEL_STYLE[m]["color"], lw=1.7,
                      linestyle=MODEL_STYLE[m]["linestyle"], label=m) for m in MODEL_ORDER]
    main.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.05, 0.995), ncol=2,
                frameon=False, handlelength=2.15, columnspacing=0.95,
                labelspacing=0.25, handletextpad=0.45)

    metric_panel(fig.add_subplot(grid[0, 1]), "b", "Stress at 15% strain", "Stress (GPa)",
                 "sigma15_GPa", summary, metrics, (7.52, 8.00), 3)
    metric_panel(fig.add_subplot(grid[1, 1]), "c", "Peak response", "Peak stress (GPa)",
                 "peak_stress_GPa", summary, metrics, (8.40, 8.95), 3, peak_labels=True)
    metric_panel(fig.add_subplot(grid[2, 1]), "d", "F50 transition", "F50 strain (%)",
                 "F50_raw_first_row_percent", summary, metrics, (17.35, 19.55), 2, suffix="%")

    outputs = {
        "pdf": HERE / "Fig2_n3.pdf", "svg": HERE / "Fig2_n3.svg",
        "png600": HERE / "Fig2_n3_600dpi.png", "tif600": HERE / "Fig2_n3_600dpi.tif",
        "preview": HERE / "Fig2_n3_preview.png",
    }
    fig.savefig(outputs["pdf"]); fig.savefig(outputs["svg"])
    fig.savefig(outputs["png600"], dpi=600)
    fig.savefig(outputs["tif600"], dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(outputs["preview"], dpi=220)
    plt.close(fig)
    print(f"font={font_path}")
    print(f"size_mm={FINAL_WIDTH_MM}x{FINAL_HEIGHT_MM}")
    print("individual_trajectories=RAW_UNSMOOTHED")
    print("mean_prepeak=centered_31_row_display_smoothing_until_15.5_percent")
    print("mean_and_sample_SD=POINTWISE_BY_COMMON_STEP_NO_INTERPOLATION")
    print("mean_postpeak=RAW_POINTWISE_MEAN")
    print("event_extraction=PER_TRAJECTORY_RAW_ROWS_ONLY")
    for key, path in outputs.items(): print(f"{key}={path}")


if __name__ == "__main__":
    build()
