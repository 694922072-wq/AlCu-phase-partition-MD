from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _supp_utils import FIG_DIR, OVERLAP_DIR, TABLE_DIR, event_label, set_matplotlib_style, write_run_status


def main() -> None:
    import matplotlib.pyplot as plt

    set_matplotlib_style()
    stats = pd.read_csv(TABLE_DIR / "layerwise_stats_selected_keyframes_x25_clean.csv")
    cols = [
        "model",
        "event_name",
        "event_order",
        "actual_strain_percent",
        "step",
        "frame_index",
        "region_name",
        "region_type",
        "source",
        "z_min",
        "z_max",
        "atom_count",
        "high_d2min_count",
        "high_pe_count",
        "overlap_count",
        "overlap_ratio_to_highD2min",
        "overlap_ratio_to_highPE",
        "overlap_ratio_to_region_atoms",
    ]
    overlap = stats[cols].copy()
    overlap.to_csv(OVERLAP_DIR / "pe_d2min_overlap_by_region.csv", index=False)

    key = overlap[(overlap["event_name"].isin(["15pct", "peak", "first_50_drop", "20pct", "25pct"])) & (overlap["region_name"] == "whole_model")]
    fig, ax = plt.subplots(figsize=(7.5, 4.4), dpi=180)
    for model, sub in key.groupby("model"):
        sub = sub.sort_values("event_order")
        ax.plot(sub["actual_strain_percent"], sub["overlap_ratio_to_region_atoms"], marker="o", lw=1.8, label=model)
    ax.set_xlabel("Engineering strain (%)")
    ax.set_ylabel("Overlap / all atoms")
    ax.set_title("PE-D2min high-value overlap by model")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "PE_D2min_overlap_by_model_event.png", dpi=300)
    plt.close(fig)

    sub = overlap[overlap["source"].isin(["layer", "interface_10A"])].copy()
    peak = sub[sub["event_name"] == "peak"].set_index(["model", "region_name"])
    drop = sub[sub["event_name"] == "first_50_drop"].set_index(["model", "region_name"])
    common = peak.index.intersection(drop.index)
    delta = drop.loc[common][["source", "overlap_ratio_to_region_atoms"]].copy()
    delta["delta_overlap_peak_to_first50"] = (
        drop.loc[common]["overlap_ratio_to_region_atoms"] - peak.loc[common]["overlap_ratio_to_region_atoms"]
    )
    delta = delta.reset_index()
    # Show regions with the largest rise in coupled PE-D2min high-value atoms.
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), dpi=180)
    for ax, (model, g) in zip(axes, delta.groupby("model")):
        top = g.sort_values("delta_overlap_peak_to_first50", ascending=False).head(8)
        ax.barh(top["region_name"], top["delta_overlap_peak_to_first50"], color="#4c78a8")
        ax.invert_yaxis()
        ax.set_title(model)
        ax.set_xlabel("Delta overlap / region atoms")
        ax.tick_params(axis="y", labelsize=7)
    fig.suptitle("Peak to first_50_drop: PE-D2min overlap increase")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "PE_D2min_overlap_peak_to_drop.png", dpi=300)
    plt.close(fig)
    write_run_status("04_pe_d2min_overlap_analysis.py: success")
    print("04: PE-D2min overlap analysis complete")


if __name__ == "__main__":
    main()
