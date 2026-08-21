from pathlib import Path
import csv
# Replot-only provenance wrapper. No analysis is performed.
SOURCE = Path(r"D:\leng\AlCu_MAIN_PPP_Zhou\cms_manuscript_finalization\09_text_references_figures_integration\stage3_figure_rebuild\05_stageB_fig6_fig7_relayout_replot\04_Fig7_stageB\Fig7_quant_replot_source.csv")
MODEL_COLORS = {"M3_SYM":"#0072B2","M4_RATIO":"#E69F00","M4_SYM":"#7B2CBF"}
MODEL_MARKERS = {"M3_SYM":"o","M4_RATIO":"s","M4_SYM":"^"}
with SOURCE.open(encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
print(f"Loaded {len(rows)} locked rows from {SOURCE}; plotting parameters are registered in the Stage B source manifest.")
