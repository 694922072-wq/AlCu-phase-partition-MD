"""Export the editable Fig5 SVG at the exact legacy raster dimensions."""

from pathlib import Path

import fitz
from PIL import Image

HERE = Path(__file__).resolve().parent
SVG = HERE / "Fig5_n3.svg"
PNG = HERE / "Fig5_n3_600dpi.png"
TIF = HERE / "Fig5_n3_600dpi.tif"
PREVIEW = HERE / "Fig5_n3_preview.png"
PDF = HERE / "Fig5_n3.pdf"
TARGET_SIZE = (4015, 2952)

doc = fitz.open(SVG)
page = doc[0]
pix = page.get_pixmap(matrix=fitz.Matrix(600 / 72, 600 / 72), alpha=False)
tmp = HERE / "_Fig5_n3_render_tmp.png"
pix.save(tmp)
with Image.open(tmp) as rendered:
    exact = rendered.resize(TARGET_SIZE, Image.Resampling.LANCZOS)
    exact.save(PNG, dpi=(600, 600))
    exact.save(TIF, dpi=(600, 600), compression="tiff_lzw")
    preview = exact.copy()
    preview.thumbnail((1600, 1600), Image.Resampling.LANCZOS)
    preview.save(PREVIEW)
tmp.unlink()
PDF.write_bytes(doc.convert_to_pdf())
print(f"svg_points={page.rect.width}x{page.rect.height}")
print(f"raster_pixels={TARGET_SIZE[0]}x{TARGET_SIZE[1]}")
print(f"png={PNG}")
print(f"pdf={PDF}")
