from pathlib import Path
import pymupdf as pdf
folder = Path(__file__).resolve().parent
source = folder / "Fig5_partial.pdf"
doc = pdf.open(source)
page = doc[0]
boxes = page.search_for("B  Program embedding in the OSCC network")
assert len(boxes) == 1
box = boxes[0] + (-1, -1, 1, 1)
assert box.y1 < 25
page.add_redact_annot(box, fill=(1, 1, 1))
page.apply_redactions(images=0, graphics=0)
page.insert_text((331.0234, 18.1), "B  Program embedding in the TCGA-HNSC network", fontname="hebo", fontsize=10, color=(.13, .14, .15))
output = folder / "Fig5_repaired.pdf"
doc.save(output, garbage=4, deflate=True)
with pdf.open(source) as original, pdf.open(output) as fixed:
    area = pdf.Rect(0, 26, original[0].rect.width, original[0].rect.height)
    assert original[0].get_pixmap(clip=area).samples == fixed[0].get_pixmap(clip=area).samples
    fixed[0].get_pixmap(dpi=300).save(folder / "Fig5_repaired.png")
