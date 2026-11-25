import os
pass
import csv
from pathlib import Path
import pymupdf as pdf
from PIL import Image

folder=Path(__file__).resolve().parent
source=next(Path((os.environ["OSCC_DATA"] + "/Integration_Manuscript/Figure5_SharedNetworkArchitecture")).glob("*FIXED.pdf"))
with (folder.parent/"11_network/all_score_pairs_by_module.tsv").open() as f:
    data=list(csv.DictReader(f,delimiter="\t"))
assert len(data)==105
pairs=[("Stromal5","Ultra7","MEturquoise"),("Proliferation2","Final8Risk","MEblack"),
       ("Proliferation2","Ultra7","MEblack")]
rows=[next(r for r in data if (r["score1"],r["score2"],r["module"])==key) for key in pairs]
original=pdf.open(source); doc=pdf.open(source); page=doc[0]
panel=pdf.Rect(0,305,395,576)
page.add_redact_annot(panel,fill=(1,1,1))
page.apply_redactions(graphics=pdf.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED)

page.draw_rect(panel,color=None,fill=(1,1,1))
dark=(.16,.17,.18); blue=(.19,.42,.64); red=(.73,.26,.24); grey=(.72,.74,.76)
def txt(x,y,s,size=8,bold=False,align="left",color=dark):
    font="hebo" if bold else "helv"
    w=pdf.get_text_length(s,fontname=font,fontsize=size)
    if align=="center": x-=w/2
    if align=="right": x-=w
    page.insert_text((x,y),s,fontname=font,fontsize=size,color=color)
def coord(r): return 128+(r+.35)*262/.95
txt(12,322,"C  Sensitivity to a shared co-expression module",10,True)
page.draw_line((coord(0),339),(coord(0),492),color=grey,width=.6,dashes="[3 3] 0")
for i,(r,y) in enumerate(zip(rows,[362,416,470])):
    a,b=[float(r[k]) for k in ("raw_rho","partial_rho")]
    label=r["score2"].replace("Final8Risk","Final8")
    txt(119,y-2,r["score1"]+" / "+label,7.7,align="right")
    txt(119,y+10,r["module"][2:]+"; n="+r["N"],7,align="right",color=(.32,.34,.36))
    page.draw_line((coord(a),y-2),(coord(b),y+2),color=grey,width=.9)
    page.draw_circle((coord(a),y-2),2.7,color=blue,fill=blue)
    x=coord(b)
    page.draw_rect(pdf.Rect(x-2.4,y-.4,x+2.4,y+4.4),color=red,fill=red)
    txt(coord(a),y-9,f"{a:.2f}",7.5,align="center",color=blue)
    txt(coord(b),y+15,f"{b:.2f}",7.5,align="center",color=red)
page.draw_line((coord(-.35),499),(coord(.60),499),color=dark,width=.65)
for v in (-.2,0,.2,.4,.6):
    x=coord(v);page.draw_line((x,499),(x,502),color=dark,width=.65)
    txt(x,513,f"{v:.1f}",7.7,align="center")
txt(252,531,"Rank correlation coefficient",8,align="center")
page.draw_circle((149,551),2.7,color=blue,fill=blue);txt(157,554,"Raw",8)
page.draw_rect(pdf.Rect(200,548.6,204.8,553.4),color=red,fill=red)
txt(211,554,"Module-adjusted",8)
dest=folder/"Fig5_partial.pdf"
doc.save(dest,garbage=4,deflate=True)
with pdf.open(dest) as check:

    for box in (pdf.Rect(0,0,734,304),pdf.Rect(396,305,734,576)):
        assert check[0].get_pixmap(clip=box).samples==original[0].get_pixmap(clip=box).samples
    for r in rows:
        assert f"{float(r['partial_rho']):.2f}" in check[0].get_text(clip=panel)
    for dpi,suffix in ((300,".png"),(100,"_preview.png")):
        check[0].get_pixmap(dpi=dpi).save(folder/("Fig5_partial"+suffix))
print("Panel C recomputed; unchanged panels passed pixel comparison.")
