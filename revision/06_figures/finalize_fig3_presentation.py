from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import argparse, csv

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--counts", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--font-regular", required=True)
parser.add_argument("--font-bold", required=True)
args = parser.parse_args()
with open(args.counts, encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
vals = [int(row["genes"]) for row in rows]
assert vals == [43, 32, 54]

def repair_image(im):
 im=im.convert('RGB');w,h=im.size;sx=w/1600;sy=h/1020
 box=tuple(round(v*(sx if i%2==0 else sy)) for i,v in enumerate((1073,20,1595,413)))
 pw,ph=box[2]-box[0],box[3]-box[1]
 patch=Image.new('RGB',(pw,ph),'white');pd=ImageDraw.Draw(patch);k=pw/522
 font=lambda size,bold=False:ImageFont.truetype(args.font_bold if bold else args.font_regular,round(size*k))
 def xy(x,y):return (round(x*k),round(y*ph/393))
 def line(coords):pd.line([xy(*p) for p in coords],fill='black',width=max(1,round(1.4*k)))
 pd.text(xy(75,12),'Shared DEG directions',font=font(20,True),fill='black')
 line([(75,45),(75,295),(495,295)])
 for tick in [0,20,40,60]:
  y=295-tick*250/65;line([(70,y),(75,y)])
  pd.text(xy(65,y),str(tick),font=font(15),fill='black',anchor='rm')
 for j,(v,label) in enumerate(zip(vals,['Down in both','Opposite','Up in both'])):
  x=145+140*j;y=295-v*250/65
  pd.rectangle([xy(x-42,y),xy(x+42,295)],fill='#B34E43')
  pd.text(xy(x,y-7),str(v),font=font(18),fill='black',anchor='mb')
  lab=Image.new('RGBA',(round(200*k),round(38*k)),(255,255,255,0));ld=ImageDraw.Draw(lab)
  ld.text((round(195*k),round(3*k)),label,font=font(15),fill='black',anchor='rt')
  lab=lab.rotate(30,expand=True,resample=Image.Resampling.BICUBIC)
  patch.paste(lab,(round((x+10)*k)-lab.width,round(304*ph/393)),lab)
 ylabel=Image.new('RGBA',(round(190*k),round(30*k)),(255,255,255,0));ImageDraw.Draw(ylabel).text((0,0),'Number of genes',font=font(16),fill='black');ylabel=ylabel.rotate(90,expand=True)
 patch.paste(ylabel,(round(6*k),round(118*ph/393)),ylabel)
 im.paste(patch,box)
 draw=ImageDraw.Draw(im)
 def label(rect,txt,size,bold=False,color='black'):
  r=tuple(round(v*(sx if i%2==0 else sy)) for i,v in enumerate(rect))
  draw.rectangle(r,fill='white')
  font=ImageFont.truetype(args.font_bold if bold else args.font_regular,round(size*sx))
  draw.text((r[0],r[1]),txt,font=font,fill=color)
 label((1033,442,1590,468),'Concordant Hallmark enrichment in resistant cells',16,True)
 label((640,62,1022,86),'Independent acquired-resistance model',16,False,'#555555')
 return im


im = Image.open(args.input)
dpi = im.info.get("dpi", (600, 600))
result = repair_image(im)
if Path(args.output).suffix.lower() in (".tif", ".tiff"):
    result.save(args.output, compression="tiff_lzw", dpi=dpi)
else:
    result.save(args.output, dpi=dpi)
