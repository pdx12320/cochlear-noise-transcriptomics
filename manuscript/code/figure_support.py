from pathlib import Path
import math, json, shutil, subprocess, html, os
import reportlab
import numpy as np
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'source_inputs'; OUT=ROOT/'figures'; DAT=ROOT/'source_data'
# Display artifact numbers follow first citation; source-table identifiers stay fixed.
_figure_map = pd.read_csv(ROOT/'SUPPLEMENTARY_FIGURE_MAP.tsv', sep='\t')
FIGURE_OUTPUT_NAMES = dict(zip(_figure_map.original_basename, _figure_map.final_basename))

# Standalone LaTeX bundle and repository manuscript layout are both supported.
BULK_ROOT = ROOT/'reanalysis/GSE160639'
if not BULK_ROOT.is_dir():
    BULK_ROOT = ROOT.parent/'analysis/GSE160639'
OUT.mkdir(parents=True,exist_ok=True); DAT.mkdir(parents=True,exist_ok=True)
# Optional FIGURE_FONT_DIR points to a directory containing Arial font files.
# Standard ReportLab Vera fonts provide a portable fallback when Arial is absent.
_font_dirs = [Path(os.environ['FIGURE_FONT_DIR'])] if os.environ.get('FIGURE_FONT_DIR') else []
_font_dirs.append(Path('/System/Library/Fonts/Supplemental'))
_font_files = ['Arial.ttf', 'Arial Bold.ttf', 'Arial Italic.ttf']
_font_dir = next((p for p in _font_dirs if all((p/n).is_file() for n in _font_files)), None)
if _font_dir is not None:
    FONT_FAMILY = 'Arial'
else:
    _font_dir = Path(reportlab.__file__).resolve().parent/'fonts'
    _font_files = ['Vera.ttf', 'VeraBd.ttf', 'VeraIt.ttf']
    FONT_FAMILY = 'Bitstream Vera Sans'
for _alias, _filename in zip(['Arial', 'Arial-Bold', 'Arial-Italic'], _font_files):
    pdfmetrics.registerFont(TTFont(_alias, str(_font_dir/_filename)))

# Poppler may be supplied through PDFTOPPM or found on PATH.
_pdftoppm_candidates = [os.environ.get('PDFTOPPM'), shutil.which('pdftoppm'),
    str(Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm')]
PDFTOPPM = next((p for p in _pdftoppm_candidates if p and Path(p).is_file()), None)
if PDFTOPPM is None:
    raise RuntimeError('Poppler pdftoppm is required for PNG export. Set PDFTOPPM to its executable path.')
BLUE='#377DA8'; RED='#CA5A49'; LIGHT='#D5DADF'; GRAY='#7B858E'; BLACK='#171717'
colors=['#4477AA','#66A699','#CCBB44','#AA4499','#EE8866','#77AADD','#AA7744','#999933','#B8BEC5']
types=['Endothelial','Epithelial / supporting-associated','Fibrocyte / mesenchymal','Glial / Schwann-like','Immune / myeloid','Melanocyte / intermediate','Neuronal / SGN','Pericyte / vascular smooth muscle','UNRESOLVED']
col=dict(zip(types,colors)); short=dict(zip(types,['Endothelial','Supporting epithelial','Fibrocyte / mesenchymal','Glial / Schwann-like','Immune / myeloid','Melanocyte / intermediate','Neuronal','Pericyte / vascular muscle','Unassigned']))

def file(name): return next(SRC.rglob(name))
def read(name): return pd.read_csv(file(name),sep='\t' if name.endswith('.tsv') else ',')
def save_table(df,name):df.to_csv(DAT/name,sep='\t',index=False)
def rgb(c):return tuple(int(c[i:i+2],16)/255 for i in (1,3,5))
def mix(a,b,t):
    av=rgb(a);bv=rgb(b);return '#'+''.join(f'{round((x+(y-x)*t)*255):02x}' for x,y in zip(av,bv))
def div(v,lim=3.3):return mix('#F5F5F3',RED if v>0 else BLUE,min(abs(v)/lim,1))
def seq(v,lim=3.1):return mix('#E6EEF3','#225C86',min(v/lim,1))
def fmtq(v):return f'{v:.1e}' if v<.001 else f'{v:.3f}'
class Fig:
 def __init__(self,name,w=510,h=500):
  self.name=name;self.output_name=FIGURE_OUTPUT_NAMES.get(name,name);self.w=w;self.h=h;self.pdf=canvas.Canvas(str(OUT/(self.output_name+'.pdf')),pagesize=(w,h));self.pdf.setTitle(name);self.pdf.setAuthor('Dengxiang Yang');self.svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w/72:.6f}in" height="{h/72:.6f}in" viewBox="0 0 {w} {h}">','<rect width="100%" height="100%" fill="white"/>']
 def text(self,x,y,s,size=8,color=BLACK,bold=False,italic=False,align='left',rotate=0):
  s=str(s); font='Arial-Bold' if bold else 'Arial-Italic' if italic else 'Arial';p=self.pdf
  p.saveState();p.setFillColorRGB(*rgb(color));p.setFont(font,size);p.translate(x,self.h-y);p.rotate(-rotate)
  (p.drawCentredString if align=='center' else p.drawRightString if align=='right' else p.drawString)(0,0,s);p.restoreState()
  anchor={'left':'start','center':'middle','right':'end'}[align]
  self.svg.append(f'<text x="{x}" y="{y}" font-family="{FONT_FAMILY}, sans-serif" font-size="{size}" font-weight="{"bold" if bold else "normal"}" font-style="{"italic" if italic else "normal"}" text-anchor="{anchor}" fill="{color}" transform="rotate({rotate} {x} {y})">{html.escape(s)}</text>')
 def line(self,x1,y1,x2,y2,color='#90989E',width=.5,dash=None):
  p=self.pdf;p.saveState();p.setStrokeColorRGB(*rgb(color));p.setLineWidth(width)
  if dash:p.setDash(dash)
  p.line(x1,self.h-y1,x2,self.h-y2);p.restoreState();self.svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}" '+(f'stroke-dasharray="{dash},{dash}"' if dash else '')+'/>')
 def rect(self,x,y,w,h,fill,stroke=None,lw=.5):
  p=self.pdf;p.saveState();p.setFillColorRGB(*rgb(fill));p.setLineWidth(lw)
  if stroke:p.setStrokeColorRGB(*rgb(stroke))
  p.rect(x,self.h-y-h,w,h,fill=1,stroke=int(stroke is not None));p.restoreState();self.svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{lw}"/>')
 def circle(self,x,y,r,fill,stroke=None,lw=.5):
  p=self.pdf;p.saveState();p.setFillColorRGB(*rgb(fill));p.setLineWidth(lw)
  if stroke:p.setStrokeColorRGB(*rgb(stroke))
  p.circle(x,self.h-y,r,fill=1,stroke=int(stroke is not None));p.restoreState();self.svg.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{lw}"/>')
 def panel(self,label,title,x,y):self.text(x,y,label.lower(),11,bold=True);self.text(x+16,y,title,9,bold=True)
 def finish(self):
  self.pdf.showPage();self.pdf.save();self.svg.append('</svg>');(OUT/(self.output_name+'.svg')).write_text('\n'.join(self.svg))
  subprocess.run([PDFTOPPM,'-png','-r','600','-singlefile',str(OUT/(self.output_name+'.pdf')),str(OUT/self.output_name)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
  print(self.output_name,flush=True)

def legend(f,items,x,y,cols=3,dx=160,dy=12,size=7):
 for i,(name,color) in enumerate(items):
  xx=x+(i%cols)*dx;yy=y+(i//cols)*dy;f.rect(xx,yy-6,6,6,color);f.text(xx+9,yy,name,size)
def h_axis(f,x0,x1,y,ticks,lo,hi,label='',label_y=None):
 f.line(x0,y,x1,y,BLACK,.6)
 for v in ticks:
  x=x0+(v-lo)/(hi-lo)*(x1-x0);f.line(x,y,x,y+3,BLACK,.5);f.text(x,y+13,f'{v:g}',7,align='center')
 if label:f.text((x0+x1)/2,label_y or y+27,label,8,align='center')
def y_axis(f,x,y0,y1,ticks,lo,hi,label=''):
 f.line(x,y0,x,y1,BLACK,.6)
 for v in ticks:
  y=y1-(v-lo)/(hi-lo)*(y1-y0);f.line(x-3,y,x,y,BLACK,.5);f.text(x-5,y+2,f'{v:g}',7,align='right')
 if label:f.text(x-29,(y0+y1)/2,label,7.5,align='center',rotate=-90)
