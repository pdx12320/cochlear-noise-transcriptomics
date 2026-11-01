from pathlib import Path
import argparse, io
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

_parser=argparse.ArgumentParser(description='Render figures from archived result tables.')
_parser.add_argument('--bulk-only',action='store_true')
_parser.add_argument('--figures',help='Comma-separated figure basenames to export; other artwork and source tables are left unchanged.')
_args=_parser.parse_args()
_final_to_source={final:original for original,final in FIGURE_OUTPUT_NAMES.items()}
_selected_figures={_final_to_source.get(n,n) for n in _args.figures.split(',')} if _args.figures else None
_known_figures={'Figure2_cell_type_pathway_responses','Figure4_noise_exposure_comparisons','Figure5_condition_expression','FigureS1_bulk_volcanoes','FigureS2_all_twenty_genes','FigureS3_localization_and_cell_state','FigureS4_cell_type_landscape'}
if _selected_figures is not None and not _selected_figures.issubset(_known_figures):
 raise ValueError('Unknown figure basename in --figures')

def file(name): return next(SRC.rglob(name))
def read(name): return pd.read_csv(file(name),sep='\t' if name.endswith('.tsv') else ',')
def save_table(df,name):
 if _selected_figures is None or any(name.startswith(n.split('_',1)[0]) for n in _selected_figures):
  df.to_csv(DAT/name,sep='\t',index=False)
def rgb(c):return tuple(int(c[i:i+2],16)/255 for i in (1,3,5))
def mix(a,b,t):
    av=rgb(a);bv=rgb(b);return '#'+''.join(f'{round((x+(y-x)*t)*255):02x}' for x,y in zip(av,bv))
def div(v,lim=3.3):return mix('#F5F5F3',RED if v>0 else BLUE,min(abs(v)/lim,1))
def seq(v,lim=3.1):return mix('#E6EEF3','#225C86',min(v/lim,1))
def fmtq(v):return f'{v:.1e}' if v<.001 else f'{v:.3f}'
class Fig:
 def __init__(self,name,w=510,h=500):
  self.name=name;self.output_name=FIGURE_OUTPUT_NAMES.get(name,name);self.w=w;self.h=h;self.selected=_selected_figures is None or name in _selected_figures;self.pdf=canvas.Canvas(str(OUT/(self.output_name+'.pdf')) if self.selected else io.BytesIO(),pagesize=(w,h));self.pdf.setTitle(name);self.pdf.setAuthor('Dengxiang Yang');self.svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w/72:.6f}in" height="{h/72:.6f}in" viewBox="0 0 {w} {h}">','<rect width="100%" height="100%" fill="white"/>']
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
  self.pdf.showPage();self.pdf.save()
  if not self.selected:return
  self.svg.append('</svg>');(OUT/(self.output_name+'.svg')).write_text('\n'.join(self.svg))
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

BULK_ONLY = _args.bulk_only
if not BULK_ONLY:
 # Supplementary Figure S4
 qc=read('FIG_QC_CELL_COUNTS_RAW_QC_DOUBLET.tsv');qc['sort']=qc.tissue.map({'SGN':0,'Lateral_wall':1});qc=qc.sort_values(['sort','GSM']).drop(columns='sort');save_table(qc,'FigureS4A_cell_counts.tsv')
 de=read('S1_DS_SUMMARY.tsv');save_table(de.drop(columns=['result_file'], errors='ignore'),'FigureS4D_differential_summary.tsv')
 f=Fig('FigureS4_cell_type_landscape',h=538)
 f.panel('A','Cell recovery across biological samples',8,15)
 x0=53;x1=493;y0=35;y1=122
 for v in [0,2,4,6]:
  y=y1-v/7*(y1-y0);f.line(x0,y,x1,y,'#E8EBED',.4)
 y_axis(f,x0,y0,y1,[0,2,4,6],0,7,'Cells (thousands)')
 for i,r in enumerate(qc.itertuples()):
  x=x0+12+i*27; vals=[r.n_cells_raw,r.n_cells_qc_pass,r.n_cells_final]
  for j,(v,c) in enumerate(zip(vals,[LIGHT,'#8AA5B6',BLUE])): f.rect(x-9+j*5,y1-v/1000/7*(y1-y0),4,v/1000/7*(y1-y0),c)
  f.text(x,y1+10,('C' if r.condition=='Naive' else 'N')+str(i%4+1),7,align='center')
 f.text(160,147,'SGN compartment',8,align='center');f.text(378,147,'Lateral wall',8,align='center');f.line(268,31,268,124,'#B7BEC3',.6)
 legend(f,[('Raw',LIGHT),('After QC','#8AA5B6'),('Final singlets',BLUE)],285,17,dx=66,size=7)
 for lab,tissue,xstart in [('B','SGN',8),('C','Lateral_wall',267)]:
  f.panel(lab,'SGN composition' if tissue=='SGN' else 'Lateral-wall composition',xstart,169)
  comp=read('FIG_'+tissue+'_BROAD_CELLTYPE_SAMPLE_COUNTS.tsv'); save_table(comp,'FigureS4'+lab+'_composition_counts.tsv');cs=comp.set_index('broad_annotation');xx=xstart+42; yy=192; hh=91;ww=187
  y_axis(f,xx,yy,yy+hh,[0,50,100],0,100,'Cells (%)')
  for i,gsm in enumerate(cs.columns):
   top=yy+hh;tot=cs[gsm].sum();x=xx+5+i*22.8
   for ct in types:
    if ct not in cs.index:continue
    v=cs.loc[ct,gsm]/tot*hh;top-=v;f.rect(x,top,16,v,col[ct])
   f.text(x+8,yy+hh+11,('C' if i<4 else 'N')+str(i%4+1),7,align='center')
  f.line(xx,yy+hh,xx+ww,yy+hh,BLACK,.5)
 legend(f,[(short[t],col[t]) for t in types],18,311,cols=3,dx=165,dy=12,size=7)
 f.panel('D','Cell-type differential expression after noise exposure',8,357)
 x0=191;x1=491;yy=377;step=9.8;maxval=2400;mid=(x0+x1)/2
 for v in [-2000,-1000,0,1000,2000]:
  x=mid+v/maxval*(x1-x0)/2;f.line(x,yy-6,x,yy+12*step+5,'#E8EBED' if v else '#9CA4AA',.5)
 for i,r in enumerate(de.itertuples()):
  y=yy+i*step; name=('SGN' if r.tissue=='SGN' else 'LW')+' | '+short[r.cell_type];f.text(x0-5,y+2,name,7,align='right')
  up=r.fdr_up_noise;down=r.fdr_down_noise;scale=(x1-x0)/2/maxval
  f.rect(mid-down*scale,y-3,down*scale,6,BLUE);f.rect(mid,y-3,up*scale,6,RED)
  f.text(mid-down*scale-3,y+2,str(down),6.5,align='right');f.text(mid+up*scale+3,y+2,str(up),6.5)
 h_axis(f,x0,x1,506,[-2000,-1000,0,1000,2000],-2400,2400)
 f.text(240,525,'Lower in noise',7,color=BLUE,align='center');f.text(438,525,'Higher in noise',7,color=RED,align='center')
 f.finish()

 # Figure 2
 focus=read('FIG_FOCUSED_PROGRAM_NES_HEATMAP.tsv');save_table(focus,'Figure2B_cell_type_pathway_results.tsv');null=read('GSE252285_MATCHED_NULL_PRIMARY.tsv');save_table(null,'FigureS3_localization_results.tsv')
 f=Fig('FigureS3_localization_and_cell_state',h=516);f.panel('A','Lateral-wall pathway responses in GSE168041',8,16)
 paths=[('HALLMARK_TNFA_SIGNALING_VIA_NFKB','TNF signaling via NF-kB','Hallmark'),('HALLMARK_IL6_JAK_STAT3_SIGNALING','IL6–JAK–STAT3 signaling','Hallmark'),('HALLMARK_ANGIOGENESIS','Angiogenesis','Hallmark'),('REACTOME_SIGNALING_BY_VEGF','Signaling by VEGF','Reactome'),('REACTOME_MAPK_FAMILY_SIGNALING_CASCADES','MAPK family signaling cascades','Reactome'),('REACTOME_ACTIVATION_OF_THE_AP_1_FAMILY_OF_TRANSCRIPTION_FACTORS','Activation of AP-1 family factors','Reactome')]
 cts=['Endothelial','Fibrocyte / mesenchymal','Immune / myeloid'];heads=['Endothelial','Fibrocyte / mesenchymal','Immune / myeloid']
 for j,h in enumerate(heads):f.text(258+j*94,40,h,7,align='center')
 for i,(pid,label,coll) in enumerate(paths):
  y=52+i*32;f.text(197,y+13,label,7.7,align='right');f.text(197,y+23,coll,6.8,color=GRAY,align='right')
  for j,ct in enumerate(cts):
   x=213+j*94;row=focus[(focus.pathway_id.str.endswith('::'+pid))&(focus.cell_type==ct)]
   if row.empty:f.rect(x,y,90,29,'#EFEFEF');f.text(x+45,y+18,'—',8,color=GRAY,align='center');continue
   r=row.iloc[0];f.rect(x,y,90,29,div(r.NES));tc='#FFFFFF' if abs(r.NES)>2.4 else BLACK
   f.text(x+45,y+12,f'{r.NES:+.2f}'+('*' if r.FDR<.05 else ''),8.5,color=tc,bold=True,align='center');f.text(x+45,y+23,'q = '+fmtq(r.FDR),7,color=tc,align='center')
 f.text(214,263,'NES',7.5,bold=True)
 for i in range(100):f.rect(240+i*1.4,254,1.5,8,div(-3.3+6.6*i/99))
 f.text(240,274,'−3.3',7,align='center');f.text(310,274,'0',7,align='center');f.text(380,274,'+3.3',7,align='center');f.text(493,267,'* q < 0.05',7,align='right')
 f.panel('B','Pathway localization and condition response',8,298)
 f.text(231,320,'GSE252285 localization',8,bold=True,align='center');f.text(415,320,'GSE168041 noise response',8,bold=True,align='center')
 pairs=[('TNFA_NFKB','Endothelial_vascular','TNF–NF-kB | Endothelial','HALLMARK_TNFA_SIGNALING_VIA_NFKB','Endothelial'),('TNFA_NFKB','Fibrocyte_mesenchymal','TNF–NF-kB | Fibrocyte','HALLMARK_TNFA_SIGNALING_VIA_NFKB','Fibrocyte / mesenchymal'),('VEGF_ANGIOGENESIS','Endothelial_vascular','Angiogenesis | Endothelial','HALLMARK_ANGIOGENESIS','Endothelial'),('IL6_JAK_STAT3','Endothelial_vascular','IL6–JAK–STAT3 | Endothelial','HALLMARK_IL6_JAK_STAT3_SIGNALING','Endothelial'),('IL6_JAK_STAT3','Immune_myeloid','IL6–JAK–STAT3 | Myeloid','HALLMARK_IL6_JAK_STAT3_SIGNALING','Immune / myeloid'),('MAPK_AP1','Endothelial_vascular','MAPK/AP-1 | Endothelial','REACTOME_MAPK_FAMILY_SIGNALING_CASCADES','Endothelial')]
 # Source program is matched flexibly only for the exact sixth row identifier.
 rows=[]
 for i,(prog,ct,label,pid,statect) in enumerate(pairs):
  nr=null[(null.program==prog)&(null.target_cell_type==ct)]
  if nr.empty and i==5:nr=null[null.program.str.contains('MAPK')]
  nr=nr.iloc[0]; sr=focus[(focus.pathway_id.str.endswith('::'+pid))&(focus.cell_type==statect)].iloc[0];y=342+i*21
  f.text(163,y+3,label,7,align='right');f.line(177,y,296,y,'#E5E8EB',.4);f.line(341,y,492,y,'#E5E8EB',.4)
  xx=177+nr.null_z/10*119;f.circle(xx,y,4,BLUE if nr.empirical_FDR<.05 else '#FFFFFF',BLUE,.8)
  sx=341+(sr.NES+3.5)/7*151;f.circle(sx,y,4,(RED if sr.NES>0 else BLUE) if sr.FDR<.05 else '#FFFFFF',RED if sr.NES>0 else BLUE,.8)
  rows.append({'program':prog,'target_cell_type':ct,'localization_null_z':nr.null_z,'localization_empirical_FDR':nr.empirical_FDR,'state_pathway':pid,'state_NES':sr.NES,'state_FDR':sr.FDR})
 f.line(416.5,334,416.5,453,'#9CA4AA',.6)
 h_axis(f,177,296,463,[0,5,10],0,10,'Matched-null z score');h_axis(f,341,492,463,[-3,0,3],-3.5,3.5,'Normalized enrichment score')
 f.circle(178,501,3.2,BLACK);f.text(185,504,'q < 0.05',7);f.circle(250,501,3.2,'#FFFFFF',BLACK);f.text(257,504,'q ≥ 0.05',7)
 save_table(pd.DataFrame(rows),'FigureS3_paired_display.tsv');f.finish()

 # Main Figure 2: full cell-type NF-kB landscape and focused lateral-wall programs
 fullfocus=read('S1_GSEA_FOCUSED_PROGRAMS.tsv')
 nf=fullfocus[fullfocus.pathway_id.str.endswith('::HALLMARK_TNFA_SIGNALING_VIA_NFKB')].copy()
 save_table(nf,'Figure2A_all_cell_type_NFkB.tsv')
 f=Fig('Figure2_cell_type_pathway_responses',h=565)
 f.panel('A','TNF signaling via NF-kB across cochlear cell types',8,16)
 x0=195;x1=427;lo=-3.5;hi=2.7;zero=x0+(0-lo)/(hi-lo)*(x1-x0)
 for v in [-3,-2,-1,0,1,2]:
  xx=x0+(v-lo)/(hi-lo)*(x1-x0);f.line(xx,34,xx,218,'#E7EAED' if v else '#A3ABB1',.5)
 f.text(477,32,'FDR',7.5,bold=True,align='center')
 for i,r in enumerate(nf.itertuples()):
  yy=47+i*13;label=('SGN' if r.tissue=='SGN' else 'LW')+' | '+short[r.cell_type];f.text(186,yy+3,label,7.5,align='right')
  xx=x0+(r.NES-lo)/(hi-lo)*(x1-x0);cc=RED if r.NES>0 else BLUE
  f.line(zero,yy,xx,yy,cc,1.3);f.circle(xx,yy,3.5,cc if r.FDR<.05 else '#FFFFFF',cc,.8)
  qstr=f'{r.FDR:.5f}' if .049<r.FDR<.051 else fmtq(r.FDR)
  f.text(477,yy+3,qstr,7,align='center')
 h_axis(f,x0,x1,225,[-3,-2,-1,0,1,2],lo,hi,'Normalized enrichment score')
 f.circle(30,256,3.2,BLACK);f.text(37,259,'q < 0.05',7);f.circle(102,256,3.2,'#FFFFFF',BLACK);f.text(109,259,'q ≥ 0.05',7)
 f.panel('B','Focused pathway responses in lateral-wall cell types',8,288)
 for j,h in enumerate(heads):f.text(258+j*94,313,h,7,align='center')
 for i,(pid,label,coll) in enumerate(paths):
  y=325+i*32;f.text(197,y+13,label,7.7,align='right');f.text(197,y+23,coll,6.8,color=GRAY,align='right')
  for j,ct in enumerate(cts):
   x=213+j*94;row=focus[(focus.pathway_id.str.endswith('::'+pid))&(focus.cell_type==ct)]
   if row.empty:f.rect(x,y,90,29,'#EFEFEF');f.text(x+45,y+18,'—',8,color=GRAY,align='center');continue
   r=row.iloc[0];f.rect(x,y,90,29,div(r.NES));tc='#FFFFFF' if abs(r.NES)>2.4 else BLACK
   f.text(x+45,y+12,f'{r.NES:+.2f}'+('*' if r.FDR<.05 else ''),8.5,color=tc,bold=True,align='center');f.text(x+45,y+23,'q = '+fmtq(r.FDR),7,color=tc,align='center')
 f.text(214,536,'NES',7.5,bold=True)
 for i in range(100):f.rect(240+i*1.4,527,1.5,8,div(-3.3+6.6*i/99))
 f.text(240,549,'−3.3',7,align='center');f.text(310,549,'0',7,align='center');f.text(380,549,'+3.3',7,align='center');f.text(493,540,'* q < 0.05',7,align='right')
 f.finish()

# Figure 4
# Current targeted results are the explicitly identified 2026-11-01 re-estimation.
gsea=pd.read_csv(BULK_ROOT/'enrichment/results/TARGETED_27_RESULTS.tsv',sep='\t').rename(columns={'contrast':'comparison'})
_bulk_axes = {
 'GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY':'inflammation_nfkb',
 'GOBP_PROTEIN_FOLDING':'proteostasis_er_stress',
 'GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS':'proteostasis_er_stress',
 'GOBP_RESPONSE_TO_OXIDATIVE_STRESS':'oxidative_stress',
 'HALLMARK_APOPTOSIS':'cell_death',
 'HALLMARK_INFLAMMATORY_RESPONSE':'inflammation_nfkb',
 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY':'oxidative_stress',
 'HALLMARK_TNFA_SIGNALING_VIA_NFKB':'inflammation_nfkb',
 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE':'proteostasis_er_stress'}
gsea['axis']=gsea.pathway.map(_bulk_axes)
assert not gsea.axis.isna().any()
assert len(gsea)==27 and not gsea[['NES','ES','P','global_27_test_FDR']].isna().any().any()
save_table(gsea[['comparison','pathway','axis','NES','ES','P','within_contrast_FDR','global_27_test_FDR']],'Figure4A_bulk_pathway_results.tsv')
summ=pd.read_csv(BULK_ROOT/'deseq2/results/NEW_DESEQ2_SUMMARY.tsv',sep='\t');save_table(summ,'Figure4B_bulk_DE_summary.tsv')
f=Fig('Figure4_noise_exposure_comparisons',h=474);f.panel('A','Whole-cochlea pathway responses in GSE160639',8,16)
pp=[('HALLMARK_TNFA_SIGNALING_VIA_NFKB','TNF signaling via NF-kB'),('HALLMARK_INFLAMMATORY_RESPONSE','Inflammatory response'),('GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY','Cytokine-mediated signaling'),('HALLMARK_UNFOLDED_PROTEIN_RESPONSE','Unfolded protein response'),('GOBP_PROTEIN_FOLDING','Protein folding'),('GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS','Response to ER stress'),('HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','Reactive oxygen species pathway'),('GOBP_RESPONSE_TO_OXIDATIVE_STRESS','Response to oxidative stress'),('HALLMARK_APOPTOSIS','Apoptosis')]
comps=['94dB_vs_70dB','105dB_vs_70dB','105dB_vs_94dB']
for j,c in enumerate(comps):f.text(262+j*94,39,c.replace('_vs_',' vs ').replace('dB',' dB'),8,align='center')
for i,(pid,label) in enumerate(pp):
 y=50+i*27;f.text(207,y+17,label,8,align='right')
 for j,c in enumerate(comps):
  r=gsea[(gsea.pathway==pid)&(gsea.comparison==c)].iloc[0];x=217+j*94;f.rect(x,y,89,25,div(r.NES));tc='#FFFFFF' if abs(r.NES)>2.4 else BLACK
  f.text(x+44.5,y+11,f'{r.NES:+.2f}'+('*' if r.global_27_test_FDR<.05 else ''),8,color=tc,bold=True,align='center');f.text(x+44.5,y+21,'q = '+fmtq(r.global_27_test_FDR),6.8,color=tc,align='center')
f.text(217,309,'NES',7.5,bold=True)
for i in range(100):f.rect(242+i*1.4,299,1.5,8,div(-3.3+6.6*i/99))
f.text(242,319,'−3.3',7,align='center');f.text(312,319,'0',7,align='center');f.text(382,319,'+3.3',7,align='center');f.text(495,309,'* Global q < 0.05',7,align='right')
f.panel('B','Differentially expressed genes across exposure comparisons',8,342)
x0=160;x1=491;mid=(x0+x1)/2;mx=3500
for i,c in enumerate(comps):
 r=summ[summ.comparison==c].iloc[0];y=366+i*24;up=r.positive;down=r.negative
 f.text(151,y+3,c.replace('_vs_',' vs ').replace('dB',' dB'),8,align='right');scale=(x1-x0)/2/mx;f.rect(mid-down*scale,y-5,down*scale,10,BLUE);f.rect(mid,y-5,up*scale,10,RED);f.text(mid-down*scale-4,y+3,int(down),7,align='right');f.text(mid+up*scale+4,y+3,int(up),7)
f.line(mid,355,mid,424,BLACK,.6);h_axis(f,x0,x1,430,[-3000,-1500,0,1500,3000],-mx,mx)
f.text(235,459,'Lower in numerator',7,color=BLUE,align='center');f.text(419,459,'Higher in numerator',7,color=RED,align='center');f.finish()

if not BULK_ONLY:
 # Figure 5 and supplement: gene detection/expression
 expr=read('TOP20_BY_BROAD_CELLTYPE_CONDITION.tsv');save_table(expr,'Figure5_and_S2_expression.tsv');st=read('TOP20_LOCALIZATION_STABILITY.tsv');save_table(st,'Figure5_gene_localization.tsv')
 ct4=['Immune_myeloid','Fibrocyte_mesenchymal','Endothelial_vascular','Glial','Supporting_epithelial','Hair_sensory','UNRESOLVED']; ct4lab=['Myeloid','Fibrocyte','Endothelial','Glial','Epithelial','Hair','Unassigned']; ct4colors=[col['Immune / myeloid'],col['Fibrocyte / mesenchymal'],col['Endothelial'],col['Glial / Schwann-like'],col['Epithelial / supporting-associated'],'#665566',col['UNRESOLVED']]
 conditions=['Control','HHL','NIHL'];stable=st[st.localization_stability=='STABLE_SUPPORTED'].gene.tolist(); allgenes=st.gene.tolist()

 def dot_fig(name,genes,full=False):
  h=674 if full else 520;f=Fig(name,h=h);f.panel('A','Cell-type expression across GSE252285 conditions',8,16)
  step=18 if full else 19;top=57
  for j,c in enumerate(conditions):
   start=93+j*137;f.text(start+51,37,c,9,bold=True,align='center')
   for k,ct in enumerate(ct4):
    xx=start+k*17;f.line(xx-7,top-8,xx-7,top+(len(genes)-1)*step+8,'#F0F2F4',.4)
    f.text(xx,top+(len(genes)-1)*step+20,ct4lab[k],7,rotate=60)
   for i,g in enumerate(genes):
    yy=top+i*step
    if j==0:f.text(78,yy+3,g,8,italic=True,align='right')
    f.line(start-8,yy,start+110,yy,'#E9ECEE',.4)
    for k,ct in enumerate(ct4):
     r=expr[(expr.condition==c)&(expr.broad_cell_type==ct)&(expr.gene==g)].iloc[0];rr=math.sqrt(r.pct_expressing/100)*6.5;xx=start+k*17
     if rr>0:f.circle(xx,yy,rr,seq(r.mean_expression),BLACK if r.n_total_cells<20 else None,.6)
  bottom=top+(len(genes)-1)*step+82
  f.text(88,bottom,'Detection (%)',7.5,bold=True)
  for k,pct in enumerate([25,50,75,100]):
   x=176+k*30;f.circle(x,bottom-3,math.sqrt(pct/100)*6.5,'#788C9B');f.text(x,bottom+14,pct,7,align='center')
  f.text(320,bottom-9,'Mean expression',7.5,bold=True)
  for i in range(80):f.rect(320+i*1.8,bottom-3,1.9,8,seq(3.1*i/79))
  for v in [0,1,2,3]:f.text(320+v/3.1*144,bottom+17,str(v),7,align='center')
  ytitle=bottom+46;f.panel('B','Cell composition by condition',8,ytitle)
  comp=expr[['broad_cell_type','condition','n_total_cells']].drop_duplicates();save_table(comp,'Figure5B_composition.tsv')
  for j,c in enumerate(conditions):
   y=ytitle+24+j*18;f.text(80,y+4,c,8,align='right');x=94;sub=comp[comp.condition==c].set_index('broad_cell_type');total=sub.n_total_cells.sum()
   for ct,cc in zip(ct4,ct4colors):
    ww=sub.loc[ct,'n_total_cells']/total*353;f.rect(x,y-5,ww,10,cc);x+=ww
   f.text(456,y+3,f'n = {total:,}',7)
   if j==2:h_axis(f,94,447,y+15,[0,25,50,75,100],0,100)
  legend(f,list(zip(ct4lab,ct4colors)),24,h-30,cols=4,dx=121,dy=13,size=7)
  f.finish()
 dot_fig('Figure5_condition_expression',stable);dot_fig('FigureS2_all_twenty_genes',allgenes,True)

# Supplement Figure S1: all tested genes with a finite adjusted P value
f=Fig('FigureS1_bulk_volcanoes',h=225); excluded={}
for j,c in enumerate(comps):
 d=pd.read_csv(BULK_ROOT/'deseq2/results'/('NEW_DESeq2_'+c+'_all_genes.csv.gz'));a=d[np.isfinite(d.log2FoldChange)&np.isfinite(d.padj)&(d.padj>0)].copy();a['minus_log10_FDR']=-np.log10(a.padj);save_table(d,'FigureS1_'+c+'.tsv');excluded[c]=len(d)-len(a)
 xleft=41+j*170;xright=xleft+121;y0=40;y1=171;lim=8; ymax=130
 f.panel(chr(65+j),c.replace('_vs_',' vs ').replace('dB',' dB'),6+j*170,17)
 for v in [0,40,80,120]:
  y=y1-v/ymax*(y1-y0);f.line(xleft,y,xright,y,'#EBEDEF',.4)
 y_axis(f,xleft,y0,y1,[0,40,80,120],0,ymax,'−log10 FDR')
 for is_sig in [False,True]:
  sub=a[(a.padj<.05)==is_sig]
  for r in sub.itertuples():
   xx=xleft+(np.clip(r.log2FoldChange,-lim,lim)+lim)/(2*lim)*(xright-xleft);yy=y1-min(r.minus_log10_FDR,ymax)/ymax*(y1-y0)
   cc=(RED if r.log2FoldChange>0 else BLUE) if is_sig else '#CBD0D5';f.circle(xx,yy,.6,cc)
 h_axis(f,xleft,xright,y1,[-8,-4,0,4,8],-8,8,'log2 fold change')
 ty=y1-(-math.log10(.05))/ymax*(y1-y0);f.line(xleft,ty,xright,ty,GRAY,.5,dash=2)
f.finish()

# Preserve reviewed captions in figure_legends.txt and figure_legends.json.
_parameter_file = OUT/'figure_generation_parameters.json'
_parameters = json.loads(_parameter_file.read_text()) if _parameter_file.exists() else {}
_parameters.update({'width_points':510,'raster_dpi':600,'font':FONT_FAMILY,
    'render_engine':'ReportLab and Poppler','python':__import__('sys').version,
    'reportlab':reportlab.Version,'numpy':np.__version__,'pandas':pd.__version__,
    'inference':'This rendering script reads the current archived statistical outputs and does not fit models or run enrichment.',
    'bulk_analysis':'2026-11-01 public-count DESeq2 reconstruction and fixed nine-pathway enrichment',
    'bulk_volcano_display_omissions':excluded})
_parameter_file.write_text(json.dumps(_parameters,indent=2))
