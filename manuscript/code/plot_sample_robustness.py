#!/usr/bin/env python3
"""Render animal-level expression and deletion sensitivity from archived result tables."""
from pathlib import Path
import argparse,json,math,hashlib
import numpy as np,pandas as pd
from figure_support import Fig,ROOT,OUT,DAT,BLUE,RED,GRAY,BLACK,div,h_axis,y_axis,rgb
INP=ROOT/'revision_inputs/sample_robustness'
TNF='Hallmark::HALLMARK_TNFA_SIGNALING_VIA_NFKB';MAPK='Reactome::REACTOME_MAPK_FAMILY_SIGNALING_CASCADES'
ORDER=[('Lateral_wall','Fibrocyte / mesenchymal',TNF),('SGN','Fibrocyte / mesenchymal',TNF),('Lateral_wall','Endothelial',TNF),('SGN','Endothelial',TNF),('Lateral_wall','Immune / myeloid',TNF),('SGN','Immune / myeloid',TNF),('Lateral_wall','Endothelial',MAPK)]
CELL={'Fibrocyte / mesenchymal':'fibrocyte/mesenchymal','Endothelial':'endothelial','Immune / myeloid':'immune/myeloid'}
def label(t,c):return ('LW' if t=='Lateral_wall' else 'SGN')+' '+CELL[c]
def diamond(f,x,y,r=3.2,color=BLACK):
 pts=[(x,y-r),(x+r,y),(x,y+r),(x-r,y)];p=f.pdf;p.saveState();p.setFillColorRGB(*rgb(color));p.setStrokeColorRGB(*rgb(color));path=p.beginPath();path.moveTo(x,f.h-y+r)
 for xx,yy in pts[1:]:path.lineTo(xx,f.h-yy)
 path.close();p.drawPath(path,fill=1,stroke=0);p.restoreState();f.svg.append('<polygon points="'+' '.join(f'{xx},{yy}' for xx,yy in pts)+f'" fill="{color}"/>')
def sample_code(s):return ('C' if '_naive_' in s else 'N')+s.rsplit('_',1)[1]
def profiles():
 d=pd.read_csv(INP/'PER_ANIMAL_REPRESENTATIVE_GENES.tsv',sep='\t')
 panels=[('Fibrocyte / mesenchymal',g) for g in ['Ccl5','Ifih1','Mcl1','Ifit2']]+[('Endothelial',g) for g in ['Map2k2','Mapk3','Dusp2']]
 f=Fig('FigureS7_sample_gene_expression',560,650);allrows=[]
 for i,(ct,gene) in enumerate(panels):
  col=i%2;row=i//2;x=20+col*277;y=18+row*156
  f.text(x,y,chr(65+i).lower(),11,bold=True);f.text(x+17,y,gene,10,italic=True)
  f.text(x+17,y+13,'Lateral-wall '+CELL[ct],7.7,color=GRAY)
  z=d[(d.tissue=='Lateral_wall')&(d.cell_type==ct)&(d.gene_symbol==gene)].copy();assert len(z)==8 and z.GSM.nunique()==8 and set(z.status)=={'AVAILABLE'}
  z['sample_code']=z.sample_id.map(sample_code);z['panel']=chr(65+i);allrows.append(z)
  x0=x+36;x1=x+235;y0=y+28;y1=y+125;hi=max(1,math.ceil(float(z.log1p_normalized_count.max())+.4));step=1 if hi<=5 else 2;ticks=np.arange(0,hi+.1,step)
  y_axis(f,x0,y0,y1,ticks,0,hi,'log1p normalized count')
  f.line(x0,y1,x1,y1,BLACK,.5)
  for cond,xc,color in [('Naive',x0+48,BLUE),('Noise',x0+143,RED)]:
   group=z[z.condition==cond].sort_values('sample_id');assert len(group)==4
   for (rr,dx) in zip(group.itertuples(),[-19,-6.5,6.5,19]):
    yy=y1-rr.log1p_normalized_count/hi*(y1-y0);f.circle(xc+dx,yy,3,color,stroke='#FFFFFF',lw=.4);f.text(xc+dx,yy-6,rr.sample_code,6.5,color=color,align='center')
   f.text(xc,y1+14,cond+' (n = 4)',7.6,align='center')
 f.finish();shown=pd.concat(allrows,ignore_index=True);shown.to_csv(DAT/'FigureS7_per_animal_gene_expression.tsv',sep='\t',index=False)
 return {'figure':f.output_name,'panels':7,'animals_per_panel':8,'points':56,'genes':[p[1] for p in panels],'SEM_or_cell_level_tests':False}
def mainfig():
 score=pd.read_csv(INP/'PER_ANIMAL_PATHWAY_SCORES.tsv',sep='\t');loo=pd.read_csv(INP/'LOO_FOCUSED_ALL_56.tsv',sep='\t');assert len(loo)==56
 f=Fig('Figure3_sample_robustness',560,568)
 f.panel('A','Pathway expression in individual animals',16,19)
 x0=220;cw=38;rh=22;y0=66
 f.text(x0+cw*2,42,'Unexposed (Naive)',8,color=BLUE,align='center');f.text(x0+cw*6,42,'Noise exposed',8,color=RED,align='center')
 f.line(x0,47,x0+4*cw-3,47,BLUE,1.5);f.line(x0+4*cw+3,47,x0+8*cw,47,RED,1.5)
 for j,code in enumerate(['C1','C2','C3','C4','N1','N2','N3','N4']):f.text(x0+(j+.5)*cw,59,code,8,align='center')
 rowsA=[]
 for i,(t,c,p) in enumerate(ORDER):
  y=y0+i*rh;z=score[(score.tissue==t)&(score.cell_type==c)&(score.pathway_id==p)].copy();z['sample_code']=z.sample_id.map(sample_code);z=z.set_index('sample_code').loc[['C1','C2','C3','C4','N1','N2','N3','N4']].reset_index();assert len(z)==8
  f.text(18,y+9,label(t,c),8.2,bold=True);f.text(18,y+19,'TNF-NF-kB' if p==TNF else 'MAPK-family signaling',7.2,color=GRAY)
  for j,r in enumerate(z.itertuples()):
   v=r.mean_gene_z_log1p_normalized;f.rect(x0+j*cw,y,cw-1,rh-1,div(v,1));f.text(x0+(j+.5)*cw-.5,y+14,(f'{v:+.2f}' if abs(v)>=.005 else '0.00'),7.8,color='#FFFFFF' if abs(v)>.6 else BLACK,align='center')
  z['figure_row']=i+1;rowsA.append(z)
 bx=272;by=237;bw=185
 for k in range(120):v=-1+2*k/119;f.rect(bx+k*bw/120,by,bw/120+.1,7,div(v,1))
 for v in [-1,0,1]:f.text(bx+(v+1)/2*bw,by+19,f'{v:g}',7,align='center')
 f.text(18,245,'Mean gene-standardized expression',7.7)
 f.panel('B','Leave-one-animal-out pathway estimates',16,278)
 lo=math.floor(min(loo.LOO_NES.min(),loo.baseline_NES.min())-.15);hi=math.ceil(max(loo.LOO_NES.max(),loo.baseline_NES.max())+.15)
 px0=207;px1=433;xmap=lambda v:px0+(v-lo)/(hi-lo)*(px1-px0)
 f.text(471,296,'Same direction',7.4,align='center');f.text(531,296,'q < 0.05',7.4,align='center')
 f.line(447,302,447,496,'#D5DADF',.5)
 for tick in np.arange(lo,hi+1):f.line(xmap(tick),304,xmap(tick),496,'#E5E8EA',.35,dash=1 if tick==0 else None)
 rowsB=[]
 for i,(t,c,p) in enumerate(ORDER):
  y=318+i*27;z=loo[(loo.tissue==t)&(loo.cell_type==c)&(loo.pathway_id==p)].copy();assert len(z)==8 and set(z.test_status)=={'TESTED'}
  z['omitted_code']=z.omitted_sample_id.map(sample_code);z=z.set_index('omitted_code').loc[['C1','C2','C3','C4','N1','N2','N3','N4']].reset_index()
  f.text(18,y-1,label(t,c),8.2,bold=True);f.text(18,y+10,'TNF-NF-kB' if p==TNF else 'MAPK-family signaling',7.1,color=GRAY)
  mn=z.LOO_NES.min();mx=z.LOO_NES.max();f.line(xmap(mn),y,xmap(mx),y,'#8A949C',1.2);f.line(xmap(mn),y-3,xmap(mn),y+3,'#8A949C',.7);f.line(xmap(mx),y-3,xmap(mx),y+3,'#8A949C',.7)
  for j,r in enumerate(z.itertuples()):
   col=BLUE if r.omitted_condition=='Naive' else RED;yy=y+[-6,-2,2,6][j%4];f.circle(xmap(r.LOO_NES),yy,2.3,col if r.LOO_FDR<.05 else '#FFFFFF',stroke=col,lw=.85)
  diamond(f,xmap(z.baseline_NES.iloc[0]),y-10,3.1)
  f.text(471,y+3,f'{int(z.same_direction_as_full.sum())}/8',9,bold=True,align='center');f.text(531,y+3,f'{int(z.FDR_below_005.sum())}/8',9,bold=True,align='center')
  z['figure_row']=i+1;rowsB.append(z)
 h_axis(f,px0,px1,504,np.arange(lo,hi+1),lo,hi,'NES (Noise vs Naive)',label_y=531)
 diamond(f,22,554,3.1);f.text(30,557,'Full 8-animal estimate',7.2)
 f.circle(178,554,2.4,BLUE,stroke=BLUE);f.text(186,557,'Omit naive',7.2);f.circle(235,554,2.4,RED,stroke=RED);f.text(243,557,'Omit noise',7.2)
 f.text(306,557,'Filled: q < 0.05; open: q >= 0.05',7.2)
 f.finish();pd.concat(rowsA,ignore_index=True).to_csv(DAT/'Figure3A_per_animal_pathway_expression.tsv',sep='\t',index=False);pd.concat(rowsB,ignore_index=True).to_csv(DAT/'Figure3B_leave_one_animal_out.tsv',sep='\t',index=False)
 return {'figure':'Figure3_sample_robustness','expression_rows':7,'animal_expression_cells':56,'LOO_points':56,'full_estimate_diamonds':7,'range_not_CI':True,'NES_axis_limits':[lo,hi],'no_independent_replicate_claim':True}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--profiles-only',action='store_true');a=p.parse_args();records=[profiles()]
 if not a.profiles_only:records.append(mainfig())
 record={'status':'RENDERED_REQUIRES_VISUAL_REVIEW','figures':records,'input_SHA256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in INP.iterdir() if f.is_file() and f.name.endswith('.tsv')}}
 (ROOT/'analysis_records/SAMPLE_ROBUSTNESS_FIGURE_BUILD.json').write_text(json.dumps(record,indent=2))
