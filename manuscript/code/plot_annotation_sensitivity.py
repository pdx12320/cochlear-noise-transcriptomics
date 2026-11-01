"""Supplementary Figure S6: separate c19 annotation sensitivity, no label overwrite."""
from pathlib import Path
import argparse
import math
import pandas as pd
from figure_support import Fig, ROOT, DAT, BLUE, RED, BLACK, GRAY, h_axis, legend, save_table

parser=argparse.ArgumentParser()
parser.add_argument('--source-dir',type=Path,default=ROOT/'source_inputs/annotation_sensitivity')
args=parser.parse_args();src=args.source_dir
assign=pd.read_csv(src/'ALL20_TOP_ASSIGNMENTS_ORIGINAL_VS_C19.tsv',sep='\t')
stable=pd.read_csv(src/'ALL20_MEAN_LOCALIZATION_STABILITY.tsv',sep='\t')
null=pd.read_csv(src/'C19_VS_ORIGINAL_MATCHED_NULL.tsv',sep='\t')
assert len(assign)==160 and len(stable)==20 and len(null)==12
assert set(null.match_mode)=={'primary','length_sensitivity'}
assert (null.iterations_original==10000).all() and (null.iterations_c19_revised==10000).all()
assert stable.original_stable_mean_top.sum()==11
assert stable.c19_fibrocyte_stable_mean_top.sum()==12
assert stable.original_main11_retained.all()
assert stable.loc[stable.newly_stable_in_sensitivity,'gene'].tolist()==['Cebpb']
for mode,group in null.groupby('match_mode'):
 assert len(group)==6 and group.test_id.nunique()==6
save_table(assign,'FigureS6A_all20_top_assignments.tsv')
save_table(stable,'FigureS6A_mean_localization_stability.tsv')
save_table(null,'FigureS6B_matched_null_before_after.tsv')

# Match the broad-lineage colors used in the main figures.
ct={
 'Endothelial_vascular':('#4477AA','Endothelial','Endothelial'),
 'Fibrocyte_mesenchymal':('#CCBB44','Fibrocyte','Fibrocyte / mesenchymal'),
 'Glial':('#AA4499','Glial','Glial'),
 'Immune_myeloid':('#EE8866','Myeloid','Immune / myeloid'),
 'Supporting_epithelial':('#66A699','Supporting','Supporting epithelial'),
 'UNRESOLVED':('#B8BEC5','Unassigned','Unassigned'),
}
conds=['Control','HHL','NIHL'];gene_order=stable.gene.tolist()
f=Fig('FigureS6_annotation_sensitivity',w=600,h=734)
f.panel('A','Effect of assigning cluster 19 to fibrocytes',12,20)
legend(f,[(v[2],v[0]) for v in ct.values()],30,43,cols=3,dx=185,dy=14,size=7.8)
lefts=[114,266,418];cw=48;gap=5;row_y=102;rh=12.4
for c,x in zip(conds,lefts):
 f.text(x+cw+gap/2,70,c,9,bold=True,align='center')
 for offset,line in enumerate(['Original','annotation']):
  f.text(x+cw/2,84+10*offset,line,7.2,align='center')
 for offset,line in enumerate(['Cluster 19','as fibrocytes']):
  f.text(x+cw+gap+cw/2,84+10*offset,line,7.2,align='center')
for i,g in enumerate(gene_order):
 y=row_y+i*rh;s=stable.loc[stable.gene==g].iloc[0]
 if s.original_stable_mean_top:f.text(24,y+8.7,'*',10,bold=True)
 elif s.newly_stable_in_sensitivity:f.text(24,y+8.7,'+',10,color=RED,bold=True)
 f.text(39,y+8.7,g,8.5,italic=True)
 for c,x in zip(conds,lefts):
  z=assign.loc[(assign.gene==g)&(assign.condition==c)&(assign.metric=='mean_log_normalized')].iloc[0]
  for j,k in enumerate(['original_top_type','c19_fibrocyte_top_type']):
   category=z[k];color,abbr,_=ct[category];xx=x+j*(cw+gap)
   changed=bool(z.top_type_changed) and j==1
   f.rect(xx,y,cw,rh-1.7,color,stroke=RED if changed else None,lw=1.8)
   f.text(xx+cw/2,y+8.4,abbr,7.3,align='center',color='#FFFFFF' if category in ['Endothelial_vascular','Glial'] else BLACK)
# The outlined tile is the changed NIHL Cebpb assignment; no mean-ranking change is hidden.
f.text(28,365,'* Primary recurrent genes',7.8)
f.text(199,365,'+ Additional recurrent gene',7.8)
f.rect(388,358,8,8,'#FFFFFF',stroke=RED,lw=1.5)
f.text(402,365,'Changed category',7.8)

f.panel('B','Matched-null pathway localization',12,390)
f.circle(30,405,2.8,BLUE);f.text(37,408,'Original annotation',8)
f.circle(148,405,2.8,RED);f.text(155,408,'Cluster 19 as fibrocytes',8)
f.text(292,408,'Filled: FDR < 0.05 and positive delta; open: otherwise',7.2)
labels=[
 ('MAPK_AP1__Endothelial_vascular','MAPK / AP-1 - endothelial'),
 ('VEGF_ANGIOGENESIS__Endothelial_vascular','VEGF / angiogenesis - endothelial'),
 ('TNFA_NFKB__Endothelial_vascular','TNF / NF-kB - endothelial'),
 ('TNFA_NFKB__Fibrocyte_mesenchymal','TNF / NF-kB - fibrocyte'),
 ('IL6_JAK_STAT3__Endothelial_vascular','IL6 / JAK-STAT3 - endothelial'),
 ('IL6_JAK_STAT3__Immune_myeloid','IL6 / JAK-STAT3 - myeloid'),
]
x0,x1=200,386;lo,hi=0,.06
assert null[['localization_delta_original','localization_delta_c19_revised']].min().min()>=lo
assert null[['localization_delta_original','localization_delta_c19_revised']].max().max()<=hi
xx=lambda value:x0+(value-lo)/(hi-lo)*(x1-x0)
def qfmt(q):return f'{q:.1e}' if q<.001 else f'{q:.4f}'
for mode,heading,hy,start in [('primary','Expression + detection matching',429,448),('length_sensitivity','Additional gene-length matching',562,581)]:
 f.text(28,hy,heading,9,bold=True)
 f.text(449,hy-4,'Original annotation',7.1,bold=True,align='center',color=BLUE)
 f.text(543,hy-4,'Cluster 19 as fibrocytes',7.1,bold=True,align='center',color=RED)
 f.text(449,hy+5,'FDR',7.1,bold=True,align='center',color=BLUE)
 f.text(543,hy+5,'FDR',7.1,bold=True,align='center',color=RED)
 for i,(test,label) in enumerate(labels):
  y=start+i*18;z=null.loc[(null.match_mode==mode)&(null.test_id==test)].iloc[0]
  if i%2==0:f.rect(23,y-8,554,17,'#F4F6F7')
  f.text(28,y+3,label,7.9)
  f.line(x0,y,x1,y,'#E2E6E9',.45)
  old=float(z.localization_delta_original);new=float(z.localization_delta_c19_revised)
  f.line(xx(old),y-2.5,xx(new),y+2.5,'#747D83',.75)
  for variant,offset,color in [('original',-2.5,BLUE),('c19_revised',2.5,RED)]:
   q=float(z['empirical_FDR_'+variant]);delta=float(z['localization_delta_'+variant]);sig=q<.05 and delta>0
   f.circle(xx(delta),y+offset,2.7,color if sig else '#FFFFFF',stroke=color,lw=1.1)
  f.text(449,y+3,qfmt(float(z.empirical_FDR_original)),8,align='center',color=BLUE)
  f.text(543,y+3,qfmt(float(z.empirical_FDR_c19_revised)),8,align='center',color=RED)
# Shared effect scale; delta is a score difference, not a fold change or a confidence interval.
h_axis(f,x0,x1,690,[0,.02,.04,.06],lo,hi,'Observed UCell localization delta',label_y=718)
f.finish()
