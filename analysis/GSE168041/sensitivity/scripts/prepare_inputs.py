#!/usr/bin/env python3
"""Freeze inputs and export original per-animal profiles; never refit a model."""
from pathlib import Path
import argparse,hashlib,json,re,shutil,sys
import numpy as np,pandas as pd,anndata as ad
from importlib.metadata import version
P=argparse.ArgumentParser();P.add_argument('--source',required=True);P.add_argument('--out',required=True);a=P.parse_args()
src=Path(a.source);out=Path(a.out);dest=out/'inputs'
if dest.exists():raise RuntimeError('Refusing existing frozen inputs')
dest.mkdir(parents=True)
assert version('pydeseq2')=='0.5.4' and sys.version_info[:3]==(3,12,13)
cts={'Fibrocyte_mesenchymal':'Fibrocyte / mesenchymal','Endothelial':'Endothelial','Immune_myeloid':'Immune / myeloid'}
termfile=src/'10_gsea/MSIGDB_2026.1.Mm_Hallmark_Reactome_TERM2GENE.tsv.gz'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(termfile)=='98dcdf8081f93b433396b8cbc1eaa6ad0e34ce2afc84ae82fbda99df2e641caf'
shutil.copy2(termfile,dest/termfile.name)
term=pd.read_csv(termfile,sep='\t')
TNF='Hallmark::HALLMARK_TNFA_SIGNALING_VIA_NFKB';MAPK='Reactome::REACTOME_MAPK_FAMILY_SIGNALING_CASCADES'
repgenes=['Ccl5','Ifih1','Mcl1','Ifit2','Map2k2','Mapk3','Dusp2']
manifest=[];models=[];jobs=[];sample_profiles=[];representative=[];pathway_counts=[];members=[];original=[]
for tissue in ['SGN','Lateral_wall']:
 for safe,ct in cts.items():
  s=src/'08_pseudobulk'/tissue/safe;d=dest/tissue/safe;d.mkdir(parents=True)
  for name in ['counts_samples_by_gene.tsv.gz','metadata.tsv','gene_metadata.tsv']:
   shutil.copy2(s/name,d/name);original.append({'file':str((s/name).relative_to(src)),'sha256':sha(s/name)})
  counts=pd.read_csv(s/'counts_samples_by_gene.tsv.gz',sep='\t',index_col=0)
  meta=pd.read_csv(s/'metadata.tsv',sep='\t',index_col='sample_id')
  baseline_path=src/'09_differential'/tissue/f'S1_DS_ALL_{safe}.tsv'
  b=pd.read_csv(baseline_path,sep='\t');original.append({'file':str(baseline_path.relative_to(src)),'sha256':sha(baseline_path)})
  assert counts.index.tolist()==meta.index.tolist() and len(meta)==8
  assert meta.condition.value_counts().to_dict()=={'Naive':4,'Noise':4}
  assert np.issubdtype(counts.values.dtype,np.integer) and (counts.values>=0).all()
  fixed=counts.columns[counts.sum(axis=0)>=10].tolist();assert b.gene_id.tolist()==fixed
  assert b.gene_id.is_unique and np.isfinite(b.stat).all()
  b.to_csv(d/'BASELINE_DE.tsv.gz',sep='\t',index=False,compression='gzip')
  pd.DataFrame({'gene_id':fixed}).to_csv(d/'FROZEN_GENE_IDS.tsv',sep='\t',index=False)
  gsfile=src/'10_gsea'/tissue/safe/'S1_GSEA_ALL.tsv';g=pd.read_csv(gsfile,sep='\t')
  assert set(g.tissue)=={tissue} and set(g.cell_type)=={ct} and g.pathway_id.is_unique
  original.append({'file':str(gsfile.relative_to(src)),'sha256':sha(gsfile)})
  g.to_csv(d/'BASELINE_GSEA.tsv.gz',sep='\t',index=False,compression='gzip')
  g[['pathway_id','size']].sort_values('pathway_id').to_csv(d/'FROZEN_TESTING_FAMILY.tsv',sep='\t',index=False)
  qafile=src/'09_differential'/tissue/f'S1_DDS_QA_{safe}.tsv';qa=pd.read_csv(qafile,sep='\t').set_index('sample_id')
  qa.to_csv(d/'BASELINE_SIZE_FACTORS.tsv',sep='\t',index=True)
  objfile=src/'09_differential'/tissue/f'S1_DDS_OBJECT_{safe}.h5ad';obj=ad.read_h5ad(objfile)
  original.extend([{'file':str(p.relative_to(src)),'sha256':sha(p)} for p in [qafile,objfile]])
  assert list(obj.obs_names)==meta.index.tolist() and list(obj.var_names)==fixed
  norm=np.asarray(obj.layers['normed_counts']);calc=counts[fixed].to_numpy()/qa.loc[meta.index,'size_factor'].to_numpy()[:,None]
  assert np.allclose(norm,calc,rtol=1e-12,atol=1e-8)
  nf=pd.DataFrame(norm,index=meta.index,columns=fixed);nf.index.name='sample_id';nf.to_csv(d/'BASELINE_NORMALIZED_COUNTS.tsv.gz',sep='\t',compression='gzip')
  z=b[['gene_id','gene_symbol','stat']].dropna().copy();z['gene_symbol']=z.gene_symbol.str.strip();z=z[z.gene_symbol!=''];z['absstat']=z.stat.abs()
  z=z.sort_values(['gene_symbol','absstat','stat','gene_id'],ascending=[True,False,False,True],kind='stable').drop_duplicates('gene_symbol').set_index('gene_symbol')
  for sample,row in meta.iterrows():
   jobs.append({'tissue':tissue,'cell_type':ct,'safe_cell_type':safe,'omitted_sample_id':sample,'omitted_GSM':row.GSM,'omitted_condition':row.condition,'retained_naive':4-int(row.condition=='Naive'),'retained_noise':4-int(row.condition=='Noise'),'baseline_genes':len(fixed),'original_family_size':len(g)})
   for gene in repgenes:
    gid=z.loc[gene,'gene_id'] if gene in z.index else None
    representative.append({'tissue':tissue,'cell_type':ct,'sample_id':sample,'GSM':row.GSM,'condition':row.condition,'n_cells':row.n_cells,'size_factor':qa.loc[sample,'size_factor'],'gene_symbol':gene,'gene_id':gid,'normalized_count':nf.loc[sample,gid] if gid else np.nan,'log1p_normalized_count':np.log1p(nf.loc[sample,gid]) if gid else np.nan,'status':'AVAILABLE' if gid else 'OUTSIDE_FROZEN_GENE_UNIVERSE'})
  targets=[TNF]+([MAPK] if tissue=='Lateral_wall' and safe=='Endothelial' else [])
  for pathway in targets:
   genes=sorted(set(term.loc[term.pathway_id==pathway,'gene_symbol'])&set(z.index))
   gids=z.loc[genes,'gene_id'].tolist();arr=np.log1p(nf[gids].values);sd=arr.std(axis=0,ddof=1);valid=sd>0
   scores=((arr[:,valid]-arr[:,valid].mean(axis=0))/sd[valid]).mean(axis=1)
   for j,sample in enumerate(meta.index):
    row=meta.loc[sample];sample_profiles.append({'tissue':tissue,'cell_type':ct,'sample_id':sample,'GSM':row.GSM,'condition':row.condition,'n_cells':row.n_cells,'size_factor':qa.loc[sample,'size_factor'],'pathway_id':pathway,'mean_gene_z_log1p_normalized':scores[j],'mapped_genes':len(genes),'nonzero_variance_genes':int(valid.sum()),'zero_variance_genes':int((~valid).sum())})
    for gene,gid in zip(genes,gids):pathway_counts.append({'tissue':tissue,'cell_type':ct,'pathway_id':pathway,'gene_symbol':gene,'gene_id':gid,'sample_id':sample,'GSM':row.GSM,'condition':row.condition,'normalized_count':nf.loc[sample,gid]})
   for gene,gid,keep in zip(genes,gids,valid):members.append({'tissue':tissue,'cell_type':ct,'pathway_id':pathway,'gene_symbol':gene,'gene_id':gid,'used_in_z_summary':bool(keep)})
  models.append({'tissue':tissue,'cell_type':ct,'safe_cell_type':safe,'samples':8,'genes':len(fixed),'original_family_size':len(g),'normalization_max_abs_error':float(np.max(np.abs(norm-calc)))})
for name,rows in [('MODEL_MANIFEST.tsv',models),('JOB_MANIFEST.tsv',jobs),('ORIGINAL_INPUT_SHA256.tsv',original)]:pd.DataFrame(rows).to_csv(dest/name,sep='\t',index=False)
profiles=out/'sample_profiles';profiles.mkdir()
for name,rows in [('PER_ANIMAL_PATHWAY_SCORES.tsv',sample_profiles),('PER_ANIMAL_REPRESENTATIVE_GENES.tsv',representative),('PATHWAY_MEMBER_SELECTION.tsv',members)]:pd.DataFrame(rows).to_csv(profiles/name,sep='\t',index=False)
pd.DataFrame(pathway_counts).to_csv(profiles/'PER_ANIMAL_PATHWAY_GENE_COUNTS.tsv.gz',sep='\t',index=False,compression='gzip')
for p in sorted(dest.rglob('*')):
 if p.is_file():manifest.append({'file':str(p.relative_to(dest)),'sha256':sha(p)})
pd.DataFrame(manifest).to_csv(out/'FROZEN_INPUT_SHA256.tsv',sep='\t',index=False)
(out/'INPUT_PREFLIGHT.json').write_text(json.dumps({'status':'PASS','models':models,'jobs':len(jobs),'original_source_files_checked':len(original),'main_analysis_refitted':False},indent=2))
print(json.dumps({'status':'INPUT_PREFLIGHT_PASS','models':models,'jobs':len(jobs)},indent=2),flush=True)
