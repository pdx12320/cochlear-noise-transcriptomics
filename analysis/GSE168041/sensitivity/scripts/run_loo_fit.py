#!/usr/bin/env python3
"""Fit one authorized animal-deletion sensitivity model; never overwrite outputs."""
from pathlib import Path
import argparse,hashlib,json,sys,time,warnings,platform,resource
from importlib.metadata import version
import numpy as np,pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--tissue',required=True);p.add_argument('--cell',required=True);p.add_argument('--omit',required=True);a=p.parse_args()
root=Path(a.root);src=root/'inputs'/a.tissue/a.cell;out=root/'fits'/a.tissue/a.cell/a.omit
if out.exists():raise RuntimeError('Refusing existing fit output')
assert version('pydeseq2')=='0.5.4' and sys.version_info[:3]==(3,12,13)
manifest=pd.read_csv(root/'FROZEN_INPUT_SHA256.tsv',sep='\t')
for row in manifest.itertuples():
 f=root/'inputs'/row.file
 assert hashlib.sha256(f.read_bytes()).hexdigest()==row.sha256,f
counts=pd.read_csv(src/'counts_samples_by_gene.tsv.gz',sep='\t',index_col=0)
meta=pd.read_csv(src/'metadata.tsv',sep='\t',index_col='sample_id')
baseline=pd.read_csv(src/'BASELINE_DE.tsv.gz',sep='\t').set_index('gene_id')
fixed=pd.read_csv(src/'FROZEN_GENE_IDS.tsv',sep='\t').gene_id.tolist()
assert counts.index.tolist()==meta.index.tolist() and a.omit in meta.index and len(meta)==8
assert baseline.index.tolist()==fixed
omit=meta.loc[a.omit].copy();keep=[s for s in meta.index if s!=a.omit]
filtered=counts.loc[keep,fixed].astype(np.int64);metadata=meta.loc[keep].copy()
metadata['condition']=pd.Categorical(metadata.condition,categories=['Naive','Noise'],ordered=True)
assert sorted(metadata.condition.value_counts().tolist())==[3,4] and filtered.shape==(7,len(fixed))
out.mkdir(parents=True);(out/'STATUS.txt').write_text('INCOMPLETE_DO_NOT_USE\n')
metadata.reset_index().to_csv(out/'RETAINED_SAMPLES.tsv',sep='\t',index=False)
np.random.seed(20260825);start=time.monotonic()
with warnings.catch_warnings(record=True) as ws:
 warnings.simplefilter('always')
 inference=DefaultInference(n_cpus=8)
 dds=DeseqDataSet(counts=filtered,metadata=metadata,design='~condition',refit_cooks=True,inference=inference,low_memory=False)
 dds.deseq2()
 ds=DeseqStats(dds,contrast=['condition','Noise','Naive'],alpha=.05,cooks_filter=True,independent_filter=True,inference=inference,quiet=False)
 ds.summary()
 result=ds.results_df.copy();result.index.name='gene_id'
 assert result.index.tolist()==fixed
 result.insert(0,'gene_symbol',baseline.loc[result.index,'gene_symbol'].to_numpy());result.insert(0,'gene_id',result.index)
 result=result[['gene_id','gene_symbol','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']]
 result.to_csv(out/'LOO_DE_ALL.tsv.gz',sep='\t',index=False,compression='gzip')
 sf=pd.DataFrame({'sample_id':dds.obs_names,'GSM':meta.loc[dds.obs_names,'GSM'].values,'condition':meta.loc[dds.obs_names,'condition'].values,'size_factor':dds.obs['size_factors'].values,'n_cells':meta.loc[dds.obs_names,'n_cells'].values})
 sf.to_csv(out/'LOO_SIZE_FACTORS.tsv',sep='\t',index=False)
 nr=np.asarray(dds.layers['normed_counts']);assert np.allclose(nr,filtered.values/sf.size_factor.values[:,None],rtol=1e-12,atol=1e-8)
 elapsed=time.monotonic()-start
 warning_text=[str(w.message) for w in ws]
(out/'WARNINGS.txt').write_text('\n'.join(warning_text)+'\n' if warning_text else '')
info={'status':'PASS','tissue':a.tissue,'safe_cell_type':a.cell,'omitted_sample_id':a.omit,'omitted_GSM':omit.GSM,'omitted_condition':omit.condition,'retained_naive':int((metadata.condition=='Naive').sum()),'retained_noise':int((metadata.condition=='Noise').sum()),'genes':len(result),'all_zero_genes_after_omission':int((filtered.sum(axis=0)==0).sum()),'finite_wald_stats':int(np.isfinite(result.stat).sum()),'nonmissing_pvalues':int(result.pvalue.notna().sum()),'nonmissing_padj':int(result.padj.notna().sum()),'significant_genes':int((result.padj<.05).sum()),'elapsed_seconds':elapsed,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':platform.python_version(),'pydeseq2':version('pydeseq2'),'numpy':version('numpy'),'pandas':version('pandas'),'scipy':version('scipy'),'n_inference_workers':8,'seed':20260825,'original_gene_mask_preserved':True,'warnings':warning_text}
(out/'FIT_RECORD.json').write_text(json.dumps(info,indent=2));(out/'STATUS.txt').write_text('COMPLETED_LOO_DE\n')
print(json.dumps(info,indent=2),flush=True)
