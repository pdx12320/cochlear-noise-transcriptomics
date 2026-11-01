from pathlib import Path
import argparse,json,hashlib,math
import pandas as pd,numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--require-all',action='store_true');a=p.parse_args();root=Path(a.root)
def bh(v,n=None):
 v=np.asarray(v,float);n=len(v) if n is None else n;o=np.argsort(v,kind='stable');q=np.minimum.accumulate((v[o]*n/np.arange(1,len(v)+1))[::-1])[::-1];z=np.empty(len(v));z[o]=np.minimum(q,1);return z
mf=pd.read_csv(root/'FROZEN_INPUT_SHA256.tsv',sep='\t')
for r in mf.itertuples():assert hashlib.sha256((root/'inputs'/r.file).read_bytes()).hexdigest()==r.sha256
term=pd.read_csv(root/'inputs/MSIGDB_2026.1.Mm_Hallmark_Reactome_TERM2GENE.tsv.gz',sep='\t');sets={k:set(v.gene_symbol) for k,v in term.groupby('pathway_id')}
jobs=pd.read_csv(root/'inputs/JOB_MANIFEST.tsv',sep='\t');checks=[]
for j in jobs.itertuples():
 d=root/'fits'/j.tissue/j.safe_cell_type/j.omitted_sample_id
 if not (d/'gsea/STATUS.txt').is_file():
  assert not a.require_all,str(d)
  continue
 assert (d/'gsea/STATUS.txt').read_text().strip()=='COMPLETED_LOO_GSEA'
 inp=root/'inputs'/j.tissue/j.safe_cell_type
 x=pd.read_csv(d/'LOO_DE_ALL.tsv.gz',sep='\t');mask=pd.read_csv(inp/'FROZEN_GENE_IDS.tsv',sep='\t').gene_id.tolist()
 assert x.gene_id.tolist()==mask and len(x)==j.baseline_genes
 m=pd.read_csv(d/'RETAINED_SAMPLES.tsv',sep='\t');fullm=pd.read_csv(inp/'metadata.tsv',sep='\t');expected=fullm[fullm.sample_id!=j.omitted_sample_id]
 assert m.sample_id.tolist()==expected.sample_id.tolist() and len(m)==7
 assert m.condition.value_counts().to_dict()=={'Naive':j.retained_naive,'Noise':j.retained_noise}
 sf=pd.read_csv(d/'LOO_SIZE_FACTORS.tsv',sep='\t');raw=pd.read_csv(inp/'counts_samples_by_gene.tsv.gz',sep='\t',index_col=0).loc[m.sample_id,mask]
 assert sf.sample_id.tolist()==m.sample_id.tolist()
 positive=(raw.values>0).all(axis=0);logs=np.log(raw.values[:,positive]);sfc=np.exp(np.median(logs-logs.mean(axis=0),axis=1))
 sferr=float(np.max(np.abs(sfc-sf.size_factor.values)));assert np.allclose(sfc,sf.size_factor,rtol=1e-10,atol=1e-10),sferr
 mean=(raw.values/sf.size_factor.values[:,None]).mean(axis=0);assert np.allclose(x.baseMean,mean,rtol=1e-10,atol=1e-7)
 pv=x.pvalue.notna();perr=float(np.max(np.abs(np.array([math.erfc(abs(t)/math.sqrt(2)) for t in x.loc[pv,'stat']])-x.loc[pv,'pvalue'].values)));assert perr<1e-10
 v=x.padj.notna();derr=float(np.max(np.abs(bh(x.loc[v,'pvalue'])-x.loc[v,'padj'])));assert derr<1e-10
 ranks=pd.read_csv(d/'gsea/RANKING_WALD_STAT.tsv.gz',sep='\t');z=x.loc[x.gene_symbol.notna()&np.isfinite(x.stat),['gene_id','gene_symbol','stat']].copy();z.gene_symbol=z.gene_symbol.str.strip();z=z[z.gene_symbol!=''];z['absstat']=z.stat.abs();z=z.sort_values(['gene_symbol','absstat','stat','gene_id'],ascending=[True,False,False,True],kind='stable').drop_duplicates('gene_symbol').sort_values(['stat','gene_symbol','gene_id'],ascending=[False,True,True],kind='stable')
 assert ranks.gene_id.tolist()==z.gene_id.tolist() and ranks.gene_symbol.tolist()==z.gene_symbol.tolist();assert np.allclose(ranks.stat,z.stat,rtol=1e-12,atol=1e-12)
 g=pd.read_csv(d/'gsea/LOO_GSEA_ALL.tsv.gz',sep='\t');fam=pd.read_csv(inp/'FROZEN_TESTING_FAMILY.tsv',sep='\t')
 assert len(g)==len(fam)==j.original_family_size and set(g.pathway_id)==set(fam.pathway_id)
 v=g.P.notna();gerr=float(np.max(np.abs(bh(g.loc[v,'P'],len(fam))-g.loc[v,'FDR'])));assert gerr<1e-10
 universe=set(ranks.gene_symbol)
 for r in g.itertuples():
  assert r.overlap_size==len(sets[r.pathway_id]&universe)
  if pd.notna(r.P):
   le=r.leadingEdge.split(',') if isinstance(r.leadingEdge,str) else []
   assert set(le)<=sets[r.pathway_id]&universe and len(le)==len(set(le))
   assert np.sign(r.NES)==np.sign(r.ES) and 10<=r.overlap_size<=500
 checks.append({'tissue':j.tissue,'cell_type':j.cell_type,'omitted_sample_id':j.omitted_sample_id,'omitted_GSM':j.omitted_GSM,'genes':len(x),'ranked_symbols':len(ranks),'original_family_size':len(fam),'tested_pathways':int(v.sum()),'size_factor_max_abs_error':sferr,'gene_P_max_abs_error':perr,'gene_BH_max_abs_error':derr,'pathway_BH_max_abs_error':gerr,'identities_rankings_memberships_pass':True})
if a.require_all:assert len(checks)==48
report={'status':'PASS','scope':'Independent arithmetic/identity checks only; no model refit','completed_fits_checked':len(checks),'frozen_input_files_checked':len(mf),'checks':checks}
(root.parent/('FINAL_INDEPENDENT_AUDIT.json' if a.require_all else 'PILOT_INDEPENDENT_AUDIT.json')).write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
