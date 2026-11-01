from pathlib import Path
import sys,json,hashlib
import anndata as ad
import pandas as pd
import numpy as np
src=Path(sys.argv[1]); out=Path(sys.argv[2])
out.mkdir(parents=True,exist_ok=False)
records=[]
for tissue,n in [('SGN',31556),('Lateral_wall',33149)]:
 p=src/(tissue+'_annotated.h5ad')
 def sha(path):
  h=hashlib.sha256()
  with path.open('rb') as f:
   for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
  return h.hexdigest()
 before=sha(p)
 a=ad.read_h5ad(p,backed='r')
 assert a.n_obs==n
 key=next((k for k in ['X_umap','X_tsne'] if k in a.obsm),None)
 assert key, list(a.obsm.keys())
 coords=np.asarray(a.obsm[key])[:,:2]
 assert np.isfinite(coords).all()
 names=['sample_id','GSM','condition','tissue','replicate','leiden','broad_annotation']
 frame=a.obs[names].copy();frame.insert(0,'cell',a.obs_names.astype(str));frame['dim1']=coords[:,0];frame['dim2']=coords[:,1]
 frame.to_csv(out/(tissue+'_embedding.tsv.gz'),sep='\t',index=False,compression={'method':'gzip','mtime':0})
 a.file.close()
 assert sha(p)==before
 records.append({'tissue':tissue,'cells':n,'source_object':p.name,'source_sha256':before,'embedding_key':key,'no_coordinates_recomputed':True,'broad_counts':frame.broad_annotation.value_counts().to_dict(),'output_sha256':sha(out/(tissue+'_embedding.tsv.gz'))})
(out/'EMBEDDING_EXPORT.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps([{k:r[k] for k in ['tissue','cells','embedding_key']} for r in records]))
