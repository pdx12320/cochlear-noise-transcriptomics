#!/usr/bin/env python3
"""Sequential bounded LOO controller. Stops on any failure; preserves completed jobs."""
from pathlib import Path
import argparse,csv,json,os,subprocess,time,hashlib
p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--python',required=True);p.add_argument('--rscript',required=True);a=p.parse_args()
root=Path(a.root);jobs=list(csv.DictReader((root/'inputs/JOB_MANIFEST.tsv').open(),delimiter='\t'))
assert len(jobs)==48
state=root/'BATCH_STATUS.json'
if state.exists():raise RuntimeError('Refusing existing batch status; inspect before restarting')
start=time.monotonic();completed=[]
env=os.environ.copy();env.update({'LC_ALL':'C','LANG':'C','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'})
for i,j in enumerate(jobs):
 if time.monotonic()-start>3600:raise RuntimeError('One-hour batch budget reached; pause for resource review')
 t,c,s=j['tissue'],j['safe_cell_type'],j['omitted_sample_id'];fit=root/'fits'/t/c/s
 done=fit/'gsea/STATUS.txt'
 if done.is_file() and done.read_text().strip()=='COMPLETED_LOO_GSEA':
  completed.append({'tissue':t,'safe_cell_type':c,'omitted_sample_id':s,'status':'PREVIOUS_PILOT_COMPLETE'})
 else:
  if fit.exists():raise RuntimeError(f'Incomplete existing fit requires inspection: {t}/{c}/{s}')
  commands=[('fit',[a.python,str(root/'scripts/run_loo_fit.py'),'--root',str(root),'--tissue',t,'--cell',c,'--omit',s]),('gsea',[a.rscript,'--vanilla',str(root/'scripts/run_loo_gsea.R'),str(root),t,c,s])]
  wall={}
  for stage,cmd in commands:
   tic=time.monotonic()
   with (root/'logs'/f'{t}_{c}_{s}_{stage}.log').open('w') as h:
    r=subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT,env=env,timeout=600)
   wall[stage+'_wall_seconds']=time.monotonic()-tic
   if r.returncode:raise RuntimeError(f'{stage} failed for {t}/{c}/{s}; see preserved log')
  completed.append({'tissue':t,'safe_cell_type':c,'omitted_sample_id':s,'status':'COMPLETED',**wall})
 status={'status':'RUNNING','completed':len(completed),'total':len(jobs),'elapsed_seconds':time.monotonic()-start,'jobs':completed}
 state.write_text(json.dumps(status,indent=2));print(json.dumps({k:status[k] for k in ['status','completed','total','elapsed_seconds']}),flush=True)
status['status']='COMPLETED_ALL_48';state.write_text(json.dumps(status,indent=2));print('COMPLETED_ALL_48',flush=True)
