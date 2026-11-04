#!/usr/bin/env python3
"""Supported benchmark workflow. Run from any directory; see docs/workflow.md."""
import argparse,subprocess,sys
from pathlib import Path
from check_config import resolve_profile
ROOT=Path(__file__).resolve().parents[1]

def run(script,*args,stdout=None):
 command=[sys.executable,str(ROOT/script),*map(str,args)]
 print('Running:',script,flush=True)
 subprocess.run(command,cwd=ROOT,stdout=stdout,check=True)

def new_output(path,source=None,resume=False):
 path=Path(path).resolve()
 if path.is_relative_to(ROOT):
  raise ValueError('Use a separate output directory; committed results are preserved.')
 if source and path.is_relative_to(Path(source).resolve()):raise ValueError('Output must be outside the historical source tree.')
 if path.exists() and not resume:raise ValueError('Output exists. Choose a new directory; --resume is supported only for summaries.')
 path.mkdir(parents=True,exist_ok=True)
 return path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('stage',choices=['check','summarize','equation13','ld','fine-mapping','figures','manhattan'])
 p.add_argument('--profile',choices=['mac','hpcc']);p.add_argument('--local-config',type=Path)
 p.add_argument('--out',type=Path);p.add_argument('--data',type=Path,default=ROOT/'derived_data')
 p.add_argument('--resume',action='store_true');a=p.parse_args()
 if a.resume and a.stage!='summarize':p.error('--resume applies only to summarize')
 if a.stage=='check':
  run('scripts/check_config.py');run('scripts/build_result_config.py','--check')
  run('tests/check_static.py')
  subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=ROOT,check=True)
  run('scripts/validate_manuscript.py','--check');run('tests/check_validation_status.py');return
 if a.out is None:p.error('--out is required for analysis and figures')
 if a.stage=='figures':
  data=a.data.resolve()
  if not data.is_dir():p.error('--data directory does not exist')
  out=new_output(a.out)
  run('scripts/13_figures/plot_aggregates.py','--data',data,'--out',out)
  subprocess.run(['Rscript',str(ROOT/'scripts/13_figures/plot_ld.R'),str(data),str(out)],check=True);return
 if a.profile is None:p.error('--profile is required for restricted-input stages')
 profile=resolve_profile(ROOT,a.profile,a.local_config,check_inputs=True)
 source=profile['benchmark_root'];inputs=profile['inputs']
 if a.stage=='equation13' and 'equation13_cached_theory' not in inputs:p.error('This profile has no verified cached-theory file; use the Mac profile.')
 if a.stage=='ld' and 'ld_pairs_10mb' not in inputs:p.error('This profile has no verified pairwise LD file; run on HPCC.')
 if a.stage=='fine-mapping' and 'susie_adj_saved_pips' not in inputs:p.error('This profile has no verified fine-mapping files; run on HPCC.')
 out=new_output(a.out,source,a.resume)
 if a.stage=='summarize':run('scripts/benchmark_core.py','summarize','--root',source,'--out',out/'derived_data','--cache',out/'.restricted_cache')
 elif a.stage=='equation13':run('scripts/benchmark_core.py','equation13','--root',source,'--theory',inputs['equation13_cached_theory'],'--out',out/'equation13_summary')
 elif a.stage=='manhattan':run('scripts/13_figures/plot_manhattan.py','--root',source,'--out',out)
 elif a.stage in ['ld','fine-mapping']:
  name='ld' if a.stage=='ld' else 'finemapping';evidence=out/(name+'_verification.json')
  with evidence.open('w') as f:
   if name=='ld':run('scripts/02_ld_analysis/verify_ld_means.py','--input',inputs['ld_pairs_10mb'],stdout=f)
   else:run('scripts/12_finemapping/verify_saved_results.py','--root',source,stdout=f)
  run('scripts/export_verification.py','--kind',name,'--json',evidence,'--out',out/(name+'_summary'))
if __name__=='__main__':main()
