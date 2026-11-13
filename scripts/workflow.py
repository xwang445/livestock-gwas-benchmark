#!/usr/bin/env python3
"""Current inclusive-threshold revision workflow; --historical selects the old audit."""
import argparse,subprocess,sys
from pathlib import Path
from check_config import resolve_profile
ROOT=Path(__file__).resolve().parents[1]

def run(script,*args,stdout=None):
 command=[sys.executable,str(ROOT/script),*map(str,args)]
 print('Running:',script,flush=True)
 subprocess.run(command,cwd=ROOT,stdout=stdout,check=True)

def new_output(path,source=None,resume=False,create=True):
 path=Path(path).resolve()
 if path.is_relative_to(ROOT):
  raise ValueError('Use a separate output directory; committed results are preserved.')
 if source and path.is_relative_to(Path(source).resolve()):raise ValueError('Output must be outside the historical source tree.')
 if path.exists() and not resume:raise ValueError('Output exists. Choose a new directory; --resume is supported only for summaries.')
 if create:path.mkdir(parents=True,exist_ok=True)
 return path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('stage',choices=['check','summarize','equation13','ld','fine-mapping','figures','submission','manhattan'])
 p.add_argument('--profile',choices=['mac','hpcc']);p.add_argument('--local-config',type=Path)
 p.add_argument('--out',type=Path);p.add_argument('--data',type=Path,default=ROOT/'derived_data')
 p.add_argument('--historical',action='store_true',help='Explicitly reproduce the preserved pre-revision nine-configuration audit')
 p.add_argument('--simulation',type=Path);p.add_argument('--genotype-prefix',type=Path)
 p.add_argument('--gcta-directory',type=Path);p.add_argument('--rss-directory',type=Path)
 p.add_argument('--allow-incomplete',action='store_true',help='Label incomplete configuration panels in diagnostic Manhattan plots')
 p.add_argument('--individual',action='store_true',help='Manhattan: also emit one PDF/PNG for each completed configuration/phenotype')
 p.add_argument('--resume',action='store_true');a=p.parse_args()
 if a.resume and a.stage!='summarize':p.error('--resume applies only to summarize')
 if a.resume and not a.historical:p.error('Current derivation uses a new output directory; --resume is only for the historical hash-checked cache')
 if a.individual and (a.stage!='manhattan' or a.historical):p.error('--individual applies only to current Manhattan plots')
 if a.stage=='check':
  run('scripts/check_config.py');run('scripts/build_result_config.py','--check')
  run('tests/check_static.py')
  subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=ROOT,check=True)
  run('scripts/validate_manuscript.py','--check');run('tests/check_validation_status.py');return
 if a.out is None:p.error('--out is required for analysis and figures')
 if a.stage=='submission':
  source=a.data.resolve()
  if (source/'submission_2026-11-12').is_dir():source=source/'submission_2026-11-12'
  if not (source/'figure08_pooled_metrics.tsv').is_file():p.error('--data must contain the submission_2026-11-12 aggregate sources')
  if a.historical:p.error('--historical does not apply to submission artwork')
  out=new_output(a.out)
  for name in ['plot_figures_1_2_7.py','plot_figures_3_5_6.py','plot_figure8.py']:
   run('scripts/submission_2026-11-12/'+name,'--source',source,'--output',out)
  run('scripts/submission_2026-11-12/make_supplementary_table_s1.py','--tables',source/'tables','--output',out/'supplement/tables')
  print('Public aggregate figures and Table S1 complete. Figure 4/S1–S10 require the documented restricted-input Manhattan command.',flush=True)
  return
 if a.stage=='figures':
  data=a.data.resolve()
  if not data.is_dir():p.error('--data directory does not exist')
  out=new_output(a.out)
  args=['--data',data,'--out',out]
  if a.historical:args.append('--historical')
  run('scripts/13_figures/plot_aggregates.py',*args)
  ld_data=data if (data/'ld_summary/ld_decay_10kb.tsv').is_file() else ROOT/'derived_data'
  subprocess.run(['Rscript',str(ROOT/'scripts/13_figures/plot_ld.R'),str(ld_data),str(out)],check=True);return
 if a.profile is None:p.error('--profile is required for restricted-input stages')
 profile=resolve_profile(ROOT,a.profile,a.local_config,check_inputs=True,stage=a.stage,historical=a.historical)
 source=profile['benchmark_root'];inputs=profile['inputs']
 if a.stage=='equation13' and 'equation13_cached_theory' not in inputs:p.error('This profile has no verified cached-theory file; use the Mac profile.')
 if a.stage=='ld' and 'ld_pairs_10mb' not in inputs:p.error('This profile has no verified pairwise LD file; run on HPCC.')
 if a.stage=='fine-mapping' and 'susie_adj_saved_pips' not in inputs:p.error('This profile has no verified fine-mapping files; run on HPCC.')
 options=profile.get('revision_inputs',{})
 simulation=a.simulation or options.get('simulation_directory')
 genotype=a.genotype_prefix or options.get('genotype_prefix')
 gcta=a.gcta_directory or options.get('gcta_directory') or Path(source)/'GCTA_mlma'
 rss=a.rss_directory or options.get('rss_directory')
 if a.stage=='summarize' and not a.historical and (not simulation or not genotype):
  p.error('Current summaries require verified --simulation and --genotype-prefix, or those keys in the private profile. No source paths are inferred.')
 if a.stage=='fine-mapping' and not a.historical and not rss:
  p.error('Current fine-mapping verification requires --rss-directory pointing to the completed ERV_TRUE run, or its private-profile key.')
 out=new_output(a.out,source,a.resume,create=not(a.stage=='manhattan' and not a.historical))
 if a.stage=='summarize':
  if a.historical:
   run('scripts/benchmark_core.py','summarize','--historical','--root',source,'--out',out/'derived_data','--cache',out/'.restricted_cache')
  else:
   run('scripts/revision_2026-11-05/association/derive_association_tables.py','--benchmark',source,'--simulation',simulation,'--genotype-prefix',genotype,'--gcta-directory',gcta,'--output',out)
   run('scripts/revision_2026-11-05/association/derive_power_sensitivity.py','--benchmark',source,'--gcta-directory',gcta,'--tables',out/'tables')
   run('scripts/revision_2026-11-05/association/derive_spatial_closed.py','--benchmark',source,'--gcta-directory',gcta,'--tables',out/'tables')
   run('scripts/revision_2026-11-05/association/derive_overlap_summary.py','--tables',out/'tables')
 elif a.stage=='equation13':
  args=['equation13','--root',source,'--theory',inputs['equation13_cached_theory'],'--out',out/'equation13_summary']
  if a.historical:args.append('--historical')
  run('scripts/benchmark_core.py',*args)
 elif a.stage=='manhattan':
  args=['--root',source,'--out',out]
  if a.historical:args.append('--historical')
  else:
   args.extend(['--gcta-dir',gcta])
   if not a.allow_incomplete:args.append('--require-complete')
   if a.individual:args.append('--individual')
  run('scripts/13_figures/plot_manhattan.py',*args)
 elif a.stage in ['ld','fine-mapping']:
  name='ld' if a.stage=='ld' else 'finemapping';evidence=out/(name+'_verification.json')
  with evidence.open('w') as f:
   if name=='ld':run('scripts/02_ld_analysis/verify_ld_means.py','--input',inputs['ld_pairs_10mb'],stdout=f)
   else:
    args=['--root',source]
    if a.historical:args.append('--historical')
    else:args.extend(['--rss-root',Path(rss)/'results'])
    run('scripts/12_finemapping/verify_saved_results.py',*args,stdout=f)
  run('scripts/export_verification.py','--kind',name,'--json',evidence,'--out',out/(name+'_summary'))
if __name__=='__main__':main()
