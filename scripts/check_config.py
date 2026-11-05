#!/usr/bin/env python3
"""Validate repository paths; optionally resolve/check a private machine profile.

This command reads metadata only. It never launches GWAS/fine-mapping jobs or SSH.
Absolute machine paths are loaded from ignored config/local_paths.yaml.
"""
import argparse,json
from pathlib import Path
import yaml

def public_paths(repo):
    repo=Path(repo);config=yaml.safe_load((repo/'config/paths.yaml').read_text())
    for group in ['repository_paths','result_paths','historical_result_paths','evidence_paths']:
        for key,value in config[group].items():
            p=Path(value)
            if p.is_absolute() or '..' in p.parts:raise ValueError('Public paths must be repository-relative: '+key)
            if not (repo/p).exists():raise ValueError('Missing public artifact: '+value)
    methods=yaml.safe_load((repo/'config/methods.yaml').read_text())['methods']
    for method,settings in methods.items():
        for profile,pattern in settings['output_patterns'].items():
            sample=Path(pattern.format(phenotype=1))
            if sample.is_absolute() or '..' in sample.parts:raise ValueError('Invalid input pattern: '+method)
    return config,methods

def resolve_profile(repo,profile,local_config=None,check_inputs=False,stage=None,historical=False):
    config,methods=public_paths(repo)
    local=Path(local_config) if local_config else Path(repo)/'config/local_paths.yaml'
    if not local.exists():raise ValueError('Private machine profile is absent. See config/README.md.')
    local=yaml.safe_load(local.read_text());settings=local['profiles'][profile]
    root=Path(settings['benchmark_root'])
    selected={key:item for key,item in config['restricted_inputs'].items() if profile in item['available_on']}
    inputs={key:str(root/item['path']) for key,item in selected.items()}
    revision_inputs={k:settings[k] for k in ['simulation_directory','genotype_prefix','gcta_directory','rss_directory'] if k in settings}
    selected_methods={m:s for m,s in methods.items() if not historical or m!='GCTA_MLMA'}
    patterns={m:s.get('historical_output_patterns',s['output_patterns'])[profile] if historical else s['output_patterns'][profile] for m,s in selected_methods.items()}
    gwas={m:[str((Path(revision_inputs['gcta_directory'])/Path(pattern.format(phenotype=i)).name) if m=='GCTA_MLMA' and 'gcta_directory' in revision_inputs else root/pattern.format(phenotype=i)) for i in range(1,11)] for m,pattern in patterns.items()}
    stage_inputs={'summarize':{'effect_windows','causal_snp_list','slemm_phenotypes'},
                  'equation13':{'equation13_cached_theory'},'ld':{'ld_pairs_10mb'},
                  'manhattan':set(),'fine-mapping':{'causal_snp_list','susie_adj_saved_pips','susie_adj_diagnostics','saved_maf_annotations'}}
    check_gwas=stage in [None,'summarize','manhattan','equation13']
    optional_missing=[]
    if check_inputs:
        if not root.is_dir():raise ValueError('Source root is not mounted on this machine; run the input check on its host: '+str(root))
        missing=[]
        for key,path in inputs.items():
            if stage is not None and key not in stage_inputs.get(stage,set()):continue
            item=selected[key]
            checks=[Path(path+s) for s in ['.bed','.bim','.fam']] if item['kind']=='plink_bed_prefix' else [Path(path)]
            missing.extend(str(p) for p in checks if not p.exists())
        if check_gwas:
            for method,paths in gwas.items():
                if stage=='equation13' and method not in ['BOLT-LMM','REGENIE_LOCO','REGENIE_WO_LOCO','SLEMM']:continue
                absent=[p for p in paths if not Path(p).is_file()]
                (missing if selected_methods[method].get('required_for_saved_input_check',True) else optional_missing).extend(absent)
        if missing:raise ValueError('Missing configured inputs: '+', '.join(missing))
    return {'profile':profile,'host':settings.get('ssh_host'),'benchmark_root':str(root),
        'repository_root':local['repository_root'],'public_results_root':str(Path(local['repository_root'])/config['repository_paths']['derived_data']),
        'inputs':inputs,'gwas_patterns':{m:str(Path(revision_inputs['gcta_directory'])/Path(pattern).name if m=='GCTA_MLMA' and 'gcta_directory' in revision_inputs else root/pattern) for m,pattern in patterns.items()},
        'revision_inputs':revision_inputs,'analysis_profile':'historical_audit' if historical else 'current_inclusive_revision',
        'gwas_files_checked':sum(Path(p).is_file() for paths in gwas.values() for p in paths) if check_inputs and check_gwas else 0,
        'missing_optional_gwas_files':optional_missing,
        'input_check':'PASSED' if check_inputs else 'NOT_RUN',
        'unavailable_on_this_profile':[k for k,v in config['restricted_inputs'].items() if profile not in v['available_on']]}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--profile',choices=['mac','hpcc']);p.add_argument('--local-config',type=Path);p.add_argument('--check-inputs',action='store_true');p.add_argument('--historical',action='store_true');a=p.parse_args()
    if a.profile:print(json.dumps(resolve_profile(a.repo,a.profile,a.local_config,a.check_inputs,historical=a.historical),indent=2))
    else:
        if a.check_inputs:p.error('--check-inputs requires --profile')
        public_paths(a.repo);print('All configured public paths exist; method output patterns are relative and valid.')
