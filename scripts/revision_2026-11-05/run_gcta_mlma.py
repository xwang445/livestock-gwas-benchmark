#!/usr/bin/env python3
"""Run verified full-GRM MLMA without replacing source data or completed fits.

All manifests, logs and per-SNP outputs are restricted and remain outside Git.
The per-output-directory lock prevents a duplicate invocation of this runner.
Check the host/scheduler for matching jobs before using another output directory.
"""
import argparse,csv,fcntl,hashlib,json,math,os,subprocess,time
from pathlib import Path

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def rows(path):
    with path.open() as f:return [line.split() for line in f if line.strip()]

def validate(prefix,bim):
    log=Path(str(prefix)+'.log');mlma=Path(str(prefix)+'.mlma')
    if not log.is_file() or not mlma.is_file():return False
    text=log.read_text()
    if not ('Analysis finished' in text and '20000 individuals are in common' in text):return False
    with mlma.open() as f:
        d=list(csv.DictReader(f,delimiter='\t'))
    if len(d)!=len(bim):return False
    for r,b in zip(d,bim):
        if [r['Chr'],r['SNP'],r['bp'],r['A1'],r['A2']]!=[b[0],b[1],b[3],b[4],b[5]]:return False
        p=float(r['p'])
        if not math.isfinite(p) or not 0<=p<=1:return False
    return True

def main():
    p=argparse.ArgumentParser();p.add_argument('--gcta',required=True,type=Path)
    p.add_argument('--bfile',required=True,type=Path);p.add_argument('--grm',required=True,type=Path)
    p.add_argument('--phenotypes',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--threads',type=int,default=16);p.add_argument('--first',type=int,default=1);p.add_argument('--last',type=int,default=10)
    a=p.parse_args();assert 1<=a.first<=a.last<=10 and a.threads>0
    a.out.mkdir(parents=True,exist_ok=True)
    with (a.out/'mlma_runner.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        fam=rows(Path(str(a.bfile)+'.fam'));bim=rows(Path(str(a.bfile)+'.bim'))
        ids=[tuple(x[:2]) for x in fam];idset=set(ids)
        grmids=[tuple(x[:2]) for x in rows(Path(str(a.grm)+'.grm.id'))]
        assert len(ids)==31301 and len(set(ids))==len(ids) and ids==grmids
        assert len(bim)==135358 and {int(x[0]) for x in bim}==set(range(1,19))
        expected_bytes=4*len(ids)*(len(ids)+1)//2
        inputs=[Path(str(a.bfile)+x) for x in ('.bed','.bim','.fam')]
        inputs += [Path(str(a.grm)+x) for x in ('.grm.bin','.grm.N.bin','.grm.id')]
        assert all(Path(str(a.grm)+x).stat().st_size==expected_bytes for x in ('.grm.bin','.grm.N.bin'))
        for ph in range(a.first,a.last+1):
            file=a.phenotypes/f'pheno.{ph}.txt';r=rows(file)
            if r[0][:2]==['FID','IID']:r=r[1:]
            assert len(r)==20000 and len({tuple(x[:2]) for x in r})==20000
            assert all(tuple(x[:2]) in idset and math.isfinite(float(x[2])) and float(x[2])!=-9 for x in r)
            inputs.append(file)
        banner=subprocess.run([str(a.gcta)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True).stdout
        assert '1.94.1' in banner,'Expected verified GCTA 1.94.1'
        (a.out/'gcta_banner.txt').write_text(banner)
        signature={'sha256':{str(x.resolve()):digest(x) for x in inputs},'threads':a.threads,'model':'--mlma; full GRM; intercept only'}
        state=a.out/'input_signature.json'
        if state.exists():assert json.loads(state.read_text())==signature,'Inputs/settings differ: use a new output directory'
        else:
            assert not list(a.out.glob('*.mlma')),'Unverified existing MLMA outputs: preserve and verify separately'
            state.write_text(json.dumps(signature,indent=2)+'\n')
        for ph in range(a.first,a.last+1):
            prefix=a.out/f'spg.even.100k.pheno{ph}';status=a.out/f'pheno{ph}.status.json'
            if validate(prefix,bim):
                print(f'Phenotype {ph}: completed output verified, skipping',flush=True);continue
            assert not Path(str(prefix)+'.mlma').exists(),'Incomplete existing output: preserve it before retrying'
            cmd=[str(a.gcta),'--bfile',str(a.bfile),'--grm',str(a.grm),'--mlma','--pheno',str(a.phenotypes/f'pheno.{ph}.txt'),'--thread-num',str(a.threads),'--out',str(prefix)]
            record={'command':cmd,'pid':os.getpid(),'start_unix':time.time(),'status':'RUNNING'}
            status.write_text(json.dumps(record,indent=2)+'\n')
            with (a.out/f'pheno{ph}.console.log').open('w') as log:
                r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
            record.update(exit_code=r.returncode,end_unix=time.time(),status='COMPLETED' if r.returncode==0 and validate(prefix,bim) else 'FAILED_VALIDATION')
            status.write_text(json.dumps(record,indent=2)+'\n')
            if record['status']!='COMPLETED':raise RuntimeError(f'Phenotype {ph} did not complete and validate')

if __name__=='__main__':main()
