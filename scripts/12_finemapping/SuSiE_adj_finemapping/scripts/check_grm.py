from pathlib import Path
import struct, json
prefix=Path('grm/spg.20k')
ids=Path(str(prefix)+'.grm.iid').read_text().split()
fam=[x.split()[1] for x in Path('genotype/spg.20k.fam').read_text().splitlines()]
assert ids==fam and len(ids)==20000 and len(set(ids))==20000
with open(str(prefix)+'.grm.bin','rb') as f:
    n,m=struct.unpack('<if',f.read(8))
    assert n==20000 and 0<m<=135358
    # MPH/LD Adjuster store lower triangle column by column, diagonal included.
    size=8+4*n*(n+1)//2
    assert Path(str(prefix)+'.grm.bin').stat().st_size==size
    diagonal=[]
    for j in (0,1,100,n-1):
        offset=j*n-j*(j-1)//2
        f.seek(8+4*offset)
        v=struct.unpack('<f',f.read(4))[0]/m
        assert 0<v<10
        diagonal.append(v)
Path('qc/grm.json').write_text(json.dumps(dict(n=n,markers=m,bytes=size,
    format='MPH int32 n; float32 m; float32 column-major lower triangle divided by m on reading',
    IID_order_matches=True,diagonal_checks=diagonal),indent=2)+'\n')
Path('status/grm.PASS').write_text('PASS\n')
