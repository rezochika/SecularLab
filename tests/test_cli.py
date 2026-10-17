#!/usr/bin/env python3
"""Standard-library integration tests. No network and no numerical Python libraries."""
import json, pathlib, subprocess, sys, tempfile
exe=pathlib.Path(sys.argv[1]).resolve();root=pathlib.Path(sys.argv[2]).resolve();matrix=root/'matrices/examples'
checks=0
with tempfile.TemporaryDirectory(prefix='secular-tests-') as td:
    d=pathlib.Path(td)
    def run(args, code=0):
        global checks
        r=subprocess.run([str(exe),*map(str,args)],text=True,capture_output=True,timeout=180)
        assert r.returncode==code,(r.returncode,code,r.stdout,r.stderr,args)
        checks+=1
        return r
    for name in ('spd2','spd3','spring10','spring20','cluster3'):
        out=d/name
        run(['solve','--matrix',matrix/(name+'.mtx'),'--digits','30','--out',out])
        result=json.loads((out/'result.json').read_text())
        assert result['certified'] and result['coefficient_identity_verified']; checks+=1
        run(['verify','--matrix',out/'input_K.mtx','--mass',out/'input_M.mtx','--cert',out/'eigenvalues.tsv'])
    run(['solve','--matrix',matrix/'dense_stiffness3.mtx','--mass',matrix/'dense_mass3.mtx','--digits','50','--out',d/'dense'])
    run(['verify','--matrix',d/'dense/input_K.mtx','--mass',d/'dense/input_M.mtx','--cert',d/'dense/eigenvalues.tsv'])
    run(['solve','--matrix',matrix/'galerkin6_K0.mtx','--mass',matrix/'galerkin6_M.mtx','--affine',matrix/'galerkin6_K1.mtx','--legendre-degree','128','--argument','1/3','--out',d/'galerkin6'])
    run(['verify','--matrix',d/'galerkin6/input_K.mtx','--mass',d/'galerkin6/input_M.mtx','--cert',d/'galerkin6/eigenvalues.tsv'])
    run(['coefficients','--matrix',matrix/'galerkin2_K0.mtx','--mass',matrix/'galerkin2_M.mtx','--affine',matrix/'galerkin2_K1.mtx','--out',d/'affine'])
    run(['coefficients','--matrix',matrix/'nonsymmetric2.mtx','--out',d/'nonsym'])
    run(['solve','--matrix',matrix/'ties2.mtx','--digits','1','--out',d/'ties'])
    run(['verify','--matrix',matrix/'ties2.mtx','--cert',d/'ties/eigenvalues.tsv'])
    run(['solve','--matrix',matrix/'repeated3.mtx','--out',d/'repeated'],1)
    run(['solve','--matrix',matrix/'indefinite2.mtx','--out',d/'indef'],2)
    run(['solve','--matrix',matrix/'nonsymmetric2.mtx','--out',d/'reject_nonsym'],2)
    run(['solve','--matrix',matrix/'spd2.mtx','--mass',matrix/'singular_mass2.mtx','--out',d/'singular'],2)
    run(['solve','--matrix',matrix/'spd2.mtx','--precision','64','--max-precision','64','--max-squarings','1','--out',d/'limited'],1)
    run(['solve','--matrix',matrix/'spd2.mtx','--out',d/'spd2'],2) # never silently overwrite
    lines=(d/'spd2/eigenvalues.tsv').read_text().splitlines();cols=lines[3].split('\t');cols[1]='9.'+'0'*30;lines[3]='\t'.join(cols)
    (d/'tampered.tsv').write_text('\n'.join(lines)+'\n')
    run(['verify','--matrix',matrix/'spd2.mtx','--cert',d/'tampered.tsv'],1)
    run(['not-a-command'],2)
    run(['solve','--matrix',matrix/'spd2.mtx','--digits','-1','--out',d/'invalid_options'],2)
    assert not (d/'invalid_options').exists(), 'invalid options created output files'; checks+=1
    run(['solve','--matrix',matrix/'spd2.mtx','--precision','32','--max-precision','16','--out',d/'invalid_precision'],2)
    assert not (d/'invalid_precision').exists(); checks+=1
print('PASS',checks,'CLI/integration checks')
