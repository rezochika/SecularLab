#!/usr/bin/env python3
"""Standard-library integration tests. No network and no numerical Python libraries."""
import json, pathlib, subprocess, sys, tempfile
from fractions import Fraction


def check_dense_known_spectrum(certificate):
    # P = [[1,2,1], [0,1,-1], [0,0,1]] is invertible. The fixtures are
    # M=P^T P and K=P^T diag(1,2,4) P, so their exact spectrum is 4,2,1.
    rows=[line.split('\t') for line in certificate.read_text().splitlines()
          if line and not line.startswith('#')]
    assert rows[0]==['index','decimal','lower','upper','radius'], 'dense certificate columns'
    rows=rows[1:]
    assert len(rows)==3, 'dense certificate must contain the complete three-root spectrum'
    count=2
    previous_lower=None
    for index,(fields,root) in enumerate(zip(rows,(Fraction(4),Fraction(2),Fraction(1)))):
        assert len(fields)==5, ('dense certificate field count',index,fields)
        assert fields[0]==str(index), ('dense root index',index,fields[0])
        assert fields[1]==str(root)+'.'+'0'*50, ('dense canonical 50-place rounding',index,fields[1])
        printed,lower,upper,radius=map(Fraction,fields[1:])
        assert lower<=root<=upper, ('dense known root not enclosed',index,root,lower,upper)
        assert previous_lower is None or upper<previous_lower, ('dense bounds overlap or are out of order',index)
        assert radius>=max(abs(printed-lower),abs(printed-upper),abs(printed-root)), ('dense radius too small',index)
        previous_lower=lower
        count+=6
    return count


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
    # Serial and parallel exact proof must produce identical portable artifacts.
    for workers in (1, 4):
        out=d/('workers'+str(workers))
        run(['solve','--matrix',matrix/'spd3.mtx','--digits','30','--workers',workers,'--certification','auto','--out',out])
        result=json.loads((out/'result.json').read_text())
        assert result['certification_method']=='polynomial' and result['polynomial_sign_evaluations']>0
        assert result['coefficient_oracle']=='modular' and result['oracle_moduli']>=1
        # Each phase clamps the requested workers to its own jobs: moduli, or rows of a stage.
        assert result['oracle_workers']==min(workers,result['oracle_moduli']) and result['recurrence_workers']==min(workers,3)
        assert result['inertia_evaluations']==2 # mandatory K/M SPD guards only
        assert result['times_seconds']['initialization']>=0 and result['times_seconds']['squarefree']>=0
        assert result['graeffe_unbounded_exponents'] is False # 64-bit binary exponents sufficed
        assert result['certification_attempts']==1 # the first stable candidates certify
        checks+=7
        run(['coefficients','--matrix',matrix/'spd3.mtx','--workers',workers,'--out',d/('coeff'+str(workers))])
    for artifact in ('coefficients.tsv','eigenvalues.tsv'):
        assert (d/'workers1'/artifact).read_bytes()==(d/'workers4'/artifact).read_bytes(); checks+=1
    assert (d/'coeff1/coefficients.tsv').read_bytes()==(d/'coeff4/coefficients.tsv').read_bytes(); checks+=1
    # The explicit n+1-point Bareiss oracle proves the same identity with unchanged artifacts.
    for workers in (1, 4):
        out=d/('bareiss'+str(workers))
        run(['solve','--matrix',matrix/'spd3.mtx','--digits','30','--workers',workers,'--oracle','bareiss','--out',out])
        result=json.loads((out/'result.json').read_text())
        assert result['coefficient_oracle']=='bareiss' and result['oracle_moduli']==0 and result['oracle_workers']==workers
        assert result['coefficient_identity_verified'] and result['certification_method']=='polynomial'
        for artifact in ('coefficients.tsv','eigenvalues.tsv'):
            assert (out/artifact).read_bytes()==(d/'workers1'/artifact).read_bytes()
        checks+=4
    run(['coefficients','--matrix',matrix/'spd3.mtx','--oracle','bareiss','--out',d/'coeff_bareiss'])
    assert (d/'coeff_bareiss/coefficients.tsv').read_bytes()==(d/'coeff1/coefficients.tsv').read_bytes(); checks+=1
    # Chinese remaindering of the oracle's residues replaces the recurrence, alone or beside it
    # (N-version), for diagonal and dense mass, with or without content removal: identical artifacts.
    for name,pencil in (('spd3',['--matrix',matrix/'spd3.mtx']),
                        ('dense',['--matrix',matrix/'dense_stiffness3.mtx','--mass',matrix/'dense_mass3.mtx'])):
        reference=None
        for method in ('recurrence','crt','both'):
            for content in ([],['--remove-content']):
                out=d/('method_'+name+'_'+method+('_content' if content else ''))
                run(['solve',*pencil,'--digits','30','--workers','2','--coefficients',method,*content,'--out',out])
                result=json.loads((out/'result.json').read_text())
                assert result['coefficient_method']==method and result['content_removed']==bool(content)
                assert result['certification_method']=='polynomial' and result['coefficient_identity_verified']
                assert result['coefficient_oracle']==('modular' if method=='recurrence' else 'crt')
                assert (result['coefficient_moduli']>=1)==(method!='recurrence')==(result['crt_workers']>=1)
                assert 1<=result['primitive_coefficient_bits_max']<=result['coefficient_bits_max']
                assert result['times_seconds']['crt']>=0 and (method=='crt')==(result['recurrence_workers']==0)
                artifacts=[(out/artifact).read_bytes() for artifact in ('coefficients.tsv','eigenvalues.tsv')]
                reference=reference or artifacts
                assert artifacts==reference; checks+=7
                replay=run(['verify',*pencil,'--cert',out/'eigenvalues.tsv','--coefficients',method,*content])
                assert replay.stdout.split()[0]=='VERIFIED' and 'certification_method=polynomial' in replay.stdout
        if name=='spd3':
            assert reference==[(d/'workers1'/artifact).read_bytes() for artifact in ('coefficients.tsv','eigenvalues.tsv')]; checks+=1
    for method,oracle,expected in (('crt',[],'crt'),('both',[],'crt'),('crt',['--oracle','bareiss'],'bareiss')):
        out=d/('coeff_'+method+'_'+expected)
        line=run(['coefficients','--matrix',matrix/'spd3.mtx','--coefficients',method,*oracle,'--out',out]).stdout
        assert 'determinant_identity_verified=yes' in line and 'coefficient_oracle='+expected in line
        assert 'coefficient_method='+method in line
        assert (out/'coefficients.tsv').read_bytes()==(d/'coeff1/coefficients.tsv').read_bytes(); checks+=3
    # det(E)=0 is proven from the skipped moduli; the request fails as the recurrence does.
    for method in ('recurrence','crt'):
        failed=run(['coefficients','--matrix',matrix/'spd2.mtx','--mass',matrix/'singular_mass2.mtx',
                    '--coefficients',method,'--out',d/('singular_coefficients_'+method)],2)
        assert 'singular mass/leading matrix' in failed.stderr; checks+=1
    for workers in (1, 4):
        run(['coefficients','--matrix',matrix/'galerkin2_K0.mtx','--mass',matrix/'galerkin2_M.mtx','--affine',matrix/'galerkin2_K1.mtx','--workers',workers,'--out',d/('affine'+str(workers))])
    assert (d/'affine1/coefficients_affine.tsv').read_bytes()==(d/'affine4/coefficients_affine.tsv').read_bytes(); checks+=1
    run(['coefficients','--matrix',matrix/'galerkin2_K0.mtx','--mass',matrix/'galerkin2_M.mtx','--affine',matrix/'galerkin2_K1.mtx','--oracle','bareiss','--out',d/'affine_bareiss'])
    assert (d/'affine_bareiss/coefficients_affine.tsv').read_bytes()==(d/'affine1/coefficients_affine.tsv').read_bytes(); checks+=1
    cert=d/'workers1/eigenvalues.tsv'
    for method,expected in (('auto','polynomial'),('inertia','exact_inertia')):
        out=d/('verify_'+method)
        out.mkdir() # existing empty directory is also accepted
        replay=run(['verify','--matrix',matrix/'spd3.mtx','--cert',cert,'--workers','4','--certification',method,'--out',out])
        result=json.loads((out/'result.json').read_text())
        assert replay.stdout.split()[0]=='VERIFIED' and expected in replay.stdout and 'total_seconds=' in replay.stdout
        assert result['status']=='VERIFIED' and result['certified'] and result['certification_method']==expected
        assert result['times_seconds']['total']>=0; checks+=3
        saved=(out/'result.json').read_bytes()
        run(['verify','--matrix',matrix/'spd3.mtx','--cert',cert,'--out',out],2)
        assert (out/'result.json').read_bytes()==saved; checks+=1
    # Replay must derive its polynomial from matrices, irrespective of adjacent saved files.
    (d/'workers1/coefficients.tsv').write_text('corrupt external coefficients\n')
    (d/'workers1/result.json').write_text('{"status":"FAILED","certified":false}\n')
    run(['verify','--matrix',matrix/'spd3.mtx','--cert',cert,'--workers','0'])
    # Both exact proof methods reject wrong matrices and altered certificate rows.
    rows=cert.read_text().splitlines()
    data_index=next(i for i,row in enumerate(rows) if row and not row.startswith('#') and row.split('\t')[0].isdigit())
    fields=rows[data_index].split('\t'); fields[1]='9.'+'0'*30; rows[data_index]='\t'.join(fields)
    (d/'altered.tsv').write_text('\n'.join(rows)+'\n')
    for method in ('auto','inertia'):
        run(['verify','--matrix',matrix/'spd3.mtx','--cert',d/'altered.tsv','--certification',method],1)
        run(['verify','--matrix',matrix/'dense_stiffness3.mtx','--cert',cert,'--certification',method],1)
        failed=d/('failed_'+method)
        replay=run(['verify','--matrix',matrix/'repeated3.mtx','--cert',cert,'--certification',method,'--out',failed],1)
        result=json.loads((failed/'result.json').read_text())
        assert replay.stdout.split()[0]=='FAILED' and not result['certified'] and result['status']=='FAILED'; checks+=1
    (d/'malformed.tsv').write_text('invalid certificate\n')
    run(['verify','--matrix',matrix/'spd3.mtx','--cert',d/'malformed.tsv'],2)
    for mode in ('solve','verify','coefficients'):
        base=[mode,'--matrix',matrix/'spd3.mtx']+(['--cert',cert] if mode=='verify' else [])
        for invalid in ('-1','65','3/2','abc'):
            out=d/(mode+'_invalid_workers_'+invalid.replace('/','_'))
            run([*base,'--workers',invalid,'--out',out],2)
            assert not out.exists(); checks+=1
        if mode!='coefficients':
            out=d/(mode+'_invalid_method')
            run([*base,'--certification','approximate','--out',out],2)
            assert not out.exists(); checks+=1
        out=d/(mode+'_invalid_oracle')
        run([*base,'--oracle','approximate','--out',out],2)
        assert not out.exists(); checks+=1
        out=d/(mode+'_invalid_coefficients')
        run([*base,'--coefficients','approximate','--out',out],2)
        assert not out.exists(); checks+=1
    # An oracle choice without a coefficient identity check is contradictory.
    for args in (['solve','--matrix',matrix/'spd3.mtx','--no-coefficient-check','--oracle','bareiss'],
                 ['coefficients','--matrix',matrix/'spd3.mtx','--no-coefficient-check','--oracle','auto'],
                 ['verify','--matrix',matrix/'spd3.mtx','--cert',cert,'--certification','inertia','--oracle','bareiss']):
        out=d/('oracle_without_check_'+args[0])
        run([*args,'--out',out],2)
        assert not out.exists(); checks+=1
    # Chinese remaindering is itself the identity proof, content removal acts on polynomial
    # stages only, and formal affine coefficients have no Chinese-remainder route.
    for index,args in enumerate((['solve','--matrix',matrix/'spd3.mtx','--no-coefficient-check','--coefficients','crt'],
                 ['coefficients','--matrix',matrix/'spd3.mtx','--no-coefficient-check','--coefficients','both'],
                 ['verify','--matrix',matrix/'spd3.mtx','--cert',cert,'--certification','inertia','--coefficients','crt'],
                 ['verify','--matrix',matrix/'spd3.mtx','--cert',cert,'--certification','inertia','--remove-content'],
                 ['coefficients','--matrix',matrix/'spd3.mtx','--remove-content'],
                 ['coefficients','--matrix',matrix/'galerkin2_K0.mtx','--mass',matrix/'galerkin2_M.mtx',
                  '--affine',matrix/'galerkin2_K1.mtx','--coefficients','crt'])):
        out=d/('contradictory_'+str(index))
        run([*args,'--out',out],2)
        assert not out.exists(); checks+=1
    run(['verify','--matrix',matrix/'spd3.mtx','--cert',cert,'--no-coefficient-check','--out',d/'weaken_verify'],2)
    assert not (d/'weaken_verify').exists(); checks+=1
    run(['solve','--matrix',matrix/'spd3.mtx','--no-coefficient-check','--out',d/'unchecked'])
    assert json.loads((d/'unchecked/result.json').read_text())['certification_method']=='exact_inertia'; checks+=1
    run(['solve','--matrix',matrix/'spd3.mtx','--workers','64','--certification','inertia','--out',d/'clamped'])
    result=json.loads((d/'clamped/result.json').read_text())
    assert result['oracle_workers']==min(64,result['oracle_moduli']) and result['recurrence_workers']==3
    assert result['certification_method']=='exact_inertia'; checks+=2
    run(['solve','--matrix',matrix/'spd3.mtx','--workers','64','--oracle','bareiss','--out',d/'clamped_bareiss'])
    result=json.loads((d/'clamped_bareiss/result.json').read_text())
    assert result['oracle_workers']==4 and result['coefficient_oracle']=='bareiss'; checks+=1 # n+1 points
    for name in ('spd2','spd3','spring10','spring20','cluster3'):
        out=d/name
        run(['solve','--matrix',matrix/(name+'.mtx'),'--digits','30','--out',out])
        result=json.loads((out/'result.json').read_text())
        assert result['certified'] and result['coefficient_identity_verified']; checks+=1
        run(['verify','--matrix',out/'input_K.mtx','--mass',out/'input_M.mtx','--cert',out/'eigenvalues.tsv'])
    # Proposal widths see candidates only through their nearest powers of ten, so every
    # working precision that certifies gives byte-identical intervals.
    stable=set()
    for precision in ('64','256','1024'):
        out=d/('spring20_work'+precision)
        run(['solve','--matrix',matrix/'spring20.mtx','--digits','1','--precision',precision,'--max-precision',precision,'--out',out])
        stable.add((out/'eigenvalues.tsv').read_bytes())
    assert len(stable)==1, 'spring20 certificate endpoints changed with working precision'; checks+=1
    run(['solve','--matrix',matrix/'dense_stiffness3.mtx','--mass',matrix/'dense_mass3.mtx','--digits','50','--out',d/'dense'])
    checks+=check_dense_known_spectrum(d/'dense/eigenvalues.tsv')
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
