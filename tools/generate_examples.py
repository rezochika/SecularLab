#!/usr/bin/env python3
"""Regenerate only the explicitly synthetic examples, using exact Fractions.
The external oil-rig matrices are never modified by this script.
"""
from fractions import Fraction as F
from pathlib import Path
import json
ROOT = Path(__file__).resolve().parents[1] / 'matrices' / 'examples'
ROOT.mkdir(parents=True, exist_ok=True)
def write(name, a, description):
    n = len(a); symmetric = all(a[i][j] == a[j][i] for i in range(n) for j in range(n))
    data = [(i+1,j+1,F(a[i][j])) for i in range(n) for j in range(i+1 if symmetric else n) if a[i][j]]
    text = ['%%MatrixMarket matrix coordinate real ' + ('symmetric' if symmetric else 'general'),
            '% Synthetic regression case: ' + description,
            '% Fraction tokens are a documented SecularLab extension to Matrix Market.', f'{n} {n} {len(data)}']
    text += [f'{i} {j} {x}' for i,j,x in data]
    (ROOT/name).write_text('\n'.join(text)+'\n', encoding='utf-8')
def diag(xs): return [[F(x) if i==j else F(0) for j in range(len(xs))] for i,x in enumerate(xs)]
def add(a,b): return [(a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(max(len(a),len(b)))]
def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def scale(a,s): return [s*x for x in a]
def integ(a): return sum((x*F(2,i+1) for i,x in enumerate(a) if i%2==0),F(0))
def derivative(a): return [i*a[i] for i in range(1,len(a))]
def galerkin(n):
    p=[[F(1)],[F(0),F(1)]]
    for i in range(1,n+1): p.append(scale(add(scale([F(0)]+p[-1],2*i+1),scale(p[-2],-i)),F(1,i+1)))
    basis=[add(p[i+1],scale(p[i-1],-1)) for i in range(1,n+1)]
    h=[[integ(mul(a,b)) for b in basis] for a in basis]
    d=[[integ(mul(derivative(a),derivative(b))) for b in basis] for a in basis]
    minus_x=[[-integ([F(0)]+mul(a,b)) for b in basis] for a in basis]
    for suffix,a in [('M',h),('K0',d),('K1',minus_x)]:write(f'galerkin{n}_{suffix}.mtx',a,'OPM Legendre basis, K(t)=K0+t*K1, M=H')
write('spd2.mtx',[[2,1],[1,2]],'eigenvalues 3,1')
write('spd3.mtx',[[4,1,1],[1,3,0],[1,0,2]],'dense SPD matrix, irrational simple eigenvalues')
for n in (10,20):write(f'spring{n}.mtx',[[2 if i==j else -1 if abs(i-j)==1 else 0 for j in range(n)] for i in range(n)],'fixed-end unit spring chain, M=I')
write('ties2.mtx',diag([F(1,4),F(9,4)]),'positive exact decimal ties 0.25,2.25')
write('cluster3.mtx',diag([1,1+F(1,10**12),1+F(2,10**12)]),'simple eigenvalues separated by 1e-12')
write('repeated3.mtx',diag([1,1,3]),'expected REPEATED_SPECTRUM')
write('indefinite2.mtx',diag([-1,2]),'coefficient mode allowed, positive-root solver rejects')
write('nonsymmetric2.mtx',[[1,2],[0,3]],'coefficient mode allowed, certified symmetric solver rejects')
write('singular_mass2.mtx',diag([0,1]),'singular leading matrix must be rejected')
p=[[F(1),F(2),F(1)],[F(0),F(1),F(-1)],[F(0),F(0),F(1)]]
m=[[sum(p[k][i]*p[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
k=[[sum(F([1,2,4][r])*p[r][i]*p[r][j] for r in range(3)) for j in range(3)] for i in range(3)]
write('dense_mass3.mtx',m,'M=P^T P, exact generalized eigenvalues 4,2,1')
write('dense_stiffness3.mtx',k,'K=P^T diag(1,2,4) P')
galerkin(2);galerkin(6)
print(f'Wrote {len(list(ROOT.glob("*.mtx")))} exact synthetic matrix files')
