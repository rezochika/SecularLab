#!/usr/bin/env python3
"""Independent NUMERICAL mpmath comparison, not an interval proof.
Requires mpmath. It is never used by the native solver or by its exact verifier.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
from reference_check import read_rational_matrix


def main():
    import mpmath as mp
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = 100
    k = read_rational_matrix(args.run/'input_K.mtx')
    m = read_rational_matrix(args.run/'input_M.mtx')
    def number(value):
        q = Fraction(value)
        return mp.mpf(q.numerator) / q.denominator
    n = len(k)
    if any(m[i][j] != 0 for i in range(n) for j in range(n) if i != j):
        raise ValueError('This optional reference script requires diagonal mass')
    scale = [mp.sqrt(number(m[i][i])) for i in range(n)]
    a = mp.matrix([[number(k[i][j])/(scale[i]*scale[j]) for j in range(n)] for i in range(n)])
    start = time.perf_counter()
    eigenvalues = list(mp.eigsy(a, eigvals_only=True))[::-1]
    seconds = time.perf_counter()-start
    rows = [line.split('\t') for line in (args.run/'eigenvalues.tsv').read_text().splitlines()
            if line and line[0].isdigit()]
    if len(rows) != n:
        raise ValueError('Certificate row count differs from matrix dimension')
    text = (args.run/'eigenvalues.tsv').read_text().splitlines()
    digits = int(next(line.split('=',1)[1] for line in text if line.startswith('# decimal_places=')))
    if digits > 70:
        raise ValueError('Use at most 70 output places for this fixed 100-digit numerical comparison')
    inside, rounded, distances = [], [], []
    for x, row in zip(eigenvalues, rows):
        lo, hi = number(row[2]), number(row[3])
        inside.append(lo <= x <= hi)
        r = int(mp.floor(x*10**digits + mp.mpf('0.5')))
        rounded.append(r == int(Fraction(row[1])*10**digits))
        distances.append(abs(x-(lo+hi)/2))
    report = {'scope': 'Independent mpmath NUMERICAL comparison at 100 decimal digits, not a rigorous enclosure',
              'mpmath_version': mp.__version__, 'working_decimal_digits': 100,
              'decimal_places_compared': digits, 'count': n,
              'all_numerical_references_inside_native_intervals': all(inside),
              'all_numerical_rounding_agrees': all(rounded),
              'max_distance_from_interval_center': mp.nstr(max(distances), 25),
              'eigensolver_seconds': seconds,
              'reference_values_descending': [mp.nstr(x, 95) for x in eigenvalues]}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'reference_values_descending'}, indent=2))
    return 0 if all(inside) and all(rounded) else 1


if __name__ == '__main__':
    raise SystemExit(main())
