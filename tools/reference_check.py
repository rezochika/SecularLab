#!/usr/bin/env python3
"""Optional SciPy binary64 comparison. This is NOT a high-precision certificate.
Usage: python tools/reference_check.py --run verification/oilrig50 --output reference.json
Requires numpy and scipy; neither is needed by the native solver or its exact replay.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import platform
import time


def read_rational_matrix(path):
    with Path(path).open(encoding='utf-8') as stream:
        header = stream.readline().split()
        if header[:4] != ['%%MatrixMarket', 'matrix', 'coordinate', 'real']:
            raise ValueError('Expected real coordinate Matrix Market file')
        data = (line.split() for line in stream if line.strip() and not line.lstrip().startswith('%'))
        n, cols, nnz = map(int, next(data))
        if n != cols or not 1 <= n <= 256:
            raise ValueError('Invalid dimension')
        matrix = [[Fraction(0) for _ in range(n)] for _ in range(n)]
        for _ in range(nnz):
            i, j, value = next(data)
            i, j = int(i)-1, int(j)-1
            q = Fraction(value.replace('D', 'e').replace('d', 'e'))
            matrix[i][j] += q
            if header[4] == 'symmetric' and i != j:
                matrix[j][i] += q
        return matrix


def main():
    import numpy as np
    import scipy
    import scipy.linalg
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    k = np.array(read_rational_matrix(args.run/'input_K.mtx'), dtype=float)
    m = np.array(read_rational_matrix(args.run/'input_M.mtx'), dtype=float)
    start = time.perf_counter()
    values, vectors = scipy.linalg.eigh(k, m)
    seconds = time.perf_counter() - start
    values, vectors = values[::-1], vectors[:, ::-1]
    lines = (args.run/'eigenvalues.tsv').read_text().splitlines()
    printed = [line.split('\t')[1] for line in lines if line and line[0].isdigit()]
    if len(printed) != len(values):
        raise ValueError('Wrong number of eigenvalue rows')
    native = np.array([float(x) for x in printed])
    relative = np.abs(native-values) / np.maximum(1.0, np.abs(native))
    knorm, mnorm = np.linalg.norm(k, 2), np.linalg.norm(m, 2)
    residuals = [np.linalg.norm(k@vectors[:,i]-values[i]*(m@vectors[:,i])) /
                 ((knorm+abs(values[i])*mnorm)*np.linalg.norm(vectors[:,i]))
                 for i in range(len(values))]
    passed = bool(np.max(relative) < 1e-8)
    report = {'scope': 'Independent SciPy binary64 comparison only; not a certificate of 50 digits',
              'passed': passed, 'eigenvalue_count': len(values),
              'max_scaled_value_difference': float(np.max(relative)),
              'max_relative_reference_residual': float(max(residuals)),
              'scipy_eigh_kernel_seconds': seconds,
              'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
              'reference_values_descending_binary64': values.tolist()}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'reference_values_descending_binary64'}, indent=2))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
