#!/usr/bin/env python3
"""Optional slow CTest: full matrix, all 66 roots at 50 places, then certificate replay."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

exe, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
with tempfile.TemporaryDirectory(prefix='secular-oilrig-') as temp:
    out = Path(temp) / 'results'
    stiffness = root / 'matrices/oilrig/bcsstk02.mtx'
    mass = root / 'matrices/oilrig/bcsstm02.mtx'
    subprocess.run([str(exe), 'solve', '--matrix', str(stiffness), '--mass', str(mass),
                    '--digits', '50', '--out', str(out)], check=True, timeout=3600)
    data = json.loads((out / 'result.json').read_text())
    assert data['certified'] and data['certified_eigenvalues'] == 66
    assert data['coefficient_identity_verified']
    subprocess.run([str(exe), 'verify', '--matrix', str(stiffness), '--mass', str(mass),
                    '--cert', str(out / 'eigenvalues.tsv')], check=True, timeout=3600)
print('PASS oilrig 66/66 at 50 decimal places and exact certificate replay')
