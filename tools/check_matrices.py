#!/usr/bin/env python3
"""Offline byte-identity test against the source Git blobs of the external matrices."""
import hashlib
from pathlib import Path

base = Path(__file__).resolve().parents[1] / 'matrices' / 'oilrig'
expected = {
    'bcsstk02.mtx': '0831271feb8da04c7792b34985efb00d68a21192',
    'bcsstm02.mtx': 'aafd612dbe95f77c90e3624de925003400801ccf',
}
for name, wanted in expected.items():
    data = (base / name).read_bytes()
    actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
    if actual != wanted:
        raise SystemExit(f'FAIL: {name}: source Git blob does not match: {actual}')
    print(f'PASS {name} source_blob={actual} sha256={hashlib.sha256(data).hexdigest()}')
