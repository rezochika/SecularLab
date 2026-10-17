# Local reproduction and reporting

## 1. Build and standard regression suite

Enable Release in CLion, or from the project root:

```bash
cmake --preset release
cmake --build --preset release --parallel
ctest --preset release --output-on-failure
```

For CLion's bundled Windows compiler, follow [Windows setup](WINDOWS_SETUP.md) and substitute
`windows-clion` for `release`. On MSYS2 UCRT64, substitute `windows-ucrt64`.
Python 3 enables the extra CLI/data tests; the native C++ tests always build with BUILD_TESTING.
Use the release's original matrix files and run the integrity script first:

```bash
python tools/check_matrices.py
```

## 2. Run the full external problem

The existing shared run configuration has all arguments. Or, on Windows:

```powershell
.\build\windows-clion\secular.exe solve --matrix matrices/oilrig/bcsstk02.mtx --mass matrices/oilrig/bcsstm02.mtx --digits 50 --out results/oilrig50-local
.\build\windows-clion\secular.exe verify --matrix matrices/oilrig/bcsstk02.mtx --mass matrices/oilrig/bcsstm02.mtx --cert results/oilrig50-local/eigenvalues.tsv
```

Expect CERTIFIED with 66 roots and a separate VERIFIED replay. Do not assert success merely
because the process exits or approximate roots look familiar. On an unexpected result, retain
stdout, stderr and the entire new result directory; do not change source and overwrite that run.

The exact encoded decimal input, coefficient polynomial and certified output can be compared
across machines. Wall times and Graeffe precision/iteration behavior may vary with backend and
compiler versions. Compare the entire `eigenvalues.tsv` or its numerical fields, not just the
shortened values in this documentation.

## 3. Optional different-formulation numerical comparisons

Install NumPy/SciPy or mpmath in a separate Python environment, not in the native build:

```bash
python tools/reference_check.py --run results/oilrig50-local --output results/reference-scipy-local.json
python tools/reference_high_precision.py --run results/oilrig50-local --output results/reference-mpmath-local.json
```

The SciPy script is a binary64 sanity comparison. The mpmath script uses 100 decimal digits and
is numerically independent of Graeffe, but is not a rigorous interval oracle. Neither can
replace the exact inertia certificate. Neither feeds roots back into the candidate algorithm.

## 4. What to return

Return the new `result.json`, `eigenvalues.tsv`, build and CTest logs, and the compiler,
GMP/Boost versions, operating system and hardware. Include the reference JSONs when run.
A ZIP of the new result directory is sufficient for the numerical records. For failures,
include both stdout and stderr and the exact command line.

No new maximum degree or additional application model is required to reproduce this delivery.
Extending to indefinite matrices, repeated roots, complex spectra, uncertain interval input or
very large sparse matrices would be a new algorithmic task and needs its own tests.
