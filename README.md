# SecularLab

**C++20 / CMake project for CLion: exact secular coefficients, Graeffe root recovery, and exact generalized-eigenvalue certificates.**

This is a separate reference implementation of the secular-equation workflow discussed in *Orthogonal-Polynomial Methods for Boundary-Value Problems and Secular Equations* by Tamaz Vashakmadze and Revaz Chikashua. It is **not** the existing JacobiRoots executable and does not compute Gaussian quadrature weights.

Given exact rational matrices `K` and `M`, it computes the zeros of

```text
                         det(lambda * M - K) = 0.
```

The full published **BCSSTK02 / BCSSTM02 66 x 66 oil-rig pair is included**, along with 18 small matrix files, automated tests, and recorded results. No data download occurs during configuration or execution.

## Start in CLion

### Windows: CLion bundled MinGW

The current Windows checkout uses CLion's bundled MinGW compiler with project-local GMP C++ and Boost dependencies. See [Windows setup](docs/WINDOWS_SETUP.md) for the reproducible setup command and the `windows-clion` Release and `windows-clion-debug` profiles. Both native configurations pass the five standard test groups. CLion's build action and the small example were also exercised successfully.

### Windows: alternative MSYS2 UCRT64 toolchain

Install MSYS2, open its **UCRT64** terminal, update packages as directed by MSYS2, and install:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gdb mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja mingw-w64-ucrt-x86_64-gmp mingw-w64-ucrt-x86_64-boost mingw-w64-ucrt-x86_64-python
```

1. Extract the ZIP. In CLion, **File > Open** the `SecularLab` directory containing `CMakeLists.txt`.
2. In **Settings > Build, Execution, Deployment > Toolchains**, add a **MinGW** toolchain named **MSYS2 UCRT64**. Select `C:\msys64\ucrt64` as the installation and its `g++.exe` and `gdb.exe`. Do not mix these libraries with a different MinGW/MSVC toolchain.
3. In **Settings > Build, Execution, Deployment > CMake**, enable the imported **windows-ucrt64** profile. CLion can import presets but leave newly imported profiles disabled. Reload CMake and build.
4. Select **01 - Small SPD example** first, then **02 - Oil rig 66 - 50 decimal places**. Shared configurations are in `.run/`.
5. Run **03 - Verify oil rig saved certificates** to recheck the resulting exact intervals without rerunning Graeffe.

The preset assumes MSYS2 is installed at `C:\msys64`. For another location, edit the compiler, Ninja, prefix and PATH entries in that preset or create a local `CMakeUserPresets.json` override.

**DLLs:** `C:\msys64\ucrt64\bin` must also be on the *run/debug process* PATH. Add it to your user PATH and restart CLion, or set `Path=C:\msys64\ucrt64\bin;$Path$` under **Run > Edit Configurations > Environment variables**. A CMake configure environment alone is not always the run environment. No DLLs are bundled here.

If shared configurations do not appear, create a **CMake Application** configuration: target `secular`, working directory `$PROJECT_DIR$`, and copy the corresponding command arguments below.

### Linux

On Debian/Ubuntu, install the development packages:

```bash
sudo apt-get install build-essential cmake ninja-build libgmp-dev libboost-dev python3
cmake --preset release
cmake --build --preset release --parallel
ctest --preset release
```

Or open the folder in CLion and enable the **release** profile using your GCC or Clang toolchain.

### macOS

```bash
brew install cmake ninja gmp boost
cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$(brew --prefix)"
cmake --build build/release --parallel
ctest --test-dir build/release --output-on-failure
```

In CLion use the same Homebrew prefix in a local CMake profile. Linux and Windows builds have been tested. The macOS and alternative MSYS2 configurations have not been executed here. The original delivery records remain under `verification/`; current repair validation is recorded separately.

## Commands

Run from the project root. The Linux/macOS executable is `build/release/secular`; the bundled CLion Windows preset builds `build/windows-clion/secular.exe` (MSYS2: `build/windows-ucrt64/secular.exe`).

### Complete oil-rig test

```bash
build/release/secular solve --matrix matrices/oilrig/bcsstk02.mtx --mass matrices/oilrig/bcsstm02.mtx --digits 50 --out results/oilrig50
```

### Recheck the saved certificate

```bash
build/release/secular verify --matrix matrices/oilrig/bcsstk02.mtx --mass matrices/oilrig/bcsstm02.mtx --cert results/oilrig50/eigenvalues.tsv --workers 8 --out results/oilrig50-replay
```

In its default `auto` mode, the `verify` subcommand derives its characteristic polynomial from the supplied matrices. It independently checks the complete determinant identity and squarefreeness before trying the polynomial proof of all saved intervals, with exact inertia as a fallback. Inertia-only replay checks the original pencil directly without constructing that polynomial. Saved JSON status, floating-point candidates, and coefficient files are never trusted. The optional `--out` saves a detailed `result.json` in a new or empty directory; replay without `--out` still works. Output begins with `VERIFIED` or `FAILED` and reports the actual proof method and elapsed time.

For `solve` and `verify`, `--certification auto` is the default. Use `--certification inertia` to request the original exact-inertia proof directly. `--no-coefficient-check` is available for solve and coefficients; a solve with this option always uses inertia certification. Verify rejects this option, so it cannot weaken the independent proof.

All three commands accept `--workers 0..64` for independent determinant-oracle checks, including the affine coefficient grid. The default `0` chooses a hardware-based count capped at 8. Explicit counts from 1 through 64 are capped by the number of oracle jobs; `--workers 1` requests serial checking. Each successful identity check still evaluates all `n+1` distinct points exactly.

### Small example / debugger

```bash
build/release/secular solve --matrix matrices/examples/spd3.mtx --digits 50 --out results/spd3
```

Omitted `--mass` means the identity matrix. Use a Debug profile for stepping through the recurrence, but use Release for timing.

### Paper's affine Galerkin example

```bash
build/release/secular solve --matrix matrices/examples/galerkin6_K0.mtx --mass matrices/examples/galerkin6_M.mtx --affine matrices/examples/galerkin6_K1.mtx --legendre-degree 128 --argument 1/3 --digits 50 --out results/galerkin6
```

Here `K(t)=K0+t*K1=D-t*X` and `t=P_128(1/3)` is evaluated **exactly as a rational**. The solve command specializes that parameter before constructing coefficients. This is not the paper's interval-evaluation experiment at an uncertain special-function argument.

### Exact coefficients as polynomials in a formal parameter

```bash
build/release/secular coefficients --matrix matrices/examples/galerkin2_K0.mtx --mass matrices/examples/galerkin2_M.mtx --affine matrices/examples/galerkin2_K1.mtx --out results/affine2
```

This uses the coefficientwise affine adjugate recurrence, retaining formal `t`. With the default checks it verifies the identity on an exact two-dimensional evaluation grid. Output is `coefficients_affine.tsv`. For a fixed rational specialization, add `--parameter 1/3`.

### Higher precision or a deliberately limited test

```bash
build/release/secular solve --matrix matrices/examples/cluster3.mtx --digits 100 --precision 512 --max-precision 2048 --max-squarings 128 --out results/cluster100
```

`--digits` means **decimal places after the decimal point**, not significant digits. `--precision` and `--max-precision` specify decimal working precision for the untrusted candidate calculation. More working precision or more squarings does not establish correctness until the exact certificate passes.

**Every output directory must be new or empty.** The program refuses to overwrite an existing run. Change `--out` before repeating a shared CLion configuration, or deliberately remove/archive your own previous run first. The bundled delivery records are under `verification/`, not `results/`, so a first local run does not collide with them.

## What is implemented

- Exact parsing of integers, finite decimals, scientific `e/E/d/D` notation, and rational `p/q` tokens.
- Denominator clearing, Faddeev-LeVerrier/adjugate recurrence, checked integer divisions and `S_n=0`.
- Direct Bareiss determinant checks at `n+1` distinct points, proving the whole coefficient polynomial identity.
- Parallel determinant checks with private exact-arithmetic scratch space and configurable worker count.
- A sufficient modular squarefree proof over fixed verified primes. A degree drop or nontrivial modular gcd is inconclusive and leads to another prime or the existing exact gcd fallback; it never alone establishes repeated roots.
- Formal affine coefficient polynomials for `det(zE-F0-tF1)`.
- Signed-logarithmic Graeffe iteration, avoiding giant repeatedly squared coefficients. No external eigensolver or Newton polishing supplies its roots.
- Exact symmetric inertia with a fraction-free path and a rational block-LDL fallback for zero diagonal pivots.
- Exact rational root intervals, index verification, correct half-away-from-zero decimal rounding, and radii about the printed values.
- Exact polynomial signs via homogeneous integer evaluation at rational endpoints. Once the matrix-derived degree-`n` characteristic identity is verified, `n` disjoint brackets with strict opposite signs or exact singleton zeros exhaust all roots. This full-spectrum argument does not certify a partial spectrum. Ambiguous endpoint cases retain the exact-inertia fallback.
- Portable TSV certificates, separate replay, CSV/TSV-compatible coefficient records and JSON execution summaries.

See `docs/ALGORITHM.md` for the correspondence with the paper, actual proof obligations and limitations.

JSON summaries add `certification_method` (`polynomial`, `exact_inertia`, or `none`), the actual `oracle_workers`, `squarefree_modular_prime` (zero means no modular witness: exact fallback or a check that was not run, as in inertia-only replay), and counts of polynomial-sign and inertia evaluations (including the two mandatory positive-definiteness checks). `times_seconds.initialization` and `times_seconds.squarefree` supplement the existing phase timings and total time. Coefficient and certificate TSV formats are unchanged, as are the requested decimal places and exact half-away-from-zero rounding; `--digits 50` still requests 50 places after the decimal point.

## Supported scope and honest limitations

**Root solving:** symmetric positive definite `K` and `M`, with a simple spectrum. Repeated roots are detected exactly and reported as `REPEATED_SPECTRUM`. Indefinite/nonsymmetric problems are not silently shifted, approximated by symmetric matrices, or passed to a different eigensolver. **Coefficient construction** permits nonsymmetric inputs with any nonsingular leading `M`.

The implementation is dense and intended as a reproducible research reference, not a replacement for fast sparse eigensolvers. Independent determinant identity checks can run in parallel; coefficient construction and Graeffe recovery remain serial. Limits are order 256 (formal affine mode: 24), requested decimal places 1..1000, candidate working precision 32..20000, at most 512 Graeffe squarings, and Legendre degree at most 10000. Some permitted inputs can exhaust time or memory or fail to certify within the configured limits.

Exact arithmetic uses **GMP**. Candidates use **Boost.Multiprecision's GMP floating-point backend**, not MPFR. Their logarithms and exponentials are untrusted approximations; acceptance is exact and independent of their error estimates. This backend choice differs from the manuscript's historical/native programs. No claim of formal C++ verification is made. Because the floating-point default precision is process-global, concurrent calls to `solve` in one process are not supported.

The certificate refers to the matrix whose published entries are interpreted as exact rationals. It does **not** imply that an oil rig's measured physical parameters are accurate to 50 decimal places.

## Tests and recorded results

```bash
ctest --preset release --output-on-failure
```

The standard test groups cover exact arithmetic/certificates, additional exact algebra, clustered and widely scaled root recovery, fast exact arithmetic, both certification methods, CLI examples, and byte integrity of both external matrices. Regressions include malformed numeric signs, mismatched pencil dimensions, continued Graeffe convergence after an early failed certificate, disjoint certificates for tiny roots, serial/parallel artifact equivalence, altered certificates, matrix mismatch, and replay independent of saved coefficient/status files. Python is optional for building the solver, but required for the CLI/data tests. All Python tests in CTest use the standard library only.

Optional slow full-matrix solve plus certificate replay:

```bash
cmake -S . -B build/long -G Ninja -DCMAKE_BUILD_TYPE=Release -DSECULAR_LONG_TESTS=ON
cmake --build build/long --parallel
ctest --test-dir build/long -L slow --output-on-failure
```

The tested external problem produced **66 certified roots at 50 decimal places, after 22 Graeffe squarings using 256 working decimal digits**. Exact coefficient construction and the `n+1`-point determinant check passed. See `verification/README.md` for the original delivery and [the repair verification record](verification/repairs-20261018/REPORT.txt) for the current Windows/Linux tests and native Windows oil-rig solve/replay. These are execution observations, not a speed comparison with another solver.

## Project layout

```text
CMakeLists.txt, CMakePresets.json   Build and CLion profiles
.run/                            Shared CLion run configurations
include/secular/                 Reusable C++ interface
src/exact.cpp                    Rational input, Bareiss, adjugate, inertia
src/solver.cpp                   Graeffe candidates and exact certificates
src/main.cpp                     Command-line interface
matrices/oilrig/                  Full, unchanged external matrix files
matrices/examples/                Synthetic and paper-derived regression inputs
tests/                           C++ and Python tests
tools/                           Matrix integrity, example generator, reference check
docs/                            Algorithm, formats, references and test guidance
verification/                    Actual delivery test records and oil-rig outputs
```

## Sources and permissions

The implementation is newly written for this delivery and follows the algorithmic formulas discussed in the uploaded manuscripts; it is not a claim of a new general pencil theorem. See `docs/REFERENCES.md` and `matrices/oilrig/PROVENANCE.md`. The bundled external matrix files retain their original collection headers. Their provenance is separate from the source-code license. Third-party dependencies are not vendored.
