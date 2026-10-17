# Delivery verification record

These records come from execution of the newly written SecularLab implementation, not from
an inherited JacobiRoots or manuscript experiment. See `environment.json` for source hashes,
the tested executable hash, compiler and backend versions.

## Native execution

The final source was built in Release mode with GCC 14.2.0, GMP 6.3.0 and Boost 1.83 on Linux.
CMake was 3.31.6. Both imported native presets and an additional Debug address/undefined
sanitizer build were exercised. The final standard CTest suite passed 4/4 entries:

- 35 native core assertions.
- 599 additional exact-algebra assertions, including random integer congruences with known
  inertia, nonsymmetric pencils and exact affine specializations.
- 35 CLI/integration checks, including full small examples, decimal ties, clusters, rejected
  input classes, deliberately insufficient limits, replay, mutation and no-overwrite behavior.
- Byte-identity checks of both bundled external matrix files.

The Windows/macOS builds and CLion GUI were not executed here. Their presets/configurations
and instructions are supplied for local use, not presented as completed cross-platform tests.

## Full oil-rig problem

The **entire 66 x 66 BCSSTK02/BCSSTM02 pair** was solved, not a selected submatrix:

| Quantity | Recorded result |
|---|---:|
| Certified eigenvalues | 66/66 |
| Correctly rounded decimal places | 50 |
| Working decimal digits | 256 |
| Graeffe squarings | 22 |
| Candidate-precision attempts | 1 |
| Coefficient recurrence exact divisions | 287562 |
| Largest coefficient bit length | 6707 |
| Exact coefficient identity check | Passed at 67 determinant evaluation points |
| Coefficient construction | 3.91 s |
| Coefficient identity oracle | 16.58 s |
| Graeffe candidates | 0.29 s |
| Exact acceptance checks | 118.00 s |
| Total solver time | 150.10 s |
| Process peak resident memory | 14968 KiB (about 14.6 MiB) |

The precise output values and rational enclosures are in `oilrig50/eigenvalues.tsv`.
The exact integer polynomial is in `oilrig50/coefficients.tsv`. `input_K.mtx` and
`input_M.mtx` are exact canonical copies using the rational-token extension; their
values match the original Matrix Market files. The original file bytes remain in
`matrices/oilrig/` with their original headers and source hashes.

Each certified interval was accepted by exact inertia of the original pencil, exact
rounding decisions and global disjointness. A separate invocation of `secular verify`
replays the certificates without constructing new Graeffe candidates. Its log is included.
The native replay shares the solver's inertia routine and is not an independently authored
proof implementation.

Times are one recorded execution, not a controlled benchmark, confidence interval or claim
of speed superiority. Total includes validation and polynomial-gcd overhead in addition
to the named phases. Do not compare the full certified pipeline's time to a double-precision
kernel as if they deliver the same guarantee.

## Additional different-formulation numerical checks

`reference_binary64.json` records a SciPy generalized symmetric eigensolve. All 66 values
agree at ordinary floating-point accuracy; its tolerance test is not a high-precision proof.

`reference_100_digits.json` records a separate mpmath symmetric eigensolve at 100 decimal
digits after diagonal-mass normalization. Every numerical reference is inside the native
rational interval, and all 50-place rounding decisions agree. The largest distance to a
stored interval center is approximately 4.93e-67. This is a numerical corroboration, not
an independently rigorous interval certificate. Neither reference supplied a Graeffe
candidate or was called by the native solver.

## Scope

No formal program verification, arbitrary complex-spectrum solve, uncertain physical-input
analysis or large sparse-matrix performance claim is made. The exact arithmetic treats the
published decimals as rational constants. It does not establish 50-place accuracy of the
underlying physical measurements. No new run of the manuscript's entire native campaign
was performed or relabeled as a SecularLab result.
