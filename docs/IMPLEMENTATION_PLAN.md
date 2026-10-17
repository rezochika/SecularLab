# SecularLab implementation plan

Goal: a portable CLion/CMake C++20 reference project for the OPM exact-coefficient,
Graeffe-recovery, independently certified symmetric definite pencil workflow.

Architecture: GMP integers/rationals for trusted arithmetic; Boost.Multiprecision
GMP floats only for untrusted logarithmic Graeffe candidates. Exact inertia of the
original pencil is the acceptance test, not floating-point agreement. Numerical
reference routines must never supply candidate roots.

Scope: arbitrary real rational square pencils for coefficient construction;
symmetric positive definite K and M with simple spectrum for full root recovery.
Affine coefficients and exact rational Legendre specialization are separate modes.
Finite precision/squaring limits may produce NOT_CERTIFIED. No physical precision
claim, formal verification claim, speedup claim, or modification of JacobiRoots.

1. Exact rational parsing, Matrix Market I/O, Bareiss determinant, adjugate,
   symmetric inertia including zero pivots. Test exact small matrices and errors.
2. Paper's coefficient recurrence and affine coefficient arrays. Check every
   division and S_n=0; independently compare determinants at n+1 points.
3. Signed-logarithmic Graeffe and exact rounding certificates. Test small known
   spectra, clusters, decimal ties, repeated roots, rejected matrix classes.
4. CLI, portable CLion run configurations, certificates/replay, toy and actual
   oil-rig data, independent comparison script, documentation.
5. Release build/CTest, full engineering run if feasible, independent certificate
   replay and approximate SciPy comparison; package sources/data/results only.

Ruling: high-precision candidate arithmetic uses Boost's GMP backend. MPFR is not
required. None of its approximate operations establishes an acceptance claim.
Ruling: logarithmic coefficient storage avoids exponent overflow in repeated
squaring; insignificant terms may be dropped ONLY in candidate arithmetic.

Completion record: the exact primitives, fixed and affine recurrence, logarithmic
Graeffe, exact certificate/replay, CLI, matrix inputs and shared configurations were
implemented. Release and sanitizer CTests passed. The complete external 66 x 66
pencil certified at 50 places. Final delivery evidence is under verification/.
