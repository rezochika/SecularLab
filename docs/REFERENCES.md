# Mathematical and software references

## Manuscripts supplied in the conversation

1. Tamaz Vashakmadze and Revaz Chikashua, *Orthogonal-Polynomial Methods for Boundary-Value
   Problems and Secular Equations*, supplied OPM manuscript, Sections 6-8, Theorem 9,
   Algorithms 6-7, equations 61-78 and Appendix A. The coefficient convention in the code
   follows this corrected formulation, not an uncorrected historical recurrence.
2. Tamaz S. Vashakmadze, Rezo Chikashua and Gela Manelidze, *To Problem of Reliable Calculation
   of Coefficients of Secular Equation by Special Functions*, AIP Conference Proceedings
   1558, 1164-1166 (2013), DOI 10.1063/1.4825715.
3. Historical attribution: R. Chikashua's 2009 master's thesis, supervised by Tamaz Vashakmadze.
   This archive does not redistribute the supplied papers or thesis and does not assert that
   the newly written code is their original implementation.

## Classical methods

The Faddeev-LeVerrier recurrence, determinant-adjugate identities, Graeffe root squaring,
Sylvester inertia principle, Bareiss fraction-free elimination and exact polynomial gcds
are classical. The OPM manuscript cites their antecedents. This project claims an executable
implementation and tested outputs, not new priority for those identities.

## Matrix provenance

- https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstk02.html
- https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstm02.html
- See `matrices/oilrig/PROVENANCE.md` for the exact data-copy commit and hashes.

## Dependencies and build guidance

- GMP: https://gmplib.org/ - exact integers/rationals and the candidate float backend.
- Boost.Multiprecision: https://www.boost.org/doc/libs/release/libs/multiprecision/doc/html/index.html
- CMake presets: https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html
- CLion presets: https://www.jetbrains.com/help/clion/cmake-presets.html
- CLion shared run configurations: https://www.jetbrains.com/help/clion/run-debug-configuration.html
- MSYS2 GMP: https://packages.msys2.org/packages/mingw-w64-ucrt-x86_64-gmp
- MSYS2 Boost: https://packages.msys2.org/packages/mingw-w64-ucrt-x86_64-boost

Official dependency/IDE instructions were consulted when assembling this archive.
No dependency binaries, font files, or proprietary IDE files are bundled.
