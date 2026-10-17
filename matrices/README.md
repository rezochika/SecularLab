# Included matrices

The 18 files in `examples/` are synthetic regressions or exact paper-derived discretizations.
The two files in `oilrig/` are externally supplied engineering data. They must not be
confused with the synthetic examples.

| Files | Purpose | Expected result |
|---|---|---|
| `spd2` | Small dense symmetric matrix | Eigenvalues 3 and 1 |
| `spd3` | Small irrational spectrum | Three certified distinct positive values |
| `spring10`, `spring20` | Fixed-end unit spring chain, M=I | Eigenvalues 2-2*cos(j*pi/(n+1)) |
| `ties2` | Exact decimal midpoint outputs | 9/4 and 1/4; at one place 2.3 and 0.3 |
| `cluster3` | Close but distinct roots | 1, 1+10^-12, 1+2*10^-12 |
| `dense_stiffness3`, `dense_mass3` | Non-diagonal mass, congruence of a diagonal pencil | Eigenvalues 4, 2, 1 |
| `galerkin2_K0`, `galerkin2_K1`, `galerkin2_M` | Paper's first two Legendre trial modes | Exact formal affine coefficients |
| `galerkin6_K0`, `galerkin6_K1`, `galerkin6_M` | Six-mode paper pencil, K(t)=D-tX | Positive simple spectrum for the documented parameter |
| `repeated3` | Repeated root | `REPEATED_SPECTRUM`, exit 1 |
| `indefinite2` | Negative eigenvalue | Root solver domain refusal, exit 2 |
| `nonsymmetric2` | General matrix | Coefficients allowed; root solver refused |
| `singular_mass2` | Singular leading matrix | Refused, exit 2 |
| `oilrig/bcsstk02`, `oilrig/bcsstm02` | Full external condensed oil-rig stiffness/mass pair | All 66 values certified in the delivery test |

`tools/generate_examples.py` regenerates only `examples/`. It forms the Galerkin entries
by exact polynomial integration of psi_i=P_(i+1)-P_(i-1), independently of hand-transcribed
matrix-entry formulas. It never modifies the external matrices. The files use rational
numeric tokens where needed; see `docs/INPUT_OUTPUT.md` for this Matrix Market extension.
