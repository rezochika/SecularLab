# Algorithm and acceptance conditions

## Relationship to the source papers

The coefficient construction follows the nonsingular integer-pencil adjugate recurrence
(Theorem 9 / Algorithm 6, equations 67-70) in the uploaded OPM manuscript. Graeffe uses
its equation 77 and the distinct-positive-root coefficient-ratio limit. This project
adds an explicit executable path, exact input parsing, source-independent determinant
identity tests, bounded retries, output files and a replay command.

This is new reference code, not the historical 2009 implementation, the 2013 authors'
modified special-function evaluator, or the existing JacobiRoots native release.
The root stage intentionally restricts its domain to symmetric positive definite pencils
with distinct eigenvalues. It does not implement general complex-root or equal-modulus
recovery. A refusal is preferable to silently using a different root algorithm.

## 1. Exact input and scale

Each textual matrix number is parsed as a rational, without passing through a double.
For example `1.23e-2` becomes `123/10000`. The original external matrices are not rounded
or altered. Symmetric coordinate storage is expanded exactly; duplicate coordinates
within the same triangle are summed. Mixed upper/lower symmetric storage is rejected.

Let L be the least common multiple of the denominators in K and M. Set E=LM and F=LK.
Then E,F are integer matrices and

    det(zE-F) = L^n det(zM-K).

The output coefficient file stores the left-hand polynomial, with L recorded. Its leading
coefficient is det(E), not 1. Eigenvalues are unchanged by this known scalar factor.

The solver proves positive definiteness of both E and F by exact inertia before root recovery.
The general coefficient command only requires nonsingular E.

## 2. Exact characteristic coefficients

Write det(zE-F) = d*z^n - sum(c_k*z^(n-k)), d=det(E). Initialize J=adj(E), S_0=J.
For k=1,...,n:

    T_k = F S_(k-1)
    c_k = trace(T_k)/k
    S_k = J(T_k-c_k I)/d.

`exact_divide` explicitly tests GMP divisibility before every stated division, including
Release builds. The last S_n must be the zero matrix. The exported vector is
`[d, -c_1, ..., -c_n]`, in descending powers. These signs matter.

For diagonal E, applying J*T/d is equivalent to dividing row i of T by E_ii. The implementation
uses this exact cancellation, with a divisibility check on every entry. It is not floating-point
inversion of the mass matrix. For general E, its determinant is computed by fraction-free
Bareiss elimination; a rational Gauss-Jordan inverse times d gives the integer adjugate.

By default the result is independently checked at n+1 distinct points z=-1,...,-(n+1).
At each point a fresh Bareiss determinant of zE-F must equal exact polynomial Horner evaluation.
Both sides have degree at most n, so agreement at n+1 distinct points proves the coefficient
identity in exact arithmetic. This is an independent formulation, not independent authorship
or formal verification of the implementations. `--no-coefficient-check` skips this additional
oracle; it never disables the recurrence's exact division and terminal-zero checks.

Before recovery, a primitive integer pseudo-remainder gcd of f and f' detects repeated roots.
They are reported as `REPEATED_SPECTRUM`, with no claim of distinct indexed enclosures.

### Formal affine mode

For F(t)=F0+tF1, the same recurrence is applied coefficientwise. Stage k stores at most k+1
matrix coefficients. The output `[k][j]` is the actual signed coefficient of z^(n-k)t^j.
With checking enabled, specialization at t=0,...,n and direct determinants at n+1 distinct
z-values prove the bivariate coefficient identity, using the degree bounds in each variable.
Formal affine mode has a separate order limit of 24.

`solve --parameter t` instead specializes exact rational K0+tK1 *before* the numeric pencil
recurrence. `--legendre-degree m --argument x` evaluates P_m(x) exactly with rational
three-term recurrence, then performs that specialization. The latter is valid for rational x
but does not reproduce the paper's interval special-function enclosure experiment.

## 3. Graeffe candidate generation

For descending coefficients a_0,...,a_n, one root-squaring step computes

    b_k = (-1)^k a_k^2 + 2 sum_{j<k} (-1)^j a_j a_(2k-j),

with indices outside 0,...,n omitted. It represents g(z^2)=(-1)^n f(z)f(-z).
The roots of g are the squares of the roots of f. After r steps, successive coefficient
ratios yield approximate roots

    exp((log|a_k^(r)| - log|a_(k-1)^(r)|)/2^r).

The positive, distinct spectrum assumption makes this a candidate-generation strategy for
all roots in decreasing order. The program never uses a numerical reference eigensolver,
a bisection eigensolver or Newton polishing to replace this stage.

To prevent the enormous exponent growth of literal repeated squaring, each coefficient
is represented by its sign and logarithmic magnitude. Signed additions implement the same
paired recurrence, with leading-log normalization. Very small relative terms may be
ignored in this *untrusted* arithmetic. Loss of the positive-root coefficient sign pattern
triggers a higher-precision restart. The logarithms/exponentials use Boost's GMP float
backend; they are not interval enclosures or a correctly-rounded special-function promise.

Candidates are considered for exact certification after every second squaring once successive
values are numerically stable to the heuristic absolute threshold 10^-(digits+16), and also
at the configured squaring limit. Stability is not proof or evidence that working precision
is exhausted. After a failed exact check, the same pass continues squaring, with subsequent
checks spaced 2, 4, 8, ... squarings apart to limit expensive inertia work. The final configured
step is always checked. Actual arithmetic precision loss, or exhausting the squaring budget,
can trigger a higher-precision pass within the fixed upper limit. No unlimited-success claim
is made.

## 4. Exact inertia and certification

For rational s=p/q with q>0, the integer symmetric matrix

    qF - pE

has the same inertia as K-sM. Because M is positive definite, its positive eigenvalue count
is exactly the number of generalized eigenvalues greater than s.

The fast inertia routine uses fraction-free symmetric elimination. Symmetric row/column
permutations are congruences and do not alter inertia. At each step the sign of the LDL
pivot equals sign(current Bareiss pivot / previous Bareiss pivot).

When all active diagonal entries are zero but an off-diagonal remains, the program falls
back to exact rational block elimination. A block [[0,b],[b,0]], b!=0, contributes one
positive and one negative eigenvalue; its Schur complement is formed exactly. An entirely
zero remainder contributes its dimension to the nullity. No tiny-pivot tolerance, guessed
sign or floating-point Cholesky decision enters this certificate.

For candidate k (0 is largest), construct exact rational endpoints ell,u around its rounded
center. Proposal widths depend on the requested decimal places, candidate magnitude, nearest
candidate gap and working precision. This lets very small or closely spaced roots have
distinct intervals even when their printed decimals are identical. Each candidate gets at
most three increasingly narrow proposals, capped by a heuristic working-precision allowance
that retains the original output-driven widths (so exact candidates can certify beyond their
working precision), and by the decimal-scale resource limit. These heuristics only propose intervals; none relaxes
the exact acceptance conditions. A nonsingular endpoint pair must satisfy

    positive_count(K-ell*M) = k+1,
    positive_count(K-u*M)   = k.

Therefore exactly the kth eigenvalue lies in (ell,u). A singleton is accepted only when
exact inertia proves nullity one and positive count k. Every accepted interval must be
disjoint from its neighbors, and the number of accepted intervals must equal n.

For requested d decimal places, use the exact half-away-from-zero map

    R_d(x) = sign(x) floor(10^d |x| + 1/2) / 10^d.

Both rational endpoints must have the same R_d value. If the interval crosses one rounding
midpoint, exact polynomial evaluation and inertia may prove that the eigenvalue equals that
midpoint, giving an exact singleton and the tie decision. Otherwise the run is not certified.
Rational endpoint and integer-size costs are real; they are not treated as constant-time work.

The stored radius is max(|printed-ell|,|printed-u|), an exact rational bound about the *printed*
value. A root at 1/4 rounded to 0.3 needs a radius at least 0.05 regardless of internal precision.

## 5. Replay and failure reporting

`verify` rereads the matrices and TSV and re-establishes positive definiteness, endpoint
inertias, root indices, disjointness, exact rounding and every radius bound. It does not trust
an earlier status string. The replay shares the native inertia implementation with the solver;
calling it an independent-authorship proof checker would overstate its independence.

- `CERTIFIED`, exit 0: all n eigenvalues accepted with stated decimal rounding and radii.
- `REPEATED_SPECTRUM`, exit 1: exact polynomial gcd detects repetition; simple-root recovery refused.
- `NOT_CERTIFIED`, exit 1: configured candidate/squaring limits did not prove every output.
- `ERROR`, exit 2: input/domain/resource/I/O exception handled by the program.
- `verify` exit 1: invalid certificate; exit 2: malformed input or an I/O error.

There is no partial-success eigenvalue file. Output directories are never silently overwritten.
GMP and the C++ runtime are trusted dependencies; operating-system termination or GMP allocation
failure cannot be promised to produce a structured report. The source is not formally verified.

## 6. Costs and interpretation

Fixed numeric dense coefficients need O(n) matrix products, hence O(n^4) ordinary integer
operations with classical multiplication, plus determinant/adjugate construction. The formal
affine recurrence needs O(n^5) coefficientwise integer operations. These are not bit-complexity
bounds. Integer coefficient length, input denominators and elimination growth can dominate.
Each Graeffe pass costs O(n^2) high-precision signed-log operations. Exact inertia is O(n^3)
arithmetic per endpoint, normally two endpoints for each root, and may dominate total time.

The result timings separate coefficients, coefficient oracle, Graeffe and certification. Total
also includes validation, denominator clearing and squarefree checking, so the named phase
sum need not equal total. The Graeffe squaring count reports the final attempted pass; phase
times include work from earlier precision attempts.

The oil-rig input consists of finite published decimal values. Exact certification applies
to those exact rational matrices, not to unknown physical-parameter uncertainty. No comparison
with a double-precision solver establishes 50-digit accuracy; the exact certificate does that.

The `recurrence_workspace_bits_max` field measures the largest entry of T_k-c_k I at
the sampled recurrence stages, not the peak bit length of every GMP temporary.
`exact_divisions` counts divisions in the coefficient recurrence, not those inside
Bareiss elimination, its oracle or the inertia routine.
