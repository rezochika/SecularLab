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
For a diagonal matrix, such as the usual lumped mass matrix, that inertia is the sign count of
its diagonal entries, in O(n) after the O(n^2) test that it is diagonal.
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
inversion of the mass matrix. For diagonal E, d is the product of the diagonal and J_ii=d/E_ii.
For general E, its determinant is computed by fraction-free Bareiss elimination, and a rational
Gauss-Jordan inverse times d gives the integer adjugate.

Matching powers of z in (zE-F)adj(zE-F) = det(zE-F)I shows that S_k is the coefficient of
z^(n-1-k) in adj(zE-F). For symmetric E and F that adjugate is symmetric, so every S_k is symmetric.
Only the entries (i,j) with j>=i are then computed, stored and exactly divided. For diagonal E,
row i of T_k is formed only on and right of its diagonal and divided straight into S_k; for general
E the full T_k is formed and only the upper triangle of J(T_k-c_k I)/d. Nonsymmetric pencils use
the full matrices. In each stage the rows are independent jobs for the `--workers` threads. Every
entry is the same exact integer for any schedule, and a failure reports the lowest failing row.
Every division that is performed is still checked, and S_n=0 is tested on every stored entry.
For symmetric pencils, symmetry makes the skipped lower-triangle divisions and tests redundant.

### Exact 2-adic and 5-adic reduction

Clearing decimal denominators multiplies every entry by the same L, and decimal inputs carry
further powers of 2 and 5. Before the recurrence, the solver and the `coefficients` command
remove that structure by an exact congruence. For each prime q in {2,5}, exponents a_i >= 0
are chosen with a_i + a_j <= v_q(E_ij) and a_i + a_j <= v_q(F_ij) for every nonzero entry; a
diagonal entry needs 2a_i <= v_q. The start a_i = floor(min_j v_q/2) over row and column i is
feasible. Each a_i is then raised, in a fixed order, to the largest value the others allow.
With P = diag(p_i), p_i = 2^(a_i) 5^(b_i), the pencil E' = P^-1 E P^-1, F' = P^-1 F P^-1 is formed
by checked exact divisions. Then alpha = content(E') and beta = content(F') are divided out,
again exactly, giving E'' and F''. The exact identities

    det(zE-F) = (prod_i p_i)^2 det(zE'-F'),
    det(z alpha E'' - beta F'') = sum_k c''_k alpha^(n-k) beta^k z^(n-k)

give c_k = c''_k alpha^(n-k) beta^k (prod_i p_i)^2, where c'' are the coefficients of
det(zE''-F''). That map uses only integer multiplication. The unchanged recurrence, with all its
divisibility and S_n=0 checks, runs on (E'',F''). P keeps zero entries, symmetry and a diagonal
E, so the same symmetric and diagonal paths apply. Any feasible exponents give the same
coefficients; P=I and alpha=beta=1 is the trivial case. The coefficient oracle below still checks
the mapped coefficients against the original (E,F), so the proof does not depend on the
reduction. On BCSSTK05 the largest workspace entry falls from 13150 to 9111 bits, and most
reduced stiffness entries fit one 64-bit limb.

By default the result is then checked independently of the recurrence by a multimodular
characteristic polynomial (`--oracle auto`). For a prime p not dividing det E, det(zE-F) =
det(E) det(zI-E^-1 F) holds over F_p. The residue of E^-1 F is reduced to upper Hessenberg form by
similarity, its characteristic polynomial follows from the Hessenberg recurrence, and all n+1
coefficients of q, including q_0=det E, must match det(E) chi modulo p.

Expanding det(zE-F) multilinearly in its rows and applying Hadamard's inequality to each mixed
determinant bounds every true coefficient by the corresponding coefficient B_m of
prod_i(ceil||E_i|| z + ceil||F_i||), indexed in descending powers like q. Primes are taken from a
fixed descending sequence in (2^61, 2^62), each proven by Miller-Rabin with the first twelve prime
bases 2,...,37. That test is deterministic below 318665857834031151167461 (about 3.18*10^23). Since
every prime exceeds 2^61, k matching primes have a product above 2^(61k). Primes are taken until
61k >= bits(max_m(|q_m|+B_m)), so the product exceeds max_m(|q_m|+B_m) >= |q_m-P_m|. Then
q_m ≡ P_m modulo that product forces q_m = P_m: the identity is proven, deterministically.
One mismatch modulo any prime refutes it.

Primes dividing det E are skipped, never trusted. B_0 = prod_i ceil||E_i|| bounds |det E|, and a
nonzero det E <= B_0 has fewer than bits(B_0)/61 prime factors above 2^61. Once that many are
skipped, det E = 0 is proven and the modular check is inconclusive; the original oracle then
decides. The same oracle also decides when a size-only estimate of limb operations says it is
clearly cheaper. That happens for tiny orders with very long entries, because the modular check
reduces every entry and coefficient once per prime. Both choices are exact proofs.

That original oracle (`--oracle bareiss`, and the auto fallback) checks n+1 distinct points
z=-1,...,-(n+1). At each point a fresh Bareiss determinant of zE-F must equal exact polynomial
Horner evaluation. Both sides have degree at most n, so agreement at n+1 distinct points proves
the coefficient identity in exact arithmetic. Both oracles are formulations independent of the
recurrence, not independent authorship or formal verification of the implementations.
`--no-coefficient-check` skips the oracle; it never disables the recurrence's exact division and
terminal-zero checks.

Before recovery, a primitive integer pseudo-remainder gcd of f and f' detects repeated roots.
They are reported as `REPEATED_SPECTRUM`, with no claim of distinct indexed enclosures.

### Formal affine mode

For F(t)=F0+tF1, the same recurrence is applied coefficientwise. Stage k stores at most k+1
matrix coefficients. The output `[k][j]` is the actual signed coefficient of z^(n-k)t^j.
With checking enabled, specialization at t=0,...,n and a complete univariate coefficient oracle
at each t prove the bivariate coefficient identity, using the degree bounds in each variable.
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

    |a_k^(r) / a_(k-1)^(r)|^(1/2^r).

The positive, distinct spectrum assumption makes this a candidate-generation strategy for
all roots in decreasing order. The program never uses a numerical reference eigensolver,
a bisection eigensolver or Newton polishing to replace this stage.

Literal repeated squaring makes the coefficients' exponents grow like 2^r. Each coefficient
is therefore stored in extended-exponent form m*2^e. The mantissa m is a signed GMP float at
the working precision with 1/2 <= |m| < 1, and e is an exact integer. The pass starts from the
exact coefficients (`mpf_set_z`, then normalization) and divides out a_0's power of two. For
each b_k, every term of the paired sum (the diagonal a_k^2 and each 2 a_j a_(2k-j)) gets an
exact exponent. The terms are aligned to the largest of these exponents, so a term costs one
mantissa product, a binary shift and an addition. The sum is then normalized. No logarithm or
exponential is evaluated.

Terms more than cut bits below the largest term of their sum are ignored in this *untrusted*
arithmetic. Here cut is the bit length of 10^(work+10), the same 10^-(work+10) relative cutoff
the earlier signed-logarithm form used. Loss of the positive-root coefficient sign pattern,
including exact cancellation to zero, triggers a higher-precision restart.

Exponents start as checked 64-bit integers. They can outgrow that range after about 50
squarings, for example near 54 for roots around 1e-100. The first overflow converts the
iterate exactly to unbounded GMP integers for the rest of the pass. The integer exponents are
the same either way, so the candidates do not depend on where the switch happens, and no step
falls back to logarithmic magnitudes. `graeffe_unbounded_exponents` in `result.json` reports
the switch for the final pass.

A candidate is recovered without transcendental functions. For q = |m_k/m_(k-1)| and
d = e_k - e_(k-1), normalized so that 1/2 <= q < 1, each of r rounds makes d even (doubling q
if needed), takes the square root of q and halves d. The result is q*2^d. Each square root
halves the relative error carried by q. All of this uses GMP floats through Boost's GMP
backend; it is not an interval enclosure or a correctly-rounded promise.

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

The fast inertia routine uses fraction-free symmetric elimination. A diagonal matrix needs no
elimination: its eigenvalues are its diagonal entries, so their signs are counted directly.
Symmetric row/column permutations are congruences and do not alter inertia. At each step the
sign of the LDL pivot equals sign(current Bareiss pivot / previous Bareiss pivot).

When all active diagonal entries are zero but an off-diagonal remains, the program falls
back to exact rational block elimination. A block [[0,b],[b,0]], b!=0, contributes one
positive and one negative eigenvalue; its Schur complement is formed exactly. An entirely
zero remainder contributes its dimension to the nullity. No tiny-pivot tolerance, guessed
sign or floating-point Cholesky decision enters this certificate.

For candidate k (0 is largest), construct exact rational endpoints ell,u around its rounded
center. Proposal widths depend on the requested decimal places, candidate magnitude, nearest
candidate gap and working precision. This lets very small or closely spaced roots have
distinct intervals even when their printed decimals are identical. The magnitude and the gap
enter only through their nearest powers of ten: the e with 10^(e-1/2) <= x < 10^(e+1/2),
decided by the exact test 10^(2e-1) <= x^2 < 10^(2e+1). That depends on the value alone, not on
the dyadic form of a converted floating candidate. Its boundaries are irrational, so candidate
noise can move it only next to some sqrt(10)*10^k. The certificates are therefore normally
identical for every working precision and candidate arithmetic that certifies. An earlier
estimate counted the decimal digits of the candidate's numerator and denominator, and varied by
one with the candidate's lowest set bit. Each candidate gets at most three increasingly narrow
proposals, capped by a heuristic working-precision allowance
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

### Polynomial proof

`--certification auto` (the default) first tries a proof that needs no shifted inertia. It is
available only after the oracle has verified q(z) = det(zE-F) for the original pencil. The same
proposals are used, but an endpoint is tested by the exact sign of q. A width is accepted when
q(ell) and q(u) are nonzero with opposite signs, so (ell,u) contains an odd number of roots.
An endpoint with q=0 is itself a root and becomes an exact singleton. The decimal-rounding gate
is the same as above. Once all n rows are formed, a final gate checks the shared row metadata and
every sign again. n disjoint intervals, each holding a root, exhaust the n roots of the degree-n
polynomial q. So each interval holds exactly one simple root, and the sorted rows carry the
correct indices. This argument proves only a complete spectrum, never a partial one.

For x=p/s in lowest terms with s>0, `polynomial_sign` returns the sign of the homogeneous
integer value sum_i q_i p^(n-i) s^i = s^n q(x). Every proposal endpoint of one width w has a
denominator dividing 10^(w+8), and the rounding midpoints have one dividing 10^(d+1). Signs at
such decimal points share a table T_i = q_i*Q^i for a power of ten Q. The table is built once,
and each sign is then the sign of sum_i T_i (xQ)^(n-i) = Q^n q(x), where xQ is an integer. Q is
always the exact key: a point uses a table only after the exact test that den(x) divides Q, and
any other point uses the general evaluation. Both give the identical sign; only the work changes.

Within one solve, the signs of the fixed coefficient vector are kept in a memo keyed by the
exact canonical rational. The final gate therefore re-reads the signs that certify_one computed
instead of evaluating them again. A missing key is still evaluated, so the gate still catches a
row whose endpoints were never tested. `verify` starts with an empty memo.
`polynomial_sign_evaluations` counts the evaluations actually performed, not memo hits.

When the polynomial attempt fails, exact inertia is tried on the same candidates only if it
could succeed. With q = det(zE-F) verified, q(x)=0 exactly when x is an eigenvalue, that is,
when qF-pE is singular. Both modes build the same proposals and differ only in the endpoint test.

- The polynomial test skips a width only when q(ell) and q(u) are nonzero with equal signs.
  Then neither endpoint is an eigenvalue, so the inertia singleton tests cannot pass. The
  interval holds an even number of roots counting multiplicity, so the inertia counts k+1 and
  k, which mean exactly one simple root, are impossible too.
- An endpoint with q=0 is a polynomial singleton, which always passes the rounding gate.
- When the polynomial test accepts a bracket, inertia either rejects it or accepts the same
  interval. It then runs the same rounding and midpoint-sign tests, plus an inertia test.
- The width sequence, its resource limit and the candidate ordering checks do not depend on the
  mode.

So if the polynomial attempt fails for some candidate, or in the ordering checks, the inertia
attempt fails too, and it is skipped. Only a failure of the final gates, after every candidate
had a row, is retried with exact inertia, because inertia may choose different widths. The same
holds for `verify`: under the verified identity, a row whose exact signs fail also fails its
inertia counts, so only a failure of the shared row checks reaches the inertia replay. Neither
change alters what is accepted; only `inertia_evaluations` and the timings change.

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
Symmetric pencils need about half of those operations, and the n rows of each stage are spread
over the workers. The modular oracle needs O(n^3) word operations per prime, plus one pass over
every entry and coefficient to form residues. The number of primes is about the bound's bit
length divided by 61. The n+1 Bareiss determinants instead cost
O(n^4) integer operations on entries as long as the coefficients.
Each Graeffe squaring costs O(n^2) exponent comparisons and one working-precision product,
shift and addition per term that is not negligible. Recovering the candidates after r
squarings costs r square roots per root. Exact inertia is O(n^3) arithmetic per endpoint,
normally two endpoints for each root, and may dominate total time. A polynomial sign costs one
Horner pass of integer products, which a shared scaled table keeps short.

The result timings separate coefficients, coefficient oracle, Graeffe and certification. Total
also includes validation, denominator clearing and squarefree checking, so the named phase
sum need not equal total. The Graeffe squaring count reports the final attempted pass; phase
times include work from earlier precision attempts.

The oil-rig input consists of finite published decimal values. Exact certification applies
to those exact rational matrices, not to unknown physical-parameter uncertainty. No comparison
with a double-precision solver establishes 50-digit accuracy; the exact certificate does that.

The `recurrence_workspace_bits_max` field measures the largest entry of T_k-c_k I at
the sampled recurrence stages, not the peak bit length of every GMP temporary. When only the
upper triangle of T_k is formed (symmetric pencil, diagonal E), each entry below the diagonal is
measured as E_jj*S_k(i,j), which equals it exactly. It is formed only when it could raise the
maximum, so the field keeps its full-matrix meaning. It measures the recurrence actually run, on
the reduced pencil (E'',F''); it equals the unreduced value only when the reduction is trivial.
`coefficient_bits_max` is measured after the coefficients are mapped back.
`exact_divisions` counts the divisions the coefficient recurrence actually performs: one trace
division per stage, plus n(n+1)/2 entry divisions per stage for symmetric pencils, or n^2
otherwise. It excludes those inside Bareiss elimination, the reduction, the oracles or the
inertia routine.
`coefficient_oracle` names the oracle that decided the identity, `oracle_moduli` the primes of a
modular decision, and `oracle_workers` and `recurrence_workers` the threads actually used.
