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

The solver proves positive definiteness of both E and F exactly before root recovery.
For a diagonal matrix, such as the usual lumped mass matrix, that is the sign count of
its diagonal entries, in O(n) after the O(n^2) test that it is diagonal. Any other matrix is
first tried with the congruence certificate of section 4 and otherwise decided by exact inertia.
The general coefficient command only requires nonsingular E.

## 2. Exact characteristic coefficients

Write det(zE-F) = d*z^n - sum(c_k*z^(n-k)), d=det(E). Initialize J=adj(E), S_0=J.
For k=1,...,n:

    T_k = F S_(k-1)
    c_k = trace(T_k)/k
    S_k = J(T_k-c_k I)/d.

`exact_divide` checks every stated division, including Release builds: one truncating GMP
division must leave remainder zero, or the run stops. The same check replaces the separate
divisibility test inside Bareiss elimination and the fraction-free inertia routine. The last S_n must be the zero matrix. The exported vector is
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
det(E) det(zI-E^-1 F) holds over F_p. The residue of E^-1 F is reduced by similarity to upper
Hessenberg form (or, for a symmetric pencil, to tridiagonal form, below), its characteristic
polynomial follows from the Hessenberg recurrence, and all n+1 coefficients of q, including
q_0=det E, must match det(E) chi modulo p.

Every residue is a field computation modulo p, so any correct algorithm gives the same residues;
the kernel is chosen for speed. It forms B = (E^-1 F)^T, which has the characteristic polynomial
of E^-1 F. A diagonal E scales the nonzero entries of F after one batch inversion of its diagonal.
A general E is factored PE = LU with partial pivoting in Crout order, and each column of F is
solved by substitution, skipping the leading zeros of L, U and the right-hand side. B is reduced to
upper Hessenberg form H by Gaussian similarity with pivoting, column by column through B L = L H
with L unit lower triangular: column j of H solves a triangular system against B l_j, and the rest
of that product gives l_(j+1). While the columns below the subdiagonal stay sparse, the classical
elimination row r -= u_r row j+1, column j+1 += u_r column r runs instead and skips zero rows and
columns; it hands over at the first column with more than a quarter of those entries nonzero.
det(zI-H) then follows from the recurrence p_k = z p_(k-1) - sum_(s<k) g_s p_s, one dot product per
coefficient.

When E and F are both symmetric, which is decided once per pencil, two routes use the symmetry.
A non-diagonal E is factored E = L D L^T without pivoting, L unit lower triangular, and the
congruence C = L^-1 F L^-T gives det(zE-F) = det(zD-C) with C symmetric. The factor, the upper
triangle of L^-1 F and the lower triangle of C take n^3/6 + n^3/6 + n^3/3 products, against
n^3/3 + n^3 for LU and the solves, again skipping the leading zeros of the rows of E. The pencil
(D, S), with D diagonal and S symmetric (S = C, or S = F when E is diagonal), then goes to Lanczos
in the inner product <x,y> = x^T D y, for which D^-1 S is self-adjoint. From a fixed
pseudo-random start, the recurrence q_(k+1) = D^-1 S q_k - alpha_k q_k - beta_k q_(k-1), with
delta_k = <q_k,q_k>, alpha_k = q_k^T S q_k/delta_k and beta_k = delta_k/delta_(k-1), keeps the q_k
pairwise D-orthogonal over any field while every delta_k is nonzero. After n steps they are a
basis in which D^-1 S is tridiagonal, and p_(k+1) = (z-alpha_k) p_k - beta_k p_(k-1) gives its
characteristic polynomial in O(n^2). Each step is one product with S and O(n) work, all in dot
products. A vanishing delta_k is the exact breakdown: an invariant subspace or an isotropic vector,
rare for 61-bit primes unless eigenvalues repeat, and common for small ones. That modulus then
reduces D^-1 S to Hessenberg form instead. A zero LDL^T pivot, possible without p dividing det E,
sends the modulus to the LU route. For a diagonal E and a sparse F, Lanczos costs about n(T+20n)
products, with T the entries of F used per product: its nonzeros when they fill less than half of
the rows' ranges, otherwise those ranges. The Hessenberg form of a banded F usually fills in, but
exact cancellation can keep it sparse, as for the band-3 rotation pencil. So unless F's rows are
mostly full, the sparse reduction runs once, for one prime, counting its products, and the pencil
keeps the cheaper route.

A banded E has a banded L but in general a dense L^-1, and so a dense C: forming C is cheap, but
Lanczos on it takes n^3 products. A third route keeps C implicit and runs Lanczos on (E, F) in
the inner product x^T E y. Each q = E^-1 u solves L y = u row by row as u is formed and then
L^T q = D^-1 y from the last row up, on the rows of L and of L^T (the columns of L, stored as rows
with their last nonzero), so a step is one product with F and one with each of L and L^T, about
n(4b+1) products for E and F of half-bandwidth b, and O(b n^2) for all n steps. The substitutions
are chains of dependent rows, so these products cost about 1.1 times the congruence's per product,
and the congruence's C need not be dense: F = L Lambda L^T, as in pencils built from a common factor
(cluster32 and the A-mid4 pencils, not cluster32lu), keeps it diagonal, and then forming it,
(n+1) sum_i (i-first_i) products, and Lanczos on it cost less. Only C shows which case holds. So unless even a dense C
would favour the congruence, as for a dense E behind a dense F, one prime forms C once per pencil,
and the pencil keeps the route with fewer products, the factored route's counted 1.1 times; the
overheads per row of a Lanczos step and per entry of C came to about the same and cancel. The
factored route falls back as the congruence does, a zero pivot to the LU route and a breakdown to
the general route.

The arithmetic is Montgomery's: products below p^2 < 2^124 are summed sixteen at a
time in 128 bits, then into 192 bits, and reduced once per dot product, and fixed multipliers use
Shoup's precomputed quotient. So nearly every product is one 64x64-bit multiplication and an
addition. The tests check every route against the plain Gauss-Jordan and Hessenberg reduction the
kernel replaced, on random and structured pencils and on small primes, where zero pivots, Lanczos
breakdowns and det E = 0 modulo p are common, and force each fallback at the large primes.

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

### Coefficients by Chinese remaindering

The oracle's work for one prime p already yields every coefficient of a pencil's determinant
modulo p, as det(E) chi with chi the characteristic polynomial of E^-1 F.
`--coefficients crt` keeps these residues and combines them instead of running the recurrence.
It takes them on the reduced pencil (E'',F'') of the 2-adic and 5-adic reduction above. With
f_k = alpha^(n-k) beta^k (prod_i p_i)^2 > 0, the identities there give q_k = f_k c''_k, where c''
are the coefficients of det(zE''-F''). Each c''_k has two bounds. B''_k is the Hadamard-type
bound of (E'',F''). Also |c''_k| = |q_k|/f_k <= B_k/f_k, and c''_k is an integer, so
|c''_k| <= floor(B_k/f_k). Let b_k be the smaller of the two; for a trivial reduction b_k = B_k.

Primes are taken from the same fixed sequence, skipping those dividing det E'', until k usable
primes satisfy 61k >= bits(max_m b_m)+1. Their product M then exceeds 2 max_m b_m. For each m,
Chinese remaindering gives the unique x_m with |x_m| <= (M-1)/2 that agrees with c''_m modulo
every prime. Since |c''_m| <= b_m < M/2, x_m = c''_m. Then q_m = f_m x_m by integer multiplication
alone. In the oracle's terms, M exceeds |x_m|+b_m <= M/2+b_m, and x agrees with det(E'') chi''
modulo each prime, so the construction is itself the multimodular proof of c''. The reduction's
identities make it a proof of q: E'' and F'' are formed by checked exact divisions, so
E = alpha P E'' P and F = beta P F'' P hold exactly. Unlike the recurrence, whose result the
oracle checks against the original (E,F), this proof uses the reduction; both are exact.
Running the modular oracle afterwards would recompute residues of the same determinant, so
`coefficient_oracle` reports `crt` and `oracle_moduli` the primes combined; `--oracle bareiss`
still adds the n+1 Bareiss determinants of the original pencil. Skipped primes are counted as in
the oracle: b_0 bounds |det E''|, and once bits(b_0)/61 primes are skipped, det E'' = 0 is
proven, so det E = f_0 det E'' = 0. The request then fails with the recurrence's error for a
singular leading matrix. Like the modular oracle, this route needs 64-bit GMP limbs and 128-bit
integers.

The reduction pays where q carries 2-adic and 5-adic structure. On the clustered A-mid4 pencil
with gap 10^-100, q has 22145-bit coefficients and needed 364 primes; c'' is bounded by 2548 bits,
and 42 primes suffice. The minimum matters: the reduced pencil's own bound B'' has 3202 bits there.
On BCSSTK05 220 primes fall to 153, on the order-256 lumped Laplacian 73 fall to 8. Pencils
without such structure, such as the dense congruences M = C^T C, keep their prime count. The
original pencil's bound is computed on a helper thread beside the reduction, so the extra serial
work is about one bound computation.

The combination uses the weights w_i = (M/p_i)^-1 mod p_i, formed once with O(k^2) word
operations. Each coefficient is then x = sum_i ((r_i w_i) mod p_i) M/p_i, summed up a product
tree whose node with children (x_L,M_L), (x_R,M_R) holds x_L M_R + x_R M_L, and reduced modulo M
into the symmetric range. Residues are computed one prime per worker, coefficients one per
worker. No recurrence runs, so `exact_divisions`, `recurrence_workspace_bits_max` and
`recurrence_workers` are zero; `coefficient_moduli` and `crt_workers` describe the combination,
`crt_bound_bits` is bits(max_m b_m), and `crt_reduced` tells whether the reduction was nontrivial.
The integers are those of the recurrence, so `coefficients.tsv` and the certificates are
byte-identical.

`--coefficients both` runs the recurrence and Chinese remaindering and requires identical
integers: an N-version check of the two constructions, which share the input pencil, the exact
reduction and GMP. The identity proof is the Chinese-remainder one. Both routes are themselves
proofs, so they reject `--no-coefficient-check`; neither applies to formal affine coefficients.

### Why the stopping point is a bound

The primes follow a proven bound on c'', not the actual length of q or of q/content(q). An
output-sensitive variant would reconstruct the monic polynomial chi = q/q_0 of E^-1 F, or the
primitive part q/content(q), by rational reconstruction, and stop once the result is stable over a
few extra primes. Stability is a heuristic, not a proof, and no cheaper proof is available here:

- Residues alone prove a rational a/b only through a priori bounds. Two fractions with numerators
  at most N and denominators at most D that agree modulo M are equal once M > 2ND. For
  chi_k = q_k/q_0 the only a priori denominator bound is |det E|, and then ND is at least as long
  as q_k itself.
- Scaling a candidate back to q~ = (det E/lc) times the primitive part and running the modular
  oracle needs 61k >= bits(max_m(|q~_m|+B_m)) >= bits(max_m B_m): the primes Chinese remaindering
  uses anyway, which would leave nothing for the early stop to save.
- Fewer primes need a proven common factor of the coefficients, which shortens what remains to be
  proven. The reduction supplies such factors f_k. The rest of the content, for example
  det(C)^2 for M = C^T C and K = C^T diag(lambda) C, would need its own certificate: determinantal
  divisors of E and F, a Smith form of the pencil, or a factorization of det E. None is cheap in
  general.

On the 16 pencils of the FLINT/Arb comparison, max_m b_m exceeds the longest reduced coefficient
by at most 206 bits, four primes, except for the dense range pencil of order 64 (13466 against
7461 bits). There a spread spectrum makes Hadamard's inequality loose; a sharper proven bound
would need eigenvalue information.

### Content removal

With decimal input, q = det(zE-F) usually carries a large content g = gcd(q_0,...,q_n): on
BCSSTK05 g has 4387 bits, and q/g has 8778-bit coefficients against 13165 for q. Every
`result.json` reports `coefficient_bits_max` (q), `primitive_coefficient_bits_max` (q/g) and
`coefficient_content_bits`. `--remove-content` runs the squarefree test, the Graeffe stage and
every polynomial sign on q/g. With g > 0, (q/g)(x) has the sign of q(x) at every rational x, so
both have the same roots, brackets and singletons, and every acceptance test is unchanged.
`coefficients.tsv` still holds q. The Graeffe mantissas of q/g round differently, so candidates
may differ in their last working bits; proposal widths see them only through powers of ten and
centers rounded at width+8 places, so certificates normally do not change. On the 16 pencils of
the FLINT/Arb comparison they were byte-identical with and without `--remove-content`.

Before recovery, squarefreeness is tested modulo the fixed primes 2, 3, 5, 7, 11, 65537, 2^31-1,
2^32-5 and 2^32-17, in that order. A prime that keeps the degree and gives gcd(f, f') = 1 modulo p
proves f squarefree. Otherwise a primitive integer pseudo-remainder gcd of f and f' decides
exactly. The primes below 2^32 make that slower fallback rare when small primes divide the
leading coefficient or the discriminant. Repeated roots are reported as `REPEATED_SPECTRUM`,
with no claim of distinct indexed enclosures.

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

Every inertia, of an endpoint matrix as of F and E in the positive-definiteness checks, is
decided by one of three exact tests, tried in this order. A diagonal matrix needs no
elimination: its eigenvalues are its diagonal entries, so their signs are counted directly.
Otherwise the congruence certificate at the end of this section is tried, in which floating
point only proposes and an exact integer test decides. When it fails, fraction-free symmetric
elimination decides. Symmetric row/column permutations are congruences and do not alter
inertia. At each step the sign of the LDL pivot equals sign(current Bareiss pivot / previous
Bareiss pivot).

When all active diagonal entries are zero but an off-diagonal remains, the elimination falls
back to exact rational block elimination. A block [[0,b],[b,0]], b!=0, contributes one
positive and one negative eigenvalue; its Schur complement is formed exactly. An entirely
zero remainder contributes its dimension to the nullity. No tiny-pivot tolerance, guessed
sign or floating-point Cholesky decision enters any of the three tests.

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
A table is not built when it would exceed 2^28 bits, or when Q is much longer than den(x), as
for den(x)=2^a; such points use the general evaluation too.

Within one certification attempt, the signs of the fixed coefficient vector are kept in a memo
keyed by the exact canonical rational. The final gate therefore re-reads the signs that
certify_one computed instead of evaluating them again. A missing key is still evaluated, so the gate still catches a
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

### Congruence certificate

Elimination is exact but costly, because Bareiss entries grow into leading minors about k times
as long as an entry. On BCSSTK05 one elimination took 0.84 s for F, and 7-8 s for an endpoint
matrix qF-pE, whose entries are longer by the bits of q. For a symmetric A that is not diagonal,
every inertia therefore first tries a certificate whose integers stay short.

- Proposal, untrusted. A GMP floating-point LDL^T of A without pivoting, at a precision P set on
  each object, gives W, approximately L^-1. The process-global default precision is neither read
  nor changed. X = trunc(2^(P-8)*W^T) is an integer matrix, upper triangular with 2^(P-8) on its
  diagonal by construction. A zero floating pivot abandons the proposal, and so does an entry of
  W above 2^(P+b), where b is the bit length of the longest entry of A; that bounds the exact work.
  One n-by-n floating array suffices: W^T overwrites the factored upper triangle, which is no
  longer needed once L and the pivots are known.
- Check, exact. X must be upper triangular with a nonzero diagonal, so det X != 0. B = X^T A X
  is formed in integers, and every row must satisfy |b_ii| > sum_(j!=i) |b_ij|. B is formed one
  row of its upper triangle at a time and not kept. That row completes the off-diagonal sum of
  row i, since the entries b(k,i), k<i, came with earlier rows, so a failing row ends the check.

When the check passes, B is symmetric and strictly diagonally dominant. B and every
diag(B)+t(B-diag(B)), 0<=t<=1, are then nonsingular (Levy-Desplanques), so no eigenvalue crosses
zero along that path, and B has the inertia of diag(B). By Sylvester's law of inertia so does A,
with nullity zero. A is positive definite exactly when every b_ii is positive. A negative b_ii
refutes positive definiteness just as exactly.

The floating-point stage decides nothing: a poor proposal can only fail the exact check. A tie in
the dominance test is a failure. A singular A always fails, because a singular B cannot be
strictly diagonally dominant. Each A gets a short, deterministic schedule of precisions, tried in
order; when none passes, the fraction-free inertia above decides, unchanged. No precision above
16384 bits is tried. The work is O(n^3) floating-point operations at a few limbs, plus about
n^3/6 integer products of entries near 2P+b bits.

- Positive definiteness: P = 64, 128 and 256. On BCSSTK05 the 64-bit proposal certifies F in
  0.02 s. The oil-rig K certifies at 64 bits too. Badly conditioned matrices need more: integer
  Hilbert matrices of order 16 need 128 bits, order 30 needs 256, and from about order 60 they
  fall back to elimination.
- Endpoints: a proposal endpoint s=p/q of width w lies within about 10^-w of an eigenvalue, and
  q divides 10^(w+8). Near an eigenvalue qF-pE is nearly singular, with a condition number that
  grows like 1/|s-lambda|, so the precision grows with q: P = 64+bits(q), rounded up to whole
  64-bit limbs, then 2P and 4P. At 50 places that is 320, 640 and 1280 bits. All 132 oil-rig
  endpoints and all 306 BCSSTK05 endpoints certify at 320 bits; an earlier measurement found
  only 74 of the 132 oil-rig endpoints at 192 bits, and all at 256. An oil-rig endpoint takes
  about 8 ms instead of 0.62 s, a BCSSTK05 endpoint about 0.03 s instead of 6.8 s. An exact
  singleton endpoint is an eigenvalue, so qF-pE is singular: its proposals always fail, and
  elimination proves the nullity. A denominator longer than about 16300 bits leaves no precision
  in the schedule, and elimination decides directly.

Certificates and `inertia_evaluations` are unchanged: the counter counts inertias, whichever test
decided them.

### Parallel inertia certification

Exact inertia certifies one root per job, on as many threads as the coefficient phases use
(`--workers`, clamped to n). In `solve` a job is certify_one for one candidate; in `verify` it is
the inertia test of one row. The rows, the counters and any error are those of the serial loop
for every worker count:

- The candidate ordering checks run first, serially. Only the roots below the first candidate
  that fails them are submitted: the roots the serial loop could reach.
- Roots are handed out in increasing order. The run stops at the lowest root that fails or
  throws, as the serial loop would, and every root below it has run and passed. Roots above it
  may have run too, but their rows, counts and exceptions are discarded, and no worker starts a
  root above a known failure.
- A job counts its own inertias and records the polynomial signs it requests; inertia mode
  needs them only at rounding midpoints. It may read the attempt's sign memo, but never changes
  it. After the run the roots up to the stop, in index order, add their inertia counts and
  replay their sign requests into the memo: a point missing from it counts one evaluation and is
  stored. So `inertia_evaluations` and `polynomial_sign_evaluations` equal the serial counts,
  and so does the memo that a later step of the same attempt reads.
- An exception is rethrown after this counting, and only when it came from the stopping root. An
  exception from a root above a failure is dropped, since the serial loop never reaches that root.
- In `solve`, each pass runs in two parts: first the prefix of roots 0..f, where f is the root at
  which the previous pass stopped (initially 0), then, only if all of them pass, the remaining
  roots. Roots above the one that decides an attempt would only compete with it for cores and be
  discarded. Precision-limited attempts fail at the largest root, which then runs alone. Prefix
  then suffix is still the serial loop, so the stop, its exception and the counters are unchanged.
- A root that is already running when a lower one stops abandons its exact work. Elimination
  polls once per pivot, the proposal once per factorization step and the dominance check once per
  row, whether the run has stopped below that root; if so they throw, and the exception is
  discarded with the rest of the root. The poll never fires for a root at or below the stop, so a
  hand-edited certificate with a cheap wrong first row and huge later rows fails as fast as in
  serial replay.

Each running root holds its own proposal and, when elimination decides, its own Bareiss matrix,
so peak memory grows with the worker count; `--workers 1` bounds it to one root. Each point of
the Bareiss coefficient oracle likewise holds its own n x n determinant matrix. The default worker
count is the hardware thread count, which cannot see CPU-rate quotas or cgroup limits; under those,
pass `--workers` equal to the quota. The polynomial
proof stays serial: its signs share one evaluator's scaled tables and memo, and its whole
certification takes about 0.05 s on BCSSTK05. With eight workers `--certification inertia` now
solves the oil rig in about 0.27 s instead of 78 s, and BCSSTK05 in about 2 s instead of 2024 s
(1.7 s and 14 s with one worker).

## 5. Replay and failure reporting

`verify` rereads the matrices and TSV and re-establishes positive definiteness, endpoint
inertias, root indices, disjointness, exact rounding and every radius bound. It does not trust
an earlier status string. The replay shares the native inertia implementation and the
congruence certificate with the solver;
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
length divided by 61. Chinese remaindering needs the same residues for about as many primes,
plus O(k^2) word operations for its weights and one product tree per coefficient, so it replaces
the O(n^4) recurrence and the oracle by roughly the oracle's cost. On the 16 pencils of the
FLINT/Arb comparison (Ryzen 9 7950X, Windows, CLion MinGW Release; medians), `--coefficients crt`
cut total solve time by 1.2-364x on 32 threads and 1.2-141x on one thread, with byte-identical
artifacts in all 277 runs. BCSSTK05 took 0.17 s instead of 0.52 s on 32 threads (1.25 s instead of
4.8 s on one), the n=100 LU-mass pencil 0.08 s instead of 12.4 s (0.57 s instead of 70 s), and the
dense mass of order 256 1.7 s instead of 602 s (16 s on one thread). Taking the residues on the
reduced pencil, with the bound b of section 2, then cut the primes where q carries 2-adic and
5-adic structure. In 192 further runs on the same pencils, again with byte-identical artifacts,
the one-thread coefficient stage fell from 0.145 to 0.021 s on the A-mid4 pencil with gap
10^-100 (364 to 42 primes), from 0.057 to 0.012 s with gap 10^-40 (155 to 26), from 0.29 to
0.038 s on the lumped Laplacian of order 256 (73 to 8) and from 1.2 to 0.85 s on BCSSTK05 (220 to
153). Pencils without that structure kept their primes and times within noise; the reduction
then costs a few milliseconds of serial work. Total solve time changed less, because Graeffe and
certification dominate there. The residue kernel of section 2 then cut the work per prime without
changing a residue or the number of primes. On one thread, a prime took 7.9 ms instead of 48 ms on
the order-256 pencil with a diagonal mass and a dense stiffness, and 0.41 ms instead of 2.2 ms on
the dense order-64 pencil with a spread spectrum. In 192 further runs, again with byte-identical
artifacts, the one-thread coefficient stage fell from 10.9 to 1.6 s and from 0.50 to 0.098 s on
those two pencils, where FLINT's exact solve and characteristic polynomial take 4.25 s and 0.30 s;
on 32 threads it fell from 0.55 to 0.13 s and from 0.032 to 0.013 s. The dense mass of order 256
went from 15.6 to 1.9 s on one thread. Sparse pencils, whose reduction already skipped zeros,
changed least: the band-3 rotation of order 256 from 0.14 to 0.10 s. The symmetric routes then
cut the work per prime again on symmetric pencils, still without changing a residue or the number
of primes. On one thread a prime took 11.9 ms instead of 17.8 ms on the dense mass of order 256,
0.38 ms instead of 1.30 ms on BCSSTK05, 0.26 ms instead of 0.43 ms on the order-64 pencil and
0.40 ms instead of 1.02 ms on the dense separated pencil of order 100, where C has short rows. The
diagonal-mass pencil of order 256 gained least, 6.8 ms instead of 7.7 ms: Lanczos's products
with a dense F are as many as the Hessenberg form and its recurrence take, only more regular. In
192 further runs, again with byte-identical artifacts, the one-thread coefficient stage fell from
0.233 to 0.075 s on BCSSTK05, from 1.91 to 1.26 s on the dense mass, from 0.101 to 0.064 s on the
order-64 pencil and from 1.70 to 1.49 s on the diagonal-mass pencil; on 32 threads from 0.029 to
0.017, 0.146 to 0.107, 0.014 to 0.013 and 0.140 to 0.119 s. The Laplacian and the rotation pencil
keep the general route; deciding that costs the rotation pencil about 0.6 ms once, 1-3% of its
one-thread stage. The factored route then took a prime of the band-1 consistent-mass pencil of
order 256 in 1.04 ms instead of 6.6 ms, of order 64 in 0.069 ms instead of 0.147 ms, and of random
bands of order 256 in 1.0-2.0 ms instead of 6.6-6.8 ms for half-bandwidths 1-8 and of 2D grids of
order 72-240 1.8-4.0 times faster, again with identical residues. Of the X5 pencils only cluster32lu
takes it: 0.023 ms instead of 0.029 ms a prime, and its one-thread coefficient stage 0.74 ms instead
of 1.06 ms, with byte-identical artifacts in 96 runs on the eight pencils with E not diagonal.
cluster32 and the A-mid4 pencils, whose C is diagonal, keep the congruence, and pay up to 2% a prime
for the decision and code placement. The n+1 Bareiss determinants instead cost
O(n^4) integer operations on entries as long as the coefficients.
Each Graeffe squaring costs O(n^2) exponent comparisons and one working-precision product,
shift and addition per term that is not negligible. Recovering the candidates after r
squarings costs r square roots per root. Exact inertia is O(n^3) arithmetic per endpoint,
normally two endpoints for each root, on short floating-point and integer entries when the
congruence certificate passes and on growing Bareiss entries when elimination decides. The roots
of the inertia proof run in parallel. A polynomial sign costs one
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
`coefficient_method` is `recurrence`, `crt` or `both`; `coefficient_moduli`, `crt_workers`,
`crt_bound_bits` and `crt_reduced` describe Chinese remaindering, and `times_seconds.crt` is its
share of `times_seconds.coefficients`.
`content_removed` tells whether the root stages used q/content(q).
