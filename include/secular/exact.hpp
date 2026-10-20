#pragma once
#include <gmpxx.h>
#include <vector>
#include <string>
#include <istream>
#include <ostream>
#include <filesystem>
#include <stdexcept>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <optional>

namespace secular {
    using Z=mpz_class;
    using Q=mpq_class;
    constexpr std::size_t max_dimension=256;
    template<class T> struct Matrix {
        std::size_t n {};
        std::vector<T> data;
        Matrix()=default;
        explicit Matrix(std::size_t dimension):n(dimension) {
            if(n==0 || n>max_dimension) throw std::invalid_argument("dimension must be 1..256");
            data.resize(n*n);
        }
        T& operator()(std::size_t i,std::size_t j) {
            return data[i*n+j];
        }
        const T& operator()(std::size_t i,std::size_t j)const {
            return data[i*n+j];
        }
    };
    using ZMatrix=Matrix<Z>;
    using QMatrix=Matrix<Q>;
    template<class T> Matrix<T> identity(std::size_t n) {
        Matrix<T> a(n);
        for(std::size_t i=0;i<n;++i)a(i,i)=1;
        return a;
    }
    template<class T> bool symmetric(const Matrix<T>& a) {
        for(std::size_t i=0;i<a.n;++i)for(std::size_t j=0;j<i;++j)if(a(i,j)!=a(j,i))return false;
        return true;
    }
    template<class T> bool diagonal(const Matrix<T>& a) {
        for(std::size_t i=0;i<a.n;++i)for(std::size_t j=0;j<a.n;++j)if(i!=j && a(i,j)!=0)return false;
        return true;
    }
    Q parse_rational(const std::string& text);
    QMatrix read_matrix_market(std::istream& input);
    QMatrix read_matrix_market(const std::filesystem::path& file);
    void write_matrix_market(std::ostream& out,const QMatrix& a);
    Z pow10(unsigned digits);
    Z round_scaled(const Q& x,unsigned digits);
    std::string format_decimal(const Z& scaled,unsigned digits);
    Q legendre_value(unsigned degree,const Q& x);
    std::size_t bits(const Z& z);
    Z exact_divide(const Z& a,const Z& b);
    ZMatrix multiply(const ZMatrix& a,const ZMatrix& b);
    Z determinant(ZMatrix a);
    ZMatrix adjugate(const ZMatrix& a);
    struct Inertia {
        std::size_t positive {},negative {},zero {};
    };
    // A cooperative cancellation check for long exact work: true once its result is no longer
    // needed. Work that finds it true throws Cancelled; an empty check never cancels.
    using CancelCheck=std::function<bool()>;
    struct Cancelled:std::runtime_error {
        Cancelled():std::runtime_error("exact work cancelled") {}
    };
    Inertia inertia(const ZMatrix& a);
    Inertia inertia(const QMatrix& a);
    // Exact congruence certificate for any proposal x. It requires x upper triangular with a
    // nonzero diagonal, so det(x)!=0, and B = x^T*a*x strictly diagonally dominant. Then
    // inertia(a) = inertia(B) is the sign count of diag(B), with nullity zero. Returns nullopt
    // when either test fails, which says nothing about a.
    std::optional<Inertia> certify_congruence(const ZMatrix& a,const ZMatrix& x);
    // Largest proposal precision congruence_inertia accepts.
    constexpr unsigned max_congruence_precision=16384;
    // certify_congruence for x = 2^(precision_bits-8)*L^-T, truncated to integers. L is an
    // untrusted floating-point LDL^T factor of a at precision_bits (16..max_congruence_precision),
    // set per object. A zero floating pivot, or a factor far beyond the scale of a, also gives nullopt.
    std::optional<Inertia> congruence_inertia(const ZMatrix& a,unsigned precision_bits);
    struct InertiaReport {
        Inertia inertia;
        // Exact test that decided: "diagonal" (sign count), "congruence" or "elimination"
        // (the fraction-free inertia).
        std::string method;
        // Proposal precision of a deciding congruence, 0 otherwise.
        unsigned precision_bits=0;
    };
    // The exact inertia of a symmetric matrix: the diagonal sign count, else congruence_inertia at
    // each proposal precision in turn, else inertia(). Every outcome is exact and the same; the
    // precisions only change the work.
    InertiaReport scheduled_inertia(const ZMatrix& a,const std::vector<unsigned>& precisions,
                                    const CancelCheck& cancelled={});
    struct DefinitenessReport {
        bool positive_definite=false;
        // Exact test that decided: "diagonal" (sign count), "congruence" or "elimination"
        // (the fraction-free inertia).
        std::string method;
        // Proposal precision of a deciding congruence, 0 otherwise.
        unsigned precision_bits=0;
    };
    // Exact positive-definiteness test of a symmetric matrix: scheduled_inertia at 64, 128 and
    // 256 bits. Every outcome is exact.
    DefinitenessReport check_positive_definite(const ZMatrix& a);
    struct IntegerPencil {
        ZMatrix E,F;
        Z scale;
    };
    IntegerPencil integer_pencil(const QMatrix& stiffness,const QMatrix& mass);
    // Proposal precisions for the endpoint matrix den(s)*F-num(s)*E: P = 64+bits(den(s)), rounded
    // up to whole 64-bit limbs, then 2P and 4P, dropping those above max_congruence_precision. An
    // endpoint near an eigenvalue makes the matrix nearly singular, and the denominator of a
    // proposal endpoint grows with how near it may be.
    std::vector<unsigned> shifted_precisions(const Q& shift);
    // scheduled_inertia of den(s)*F-num(s)*E at shifted_precisions(s), which is the inertia of F-sE.
    InertiaReport shifted_inertia_report(const ZMatrix& f,const ZMatrix& e,const Q& shift,
                                         const CancelCheck& cancelled={});
    Inertia shifted_inertia(const ZMatrix& f,const ZMatrix& e,const Q& shift,const CancelCheck& cancelled={});
    struct CoefficientResult {
        std::vector<Z> coefficients;
        // descending powers, actual signs: det(zE-F)
        std::size_t exact_divisions {},max_coefficient_bits {},max_workspace_bits {};
        bool terminal_zero {};
        // Threads that computed the rows of each recurrence stage, or the residues of crt_characteristic.
        unsigned workers {};
        // crt_characteristic only: primes whose residues were combined, and primes skipped
        // because they divide det(E'); the bits of the coefficient bound that set their number,
        // and whether the pencil was reduced (E' != E).
        std::size_t moduli {},skipped_moduli {},bound_bits {};
        bool reduced {};
    };
    // Workers as for oracle_worker_count, clamped to the n rows of a stage.
    // Every worker count gives identical results.
    CoefficientResult characteristic(const ZMatrix& e,const ZMatrix& f,unsigned workers=1);
    // det(zE-F) by Chinese remaindering the residues that the modular oracle computes, taken on the
    // reduce_pencil() reduction (E',F'): det(E') det(zI-E'^-1 F') modulo oracle primes not dividing
    // det(E'), until their product M exceeds twice every bound b_m on the reduced coefficients c'_m.
    // b_m is the smaller of the reduced pencil's coefficient bound and floor(B_m/f_m), with B_m the
    // bound of the original pencil and q_m = f_m c'_m. The symmetric residue modulo M is then c'_m
    // itself, so the result is proven as the oracle proves a given vector, and q_m = f_m c'_m by
    // exact multiplication. The same integers as characteristic(); no recurrence statistic is
    // filled. Throws std::invalid_argument for a singular E, as characteristic() does, once the
    // skipped primes prove det(E')=0. Workers as for oracle_worker_count, clamped to the first
    // batch of primes. Requires the modular oracle (64-bit GMP limbs and 128-bit integers);
    // otherwise throws std::logic_error.
    CoefficientResult crt_characteristic(const ZMatrix& e,const ZMatrix& f,unsigned workers=1);
    // The positive gcd of the coefficients, 1 for a zero vector. q/content has the same roots,
    // and the same signs everywhere.
    Z polynomial_content(const std::vector<Z>& coefficients);
    // An exact diagonal congruence and content reduction: E = alpha*P*E'*P and F = beta*P*F'*P,
    // with P = diag(scale) and every scale entry a product of powers of 2 and 5. Then
    // det(zE-F) = prod(scale)^2 * sum_k c'_k alpha^(n-k) beta^k z^(n-k), where c' are the
    // coefficients of det(zE'-F'). Zero entries, symmetry and a diagonal E are preserved.
    struct PencilReduction {
        ZMatrix E,F;
        std::vector<Z> scale;
        Z e_content,f_content;
        // True when P=I and alpha=beta=1, so E',F' equal the input.
        bool trivial=true;
    };
    PencilReduction reduce_pencil(const ZMatrix& e,const ZMatrix& f);
    // characteristic() of the reduced pencil, mapped back exactly to det(zE-F). The coefficients
    // are identical to characteristic(e,f); the workspace statistic measures the reduced pencil.
    CoefficientResult reduced_characteristic(const ZMatrix& e,const ZMatrix& f,unsigned workers=1);
    // Output [k][j] is coefficient of z^(n-k) t^j in det(zE-F0-tF1).
    std::vector<std::vector<Z>> affine_characteristic(const ZMatrix& e,const ZMatrix& f0,const ZMatrix& f1);
    Z evaluate(const std::vector<Z>& coefficients,const Z& x);
    Q evaluate(const std::vector<Z>& coefficients,const Q& x);
    // Exact sign at a rational point, using only integer arithmetic.
    int polynomial_sign(const std::vector<Z>& coefficients,const Q& x);
    // Exact signs at rationals whose canonical denominator divides a fixed q>0. The table
    // T_i=c_i*q^i is formed once; each sign is then a Horner pass over the integer x*q, since
    // sum_i T_i*(x*q)^(n-i) = q^n*P(x). The sign always equals polynomial_sign.
    class ScaledPolynomial {
    public:
        ScaledPolynomial(const std::vector<Z>& coefficients,const Z& denominator);
        const Z& denominator()const {
            return q;
        }
        // True when the denominator of x is positive and divides q.
        bool covers(const Q& x)const;
        // Throws std::invalid_argument unless covers(x).
        int sign(const Q& x)const;
    private:
        Z q;
        std::vector<Z> table;
    };
    struct SquarefreeResult {
        bool repeated=false;
        // Nonzero only for a degree-preserving modular squarefree proof.
        unsigned long modular_prime=0;
    };
    SquarefreeResult check_squarefree(const std::vector<Z>& coefficients);
    // 0 selects the hardware thread count (at most 64); explicit counts are 1..64.
    unsigned oracle_worker_count(unsigned requested,std::size_t jobs);
    struct OrderedRun {
        // The first index whose job returned false or threw; the job count when none did.
        std::size_t stop=0;
        // The exception of job `stop`, when it threw.
        std::exception_ptr error;
    };
    // Runs job(0), job(1), ... on up to `workers` threads, the caller included, handing out indices
    // in increasing order. The outcome is that of the serial loop that stops at the first job
    // returning false or throwing: every job below `stop` ran once and returned true, and job `stop`
    // ran once. Later jobs may run or be skipped, so their effects must be discarded. Jobs must be
    // safe to run concurrently. A job's exception is returned, never rethrown here.
    // job(i,superseded) may stop early once superseded() is true: a lower job has then stopped the
    // run, so job i's result is discarded whatever it is. It never turns true for jobs up to `stop`.
    OrderedRun run_ordered(std::size_t jobs,unsigned workers,
                           const std::function<bool(std::size_t,const CancelCheck&)>& job);
    enum class CoefficientOracle {
        // Multimodular characteristic polynomial; n+1 Bareiss determinants if inconclusive.
        automatic,
        // Always the n+1-point Bareiss determinant check.
        bareiss
    };
    struct OracleReport {
        bool verified=false;
        // Check that decided: "modular", "bareiss", or "none" for rejected storage.
        std::string method="none";
        // Primes whose residues decided a modular proof or refutation.
        std::size_t moduli=0;
        // Threads used by the deciding check.
        unsigned workers=0;
    };
    // Exact proof or refutation of det(zE-F) = sum_k coefficients[k] z^(n-k).
    OracleReport check_coefficients(const ZMatrix& e,const ZMatrix& f,const std::vector<Z>& coefficients,
                                    unsigned workers,CoefficientOracle oracle=CoefficientOracle::automatic);
    bool verify_coefficients(const ZMatrix& e,const ZMatrix& f,const std::vector<Z>& coefficients,unsigned workers=1);
    // The index-th prime of the fixed descending sequence of proven primes below 2^62
    // used by the modular oracle (index 0 is the largest).
    std::uint64_t oracle_prime(std::size_t index);
    // Implementations of the residue kernel behind the modular oracle and crt_characteristic.
    // Every residue is a field computation modulo p, so both give the same residues; `reference`
    // is the plain Gauss-Jordan and Hessenberg reduction kept to test `fast` against.
    enum class ResidueKernel {
        fast,reference
    };
    // det(zE-F) modulo each odd prime p < 2^62 in descending powers, as residues below p, or nullopt
    // when p divides det(E). One workspace serves every prime, as in the oracle. Exposed for tests and
    // benchmarks. Requires the modular oracle.
    std::vector<std::optional<std::vector<std::uint64_t>>> pencil_residues(const ZMatrix& e,const ZMatrix& f,
                                                                           const std::vector<std::uint64_t>& primes,
                                                                           ResidueKernel kernel=ResidueKernel::fast);
    bool has_repeated_roots(const std::vector<Z>& coefficients);
}
