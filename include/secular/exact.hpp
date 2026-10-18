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
    Inertia inertia(const ZMatrix& a);
    Inertia inertia(const QMatrix& a);
    struct IntegerPencil {
        ZMatrix E,F;
        Z scale;
    };
    IntegerPencil integer_pencil(const QMatrix& stiffness,const QMatrix& mass);
    Inertia shifted_inertia(const ZMatrix& f,const ZMatrix& e,const Q& shift);
    struct CoefficientResult {
        std::vector<Z> coefficients;
        // descending powers, actual signs: det(zE-F)
        std::size_t exact_divisions {},max_coefficient_bits {},max_workspace_bits {};
        bool terminal_zero {};
        // Threads that computed the rows of each recurrence stage.
        unsigned workers {};
    };
    // Workers as for oracle_worker_count, clamped to the n rows of a stage.
    // Every worker count gives identical results.
    CoefficientResult characteristic(const ZMatrix& e,const ZMatrix& f,unsigned workers=1);
    // Output [k][j] is coefficient of z^(n-k) t^j in det(zE-F0-tF1).
    std::vector<std::vector<Z>> affine_characteristic(const ZMatrix& e,const ZMatrix& f0,const ZMatrix& f1);
    Z evaluate(const std::vector<Z>& coefficients,const Z& x);
    Q evaluate(const std::vector<Z>& coefficients,const Q& x);
    // Exact sign at a rational point, using only integer arithmetic.
    int polynomial_sign(const std::vector<Z>& coefficients,const Q& x);
    struct SquarefreeResult {
        bool repeated=false;
        // Nonzero only for a degree-preserving modular squarefree proof.
        unsigned long modular_prime=0;
    };
    SquarefreeResult check_squarefree(const std::vector<Z>& coefficients);
    // 0 selects hardware concurrency, capped at eight; explicit counts are 1..64.
    unsigned oracle_worker_count(unsigned requested,std::size_t jobs);
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
    bool has_repeated_roots(const std::vector<Z>& coefficients);
}
