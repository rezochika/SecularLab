#pragma once
#include "secular/exact.hpp"
#include <functional>

namespace secular {
    enum class CoefficientMethod {
        // The adjugate recurrence (on the reduced pencil), then the coefficient oracle.
        recurrence,
        // crt_characteristic: Chinese remaindering of the oracle's residues, a proof by itself.
        crt,
        // Both, which must give the same integers (an N-version check); the proof is the CRT one.
        both
    };
    struct Options {
        unsigned digits=50;
        unsigned precision=256;
        // decimal working precision, candidate stage only
        unsigned max_precision=2048;
        unsigned max_squarings=128;
        bool full_coefficient_check=true;
        bool check_squarefree=true;
        // Exact parallel work: coefficient recurrence rows and oracle jobs.
        unsigned oracle_workers=0;
        CoefficientOracle coefficient_oracle=CoefficientOracle::automatic;
        // crt and both are proofs in themselves, so they require full_coefficient_check; an
        // explicit bareiss oracle still checks their result as well.
        CoefficientMethod coefficient_method=CoefficientMethod::recurrence;
        // Graeffe, squarefree and sign stages use q/content(q): the same roots and, for a positive
        // content, the same signs. coefficients.tsv still holds q.
        bool remove_content=false;
        bool polynomial_certification=true;
        std::function<void(const std::string&)> progress;
    };
    struct RootCertificate {
        std::size_t index {};
        // 0 = largest eigenvalue
        Q lower,upper,radius;
        std::string decimal;
    };
    struct SolveResult {
        bool certified=false;
        std::string status="NOT_CERTIFIED",message;
        std::vector<RootCertificate> roots;
        std::vector<Z> coefficients;
        Z determinant_scale;
        // det(zE-F)=scale^n det(zM-K)
        std::size_t order {},exact_divisions {},max_coefficient_bits {},max_workspace_bits {};
        unsigned working_digits {},graeffe_steps {},attempts {};
        // Times the candidates of some pass were submitted to exact certification.
        unsigned certification_attempts {};
        // The final Graeffe pass outgrew 64-bit binary exponents and continued with unbounded ones.
        bool graeffe_unbounded_exponents=false;
        bool coefficient_identity_verified=false;
        double initialization_seconds=0,squarefree_seconds=0;
        unsigned oracle_workers=0,recurrence_workers=0;
        // "modular" or "bareiss" when the coefficient identity was checked, "crt" when Chinese
        // remaindering proved it, else "none".
        std::string coefficient_oracle="none";
        std::size_t oracle_moduli=0;
        // "recurrence", "crt" or "both"; the primes Chinese remaindering combined and its threads.
        std::string coefficient_method="recurrence";
        std::size_t coefficient_moduli=0;
        unsigned crt_workers=0;
        // The positive content of q and the longest coefficient of q/content.
        Z coefficient_content=1;
        std::size_t primitive_coefficient_bits_max=0;
        // The Graeffe, squarefree and sign stages used q/content.
        bool content_removed=false;
        // Part of coefficients_seconds spent in Chinese remaindering.
        double crt_seconds=0;
        unsigned long squarefree_modular_prime=0;
        std::string certification_method="none";
        // Includes the two positive-definiteness guards as well as shifted checks.
        std::size_t polynomial_sign_evaluations=0,inertia_evaluations=0;
        double coefficients_seconds {},oracle_seconds {},graeffe_seconds {},certification_seconds {},total_seconds {};
    };
    struct GraeffeCandidates {
        // Exact values of the floating candidates, largest first.
        std::vector<Q> roots;
        // The pass outgrew its narrow binary exponents and continued with unbounded ones.
        bool unbounded_exponents=false;
    };
    // The UNTRUSTED candidate stage alone, for tests and diagnostics: `squarings` root-squaring
    // steps of descending coefficients with the positive-root sign pattern, at `work` decimal
    // digits. It sets the process-global floating default precision, as solve does, and throws
    // std::runtime_error when the pass loses the sign pattern. narrow_exponent_bits below 63
    // only moves the switch to unbounded exponents earlier; the candidates never depend on it.
    GraeffeCandidates graeffe_candidates(const std::vector<Z>& coefficients,unsigned work,unsigned squarings,
                                         unsigned narrow_exponent_bits=63);
    void validate_options(const Options& options);
    // The coefficient phase of solve alone: q = det(zE-F) by options.coefficient_method and, with
    // `check`, its identity proof. Fills the coefficient, oracle and content fields of the result;
    // throws when a check fails.
    SolveResult exact_coefficients(const IntegerPencil& pencil,const Options& options,bool check);
    SolveResult solve(const QMatrix& stiffness,const QMatrix& mass,const Options& options);
    SolveResult verify_certificates(const QMatrix& stiffness,const QMatrix& mass,
        const std::vector<RootCertificate>& roots,unsigned digits,const Options& options=Options{});
    bool verify_root_certificates(const QMatrix& stiffness,const QMatrix& mass,const std::vector<RootCertificate>& roots,unsigned digits);
    void write_certificates(std::ostream& out,const std::vector<RootCertificate>& roots,unsigned digits);
    std::vector<RootCertificate> read_certificates(std::istream& in,unsigned& digits);
    void write_result_json(std::ostream& out,const SolveResult& r,unsigned digits);
}
