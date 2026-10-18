#pragma once
#include "secular/exact.hpp"
#include <functional>

namespace secular {
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
        bool coefficient_identity_verified=false;
        double initialization_seconds=0,squarefree_seconds=0;
        unsigned oracle_workers=0,recurrence_workers=0;
        // "modular" or "bareiss" when the coefficient identity was checked, else "none".
        std::string coefficient_oracle="none";
        std::size_t oracle_moduli=0;
        unsigned long squarefree_modular_prime=0;
        std::string certification_method="none";
        // Includes the two positive-definiteness guards as well as shifted checks.
        std::size_t polynomial_sign_evaluations=0,inertia_evaluations=0;
        double coefficients_seconds {},oracle_seconds {},graeffe_seconds {},certification_seconds {},total_seconds {};
    };
    void validate_options(const Options& options);
    SolveResult solve(const QMatrix& stiffness,const QMatrix& mass,const Options& options);
    SolveResult verify_certificates(const QMatrix& stiffness,const QMatrix& mass,
        const std::vector<RootCertificate>& roots,unsigned digits,const Options& options=Options{});
    bool verify_root_certificates(const QMatrix& stiffness,const QMatrix& mass,const std::vector<RootCertificate>& roots,unsigned digits);
    void write_certificates(std::ostream& out,const std::vector<RootCertificate>& roots,unsigned digits);
    std::vector<RootCertificate> read_certificates(std::istream& in,unsigned& digits);
    void write_result_json(std::ostream& out,const SolveResult& r,unsigned digits);
}
