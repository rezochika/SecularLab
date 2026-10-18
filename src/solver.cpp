#include "secular/solver.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <optional>
#include <limits>

namespace secular {
    using Real=boost::multiprecision::mpf_float;
    using Clock=std::chrono::steady_clock;
    namespace {
        double elapsed(Clock::time_point a) {
            return std::chrono::duration<double>(Clock::now()-a).count();
        }
        struct SignedLog {
            int sign=0;
            Real value=0;
        };
        struct PrecisionLoss:std::runtime_error {
            using std::runtime_error::runtime_error;
        };
        Q rational(const Real& x) {
            Q q;
            mpq_set_f(q.get_mpq_t(),x.backend().data());
            return q;
        }
        Q qabs(Q q) {
            return q<0?-q:q;
        }
        SignedLog plus(const SignedLog& x,const SignedLog& y,const Real& cutoff) {
            if(x.sign==0)return y;
            if(y.sign==0)return x;
            const auto* hi=&x;
            const auto* lo=&y;
            if(x.value<y.value)std::swap(hi,lo);
            Real diff=lo->value-hi->value;
            // Only the UNTRUSTED candidate stage drops negligible terms. Acceptance
            // later proves the complete rootset with verified exact coefficients
            // or recomputes exact inertia of the original rational pencil.
            if(diff<cutoff)return *hi;
            Real t=exp(diff);
            if(hi->sign==lo->sign)return {
                hi->sign,Real(hi->value+log(Real(1)+t))
            };
            if(diff==0)return {};
            Real remainder=Real(1)-t;
            if(remainder<=0)throw PrecisionLoss("cancellation exhausted Graeffe working precision");
            return {
                hi->sign,Real(hi->value+log(remainder))
            };
        }
        std::vector<SignedLog> square_roots(const std::vector<SignedLog>& a,const Real& cutoff) {
            const std::size_t n=a.size()-1;
            std::vector<SignedLog> b(n+1);
            Real ln2=log(Real(2));
            for(std::size_t k=0;k<=n;++k) {
                SignedLog sum {
                    (k%2)?-1:1,Real(2*a[k].value)
                };
                std::size_t start=2*k>n?2*k-n:0;
                // Paired terms of (-1)^n f(z)f(-z), equation (77) of OPM.
                for(std::size_t j=start;j<k;++j) {
                    auto h=2*k-j;
                    SignedLog term {
                        ((j%2)?-1:1)*a[j].sign*a[h].sign,Real(a[j].value+a[h].value+ln2)
                    };
                    sum=plus(sum,term,cutoff);
                }
                const int expected=(k%2)?-1:1;
                if(sum.sign!=expected)throw PrecisionLoss("Graeffe coefficient lost the positive-root sign pattern");
                b[k]=std::move(sum);
            }
            Real shift=b[0].value;
            for(auto& c:b)c.value-=shift;
            return b;
        }
        std::vector<Real> recover(const std::vector<SignedLog>& a,unsigned steps) {
            std::vector<Real> roots;
            roots.reserve(a.size()-1);
            for(std::size_t i=1;i<a.size();++i) {
                Real v=a[i].value-a[i-1].value;
                mpf_div_2exp(v.backend().data(),v.backend().data(),steps);
                roots.emplace_back(exp(v));
            }
            return roots;
        }
        long decimal_scale(const Q& positive) {
            // If this returns e, 10^(e-1) < positive < 10^(e+1).
            // The estimate only sets UNTRUSTED proposal widths, never acceptance.
            return static_cast<long>(positive.get_num().get_str().size())-
                   static_cast<long>(positive.get_den().get_str().size());
        }
        int counted_sign(const std::vector<Z>& coefficients,const Q& x,SolveResult& result) {
            ++result.polynomial_sign_evaluations;
            return polynomial_sign(coefficients,x);
        }
        Inertia counted_inertia(const IntegerPencil& p,const Q& x,SolveResult& result) {
            ++result.inertia_evaluations;
            return shifted_inertia(p.F,p.E,x);
        }
        // All proof paths share the complete metadata and decimal-rounding gate.
        bool valid_rows(const std::vector<RootCertificate>& roots,std::size_t n,unsigned digits) {
            if(n==0||roots.size()!=n||digits<1||digits>1000)return false;
            for(std::size_t i=0;i<n;++i) {
                const auto& r=roots[i];
                if(r.index!=i||r.lower>r.upper||r.radius<0||(i&&r.upper>=roots[i-1].lower))return false;
                Q printed=parse_rational(r.decimal);
                if(format_decimal(round_scaled(printed,digits),digits)!=r.decimal)return false;
                if(format_decimal(round_scaled(r.lower,digits),digits)!=r.decimal||format_decimal(round_scaled(r.upper,digits),digits)!=r.decimal)return false;
                if(r.radius<qabs(printed-r.lower)||r.radius<qabs(printed-r.upper))return false;
            }
            return true;
        }
        bool polynomial_rootset(const std::vector<RootCertificate>& roots,const std::vector<Z>& co,
                                std::size_t n,unsigned digits,SolveResult& result) {
            // n disjoint brackets/singletons each contain a distinct real zero.
            // Only the verified original degree-n characteristic polynomial
            // makes this a proof of the entire spectrum and its sorted indices.
            if(!result.coefficient_identity_verified||co.size()!=n+1||co.front()==0||!valid_rows(roots,n,digits))return false;
            for(const auto& r:roots) {
                int lo=counted_sign(co,r.lower,result);
                if(r.lower==r.upper) {
                    if(lo!=0)return false;
                }
                else {
                    int hi=counted_sign(co,r.upper,result);
                    if(lo==0||hi==0||lo==hi)return false;
                }
            }
            return true;
        }
        bool inertia_rootset(const IntegerPencil& p,const std::vector<RootCertificate>& roots,
                             unsigned digits,SolveResult& result) {
            if(!valid_rows(roots,p.E.n,digits))return false;
            for(std::size_t i=0;i<roots.size();++i) {
                const auto& r=roots[i];
                auto l=counted_inertia(p,r.lower,result);
                if(r.lower==r.upper) {
                    if(l.zero!=1||l.positive!=i)return false;
                }
                else {
                    auto u=counted_inertia(p,r.upper,result);
                    if(l.zero||u.zero||l.positive!=i+1||u.positive!=i)return false;
                }
            }
            return true;
        }
        bool valid_matrix_shapes(const QMatrix& k,const QMatrix& m) {
            return k.n>0&&k.n<=max_dimension&&k.n==m.n&&k.data.size()==k.n*k.n&&m.data.size()==m.n*m.n;
        }
        void compute_coefficients(const IntegerPencil& p,const Options& o,SolveResult& out,bool verify) {
            if(o.progress)o.progress("Computing exact characteristic coefficients...");
            auto timer=Clock::now();
            auto c=characteristic(p.E,p.F,o.oracle_workers);
            out.coefficients_seconds=elapsed(timer);
            out.coefficients=std::move(c.coefficients);
            out.exact_divisions=c.exact_divisions;
            out.max_coefficient_bits=c.max_coefficient_bits;
            out.max_workspace_bits=c.max_workspace_bits;
            out.recurrence_workers=c.workers;
            if(verify) {
                if(o.progress)o.progress(o.coefficient_oracle==CoefficientOracle::bareiss?
                    "Checking polynomial at n+1 independent Bareiss determinants...":
                    "Checking polynomial against det(zE-F) modulo proven primes...");
                timer=Clock::now();
                const auto report=check_coefficients(p.E,p.F,out.coefficients,o.oracle_workers,o.coefficient_oracle);
                out.oracle_seconds=elapsed(timer);
                out.oracle_workers=report.workers;
                out.coefficient_oracle=report.method;
                out.oracle_moduli=report.moduli;
                out.coefficient_identity_verified=report.verified;
                if(!out.coefficient_identity_verified)throw std::logic_error("independent determinant verification failed");
            }
        }
        bool repeated_spectrum(SolveResult& out) {
            auto timer=Clock::now();
            auto squarefree=check_squarefree(out.coefficients);
            out.squarefree_seconds=elapsed(timer);
            out.squarefree_modular_prime=squarefree.modular_prime;
            return squarefree.repeated;
        }
        std::optional<RootCertificate> certify_one(const IntegerPencil& pencil,const std::vector<Z>& coefficients,const Q& approx,const Q& separation,std::size_t index,unsigned digits,unsigned work,bool polynomial,SolveResult& result) {
            if(approx<=0)return std::nullopt;
            const long scale=decimal_scale(approx);
            long places=std::max(static_cast<long>(digits)+8,9-scale);
            if(separation>0)places=std::max(places,2-decimal_scale(separation));
            // Resolve small roots and close neighbors independently of output
            // decimal places. Working precision bounds adaptive refinement, but
            // retain the original output-driven widths: an exact candidate can
            // certify more decimal places than its floating working precision.
            const long available=std::max(static_cast<long>(digits)+32,
                                          static_cast<long>(work)-scale-16);
            places=std::min(places,available);
            long last=-1;
            for(unsigned extra:{0U,8U,24U}) {
                const long width_places=std::min(places+static_cast<long>(extra),available);
                if(width_places==last)continue;
                last=width_places;
                // pow10 has an explicit resource limit; failure to propose an
                // interval within it must remain an uncertified candidate.
                if(width_places>99992)return std::nullopt;
                const auto width=static_cast<unsigned>(width_places);
                Q center(round_scaled(approx,width+8),pow10(width+8));
                center.canonicalize();
                Q radius(1,pow10(width));
                Q lower=center-radius,upper=center+radius;
                if(polynomial) {
                    int lo=counted_sign(coefficients,lower,result),hi=counted_sign(coefficients,upper,result);
                    if(lo==0)upper=lower;
                    else if(hi==0)lower=upper;
                    else if(lo==hi)continue;
                }
                else {
                    auto li=counted_inertia(pencil,lower,result),ui=counted_inertia(pencil,upper,result);
                    if(li.zero==1&&li.positive==index)upper=lower;
                    else if(ui.zero==1&&ui.positive==index)lower=upper;
                    else if(li.zero||ui.zero||li.positive!=index+1||ui.positive!=index)continue;
                }
                auto lo=round_scaled(lower,digits),hi=round_scaled(upper,digits);
                if(lo!=hi) {
                    if(hi-lo!=1)continue;
                    Q midpoint(lo+hi,2*pow10(digits));
                    midpoint.canonicalize();
                    if(counted_sign(coefficients,midpoint,result)!=0)continue;
                    if(!polynomial) {
                        auto mi=counted_inertia(pencil,midpoint,result);
                        if(mi.zero!=1||mi.positive!=index)continue;
                    }
                    lower=upper=midpoint;
                    lo=round_scaled(midpoint,digits);
                }
                const auto decimal=format_decimal(lo,digits);
                Q printed=parse_rational(decimal);
                Q rad=std::max(qabs(printed-lower),qabs(printed-upper));
                return RootCertificate {
                    index,lower,upper,rad,decimal
                };
            }
            return std::nullopt;
        }
        std::vector<RootCertificate> certify_all(const IntegerPencil& p,const std::vector<Z>& co,const std::vector<Real>& candidates,unsigned digits,unsigned work,bool polynomial,SolveResult& result) {
            std::vector<Q> approximations;
            for(const auto& candidate:candidates)approximations.push_back(rational(candidate));
            std::vector<RootCertificate> out;
            for(std::size_t i=0;i<candidates.size();++i) {
                if(approximations[i]<=0||(i&&approximations[i]>=approximations[i-1]))return {};
                Q separation=0;
                if(i)separation=approximations[i-1]-approximations[i];
                if(i+1<approximations.size()) {
                    Q gap=approximations[i]-approximations[i+1];
                    if(gap<=0)return {};
                    if(separation==0||gap<separation)separation=gap;
                }
                auto r=certify_one(p,co,approximations[i],separation,i,digits,work,polynomial,result);
                if(!r)return {};
                out.push_back(std::move(*r));
            }
            if(!valid_rows(out,p.E.n,digits))return {};
            if(polynomial&&!polynomial_rootset(out,co,p.E.n,digits,result))return {};
            return out;
        }
        std::string escape(const std::string& s) {
            std::string o;
            for(unsigned char c:s) {
                switch(c) {
                    case '\\':o+="\\\\";
                    break;
                    case '"':o+="\\\"";
                    break;
                    case '\n':o+="\\n";
                    break;
                    case '\r':o+="\\r";
                    break;
                    case '\t':o+="\\t";
                    break;
                    default:if(c<32)o+='?';
                    else o+=static_cast<char>(c);
                }
            }
            return o;
        }
    }
    void validate_options(const Options& o) {
        if(o.digits<1||o.digits>1000||o.precision<32||o.max_precision<o.precision||o.max_precision>20000||o.max_squarings<1||o.max_squarings>512)throw std::invalid_argument("invalid resource limits (digits 1..1000, precision 32..20000, squarings 1..512)");
        if(o.oracle_workers>64)throw std::invalid_argument("oracle workers must be 0..64");
    }
    SolveResult solve(const QMatrix& k,const QMatrix& m,const Options& o) {
        auto begin=Clock::now();
        if(!valid_matrix_shapes(k,m))throw std::invalid_argument("invalid or mismatched K and M dimensions");
        validate_options(o);
        if(!symmetric(k)||!symmetric(m))throw std::invalid_argument("root recovery requires symmetric K and M");
        auto p=integer_pencil(k,m);
        SolveResult out;
        out.order=k.n;
        out.determinant_scale=p.scale;
        ++out.inertia_evaluations;
        if(inertia(p.E).positive!=k.n)throw std::invalid_argument("root recovery requires positive definite M");
        ++out.inertia_evaluations;
        if(inertia(p.F).positive!=k.n)throw std::invalid_argument("root recovery currently requires positive definite K (positive eigenvalues)");
        out.initialization_seconds=elapsed(begin);
        compute_coefficients(p,o,out,o.full_coefficient_check);
        if(o.check_squarefree&&repeated_spectrum(out)) {
            out.status="REPEATED_SPECTRUM";
            out.message="Exact polynomial gcd detects a repeated root; simple-root Graeffe certification is not applicable.";
            out.total_seconds=elapsed(begin);
            return out;
        }
        for(unsigned work=o.precision;;) {
            ++out.attempts;
            out.working_digits=work;
            Real::default_precision(work);
            if(o.progress)o.progress("Graeffe candidate pass at "+std::to_string(work)+" decimal working digits...");
            auto gt=Clock::now();
            try {
                std::vector<SignedLog> coeff;
                Real lead=log(Real(out.coefficients[0].get_str()));
                for(std::size_t i=0;i<out.coefficients.size();++i) {
                    Z z=out.coefficients[i];
                    int s=mpz_sgn(z.get_mpz_t());
                    if(s!=((i%2)?-1:1))throw std::logic_error("exact positive-root coefficient sign invariant failed");
                    if(z<0)z=-z;
                    coeff.push_back( {
                        s,Real(log(Real(z.get_str()))-lead)
                    });
                }
                const Real cutoff=-Real(work+10)*log(Real(10));
                const Real threshold("1e-"+std::to_string(o.digits+16));
                std::vector<Real> previous;
                unsigned next_certification=0,certification_delay=2;
                for(unsigned step=1;step<=o.max_squarings;++step) {
                    coeff=square_roots(coeff,cutoff);
                    out.graeffe_steps=step;
                    if(step%2&&step!=o.max_squarings)continue;
                    auto candidates=recover(coeff,step);
                    bool stable=previous.size()==candidates.size();
                    for(std::size_t i=0;i<candidates.size();++i) {
                        if(candidates[i]<=0||(i&&candidates[i]>=candidates[i-1]))stable=false;
                        if(stable&&abs(candidates[i]-previous[i])>threshold)stable=false;
                    }
                    previous=candidates;
                    if((stable&&step>=next_certification)||step==o.max_squarings) {
                        out.graeffe_seconds+=elapsed(gt);
                        const bool try_polynomial=o.polynomial_certification&&out.coefficient_identity_verified;
                        if(o.progress)o.progress(std::string("Testing ")+(try_polynomial?"exact full-spectrum polynomial proof":"exact inertia")+" and decimal rounding after "+std::to_string(step)+" squarings...");
                        auto ct=Clock::now();
                        std::vector<RootCertificate> roots;
                        std::string method="exact_inertia";
                        if(try_polynomial) {
                            roots=certify_all(p,out.coefficients,candidates,o.digits,work,true,out);
                            if(roots.size()==k.n)method="polynomial";
                        }
                        if(roots.size()!=k.n)roots=certify_all(p,out.coefficients,candidates,o.digits,work,false,out);
                        out.certification_seconds+=elapsed(ct);
                        if(roots.size()==k.n) {
                            out.roots=std::move(roots);
                            out.certified=true;
                            out.status="CERTIFIED";
                            out.certification_method=method;
                            out.message="All eigenvalues isolated, indexed and correctly rounded by "+
                                std::string(method=="polynomial"?"an exact full-spectrum polynomial proof.":"exact inertia.");
                            out.total_seconds=elapsed(begin);
                            return out;
                        }
                        gt=Clock::now();
                        // Absolute candidate stability does not imply Graeffe
                        // separation is complete (especially for tiny/clustered
                        // roots). Continue squaring; geometrically space failed
                        // certification attempts to bound costly exact checks.
                        next_certification=step+certification_delay;
                        certification_delay*=2;
                    }
                }
                out.graeffe_seconds+=elapsed(gt);
            }
            catch(const PrecisionLoss& ex) {
                out.graeffe_seconds+=elapsed(gt);
                out.message=ex.what();
                if(o.progress)o.progress(out.message);
            }
            if(work>=o.max_precision)break;
            work=std::min(o.max_precision,work*2);
        }
        out.status="NOT_CERTIFIED";
        if(out.message.empty())out.message="Configured Graeffe/precision limits did not establish all requested output digits.";
        out.total_seconds=elapsed(begin);
        return out;
    }
    SolveResult verify_certificates(const QMatrix& k,const QMatrix& m,const std::vector<RootCertificate>& roots,unsigned digits,const Options& o) {
        const auto begin=Clock::now();
        SolveResult out;
        out.status="FAILED";
        out.order=k.n;
        try {
            validate_options(o);
            if(digits<1||digits>1000||!valid_matrix_shapes(k,m)||!symmetric(k)||!symmetric(m))
                throw std::invalid_argument("invalid certificate precision or matrix dimensions/symmetry");
            auto p=integer_pencil(k,m);
            out.determinant_scale=p.scale;
            ++out.inertia_evaluations;
            if(inertia(p.E).positive!=k.n)throw std::invalid_argument("verification requires positive definite M");
            ++out.inertia_evaluations;
            if(inertia(p.F).positive!=k.n)throw std::invalid_argument("verification requires positive definite K");
            out.initialization_seconds=elapsed(begin);
            if(o.polynomial_certification) {
                // Verification starts with the original matrices; neither saved
                // coefficient files nor a previous status can authorize a proof.
                compute_coefficients(p,o,out,true);
                if(repeated_spectrum(out))throw std::invalid_argument("repeated spectrum cannot have a complete distinct-root certificate");
            }
            auto timer=Clock::now();
            bool accepted=false;
            try {
                if(o.polynomial_certification&&polynomial_rootset(roots,out.coefficients,k.n,digits,out)) {
                    accepted=true;
                    out.certification_method="polynomial";
                }
                else if(inertia_rootset(p,roots,digits,out)) {
                    accepted=true;
                    out.certification_method="exact_inertia";
                }
            }
            catch(...) {
                out.certification_seconds=elapsed(timer);
                throw;
            }
            out.certification_seconds=elapsed(timer);
            if(accepted) {
                out.certified=true;
                out.status="VERIFIED";
                out.roots=roots;
                out.message="All certificate rows verified against the original matrices by "+
                    std::string(out.certification_method=="polynomial"?"an exact full-spectrum polynomial proof.":"exact inertia.");
            }
            else out.message="Certificate rows failed exact root, ordering, radius or decimal checks.";
        }
        catch(const std::exception& ex) {
            out.message=ex.what();
            if(out.initialization_seconds==0)out.initialization_seconds=elapsed(begin);
        }
        out.total_seconds=elapsed(begin);
        return out;
    }
    bool verify_root_certificates(const QMatrix& k,const QMatrix& m,const std::vector<RootCertificate>& roots,unsigned digits) {
        return verify_certificates(k,m,roots,digits).certified;
    }
    void write_certificates(std::ostream& out,const std::vector<RootCertificate>& roots,unsigned digits) {
        out<<"# SecularLab certificates v1\n# decimal_places="<<digits<<"\nindex\tdecimal\tlower\tupper\tradius\n";
        for(const auto& r:roots)out<<r.index<<'\t'<<r.decimal<<'\t'<<r.lower<<'\t'<<r.upper<<'\t'<<r.radius<<'\n';
        if(!out)throw std::runtime_error("certificate output failed");
    }
    std::vector<RootCertificate> read_certificates(std::istream& in,unsigned& digits) {
        std::string line;
        digits=0;
        std::vector<RootCertificate> r;
        bool header=false;
        while(std::getline(in,line)) {
            if(!line.empty()&&line.back()=='\r')line.pop_back();
            if(line.rfind("# decimal_places=",0)==0) {
                auto x=parse_rational(line.substr(17));
                if(x.get_den()!=1||x<1||x>1000)throw std::invalid_argument("invalid certificate precision");
                digits=static_cast<unsigned>(x.get_num().get_ui());
                continue;
            }
            if(line.empty()||line[0]=='#')continue;
            if(!header) {
                if(line!="index\tdecimal\tlower\tupper\tradius")throw std::invalid_argument("invalid certificate columns");
                header=true;
                continue;
            }
            std::istringstream row(line);
            std::vector<std::string> fields;
            std::string s;
            while(std::getline(row,s,'\t'))fields.push_back(s);
            if(fields.size()!=5)throw std::invalid_argument("invalid certificate row");
            auto ix=parse_rational(fields[0]);
            if(ix.get_den()!=1||ix<0||ix>static_cast<unsigned long>(max_dimension))throw std::invalid_argument("invalid certificate index");
            if(r.size()>=max_dimension)throw std::invalid_argument("too many certificate rows");
            r.push_back( {
                static_cast<std::size_t>(ix.get_num().get_ui()),parse_rational(fields[2]),parse_rational(fields[3]),parse_rational(fields[4]),fields[1]
            });
        }
        if(!header||!digits)throw std::invalid_argument("missing certificate metadata");
        return r;
    }
    void write_result_json(std::ostream& out,const SolveResult& r,unsigned digits) {
        out<<std::setprecision(12)<<"{\n  \"schema\": \"SecularLab-1\",\n"
           <<"  \"status\": \""<<r.status<<"\",\n"
           <<"  \"message\": \""<<escape(r.message)<<"\",\n"
           <<"  \"certified\": "<<(r.certified?"true":"false")<<",\n"
           <<"  \"certification_method\": \""<<r.certification_method<<"\",\n"
           <<"  \"order\": "<<r.order<<",\n"
           <<"  \"decimal_places\": "<<digits<<",\n"
           <<"  \"working_decimal_digits\": "<<r.working_digits<<",\n"
           <<"  \"graeffe_squarings\": "<<r.graeffe_steps<<",\n"
           <<"  \"precision_attempts\": "<<r.attempts<<",\n"
           <<"  \"exact_divisions\": "<<r.exact_divisions<<",\n"
           <<"  \"coefficient_bits_max\": "<<r.max_coefficient_bits<<",\n"
           <<"  \"recurrence_workspace_bits_max\": "<<r.max_workspace_bits<<",\n"
           <<"  \"coefficient_identity_verified\": "<<(r.coefficient_identity_verified?"true":"false")<<",\n"
           <<"  \"oracle_workers\": "<<r.oracle_workers<<",\n"
           <<"  \"recurrence_workers\": "<<r.recurrence_workers<<",\n"
           <<"  \"coefficient_oracle\": \""<<r.coefficient_oracle<<"\",\n"
           <<"  \"oracle_moduli\": "<<r.oracle_moduli<<",\n"
           <<"  \"squarefree_modular_prime\": "<<r.squarefree_modular_prime<<",\n"
           <<"  \"polynomial_sign_evaluations\": "<<r.polynomial_sign_evaluations<<",\n"
           <<"  \"inertia_evaluations\": "<<r.inertia_evaluations<<",\n"
           <<"  \"certified_eigenvalues\": "<<r.roots.size()<<",\n"
           <<"  \"times_seconds\": {\"initialization\": "<<r.initialization_seconds
           <<", \"coefficients\": "<<r.coefficients_seconds
           <<", \"coefficient_oracle\": "<<r.oracle_seconds
           <<", \"squarefree\": "<<r.squarefree_seconds
           <<", \"graeffe\": "<<r.graeffe_seconds
           <<", \"certification\": "<<r.certification_seconds
           <<", \"total\": "<<r.total_seconds<<"},\n"
           <<"  \"denominator_clearing_multiplier\": \""<<r.determinant_scale<<"\",\n"
           <<"  \"gmp_version\": \""<<gmp_version<<"\",\n"
           <<"  \"candidate_backend\": \"Boost.Multiprecision GMP float (untrusted)\"\n}\n";
        if(!out)throw std::runtime_error("JSON output failed");
    }
}
