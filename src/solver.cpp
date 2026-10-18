#include "secular/solver.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <chrono>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <map>
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
        // Graeffe candidates in extended-exponent form (UNTRUSTED stage). A coefficient is
        // m*2^e: a working-precision mantissa with 1/2<=|m|<1 and an exact integer exponent, so
        // squaring needs no logarithms or exponentials. Exponents start as checked int64 values.
        // The first overflow converts the iterate to unbounded GMP integers for the rest of the
        // pass, so every value is the one unbounded exponents give, wherever the switch happens.
        struct ExponentOverflow {};
        using Narrow=std::int64_t;
        // Narrow exponents stay within [-limit,limit]; a result outside it throws instead.
        struct NarrowExponents {
            Narrow limit;
            void sum(Narrow& r,Narrow a,Narrow b)const {
                if((b>0&&a>limit-b)||(b<0&&a< -limit-b))throw ExponentOverflow {};
                r=a+b;
            }
            void difference(Narrow& r,Narrow a,Narrow b)const {
                if((b<0&&a>limit+b)||(b>0&&a< -limit+b))throw ExponentOverflow {};
                r=a-b;
            }
            void add(Narrow& r,long v)const {
                sum(r,r,v);
            }
            // The alignment shift if it is at most cut, else -1: the term is negligible.
            long shift(Narrow v,long cut)const {
                return v>cut?-1:static_cast<long>(v);
            }
            bool odd(Narrow v)const {
                return v&1;
            }
            void halve_even(Narrow& v)const {
                v/=2;
            }
            long small(Narrow v)const {
                if(v>std::numeric_limits<long>::max()||v<std::numeric_limits<long>::min())throw ExponentOverflow {};
                return static_cast<long>(v);
            }
        };
        struct WideExponents {
            void sum(Z& r,const Z& a,const Z& b)const {
                mpz_add(r.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
            }
            void difference(Z& r,const Z& a,const Z& b)const {
                mpz_sub(r.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
            }
            void add(Z& r,long v)const {
                if(v>=0)mpz_add_ui(r.get_mpz_t(),r.get_mpz_t(),static_cast<unsigned long>(v));
                else mpz_sub_ui(r.get_mpz_t(),r.get_mpz_t(),0UL-static_cast<unsigned long>(v));
            }
            long shift(const Z& v,long cut)const {
                return mpz_cmp_si(v.get_mpz_t(),cut)>0?-1:mpz_get_si(v.get_mpz_t());
            }
            bool odd(const Z& v)const {
                return mpz_odd_p(v.get_mpz_t());
            }
            void halve_even(Z& v)const {
                mpz_tdiv_q_2exp(v.get_mpz_t(),v.get_mpz_t(),1);
            }
            long small(const Z& v)const {
                if(!mpz_fits_slong_p(v.get_mpz_t()))throw PrecisionLoss("Graeffe candidate exceeds the floating exponent range");
                return mpz_get_si(v.get_mpz_t());
            }
        };
        Z widen(Narrow v) {
            const std::uint64_t magnitude=v<0?0-static_cast<std::uint64_t>(v):static_cast<std::uint64_t>(v);
            Z z;
            mpz_import(z.get_mpz_t(),1,-1,sizeof magnitude,0,0,&magnitude);
            if(v<0)mpz_neg(z.get_mpz_t(),z.get_mpz_t());
            return z;
        }
        // Scales x into 1/2<=|x|<1 and returns the binary exponent removed (0 for x=0).
        long normalize(mpf_t x) {
            if(mpf_sgn(x)==0)return 0;
            long ex=0;
            mpf_get_d_2exp(&ex,x);
            if(ex>0)mpf_div_2exp(x,x,static_cast<mp_bitcnt_t>(ex));
            else if(ex<0)mpf_mul_2exp(x,x,static_cast<mp_bitcnt_t>(-ex));
            return ex;
        }
        void scale_2exp(mpf_t x,long e) {
            if(e>0)mpf_mul_2exp(x,x,static_cast<mp_bitcnt_t>(e));
            else if(e<0)mpf_div_2exp(x,x,static_cast<mp_bitcnt_t>(-e));
        }
        template<class Exponent> struct Iterate {
            std::vector<Real> mantissa;
            std::vector<Exponent> exponent;
        };
        class Graeffe {
        public:
            // Mantissas take Boost's current default precision; cut is the drop threshold in bits.
            // Narrow exponents have at most narrow_bits magnitude bits (63 fills an int64).
            Graeffe(const std::vector<Z>& coefficients,long cut_bits,unsigned narrow_bits=63)
                :cut(cut_bits),narrow_ops {narrow_bits>=63?std::numeric_limits<Narrow>::max():(Narrow(1)<<narrow_bits)-1} {
                const auto size=coefficients.size();
                for(auto* x: {&narrow,&narrow_next}) {
                    x->mantissa.resize(size);
                    x->exponent.resize(size);
                }
                terms.resize(size);
                for(std::size_t i=0;i<size;++i)
                    if(mpz_sgn(coefficients[i].get_mpz_t())!=((i%2)?-1:1))
                        throw std::logic_error("exact positive-root coefficient sign invariant failed");
                try {
                    initialize(narrow,coefficients,narrow_ops);
                }
                catch(const ExponentOverflow&) {
                    // Only a deliberately small narrow limit overflows on the input itself.
                    widen_iterate(false);
                    initialize(broad,coefficients,wide_ops);
                }
            }
            bool wide_exponents()const {
                return wide;
            }
            void square() {
                if(!wide) {
                    try {
                        step(narrow,narrow_next,terms,narrow_scratch,narrow_ops);
                        std::swap(narrow,narrow_next);
                        return;
                    }
                    catch(const ExponentOverflow&) {
                        widen_iterate(true);
                    }
                }
                step(broad,broad_next,wide_terms,wide_scratch,wide_ops);
                std::swap(broad,broad_next);
            }
            // root_i = |a_i/a_(i-1)|^(1/2^steps), largest first, by repeated square roots.
            std::vector<Real> roots(unsigned steps) {
                if(!wide) {
                    try {
                        return recover(narrow,steps,narrow_scratch,narrow_ops);
                    }
                    catch(const ExponentOverflow&) {
                        widen_iterate(true);
                    }
                }
                return recover(broad,steps,wide_scratch,wide_ops);
            }
        private:
            template<class Exponent,class Ops> void initialize(Iterate<Exponent>& a,const std::vector<Z>& c,const Ops& ops) {
                for(std::size_t i=0;i<c.size();++i) {
                    auto* m=a.mantissa[i].backend().data();
                    mpf_set_z(m,c[i].get_mpz_t());
                    a.exponent[i]=0;
                    ops.add(a.exponent[i],normalize(m));
                }
                const Exponent lead=a.exponent[0];
                for(auto& e:a.exponent)ops.difference(e,e,lead);
            }
            void widen_iterate(bool copy) {
                for(auto* x: {&broad,&broad_next}) {
                    x->mantissa.resize(narrow.mantissa.size());
                    x->exponent.resize(narrow.exponent.size());
                }
                wide_terms.resize(terms.size());
                if(copy)for(std::size_t i=0;i<narrow.mantissa.size();++i) {
                    mpf_set(broad.mantissa[i].backend().data(),narrow.mantissa[i].backend().data());
                    broad.exponent[i]=widen(narrow.exponent[i]);
                }
                wide=true;
            }
            // b_k = (-1)^k a_k^2 + 2 sum_{j<k} (-1)^j a_j a_(2k-j), equation 77 of OPM. Every term
            // is aligned to the largest term exponent of its sum, and terms more than cut bits
            // below it are dropped, which only the untrusted candidate stage may do.
            template<class Exponent,class Ops> void step(const Iterate<Exponent>& a,Iterate<Exponent>& b,
                                                         std::vector<Exponent>& term,Exponent& shift,const Ops& ops) {
                const std::size_t n=a.mantissa.size()-1;
                auto* sum=accumulator.backend().data();
                auto* t=product.backend().data();
                for(std::size_t k=0;k<=n;++k) {
                    const std::size_t start=2*k>n?2*k-n:0;
                    // term[k] holds the diagonal exponent, term[j] that of 2 a_j a_(2k-j).
                    ops.sum(term[k],a.exponent[k],a.exponent[k]);
                    std::size_t top=k;
                    for(std::size_t j=start;j<k;++j) {
                        ops.sum(term[j],a.exponent[j],a.exponent[2*k-j]);
                        ops.add(term[j],1);
                        if(term[j]>term[top])top=j;
                    }
                    const Exponent& reference=term[top];
                    mpf_set_ui(sum,0);
                    for(std::size_t j=start;j<=k;++j) {
                        ops.difference(shift,reference,term[j]);
                        const long s=ops.shift(shift,cut);
                        if(s<0)continue;
                        const auto h=2*k-j;
                        mpf_mul(t,a.mantissa[j].backend().data(),a.mantissa[h].backend().data());
                        mpf_div_2exp(t,t,static_cast<mp_bitcnt_t>(s));
                        if(j%2)mpf_sub(sum,sum,t);
                        else mpf_add(sum,sum,t);
                    }
                    if(mpf_sgn(sum)!=((k%2)?-1:1))throw PrecisionLoss("Graeffe coefficient lost the positive-root sign pattern");
                    auto* m=b.mantissa[k].backend().data();
                    mpf_set(m,sum);
                    b.exponent[k]=reference;
                    ops.add(b.exponent[k],normalize(m));
                }
                for(std::size_t k=n+1;k-->0;)ops.difference(b.exponent[k],b.exponent[k],b.exponent[0]);
            }
            template<class Exponent,class Ops> std::vector<Real> recover(const Iterate<Exponent>& a,unsigned steps,
                                                                         Exponent& d,const Ops& ops) {
                std::vector<Real> out(a.mantissa.size()-1);
                for(std::size_t i=1;i<a.mantissa.size();++i) {
                    auto* q=out[i-1].backend().data();
                    mpf_div(q,a.mantissa[i].backend().data(),a.mantissa[i-1].backend().data());
                    mpf_abs(q,q);
                    ops.difference(d,a.exponent[i],a.exponent[i-1]);
                    ops.add(d,normalize(q));
                    // Each round replaces q*2^d by its square root, making d even before halving
                    // it. Every square root halves the relative error carried by q.
                    for(unsigned r=0;r<steps;++r) {
                        if(ops.odd(d)) {
                            mpf_mul_2exp(q,q,1);
                            ops.add(d,-1);
                        }
                        mpf_sqrt(q,q);
                        ops.halve_even(d);
                    }
                    scale_2exp(q,ops.small(d));
                }
                return out;
            }
            long cut;
            NarrowExponents narrow_ops;
            WideExponents wide_ops;
            bool wide=false;
            Iterate<Narrow> narrow,narrow_next;
            Iterate<Z> broad,broad_next;
            std::vector<Narrow> terms;
            std::vector<Z> wide_terms;
            Narrow narrow_scratch=0;
            Z wide_scratch;
            Real accumulator,product;
        };
        // Terms below 10^-(work+10) of their sum's largest term are dropped: the cut in bits is
        // the bit length of 10^(work+10).
        long graeffe_cut(unsigned work) {
            return static_cast<long>(bits(pow10(work+10)));
        }
        // The nearest power of ten: e with 10^(e-1/2) <= positive < 10^(e+1/2), so that
        // 10^(e-1) < positive < 10^(e+1). The estimate only sets UNTRUSTED proposal widths, never
        // acceptance. It is decided by the exact test 10^(2e-1) <= positive^2 < 10^(2e+1), so it
        // depends on the value alone, not on the dyadic form of a converted floating candidate.
        // Its boundaries are irrational, so a rational candidate is never on one, and low-order
        // candidate noise can move it only next to some sqrt(10)*10^k.
        long decimal_scale(const Q& positive) {
            // The digit-count estimate is floor(log10(positive)) or one more; so is the result.
            const long digits=static_cast<long>(positive.get_num().get_str().size())-
                              static_cast<long>(positive.get_den().get_str().size());
            // Far beyond the 10^5 decimal-scale resource limit, the estimate is kept.
            if(digits>1000000||digits< -1000000)return digits;
            const Q square=positive*positive;
            for(long e=digits+1;e>digits-1;--e) {
                const long exponent=2*e-1;
                Z power;
                mpz_ui_pow_ui(power.get_mpz_t(),10,static_cast<unsigned long>(exponent<0?-exponent:exponent));
                Q bound=exponent<0?Q(1,power):Q(power);
                bound.canonicalize();
                if(square>=bound)return e;
            }
            return digits-1;
        }
        // The least m with den | 10^m, when den>0 has no prime factors besides 2 and 5.
        std::optional<unsigned long> decimal_exponent(const Z& den) {
            if(den<=0)return std::nullopt;
            Z rest;
            const auto twos=mpz_scan1(den.get_mpz_t(),0);
            mpz_tdiv_q_2exp(rest.get_mpz_t(),den.get_mpz_t(),twos);
            const Z five=5;
            const auto fives=mpz_remove(rest.get_mpz_t(),rest.get_mpz_t(),five.get_mpz_t());
            if(rest!=1)return std::nullopt;
            return std::max<unsigned long>(twos,fives);
        }
        // Exact signs of one fixed coefficient vector. A repeated point is answered from a memo
        // keyed by its canonical rational, so the final rootset gate reuses the signs certify_one
        // computed. New points with decimal denominators, as every proposal endpoint and rounding
        // midpoint has, share scaled tables at powers of ten; other points use the general
        // evaluator. Neither changes a sign, only the work: the counter records evaluations.
        class SignEvaluator {
        public:
            explicit SignEvaluator(const std::vector<Z>& coefficients):co(coefficients) {
                for(const auto& c:co)coefficient_bits+=bits(c);
            }
            const std::vector<Z>& coefficients()const {
                return co;
            }
            int operator()(const Q& x,SolveResult& result) {
                if(const auto hit=memo.find(x);hit!=memo.end())return hit->second;
                ++result.polynomial_sign_evaluations;
                const int sign=evaluate(x);
                memo.emplace(x,sign);
                return sign;
            }
        private:
            struct Table {
                unsigned long power;
                ScaledPolynomial values;
            };
            int evaluate(const Q& x) {
                const auto power=decimal_exponent(x.get_den());
                if(!power||*power>max_power)return polynomial_sign(co,x);
                // The smallest cached 10^m with den(x) | 10^m, unless it is much longer.
                const Table* best=nullptr;
                for(const auto& t:tables)
                    if(t.power>=*power&&t.power-*power<=slack&&(!best||t.power<best->power))best=&t;
                if(!best) {
                    // A table holds about n^2/2 factors 10^m in all, while one general evaluation
                    // holds two numbers of n such factors. Very long denominators, as a replayed
                    // file may contain, are therefore evaluated directly.
                    // n <= 257 and m <= 10^5 keep this estimate far inside 64 bits.
                    const std::uint64_t n=co.size(),power_bits=*power*3322/1000+1;
                    if(coefficient_bits+n*(n+1)/2*power_bits>table_bits)return polynomial_sign(co,x);
                    if(tables.size()==capacity)tables.erase(tables.begin());
                    tables.push_back({*power,ScaledPolynomial(co,pow10(static_cast<unsigned>(*power)))});
                    best=&tables.back();
                }
                return best->values.sign(x);
            }
            // Proposal widths of one candidate differ by 8 and 24 places; rounding midpoints of
            // the default widths need 15 places fewer than their endpoints. A table may hold
            // about 2^28 bits (32 MiB); order 256 at 1000 places needs under half of that.
            static constexpr unsigned long max_power=100000,slack=16;
            static constexpr std::size_t capacity=4;
            static constexpr std::uint64_t table_bits=std::uint64_t(1)<<28;
            const std::vector<Z>& co;
            std::uint64_t coefficient_bits=0;
            std::map<Q,int> memo;
            std::vector<Table> tables;
        };
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
        // On failure, row_failure tells whether the shared checks passed and an exact sign of
        // some row failed. With q = det(zE-F) verified, the row would then fail exact inertia too.
        bool polynomial_rootset(const std::vector<RootCertificate>& roots,SignEvaluator& signs,
                                std::size_t n,unsigned digits,SolveResult& result,bool* row_failure=nullptr) {
            // n disjoint brackets/singletons each contain a distinct real zero.
            // Only the verified original degree-n characteristic polynomial
            // makes this a proof of the entire spectrum and its sorted indices.
            if(row_failure)*row_failure=false;
            const auto& co=signs.coefficients();
            if(!result.coefficient_identity_verified||co.size()!=n+1||co.front()==0||!valid_rows(roots,n,digits))return false;
            if(row_failure)*row_failure=true;
            for(const auto& r:roots) {
                int lo=signs(r.lower,result);
                if(r.lower==r.upper) {
                    if(lo!=0)return false;
                }
                else {
                    int hi=signs(r.upper,result);
                    if(lo==0||hi==0||lo==hi)return false;
                }
            }
            if(row_failure)*row_failure=false;
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
            // The recurrence runs on an exactly reduced pencil and maps its coefficients back;
            // the oracle below still checks them against the original pencil.
            auto c=reduced_characteristic(p.E,p.F,o.oracle_workers);
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
        std::optional<RootCertificate> certify_one(const IntegerPencil& pencil,SignEvaluator& signs,const Q& approx,const Q& separation,std::size_t index,unsigned digits,unsigned work,bool polynomial,SolveResult& result) {
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
                    int lo=signs(lower,result),hi=signs(upper,result);
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
                    if(signs(midpoint,result)!=0)continue;
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
        // On failure, per_root_failure tells whether it happened before every candidate had a
        // row: in the mode-independent ordering checks or inside certify_one.
        std::vector<RootCertificate> certify_all(const IntegerPencil& p,SignEvaluator& signs,const std::vector<Real>& candidates,unsigned digits,unsigned work,bool polynomial,SolveResult& result,bool& per_root_failure) {
            per_root_failure=true;
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
                auto r=certify_one(p,signs,approximations[i],separation,i,digits,work,polynomial,result);
                if(!r)return {};
                out.push_back(std::move(*r));
            }
            per_root_failure=false;
            if(!valid_rows(out,p.E.n,digits))return {};
            if(polynomial&&!polynomial_rootset(out,signs,p.E.n,digits,result))return {};
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
    GraeffeCandidates graeffe_candidates(const std::vector<Z>& c,unsigned work,unsigned squarings,unsigned narrow_bits) {
        if(c.size()<2||work<32||work>20000||squarings<1||squarings>512||narrow_bits<1)
            throw std::invalid_argument("invalid Graeffe candidate request");
        Real::default_precision(work);
        Graeffe graeffe(c,graeffe_cut(work),narrow_bits);
        for(unsigned step=0;step<squarings;++step)graeffe.square();
        const auto roots=graeffe.roots(squarings);
        GraeffeCandidates out;
        out.unbounded_exponents=graeffe.wide_exponents();
        for(const auto& r:roots)out.roots.push_back(rational(r));
        return out;
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
        // Signs of the fixed coefficient vector, shared by every certification attempt.
        SignEvaluator signs(out.coefficients);
        for(unsigned work=o.precision;;) {
            ++out.attempts;
            out.working_digits=work;
            Real::default_precision(work);
            if(o.progress)o.progress("Graeffe candidate pass at "+std::to_string(work)+" decimal working digits...");
            auto gt=Clock::now();
            try {
                Graeffe graeffe(out.coefficients,graeffe_cut(work));
                out.graeffe_unbounded_exponents=false;
                const Real threshold("1e-"+std::to_string(o.digits+16));
                std::vector<Real> previous;
                unsigned next_certification=0,certification_delay=2;
                for(unsigned step=1;step<=o.max_squarings;++step) {
                    graeffe.square();
                    out.graeffe_steps=step;
                    out.graeffe_unbounded_exponents=graeffe.wide_exponents();
                    if(step%2&&step!=o.max_squarings)continue;
                    auto candidates=graeffe.roots(step);
                    out.graeffe_unbounded_exponents=graeffe.wide_exponents();
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
                        bool per_root_failure=false;
                        if(try_polynomial) {
                            roots=certify_all(p,signs,candidates,o.digits,work,true,out,per_root_failure);
                            if(roots.size()==k.n)method="polynomial";
                        }
                        // With q = det(zE-F) verified, a per-root polynomial failure proves that
                        // exact inertia fails on the same proposals (docs/ALGORITHM.md, section
                        // 4). Only a failure of the final gates is retried with exact inertia.
                        if(roots.size()!=k.n&&!(try_polynomial&&per_root_failure))
                            roots=certify_all(p,signs,candidates,o.digits,work,false,out,per_root_failure);
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
                // A fresh evaluator: no sign is carried over from an earlier computation.
                SignEvaluator signs(out.coefficients);
                bool row_failure=false;
                if(o.polynomial_certification&&polynomial_rootset(roots,signs,k.n,digits,out,&row_failure)) {
                    accepted=true;
                    out.certification_method="polynomial";
                }
                // A row that fails its exact signs under the verified identity fails exact
                // inertia too, so only failures of the shared row checks reach inertia.
                else if(!row_failure&&inertia_rootset(p,roots,digits,out)) {
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
           <<"  \"graeffe_unbounded_exponents\": "<<(r.graeffe_unbounded_exponents?"true":"false")<<",\n"
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
