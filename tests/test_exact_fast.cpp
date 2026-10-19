#include "secular/exact.hpp"
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>

using namespace secular;

int main() {
    unsigned checks=0;
    auto check=[&](bool value,const char* message) {
        ++checks;
        if(!value)throw std::runtime_error(message);
    };
    auto invalid=[&](auto action,const char* message) {
        bool caught=false;
        try { action(); } catch(const std::invalid_argument&) { caught=true; }
        check(caught,message);
    };
    try {
        // Integer homogeneous Horner must retain the sign of negative content.
        const std::vector<Z> quadratic{4,-4,-3}; // (2x-3)(2x+1)
        check(polynomial_sign(quadratic,Q(-1))==1,"negative x sign");
        check(polynomial_sign(quadratic,Q(-1,2))==0,"negative rational zero");
        check(polynomial_sign(quadratic,Q(0))==-1,"zero x sign");
        check(polynomial_sign(quadratic,Q(1,3))==-1,"rational interior sign");
        check(polynomial_sign(quadratic,Q(3,2))==0,"positive rational zero");
        check(polynomial_sign({-4,4,3},Q(2))==-1,"negative content sign");
        check(polynomial_sign({0,0},Q(5,7))==0,"zero polynomial sign");
        check(polynomial_sign({0,0,-9},Q(-7,3))==-1,"leading zero constant sign");
        const Z huge=pow10(1200);
        check(polynomial_sign({huge,-huge},Q(huge+1,huge))==1,"huge nearby rational above root");
        check(polynomial_sign({-huge,huge},Q(huge-1,huge))==1,"huge negative content below root");
        invalid([&] { polynomial_sign({},Q(0)); },"empty sign polynomial accepted");

        // A shared scaled table must reproduce the general sign at every covered point.
        const std::vector<Z> cubic{8,-36,46,-15}; // (2x-1)(2x-3)(2x-5)
        const ScaledPolynomial scaled(cubic,pow10(6));
        check(scaled.denominator()==pow10(6),"scaled table denominator");
        for(const Q& root:{Q(1,2),Q(3,2),Q(5,2)})check(scaled.covers(root)&&scaled.sign(root)==0,"scaled table root sign");
        for(long numerator=-3000000;numerator<=3000000;numerator+=12347) {
            for(const Q x:{Q(numerator,1000000),Q(numerator,8),Q(numerator,625),Q(numerator,1)}) {
                Q canonical=x;
                canonical.canonicalize();
                check(scaled.covers(canonical)&&scaled.sign(canonical)==polynomial_sign(cubic,canonical),
                      "scaled table sign differs from general sign");
            }
        }
        check(!scaled.covers(Q(1,3))&&!scaled.covers(Q(1,10000000)),"scaled table covers a non-dividing denominator");
        invalid([&] { (void)scaled.sign(Q(1,3)); },"scaled table evaluated an uncovered point");
        invalid([&] { ScaledPolynomial({},Z(10)); },"empty scaled polynomial accepted");
        invalid([&] { ScaledPolynomial(cubic,Z(0)); },"zero scaled denominator accepted");
        invalid([&] { ScaledPolynomial(cubic,Z(-10)); },"negative scaled denominator accepted");
        check(ScaledPolynomial({0,0},Z(7)).sign(Q(5,7))==0,"scaled zero polynomial sign");
        check(ScaledPolynomial({-4,4,3},Z(1)).sign(Q(2))==-1,"scaled negative content sign");
        {
            // Long random coefficients and nearby decimal points, as certification uses them.
            std::mt19937_64 generator(20261018);
            std::vector<Z> wide(41);
            for(auto& c:wide) {
                c=0;
                for(int limb=0;limb<6;++limb)c=(c<<32)+static_cast<unsigned long>(generator()&0xffffffffu);
                if(generator()%2)c=-c;
            }
            const ScaledPolynomial decimal(wide,pow10(40));
            for(int sample=0;sample<200;++sample) {
                const unsigned places=static_cast<unsigned>(generator()%41);
                Z numerator=static_cast<unsigned long>(generator()&0xffffffffu);
                numerator=(numerator<<32)+static_cast<unsigned long>(generator()&0xffffffffu);
                if(generator()%2)numerator=-numerator;
                Q x(numerator,pow10(places));
                x.canonicalize();
                check(decimal.sign(x)==polynomial_sign(wide,x),"scaled table sign differs on long coefficients");
            }
        }

        auto squarefree=check_squarefree({1,0,-2});
        check(!squarefree.repeated&&squarefree.modular_prime!=0,"squarefree witness missing");
        auto repeated=check_squarefree({1,-2,1});
        check(repeated.repeated&&repeated.modular_prime==0,"repeated root not proved exactly");
        check(has_repeated_roots({1,-2,1}),"compatibility repeated wrapper");
        check(!has_repeated_roots({1,0,-2}),"compatibility squarefree wrapper");
        auto degree_drop=check_squarefree({2,1,1});
        check(!degree_drop.repeated&&degree_drop.modular_prime!=0&&degree_drop.modular_prime!=2,
              "degree-dropping prime accepted");
        auto collision=check_squarefree({1,-2,0});
        check(!collision.repeated&&collision.modular_prime!=0&&collision.modular_prime!=2,
              "modular collision treated as repeated");
        auto last_prime=check_squarefree({1,-2310,0});
        check(!last_prime.repeated&&last_prime.modular_prime==65537,"65537 witness unavailable");
        // Primitive leading coefficient 32768 makes the derivative leading
        // residue -1 modulo 65537. Inverting it squares 65536, beyond uint32.
        auto wide_product=check_squarefree({65536,0,2310});
        check(!wide_product.repeated&&wide_product.modular_prime==65537,"modular product overflow");
        // x(x-a) is squarefree modulo p exactly when p does not divide a; the primes below
        // 2^32 follow 65537 in a fixed order, and only all of them failing reaches the gcd.
        Z small_primes=Z(2310)*65537;
        auto mersenne=check_squarefree({1,-small_primes,0});
        check(!mersenne.repeated&&mersenne.modular_prime==2147483647UL,"2^31-1 witness unavailable");
        small_primes*=2147483647UL;
        auto below_2_32=check_squarefree({1,-small_primes,0});
        check(!below_2_32.repeated&&below_2_32.modular_prime==4294967291UL,"2^32-5 witness unavailable");
        small_primes*=4294967291UL;
        auto last_32_bit=check_squarefree({1,-small_primes,0});
        check(!last_32_bit.repeated&&last_32_bit.modular_prime==4294967279UL,"2^32-17 witness unavailable");
        const Z all_primes=small_primes*4294967279UL;
        auto fallback=check_squarefree({1,-all_primes,0});
        check(!fallback.repeated&&fallback.modular_prime==0,"inconclusive reductions must use exact gcd");
        auto fallback_repeated=check_squarefree({1,-2*all_primes,all_primes*all_primes});
        check(fallback_repeated.repeated&&fallback_repeated.modular_prime==0,"exact gcd missed (x-a)^2");
        auto content=check_squarefree({-6,0,12});
        check(!content.repeated&&content.modular_prime!=0,"negative content obscured witness");
        auto content_repeated=check_squarefree({-18,36,-18});
        check(content_repeated.repeated&&content_repeated.modular_prime==0,"negative repeated content");
        check(!check_squarefree({7,-11}).repeated,"linear polynomial marked repeated");
        invalid([&] { check_squarefree({}); },"empty squarefree polynomial accepted");
        invalid([&] { check_squarefree({1}); },"constant squarefree polynomial accepted");
        invalid([&] { check_squarefree({0,1}); },"missing leading coefficient accepted");

        check(oracle_worker_count(1,10)==1,"serial worker count");
        check(oracle_worker_count(64,3)==3,"workers not clamped to jobs");
        check(oracle_worker_count(0,0)==0,"zero jobs must use no workers");
        check(oracle_worker_count(0,1)==1,"auto singleton worker count");
        // Automatic uses every hardware thread, at most the 64 that may be requested.
        const unsigned hardware=std::min(64u,std::max(1u,std::thread::hardware_concurrency()));
        check(oracle_worker_count(0,100)==std::min(hardware,100u),"auto worker count");
        check(oracle_worker_count(0,2)==std::min(hardware,2u),"auto workers not clamped to jobs");

        // Checked exact division: one truncating division, accepted only at remainder zero.
        const Z wide_value=pow10(400)+7;
        const std::vector<Z> divisors{Z(3),Z(-3),Z(Z(1)+pow10(30)),Z(-(Z(1)+pow10(30))),Z(pow10(200)+1)};
        for(const Z& divisor:divisors) {
            check(exact_divide(wide_value*divisor,divisor)==wide_value,"exact quotient");
            check(exact_divide(-wide_value*divisor,divisor)==-wide_value,"negative exact quotient");
            bool caught=false;
            try { exact_divide(wide_value*divisor+1,divisor); }
            catch(const std::logic_error&) { caught=true; }
            check(caught||divisor==1,"non-exact division accepted");
        }
        check(exact_divide(Z(0),Z(5))==0,"zero dividend");
        bool zero_divisor=false;
        try { exact_divide(Z(0),Z(0)); }
        catch(const std::logic_error&) { zero_divisor=true; }
        check(zero_divisor,"zero divisor accepted");
        invalid([&] { oracle_worker_count(65,10); },"excess workers accepted");
        invalid([&] { oracle_worker_count(65,0); },"zero jobs bypassed worker validation");

        auto e=identity<Z>(3);
        ZMatrix f(3);
        f(0,0)=1; f(1,1)=2; f(2,2)=4;
        const std::vector<Z> correct{1,-7,14,-8};
        for(unsigned workers:{0u,1u,2u,8u,64u})
            check(verify_coefficients(e,f,correct,workers),"serial/parallel correct identity rejected");
        check(verify_coefficients(e,f,correct),"default oracle compatibility");
        for(std::size_t missed=0;missed<4;++missed) {
            // Add the cubic product over all other oracle points. The corrupt
            // polynomial agrees at three points and differs only at `missed`.
            std::vector<Z> perturbation{1};
            for(std::size_t j=0;j<4;++j)if(j!=missed) {
                std::vector<Z> next(perturbation.size()+1);
                for(std::size_t k=0;k<perturbation.size();++k) {
                    next[k]+=perturbation[k];
                    next[k+1]+=perturbation[k]*static_cast<unsigned long>(j+1);
                }
                perturbation=std::move(next);
            }
            auto corrupt=correct;
            for(std::size_t k=0;k<corrupt.size();++k)corrupt[k]+=perturbation[k];
            check(!verify_coefficients(e,f,corrupt,1),"serial oracle skipped a point");
            check(!verify_coefficients(e,f,corrupt,4),"parallel oracle skipped a point");
            for(unsigned workers:{1u,4u})
                check(!check_coefficients(e,f,corrupt,workers,CoefficientOracle::bareiss).verified,
                      "Bareiss oracle skipped a point");
        }
        auto nonsymmetric=e;
        nonsymmetric(0,1)=2; nonsymmetric(1,2)=-3;
        check(verify_coefficients(nonsymmetric,f,correct,4),"general triangular pencil rejected");
        check(!verify_coefficients(e,ZMatrix(2),correct,2),"dimension mismatch accepted");
        check(!verify_coefficients(e,f,{1,2},2),"coefficient count mismatch accepted");
        auto malformed=e;
        malformed.data.pop_back();
        check(!verify_coefficients(malformed,f,correct,2),"malformed matrix storage accepted");
        check(!verify_coefficients(ZMatrix{},ZMatrix{},{1},2),"empty matrix accepted");
        invalid([&] { verify_coefficients(e,f,correct,65); },"oracle accepted excess workers");
        std::cout<<checks<<" exact fast-path checks passed\n";
    } catch(const std::exception& ex) {
        std::cerr<<"exact fast-path failure: "<<ex.what()<<"\n";
        return 1;
    }
}
