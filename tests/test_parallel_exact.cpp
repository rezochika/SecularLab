#include "secular/exact.hpp"
#include <iostream>
#include <random>
#include <stdexcept>

using namespace secular;

namespace {
    unsigned checks=0;
    void check(bool ok,const char* message) {
        ++checks;
        if(!ok)throw std::runtime_error(message);
    }
    template<class Action> void invalid(Action action,const char* message) {
        bool caught=false;
        try {
            action();
        }
        catch(const std::invalid_argument&) {
            caught=true;
        }
        check(caught,message);
    }
    // The original serial recurrence on full matrices, measuring every entry of T_k-c_k I.
    CoefficientResult reference_characteristic(const ZMatrix& e,const ZMatrix& f) {
        const auto n=e.n;
        const Z d=determinant(e);
        const ZMatrix j=adjugate(e);
        ZMatrix s=j;
        CoefficientResult r;
        r.coefficients.push_back(d);
        for(std::size_t k=1;k<=n;++k) {
            ZMatrix t=multiply(f,s);
            Z trace=0;
            for(std::size_t i=0;i<n;++i)trace+=t(i,i);
            const Z c=exact_divide(trace,Z(static_cast<unsigned long>(k)));
            ++r.exact_divisions;
            r.coefficients.push_back(-c);
            for(std::size_t i=0;i<n;++i)t(i,i)-=c;
            for(const auto& v:t.data)r.max_workspace_bits=std::max(r.max_workspace_bits,bits(v));
            if(diagonal(e))for(std::size_t i=0;i<n;++i)for(std::size_t c2=0;c2<n;++c2)s(i,c2)=exact_divide(t(i,c2),e(i,i));
            else {
                s=multiply(j,t);
                for(auto& v:s.data)v=exact_divide(v,d);
            }
            r.exact_divisions+=n*n;
        }
        r.terminal_zero=true;
        for(const auto& v:s.data)if(v!=0)r.terminal_zero=false;
        for(const auto& c:r.coefficients)r.max_coefficient_bits=std::max(r.max_coefficient_bits,bits(c));
        return r;
    }
    struct Pencil {
        ZMatrix e,f;
    };
    std::mt19937_64 generator(20261018);
    Z random_integer(unsigned magnitude_bits) {
        Z v=0;
        for(unsigned b=0;b<magnitude_bits;b+=32)v=(v<<32)+static_cast<unsigned long>(generator()&0xffffffffu);
        v>>=static_cast<unsigned long>((magnitude_bits+31)/32*32-magnitude_bits);
        return generator()%2?v:Z(-v);
    }
    // Diagonal or dense SPD E, symmetric or general F, sparse or dense, small or wide entries.
    Pencil random_pencil(std::size_t n,bool diagonal_e,bool symmetric_f,bool symmetric_e,unsigned width,unsigned density) {
        Pencil p {
            ZMatrix(n),ZMatrix(n)
        };
        if(diagonal_e)for(std::size_t i=0;i<n;++i)p.e(i,i)=abs(random_integer(width))+1;
        else {
            ZMatrix c=identity<Z>(n);
            for(std::size_t i=0;i<n;++i)for(std::size_t j=i;j<n;++j)c(i,j)=i==j?Z(abs(random_integer(width))+1):random_integer(width);
            if(symmetric_e) {
                ZMatrix ct(n);
                for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)ct(i,j)=c(j,i);
                p.e=multiply(ct,c);
            }
            else p.e=c;
        }
        for(std::size_t i=0;i<n;++i)for(std::size_t j=symmetric_f?i:0;j<n;++j) {
            if(i!=j&&generator()%100>=density)continue;
            p.f(i,j)=random_integer(width);
            if(symmetric_f)p.f(j,i)=p.f(i,j);
        }
        return p;
    }
    void recurrence_matches_reference(const Pencil& p) {
        const auto n=p.e.n;
        const auto reference=reference_characteristic(p.e,p.f);
        check(reference.terminal_zero,"reference recurrence did not terminate");
        // Small random "general" matrices can be symmetric; the solver decides from the entries.
        const std::size_t per_stage=symmetric(p.e)&&symmetric(p.f)?n*(n+1)/2:n*n;
        for(unsigned workers:{1u,2u,3u,8u,0u}) {
            const auto actual=characteristic(p.e,p.f,workers);
            check(actual.coefficients==reference.coefficients,"recurrence coefficients differ from reference");
            check(actual.terminal_zero,"recurrence did not terminate at S_n=0");
            check(actual.max_workspace_bits==reference.max_workspace_bits,"workspace statistic differs from full T_k-c_k I");
            check(actual.max_coefficient_bits==reference.max_coefficient_bits,"coefficient bits differ");
            check(actual.exact_divisions==n+n*per_stage,"recurrence division count");
            check(workers==0?actual.workers>=1&&actual.workers<=64:actual.workers==std::min<std::size_t>(workers,n),
                  "recurrence workers not clamped to stage rows");
        }
    }
    // L*D*A*D with D=diag(2^x_i 5^y_i): the decimal structure that denominator clearing creates.
    Pencil decimal_pencil(const Pencil& p,unsigned long e_factor,unsigned long f_factor) {
        const auto n=p.e.n;
        std::vector<Z> d(n);
        for(auto& v:d) {
            Z two,five;
            mpz_ui_pow_ui(two.get_mpz_t(),2,static_cast<unsigned long>(generator()%9));
            mpz_ui_pow_ui(five.get_mpz_t(),5,static_cast<unsigned long>(generator()%7));
            v=two*five;
        }
        Pencil out=p;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
            out.e(i,j)*=d[i]*d[j]*e_factor;
            out.f(i,j)*=d[i]*d[j]*f_factor;
        }
        return out;
    }
    bool decimal_power(Z v) {
        mpz_remove(v.get_mpz_t(),v.get_mpz_t(),Z(2).get_mpz_t());
        mpz_remove(v.get_mpz_t(),v.get_mpz_t(),Z(5).get_mpz_t());
        return v==1;
    }
    PencilReduction reduction_matches_recurrence(const Pencil& p) {
        const auto n=p.e.n;
        const auto r=reduce_pencil(p.e,p.f);
        check(r.scale.size()==n&&r.e_content>0&&r.f_content>0,"reduction shape");
        for(std::size_t i=0;i<n;++i) {
            check(r.scale[i]>0&&decimal_power(r.scale[i]),"reduction scale is not a power product of 2 and 5");
            for(std::size_t j=0;j<n;++j) {
                const Z factor=r.scale[i]*r.scale[j];
                check(p.e(i,j)==r.e_content*factor*r.E(i,j)&&p.f(i,j)==r.f_content*factor*r.F(i,j),
                      "reduced pencil is not an exact congruence of the input");
            }
        }
        check(symmetric(r.E)==symmetric(p.e)&&symmetric(r.F)==symmetric(p.f)&&diagonal(r.E)==diagonal(p.e),
              "reduction changed the pencil structure");
        const auto direct=characteristic(p.e,p.f,1);
        for(unsigned workers:{1u,3u,0u}) {
            const auto reduced=reduced_characteristic(p.e,p.f,workers);
            check(reduced.coefficients==direct.coefficients,"reduced recurrence coefficients differ");
            check(reduced.terminal_zero&&reduced.exact_divisions==direct.exact_divisions,"reduced recurrence checks");
            check(reduced.max_coefficient_bits==direct.max_coefficient_bits,"reduced coefficient bits");
            check(reduced.max_workspace_bits<=direct.max_workspace_bits,"reduced workspace is longer");
            check(!r.trivial||reduced.max_workspace_bits==direct.max_workspace_bits,"trivial reduction changed the workspace");
        }
        // The oracle proves the mapped coefficients against the original pencil.
        check(check_coefficients(p.e,p.f,direct.coefficients,2).verified,"oracle rejected reduced coefficients");
        return r;
    }
    void oracle_matches_bareiss(const Pencil& p) {
        const auto n=p.e.n;
        const auto c=characteristic(p.e,p.f).coefficients;
        for(unsigned workers:{1u,4u}) {
            const auto modular=check_coefficients(p.e,p.f,c,workers);
            check(modular.verified&&modular.method=="modular"&&modular.moduli>=1,"modular oracle rejected true coefficients");
            check(modular.workers==std::min<std::size_t>(workers,modular.moduli),"modular workers not clamped to moduli");
            const auto bareiss=check_coefficients(p.e,p.f,c,workers,CoefficientOracle::bareiss);
            check(bareiss.verified&&bareiss.method=="bareiss"&&bareiss.moduli==0,"Bareiss oracle rejected true coefficients");
            check(bareiss.workers==std::min<std::size_t>(workers,n+1),"Bareiss workers not clamped to points");
        }
        for(std::size_t k=0;k<=n;++k)for(int delta:{1,-1}) {
            auto corrupt=c;
            corrupt[k]+=delta;
            const auto modular=check_coefficients(p.e,p.f,corrupt,3);
            check(!modular.verified&&modular.method=="modular","modular oracle accepted a corrupt coefficient");
            check(!check_coefficients(p.e,p.f,corrupt,2,CoefficientOracle::bareiss).verified,
                  "Bareiss oracle accepted a corrupt coefficient");
        }
    }
}

int main() {
    try {
        // The fixed modulus sequence: the largest primes below 2^62, in descending order.
        const std::uint64_t top=std::uint64_t(1)<<62;
        const unsigned offsets[]= {
            57,87,117,143,153,167,171,195,203,273
        };
        for(std::size_t i=0;i<10;++i)check(oracle_prime(i)==top-offsets[i],"oracle prime sequence");

        for(std::size_t n:{1u,2u,3u,7u,12u}) {
            for(unsigned width:{3u,90u}) {
                recurrence_matches_reference(random_pencil(n,true,true,true,width,35));
                recurrence_matches_reference(random_pencil(n,true,true,true,width,100));
                recurrence_matches_reference(random_pencil(n,false,true,true,width,60));
                recurrence_matches_reference(random_pencil(n,true,false,true,width,60));
                recurrence_matches_reference(random_pencil(n,false,false,false,width,60));
                oracle_matches_bareiss(random_pencil(n,true,true,true,width,40));
                oracle_matches_bareiss(random_pencil(n,false,false,false,width,60));
            }
        }
        // 2-adic/5-adic congruence and content reduction before the recurrence.
        for(std::size_t n:{1u,2u,3u,7u,12u}) {
            for(unsigned width:{3u,90u}) {
                reduction_matches_recurrence(random_pencil(n,true,true,true,width,35));
                const auto lumped=reduction_matches_recurrence(decimal_pencil(random_pencil(n,true,true,true,width,35),1000,40));
                check(!lumped.trivial,"decimal pencil was not reduced");
                reduction_matches_recurrence(decimal_pencil(random_pencil(n,false,true,true,width,60),50,8));
                reduction_matches_recurrence(decimal_pencil(random_pencil(n,true,false,true,width,60),25,1));
                reduction_matches_recurrence(decimal_pencil(random_pencil(n,false,false,false,width,60),1,1));
            }
        }
        {
            // Odd entries with a unit gcd leave nothing to remove.
            ZMatrix e=identity<Z>(2),f(2);
            e(1,1)=3;
            f(0,0)=7; f(0,1)=f(1,0)=1; f(1,1)=5;
            const auto r=reduce_pencil(e,f);
            check(r.trivial&&r.E.data==e.data&&r.F.data==f.data,"irreducible pencil was changed");
            // A zero F has content zero, which the reduction treats as one.
            const auto zero=reduction_matches_recurrence({identity<Z>(3),ZMatrix(3)});
            check(zero.f_content==1,"zero stiffness content");
        }
        invalid([] {
            reduce_pencil(identity<Z>(2),identity<Z>(3));
        },"reduction accepted mismatched dimensions");
        invalid([] {
            reduced_characteristic(ZMatrix(2),identity<Z>(2),2);
        },"reduced recurrence accepted a singular leading matrix");
        invalid([] {
            reduced_characteristic(identity<Z>(2),identity<Z>(2),65);
        },"reduced recurrence accepted excess workers");
        invalid([] {
            characteristic(identity<Z>(2),identity<Z>(2),65);
        },"recurrence accepted excess workers");
        invalid([] {
            characteristic(ZMatrix(2),identity<Z>(2),2);
        },"singular leading matrix accepted");
        invalid([] {
            characteristic(ZMatrix{},ZMatrix{},2);
        },"empty pencil accepted");

        // Many moduli: wide entries need several CRT primes, scheduled across workers.
        const auto wide=random_pencil(6,true,true,true,240,100);
        const auto wide_c=characteristic(wide.e,wide.f).coefficients;
        const auto wide_report=check_coefficients(wide.e,wide.f,wide_c,4);
        check(wide_report.verified&&wide_report.moduli>8&&wide_report.workers==4,"multi-prime parallel proof");
        // No prime above 2^61 divides this det(E) < 2^61, so the proof used the first moduli.
        // Adding their product leaves every residue unchanged; the self-referential bound must
        // then demand further primes, which detect the change.
        const auto small=random_pencil(5,true,true,true,8,100);
        auto small_c=characteristic(small.e,small.f).coefficients;
        const auto small_report=check_coefficients(small.e,small.f,small_c,2);
        check(small_report.verified&&small_report.moduli>=1,"small pencil proof");
        Z product=1;
        for(std::size_t i=0;i<small_report.moduli;++i) {
            const auto p=oracle_prime(i);
            Z word;
            mpz_import(word.get_mpz_t(),1,-1,sizeof p,0,0,&p);
            product*=word;
        }
        for(std::size_t k=0;k<small_c.size();++k) {
            auto aliased=small_c;
            aliased[k]+=product;
            const auto report=check_coefficients(small.e,small.f,aliased,2);
            check(!report.verified&&report.method=="modular"&&report.moduli>small_report.moduli,
                  "modular oracle accepted a residue-preserving corruption");
        }
        // Moduli dividing det(E) are skipped, never trusted.
        Pencil divisible=random_pencil(4,true,true,true,20,100);
        const auto p0=oracle_prime(0),p1=oracle_prime(1);
        Z factor0,factor1;
        mpz_import(factor0.get_mpz_t(),1,-1,sizeof p0,0,0,&p0);
        mpz_import(factor1.get_mpz_t(),1,-1,sizeof p1,0,0,&p1);
        divisible.e(0,0)=factor0*factor1;
        const auto divisible_c=characteristic(divisible.e,divisible.f).coefficients;
        const auto divisible_report=check_coefficients(divisible.e,divisible.f,divisible_c,3);
        check(divisible_report.verified&&divisible_report.method=="modular","modular proof with skipped moduli");
        auto divisible_bad=divisible_c;
        divisible_bad.back()+=1;
        check(!check_coefficients(divisible.e,divisible.f,divisible_bad,3).verified,"skipped moduli hid a corruption");
        // Singular E: every modulus divides det E, which the bound proves is zero. The
        // inconclusive modular check then defers to the n+1 Bareiss determinants.
        ZMatrix singular_e(2),singular_f(2);
        singular_e(0,0)=singular_e(0,1)=singular_e(1,0)=singular_e(1,1)=1;
        singular_f(0,0)=2; singular_f(0,1)=3; singular_f(1,0)=-1; singular_f(1,1)=5;
        // det(zE-F)=(z-2)(z-5)-(z-3)(z+1) = -5z+13.
        const std::vector<Z> singular_c {
            0,-5,13
        };
        const auto fallback=check_coefficients(singular_e,singular_f,singular_c,2);
        check(fallback.verified&&fallback.method=="bareiss"&&fallback.moduli==0,"singular E fallback");
        check(!check_coefficients(singular_e,singular_f,{0,-5,12},2).verified,"singular E fallback accepted corruption");
        ZMatrix diagonal_singular=identity<Z>(3);
        diagonal_singular(1,1)=0;
        const auto diagonal_singular_f=identity<Z>(3);
        // det(z*diag(1,0,1)-I) = -(z-1)^2.
        check(check_coefficients(diagonal_singular,diagonal_singular_f,{0,-1,2,-1},4).method=="bareiss",
              "diagonal singular E fallback");
        // A tiny order with very long entries: residues of every entry and coefficient per
        // modulus would cost far more than three Bareiss determinants, which decide instead.
        const auto long_entries=random_pencil(2,true,true,true,200000,100);
        const auto long_c=characteristic(long_entries.e,long_entries.f).coefficients;
        const auto long_report=check_coefficients(long_entries.e,long_entries.f,long_c,2);
        check(long_report.verified&&long_report.method=="bareiss"&&long_report.moduli==0,"size estimate prefers Bareiss");
        auto long_bad=long_c;
        long_bad[1]+=1;
        check(!check_coefficients(long_entries.e,long_entries.f,long_bad,2).verified,"long-entry corruption accepted");
        check(!verify_coefficients(identity<Z>(2),identity<Z>(2),{1,-2},2),"coefficient count mismatch accepted");
        invalid([] {
            check_coefficients(identity<Z>(2),identity<Z>(2),{1,-2,1},65);
        },"oracle accepted excess workers");
        std::cout<<"PASS "<<checks<<" parallel recurrence and modular oracle checks\n";
        return 0;
    }
    catch(const std::exception& ex) {
        std::cerr<<"FAIL after "<<checks<<" checks: "<<ex.what()<<'\n';
        return 1;
    }
}
