#include "secular/exact.hpp"
#include "secular/solver.hpp"
#include <iostream>
#include <random>
#include <stdexcept>
#include <tuple>

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
    struct Pencil {
        ZMatrix e,f;
    };
    std::mt19937_64 generator(20261020);
    Z random_integer(unsigned magnitude_bits) {
        Z v=0;
        for(unsigned b=0;b<magnitude_bits;b+=32)v=(v<<32)+static_cast<unsigned long>(generator()&0xffffffffu);
        v>>=static_cast<unsigned long>((magnitude_bits+31)/32*32-magnitude_bits);
        return generator()%2?v:Z(-v);
    }
    ZMatrix transpose(const ZMatrix& a) {
        ZMatrix t(a.n);
        for(std::size_t i=0;i<a.n;++i)for(std::size_t j=0;j<a.n;++j)t(i,j)=a(j,i);
        return t;
    }
    // Diagonal, dense symmetric or general E with nonzero determinant; symmetric or general F.
    Pencil random_pencil(std::size_t n,bool diagonal_e,bool symmetric_e,bool symmetric_f,unsigned width,unsigned density) {
        Pencil p {
            ZMatrix(n),ZMatrix(n)
        };
        if(diagonal_e)for(std::size_t i=0;i<n;++i)p.e(i,i)=abs(random_integer(width))+1;
        else {
            ZMatrix c(n);
            for(std::size_t i=0;i<n;++i)for(std::size_t j=i;j<n;++j)c(i,j)=i==j?Z(abs(random_integer(width))+1):random_integer(width);
            p.e=symmetric_e?multiply(transpose(c),c):c;
        }
        for(std::size_t i=0;i<n;++i)for(std::size_t j=symmetric_f?i:0;j<n;++j) {
            if(i!=j&&generator()%100>=density)continue;
            p.f(i,j)=random_integer(width);
            if(symmetric_f)p.f(j,i)=p.f(i,j);
        }
        return p;
    }
    Z word(std::uint64_t w) {
        Z z;
        mpz_import(z.get_mpz_t(),1,-1,sizeof w,0,0,&w);
        return z;
    }
    // The recurrence and Chinese remaindering give the same integers for every worker count.
    CoefficientResult crt_matches_recurrence(const Pencil& p) {
        const auto reference=characteristic(p.e,p.f);
        CoefficientResult first;
        for(unsigned workers:{1u,2u,3u,0u}) {
            const auto crt=crt_characteristic(p.e,p.f,workers);
            check(crt.coefficients==reference.coefficients,"Chinese remaindering differs from the recurrence");
            check(crt.max_coefficient_bits==reference.max_coefficient_bits,"Chinese remaindering coefficient bits");
            check(crt.moduli>=1&&crt.exact_divisions==0&&crt.max_workspace_bits==0&&!crt.terminal_zero,
                  "Chinese remaindering statistics");
            check(workers==0?crt.workers>=1&&crt.workers<=64:crt.workers>=1&&crt.workers<=workers,
                  "Chinese remaindering workers not clamped");
            if(workers==1)first=crt;
            else check(crt.moduli==first.moduli&&crt.skipped_moduli==first.skipped_moduli,
                       "moduli depend on the worker count");
        }
        // The oracle accepts the result. Unreduced, the combined residues are a proof it accepts with
        // the same primes; reduced, the original pencil's bound needs at least as many, less one.
        const auto report=check_coefficients(p.e,p.f,first.coefficients,2);
        check(report.verified,"oracle rejected Chinese remaindering");
        check(report.method=="bareiss"||(first.reduced?first.moduli<=report.moduli+1:
              report.moduli<=first.moduli+first.skipped_moduli),"oracle and Chinese remaindering moduli");
        check(first.bound_bits+1<=61*first.moduli&&61*(first.moduli-1)<first.bound_bits+1,"bound sets the moduli");
        check(first.reduced==!reduce_pencil(p.e,p.f).trivial,"reduction flag");
        return first;
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
    // Chinese remaindering of a reducible pencil: the moduli follow the reduced bound, at most the
    // original pencil's oracle needs, and never more than the reduced pencil's own bound asks for.
    CoefficientResult reduced_route(const Pencil& p) {
        const auto crt=crt_matches_recurrence(p);
        const auto reduction=reduce_pencil(p.e,p.f);
        check(crt.reduced&&!reduction.trivial,"decimal pencil was not reduced");
        const auto oracle=check_coefficients(p.e,p.f,crt.coefficients,1);
        check(oracle.verified&&oracle.method=="modular"&&crt.moduli<=oracle.moduli,"reduction added moduli");
        // b_k <= B'_k: the reduced pencil passed directly, when irreducible, uses B'_k alone.
        if(reduce_pencil(reduction.E,reduction.F).trivial) {
            const auto direct=crt_characteristic(reduction.E,reduction.F,2);
            check(!direct.reduced&&crt.bound_bits<=direct.bound_bits&&crt.moduli<=direct.moduli,
                  "the per-coefficient bound exceeds the reduced pencil's bound");
        }
        return crt;
    }
    // The residue kernels against the reference kernel (Gauss-Jordan on [E | F], then the right-looking
    // Hessenberg reduction) and against the exact coefficients modulo p. Small primes make zero pivots,
    // pivoting swaps, zero subdiagonals, Lanczos breakdowns and det(E) = 0 modulo p common; nullopt must
    // mean p | det(E). The symmetric and factored kernels take the symmetric routes on every symmetric
    // pencil (the congruence and the factored route for E not diagonal), so each modulus ends in one of
    // them or in one of their fallbacks.
    unsigned singular_moduli=0;
    ResidueRoutes fast_routes,symmetric_routes,factored_routes;
    void add(ResidueRoutes& total,const ResidueRoutes& r) {
        total.lanczos+=r.lanczos;
        total.congruence+=r.congruence;
        total.factored+=r.factored;
        total.pivot_fallbacks+=r.pivot_fallbacks;
        total.breakdown_fallbacks+=r.breakdown_fallbacks;
    }
    std::size_t routed(const ResidueRoutes& r) {
        return r.lanczos+r.congruence+r.factored+r.pivot_fallbacks+r.breakdown_fallbacks;
    }
    ResidueRoutes kernels_agree(const ZMatrix& e,const ZMatrix& f,const std::vector<std::uint64_t>& primes) {
        ResidueRoutes fast_route,symmetric_route,factored_route;
        const auto fast=pencil_residues(e,f,primes,ResidueKernel::fast,&fast_route);
        const auto reference=pencil_residues(e,f,primes,ResidueKernel::reference);
        check(fast==reference,"residue kernel differs from the reference kernel");
        check(pencil_residues(e,f,primes,ResidueKernel::general)==reference,"general kernel differs from the reference kernel");
        check(pencil_residues(e,f,primes,ResidueKernel::symmetric,&symmetric_route)==reference,
              "symmetric kernel differs from the reference kernel");
        check(pencil_residues(e,f,primes,ResidueKernel::factored,&factored_route)==reference,
              "factored kernel differs from the reference kernel");
        std::size_t usable=0;
        for(const auto& r:reference)usable+=r.has_value();
        if(!symmetric(e)||!symmetric(f))check(routed(fast_route)==0&&routed(symmetric_route)==0&&routed(factored_route)==0,
                                               "symmetric route on an unsymmetric pencil");
        else if(diagonal(e))check(symmetric_route.lanczos+symmetric_route.breakdown_fallbacks==usable&&
                                  symmetric_route.congruence==0&&symmetric_route.factored==0&&symmetric_route.pivot_fallbacks==0&&
                                  factored_route.lanczos==symmetric_route.lanczos&&routed(factored_route)==routed(symmetric_route),
                                  "Lanczos routes");
        else {
            check(routed(symmetric_route)==primes.size()&&symmetric_route.lanczos==0&&symmetric_route.factored==0,"congruence routes");
            check(routed(factored_route)==primes.size()&&factored_route.lanczos==0&&factored_route.congruence==0,"factored routes");
            check(routed(fast_route)==primes.size()&&(fast_route.congruence==0||fast_route.factored==0),"one route per pencil");
        }
        add(fast_routes,fast_route);
        add(symmetric_routes,symmetric_route);
        add(factored_routes,factored_route);
        const Z det=determinant(e);
        const auto exact=det!=0?characteristic(e,f).coefficients:std::vector<Z>{};
        for(std::size_t i=0;i<primes.size();++i) {
            const Z p=word(primes[i]);
            Z r;
            mpz_mod(r.get_mpz_t(),det.get_mpz_t(),p.get_mpz_t());
            check(!fast[i]==(r==0),"residue kernel misjudged p | det(E)");
            if(!fast[i]) {
                ++singular_moduli;
                continue;
            }
            check(fast[i]->size()==e.n+1,"residue count");
            for(std::size_t k=0;k<exact.size();++k) {
                mpz_mod(r.get_mpz_t(),exact[k].get_mpz_t(),p.get_mpz_t());
                check(word((*fast[i])[k])==r,"residue differs from the exact coefficient");
            }
        }
        return fast_route;
    }
    // Dense general E, singular for some small primes and possibly over the integers.
    ZMatrix general_matrix(std::size_t n,unsigned width) {
        ZMatrix a(n);
        for(auto& v:a.data)v=random_integer(width);
        return a;
    }
    // Symmetric F supported on a band, on diagonal blocks, or on a tridiagonal (already Hessenberg).
    ZMatrix structured_matrix(std::size_t n,std::size_t band,std::size_t block,unsigned width,unsigned density) {
        ZMatrix a(n);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=i;j<n;++j) {
            if(j-i>band||i/block!=j/block||(i!=j&&generator()%100>=density))continue;
            a(i,j)=a(j,i)=random_integer(width);
        }
        return a;
    }
    QMatrix rational(const ZMatrix& a,const Q& scale) {
        QMatrix q(a.n);
        for(std::size_t i=0;i<a.data.size();++i)q.data[i]=Q(a.data[i])*scale;
        return q;
    }
}

int main() {
    try {
        // Hand-checked small pencils: det(3z-5); det(z*diag(2,3)-[[1,4],[4,-2]]) = 6z^2+z-18.
        {
            ZMatrix e(1),f(1);
            e(0,0)=3;
            f(0,0)=5;
            check(crt_characteristic(e,f).coefficients==std::vector<Z>({3,-5}),"order one");
            ZMatrix e2(2),f2(2);
            e2(0,0)=2; e2(1,1)=3;
            f2(0,0)=1; f2(0,1)=f2(1,0)=4; f2(1,1)=-2;
            check(crt_characteristic(e2,f2,2).coefficients==std::vector<Z>({6,1,-18}),"order two");
            // A general E with det(E)=-1: det(z[[0,1],[1,0]]-I) = -(z^2-1) = -z^2+1.
            ZMatrix swap(2);
            swap(0,1)=swap(1,0)=1;
            check(crt_characteristic(swap,identity<Z>(2)).coefficients==std::vector<Z>({-1,0,1}),"pivoting E");
        }
        for(std::size_t n:{1u,2u,3u,7u,12u}) {
            for(unsigned width:{3u,90u}) {
                crt_matches_recurrence(random_pencil(n,true,true,true,width,35));
                crt_matches_recurrence(random_pencil(n,true,true,true,width,100));
                crt_matches_recurrence(random_pencil(n,false,true,true,width,60));
                crt_matches_recurrence(random_pencil(n,true,true,false,width,60));
                crt_matches_recurrence(random_pencil(n,false,false,false,width,60));
            }
        }
        // Wide entries need many moduli, spread across the workers.
        {
            const auto wide=random_pencil(6,false,true,true,240,100);
            const auto crt=crt_matches_recurrence(wide);
            check(crt.moduli>8&&crt.skipped_moduli==0,"wide pencil moduli");
            check(crt_characteristic(wide.e,wide.f,4).workers==4,"wide pencil workers");
            // A zero F leaves only det(E) z^n; a very long entry in a tiny pencil still combines.
            const auto zero=crt_matches_recurrence({identity<Z>(3),ZMatrix(3)});
            check(zero.coefficients==std::vector<Z>({1,0,0,0}),"zero stiffness");
            crt_matches_recurrence(random_pencil(2,true,true,true,20000,100));
        }
        // Moduli dividing det(E) are skipped and the result is still exact: diagonal E with
        // det(E) = p0*p1*(short factors), and a general E with det(E) = p0 by a unimodular congruence.
        {
            const auto p0=word(oracle_prime(0)),p1=word(oracle_prime(1));
            auto divisible=random_pencil(4,true,true,true,20,100);
            divisible.e(0,0)=p0*p1;
            const auto diagonal_skip=crt_matches_recurrence(divisible);
            check(diagonal_skip.skipped_moduli==2,"diagonal E skipped moduli");
            ZMatrix u=identity<Z>(3),d=identity<Z>(3);
            u(0,1)=2; u(0,2)=-1; u(1,2)=3;
            d(0,0)=p0;
            Pencil general {
                multiply(transpose(u),multiply(d,u)),ZMatrix(3)
            };
            for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)general.f(i,j)=static_cast<long>(i*3+j)-4;
            check(determinant(general.e)==p0,"constructed det(E)");
            const auto general_skip=crt_matches_recurrence(general);
            check(general_skip.skipped_moduli==1,"general E skipped modulus");
        }
        // Decimal-shaped pencils run on the reduced pencil E' = P^-1 E P^-1/alpha, F' = P^-1 F P^-1/beta
        // and map back by q_k = f_k c'_k, for every structure and worker count.
        for(std::size_t n:{1u,2u,3u,7u,12u}) {
            for(unsigned width:{3u,90u}) {
                reduced_route(decimal_pencil(random_pencil(n,true,true,true,width,35),1000,40));
                reduced_route(decimal_pencil(random_pencil(n,false,true,true,width,60),50,8));
                reduced_route(decimal_pencil(random_pencil(n,true,false,true,width,60),25,1));
                reduced_route(decimal_pencil(random_pencil(n,false,false,false,width,60),1,10));
            }
        }
        {
            // Long 2-adic and 5-adic factors: far fewer moduli than the original pencil's oracle.
            const auto big=decimal_pencil(random_pencil(12,true,true,true,8,100),1000000000,10000000);
            const auto crt=reduced_route(big);
            check(3*crt.moduli<check_coefficients(big.e,big.f,crt.coefficients,1).moduli,"reduction saved few moduli");
            // Moduli dividing det(E') are skipped: det(E) = p0 * (powers of 2 and 5) * (short factors).
            auto divisible=decimal_pencil(random_pencil(4,true,true,true,20,100),100,1);
            divisible.e(0,0)*=word(oracle_prime(0));
            check(reduced_route(divisible).skipped_moduli==1,"reduced pencil skipped moduli");
            // A singular E stays singular after the reduction and is rejected.
            ZMatrix singular(2),f(2);
            singular(0,0)=singular(0,1)=singular(1,0)=singular(1,1)=400;
            f(0,0)=200; f(0,1)=f(1,0)=-100; f(1,1)=500;
            check(!reduce_pencil(singular,f).trivial,"singular decimal pencil was not reduced");
            invalid([&] {
                crt_characteristic(singular,f,2);
            },"reduced singular E accepted");
        }
        // A singular E over the integers: every modulus is skipped until the bound proves
        // det(E)=0, and the request fails as the recurrence does.
        {
            ZMatrix singular(2),f(2);
            singular(0,0)=singular(0,1)=singular(1,0)=singular(1,1)=1;
            f(0,0)=2; f(0,1)=3; f(1,0)=-1; f(1,1)=5;
            invalid([&] {
                crt_characteristic(singular,f,2);
            },"dense singular E accepted");
            ZMatrix diagonal_singular=identity<Z>(3);
            diagonal_singular(1,1)=0;
            invalid([&] {
                crt_characteristic(diagonal_singular,identity<Z>(3));
            },"diagonal singular E accepted");
            invalid([] {
                crt_characteristic(ZMatrix(2),identity<Z>(2));
            },"zero E accepted");
        }
        invalid([] {
            crt_characteristic(identity<Z>(2),identity<Z>(3));
        },"mismatched dimensions accepted");
        invalid([] {
            crt_characteristic(ZMatrix{},ZMatrix{});
        },"empty pencil accepted");
        invalid([] {
            crt_characteristic(identity<Z>(2),identity<Z>(2),65);
        },"excess workers accepted");

        // The residue kernel gives the reference kernel's residues on every structure, for oracle,
        // mid-size and small primes alike.
        {
            const std::vector<std::uint64_t> primes {3,5,7,11,13,101,65537,2147483647ULL,4294967291ULL,
                                                     2305843009213693951ULL,oracle_prime(0),oracle_prime(1),oracle_prime(7)};
            for(std::size_t n:{1u,2u,3u,4u,6u,9u,16u,24u}) {
                for(unsigned width:{2u,40u}) {
                    for(unsigned density:{100u,30u,8u}) {
                        for(const auto& [diagonal_e,symmetric_e,symmetric_f]:std::vector<std::tuple<bool,bool,bool>>{
                            {true,true,true},{false,true,true},{true,true,false},{false,false,false}}) {
                            const auto p=random_pencil(n,diagonal_e,symmetric_e,symmetric_f,width,density);
                            kernels_agree(p.e,p.f,primes);
                        }
                        kernels_agree(general_matrix(n,width),random_pencil(n,true,true,false,width,density).f,primes);
                    }
                    // Tridiagonal (already Hessenberg), 3 x 3 blocks (zero subdiagonals), a band that fills in,
                    // and a banded F behind a dense E.
                    const auto e=random_pencil(n,true,true,true,width,100).e;
                    kernels_agree(e,structured_matrix(n,1,n,width,100),primes);
                    kernels_agree(e,structured_matrix(n,n,3,width,100),primes);
                    kernels_agree(e,structured_matrix(n,3,n,width,80),primes);
                    kernels_agree(random_pencil(n,false,true,true,width,100).e,structured_matrix(n,2,n,width,100),primes);
                }
            }
            // Larger sparse pencils: columns with several nonzeros below the subdiagonal, fill-in that
            // ends the sparse reduction part way, and blocks; then dense F needing a pivoting swap.
            for(unsigned density:{6u,15u,40u}) {
                const auto p=random_pencil(48,true,true,true,2,density);
                kernels_agree(p.e,p.f,primes);
            }
            const auto e48=random_pencil(48,true,true,true,2,100).e;
            kernels_agree(e48,structured_matrix(48,4,48,2,100),primes);
            kernels_agree(e48,structured_matrix(48,48,6,2,100),primes);
            auto swap_needed=random_pencil(40,true,true,true,30,100);
            swap_needed.f(0,1)=swap_needed.f(1,0)=0;
            swap_needed.f(0,2)=swap_needed.f(2,0)=0;
            kernels_agree(swap_needed.e,swap_needed.f,primes);
            // Oracle primes dividing det(E): a diagonal E with E_00 = p0 p1, and a dense E = U^T D U with
            // det(E) = p0.
            const auto p0=oracle_prime(0),p1=oracle_prime(1);
            auto divisible=random_pencil(5,true,true,true,20,100);
            divisible.e(0,0)=word(p0)*word(p1);
            kernels_agree(divisible.e,divisible.f,primes);
            check(!pencil_residues(divisible.e,divisible.f,{p0})[0]&&!pencil_residues(divisible.e,divisible.f,{p1})[0]&&
                  pencil_residues(divisible.e,divisible.f,{oracle_prime(2)})[0],"diagonal E divisible by oracle primes");
            ZMatrix u=identity<Z>(6),d=identity<Z>(6);
            for(std::size_t i=0;i<6;++i)for(std::size_t j=i+1;j<6;++j)u(i,j)=random_integer(12);
            d(0,0)=word(p0);
            const auto dense_divisible=multiply(transpose(u),multiply(d,u));
            check(determinant(dense_divisible)==word(p0),"constructed dense det(E)");
            kernels_agree(dense_divisible,random_pencil(6,true,true,true,30,100).f,primes);
            check(!pencil_residues(dense_divisible,identity<Z>(6),{p0},ResidueKernel::reference)[0]&&
                  !pencil_residues(dense_divisible,identity<Z>(6),{p0})[0],"dense E divisible by an oracle prime");
            // E singular over the integers leaves no modulus usable.
            auto singular=general_matrix(7,20);
            for(std::size_t j=0;j<7;++j)singular(3,j)=singular(1,j)*2-singular(5,j);
            kernels_agree(singular,general_matrix(7,20),primes);
            check(singular_moduli>=primes.size()+3+20,"small primes rarely divided det(E)");
            // Fallbacks at oracle primes. A zero leading pivot of a nonsingular E: LDL^T fails, LU pivots.
            {
                auto pivot=random_pencil(6,false,true,true,20,100);
                pivot.e(0,0)=word(p0);
                pivot.e(1,1)=word(p0)*3+1;
                check(determinant(pivot.e)%word(p0)!=0,"constructed pivot-free failure");
                const auto routes=kernels_agree(pivot.e,pivot.f,{p0,p1});
                check(routes.pivot_fallbacks==1&&routes.congruence==1,"zero LDL^T pivot did not fall back");
                ResidueRoutes factored;
                pencil_residues(pivot.e,pivot.f,{p0,p1},ResidueKernel::factored,&factored);
                check(factored.pivot_fallbacks==1&&factored.factored==1,"zero LDL^T pivot: factored route");
            }
            // Every eigenvalue twice: the Krylov space of any vector is a proper invariant subspace, so
            // Lanczos breaks down at every modulus, after E diagonal and after the congruence alike.
            {
                const auto g=random_pencil(4,true,true,true,20,100);
                ZMatrix e(8),f(8);
                for(std::size_t i=0;i<8;++i)for(std::size_t j=0;j<8;++j)if(i/4==j/4) {
                    e(i,j)=g.e(i%4,j%4);
                    f(i,j)=g.f(i%4,j%4);
                }
                const std::vector<std::uint64_t> large {p0,p1,oracle_prime(2)};
                ResidueRoutes routes;
                kernels_agree(e,f,large);
                pencil_residues(e,f,large,ResidueKernel::symmetric,&routes);
                check(routes.breakdown_fallbacks==3&&routes.lanczos==0,"repeated eigenvalues: diagonal E");
                ZMatrix u=identity<Z>(8);
                for(std::size_t i=0;i<8;++i)for(std::size_t j=i+1;j<8;++j)u(i,j)=random_integer(6);
                const auto dense_e=multiply(transpose(u),multiply(e,u)),dense_f=multiply(transpose(u),multiply(f,u));
                kernels_agree(dense_e,dense_f,large);
                pencil_residues(dense_e,dense_f,large,ResidueKernel::symmetric,&routes);
                check(routes.breakdown_fallbacks==3&&routes.congruence==0,"repeated eigenvalues: dense E");
                pencil_residues(dense_e,dense_f,large,ResidueKernel::factored,&routes);
                check(routes.breakdown_fallbacks==3&&routes.factored==0,"repeated eigenvalues: dense E, factored route");
                // Two copies of an 8 x 8 pencil with tridiagonal E and F: E banded, every eigenvalue twice.
                ZMatrix band_e(16),band_f(16);
                const auto g8e=structured_matrix(8,1,8,20,100),g8f=structured_matrix(8,1,8,20,100);
                for(std::size_t i=0;i<16;++i)for(std::size_t j=0;j<16;++j)if(i/8==j/8) {
                    band_e(i,j)=g8e(i%8,j%8);
                    band_f(i,j)=g8f(i%8,j%8);
                }
                kernels_agree(band_e,band_f,large);
                pencil_residues(band_e,band_f,large,ResidueKernel::factored,&routes);
                check(routes.breakdown_fallbacks==3&&routes.factored==0,"repeated eigenvalues: banded E, factored route");
            }
            // The fast kernel decides the diagonal-E route once per pencil: dense rows take Lanczos, a
            // tridiagonal F (already Hessenberg) and 3 x 3 blocks (zero subdiagonals) the general route,
            // a band that fills the Hessenberg form Lanczos; the scattered sparse F runs on its nonzeros.
            {
                const std::vector<std::uint64_t> large {p0,p1};
                const auto d=random_pencil(64,true,true,true,20,100).e;
                check(kernels_agree(d,random_pencil(64,true,true,true,20,100).f,large).lanczos==2,"dense F: Lanczos");
                check(kernels_agree(d,structured_matrix(64,1,64,20,100),large).lanczos==0,"tridiagonal F: general");
                check(kernels_agree(d,structured_matrix(64,64,3,20,100),large).lanczos==0,"block F: general");
                check(kernels_agree(d,structured_matrix(64,3,64,20,100),large).lanczos==2,"banded F: Lanczos");
                check(kernels_agree(d,random_pencil(64,true,true,true,20,4).f,large).lanczos==2,"sparse F: Lanczos");
                // A block-diagonal E keeps C = L^-1 F L^-T banded, so Lanczos uses its rows' short ranges.
                const auto blocks=structured_matrix(64,64,4,20,100);
                kernels_agree(multiply(blocks,blocks),structured_matrix(64,6,64,20,100),large);
            }
            // For E not diagonal the fast kernel forms C at one modulus when the bounds leave the choice
            // open. A generic band fills C, so Lanczos through the factor of E is cheaper: pentadiagonal E
            // and F, E with a sparse F on wide ranges (products on F's nonzeros), and an E whose profile
            // has gaps (L^T rows with zeros in their range). F = L Lambda L^T with E = L D L^T keeps C
            // diagonal, and a dense E and F never try the factored route.
            {
                const std::vector<std::uint64_t> large {p0,p1};
                const std::size_t n=32;
                const auto band_e=structured_matrix(n,2,n,12,100);
                check(kernels_agree(band_e,structured_matrix(n,2,n,12,100),large).factored==2,"banded E and F: factored route");
                ZMatrix wide=identity<Z>(n);
                for(std::size_t i=0;i<n;++i)wide(i,i)=random_integer(12);
                for(std::size_t i=8;i<n;i+=8)wide(i,0)=wide(0,i)=wide(i,i-8)=wide(i-8,i)=random_integer(12);
                check(kernels_agree(band_e,wide,large).factored==2,"banded E, sparse F: factored route");
                auto gaps=structured_matrix(n,1,n,12,100);
                for(std::size_t i=12;i<n;i+=9)gaps(i,i-9)=gaps(i-9,i)=random_integer(12);
                check(kernels_agree(gaps,structured_matrix(n,3,n,12,100),large).factored==2,"E with profile gaps: factored route");
                ZMatrix l=identity<Z>(n),d(n),lambda(n);
                for(std::size_t i=0;i<n;++i) {
                    d(i,i)=abs(random_integer(12))+1;
                    lambda(i,i)=random_integer(12);
                    for(std::size_t j=i>2?i-2:0;j<i;++j)l(i,j)=random_integer(6);
                }
                const auto congruent=kernels_agree(multiply(l,multiply(d,transpose(l))),multiply(l,multiply(lambda,transpose(l))),large);
                check(congruent.congruence==2&&congruent.factored==0,"diagonal C: congruence");
                const auto dense=random_pencil(n,false,true,true,12,100);
                check(kernels_agree(dense.e,dense.f,large).congruence==2,"dense E and F: congruence");
            }
            // Mostly small primes: 1226 Lanczos, 522 congruence, 109 pivot and 449 breakdown fallbacks for
            // the symmetric kernel when written, 518 factored, 109 pivot and 453 breakdown fallbacks for the
            // factored kernel. The fast kernel keeps tridiagonal F (all of n <= 2) and some banded F on the
            // general route and takes the factored route where C fills (dense E behind a sparse F): 329
            // Lanczos, 305 congruence, 220 factored and 196 breakdowns.
            check(symmetric_routes.lanczos>1000&&symmetric_routes.congruence>400&&symmetric_routes.pivot_fallbacks>80&&
                  symmetric_routes.breakdown_fallbacks>300,"symmetric routes or fallbacks rarely taken");
            check(factored_routes.factored>400&&factored_routes.congruence==0&&factored_routes.pivot_fallbacks>80&&
                  factored_routes.breakdown_fallbacks>300,"factored route or its fallbacks rarely taken");
            check(fast_routes.lanczos>250&&fast_routes.congruence>250&&fast_routes.factored>150&&fast_routes.pivot_fallbacks>80&&
                  fast_routes.breakdown_fallbacks>150,"fast kernel rarely took the symmetric routes");
            invalid([] {
                pencil_residues(identity<Z>(2),identity<Z>(2),{9});
            },"composite modulus accepted");
            invalid([] {
                pencil_residues(identity<Z>(2),identity<Z>(2),{2});
            },"even modulus accepted");
            invalid([] {
                pencil_residues(identity<Z>(2),identity<Z>(2),{(std::uint64_t(1)<<62)+135});
            },"modulus above 2^62 accepted");
        }

        // Content: the positive gcd, one for a zero vector; a pencil 2E,6F has content >= 2^n.
        check(polynomial_content({6,-4,2})==2&&polynomial_content({-6})==6&&polynomial_content({0,0})==1,"content");
        {
            auto p=random_pencil(5,false,true,true,8,100);
            for(auto& v:p.e.data)v*=2;
            for(auto& v:p.f.data)v*=6;
            const auto q=crt_characteristic(p.e,p.f,3).coefficients;
            check(mpz_divisible_2exp_p(polynomial_content(q).get_mpz_t(),5)!=0,"scaled pencil content");
        }

        // solve and verify: every method, with and without content removal, gives the same
        // certificates. M = C^T C and K = C^T diag(18,10,6,4) C / 10 have the spectrum 1.8,1,0.6,0.4,
        // and q = (5z-9)(5z-5)(5z-3)(5z-2) has the content 5.
        {
            ZMatrix c=identity<Z>(4),d(4);
            c(0,1)=1; c(0,3)=-2; c(1,2)=3; c(2,3)=1;
            d(0,0)=18; d(1,1)=10; d(2,2)=6; d(3,3)=4;
            const auto mass=rational(multiply(transpose(c),c),1);
            const auto stiffness=rational(multiply(transpose(c),multiply(d,c)),Q(1,10));
            Options base;
            base.digits=30;
            const auto reference=solve(stiffness,mass,base);
            check(reference.certified&&reference.coefficient_method=="recurrence"&&reference.coefficient_oracle=="modular",
                  "recurrence solve");
            check(reference.roots.size()==4&&reference.roots[0].decimal=="1.8"+std::string(29,'0'),"known spectrum");
            check(reference.coefficient_content==5&&reference.primitive_coefficient_bits_max<reference.max_coefficient_bits,
                  "decimal pencil content");
            for(auto method:{CoefficientMethod::recurrence,CoefficientMethod::crt,CoefficientMethod::both}) {
                for(bool remove:{false,true}) {
                    auto o=base;
                    o.coefficient_method=method;
                    o.remove_content=remove;
                    const auto r=solve(stiffness,mass,o);
                    check(r.certified&&r.certification_method=="polynomial"&&r.content_removed==remove,"method solve");
                    check(r.coefficients==reference.coefficients&&r.coefficient_content==reference.coefficient_content&&
                          r.primitive_coefficient_bits_max==reference.primitive_coefficient_bits_max,"method coefficients");
                    check(r.roots.size()==reference.roots.size(),"method root count");
                    for(std::size_t i=0;i<r.roots.size();++i)
                        check(r.roots[i].decimal==reference.roots[i].decimal&&r.roots[i].lower==reference.roots[i].lower&&
                              r.roots[i].upper==reference.roots[i].upper&&r.roots[i].radius==reference.roots[i].radius,
                              "method certificates differ");
                    if(method!=CoefficientMethod::recurrence)
                        check(r.coefficient_oracle=="crt"&&r.coefficient_moduli>=1&&r.oracle_moduli==r.coefficient_moduli&&
                              r.crt_workers>=1&&r.oracle_seconds==0,"Chinese remaindering proof fields");
                    else check(r.coefficient_moduli==0&&r.crt_workers==0&&r.crt_seconds==0,"recurrence fields");
                    check(method!=CoefficientMethod::crt||(r.recurrence_workers==0&&r.exact_divisions==0),"no recurrence ran");
                    check(method==CoefficientMethod::crt||r.exact_divisions==reference.exact_divisions,"recurrence ran");
                    const auto replay=verify_certificates(stiffness,mass,r.roots,30,o);
                    check(replay.certified&&replay.certification_method=="polynomial"&&replay.content_removed==remove,
                          "method verify");
                }
            }
            // An explicit Bareiss request still checks the Chinese-remaindered coefficients.
            auto bareiss=base;
            bareiss.coefficient_method=CoefficientMethod::crt;
            bareiss.coefficient_oracle=CoefficientOracle::bareiss;
            const auto checked=solve(stiffness,mass,bareiss);
            check(checked.certified&&checked.coefficient_oracle=="bareiss"&&checked.coefficient_moduli>=1,"crt with Bareiss");
            // Unchecked Chinese remaindering is contradictory; the unchecked recurrence uses inertia.
            auto unchecked=base;
            unchecked.full_coefficient_check=false;
            unchecked.remove_content=true;
            const auto inertia=solve(stiffness,mass,unchecked);
            check(inertia.certified&&inertia.certification_method=="exact_inertia"&&inertia.content_removed,"unchecked solve");
            for(auto method:{CoefficientMethod::crt,CoefficientMethod::both}) {
                unchecked.coefficient_method=method;
                invalid([&] {
                    validate_options(unchecked);
                },"unchecked Chinese remaindering accepted");
                invalid([&] {
                    solve(stiffness,mass,unchecked);
                },"unchecked Chinese remaindering solve accepted");
            }
            // The coefficient phase alone, as the coefficients command runs it.
            auto alone=base;
            alone.coefficient_method=CoefficientMethod::crt;
            const auto pencil=integer_pencil(stiffness,mass);
            const auto phase=exact_coefficients(pencil,alone,true);
            check(phase.coefficients==reference.coefficients&&phase.coefficient_identity_verified&&
                  phase.determinant_scale==pencil.scale,"coefficient phase");
        }
        std::cout<<"PASS "<<checks<<" Chinese remaindering and content checks\n";
        return 0;
    }
    catch(const std::exception& ex) {
        std::cerr<<"FAIL after "<<checks<<" checks: "<<ex.what()<<'\n';
        return 1;
    }
}
