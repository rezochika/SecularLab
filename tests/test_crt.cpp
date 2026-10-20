#include "secular/exact.hpp"
#include "secular/solver.hpp"
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
        // The combined residues are a proof the oracle accepts with the same primes.
        const auto report=check_coefficients(p.e,p.f,first.coefficients,2);
        check(report.verified&&(report.method=="bareiss"||report.moduli<=first.moduli+first.skipped_moduli),
              "oracle rejected Chinese remaindering");
        return first;
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
