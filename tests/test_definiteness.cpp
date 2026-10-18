#include "secular/exact.hpp"
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

using namespace secular;

namespace {
    bool same(const Inertia& a,const Inertia& b) {
        return a.positive==b.positive&&a.negative==b.negative&&a.zero==b.zero;
    }
    ZMatrix transpose(const ZMatrix& a) {
        ZMatrix t(a.n);
        for(std::size_t i=0;i<a.n;++i)for(std::size_t j=0;j<a.n;++j)t(i,j)=a(j,i);
        return t;
    }
    ZMatrix from_rows(std::initializer_list<std::initializer_list<long>> rows) {
        ZMatrix a(rows.size());
        std::size_t i=0;
        for(const auto& row:rows) {
            std::size_t j=0;
            for(long v:row)a(i,j++)=v;
            ++i;
        }
        return a;
    }
    // lcm(1..2n-1)/(i+j+1): a positive definite integer Hilbert matrix, badly conditioned.
    ZMatrix hilbert(std::size_t n) {
        Z l=1;
        for(unsigned long k=1;k<2*n;++k)mpz_lcm_ui(l.get_mpz_t(),l.get_mpz_t(),k);
        ZMatrix h(n);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)h(i,j)=l/Z(static_cast<unsigned long>(i+j+1));
        return h;
    }
    // [[1,N],[N,N^2+1]] with N=2^bits: determinant 1, positive definite. Its floating LDL^T
    // pivot N^2+1-N^2 is lost below 2*bits+1 bits of precision.
    ZMatrix cancelling(unsigned long bits) {
        Z big;
        mpz_setbit(big.get_mpz_t(),bits);
        ZMatrix a(2);
        a(0,0)=1;
        a(0,1)=a(1,0)=big;
        a(1,1)=big*big+1;
        return a;
    }
}

int main(int argc,char** argv) {
    unsigned checks=0;
    auto check=[&](bool value,const std::string& message) {
        ++checks;
        if(!value)throw std::runtime_error(message);
    };
    auto invalid=[&](auto action,const char* message) {
        bool caught=false;
        try { action(); } catch(const std::invalid_argument&) { caught=true; }
        check(caught,message);
    };
    try {
        // The exact checker alone, on hand-made proposals. [[1,2],[2,1]] has eigenvalues 3, -1.
        const auto coupled=from_rows({{1,2},{2,1}});
        auto certified=certify_congruence(coupled,from_rows({{1,-2},{0,1}}));
        check(certified&&same(*certified,{1,1,0}),"exact L^-T congruence certifies the inertia");
        check(!certify_congruence(coupled,identity<Z>(2)),"a non-dominant B is rejected");
        check(!certify_congruence(coupled,from_rows({{1,0},{-2,1}})),"a lower-triangular proposal is rejected");
        check(!certify_congruence(coupled,from_rows({{1,-2},{0,0}})),"a singular proposal is rejected");
        // L^-T*diag(3,5) gives B = diag(9,-75).
        certified=certify_congruence(coupled,from_rows({{3,-10},{0,5}}));
        check(certified&&same(*certified,{1,1,0}),"any nonsingular triangular scaling is accepted");
        // Dominance must be strict: [[2,1],[1,1]] is positive definite, but its second row ties.
        check(!certify_congruence(from_rows({{2,1},{1,1}}),identity<Z>(2)),"a dominance tie is rejected");
        certified=certify_congruence(from_rows({{2,1},{1,1}}),from_rows({{1,-1},{0,2}}));
        check(certified&&same(*certified,{2,0,0}),"a dominant congruence of the tie certifies");
        // A singular matrix can never pass: B is singular too, so some row cannot dominate.
        check(!certify_congruence(from_rows({{1,1},{1,1}}),from_rows({{1,-1},{0,1}})),"a singular matrix is never certified");
        invalid([&] { (void)certify_congruence(from_rows({{1,2},{3,1}}),identity<Z>(2)); },"nonsymmetric matrix accepted");
        invalid([&] { (void)certify_congruence(coupled,identity<Z>(3)); },"proposal dimension mismatch accepted");
        invalid([&] { (void)congruence_inertia(from_rows({{1,2},{3,1}}),64); },"nonsymmetric matrix proposed");
        invalid([&] { (void)congruence_inertia(coupled,8); },"too little proposal precision accepted");
        invalid([&] { (void)check_positive_definite(from_rows({{1,2},{3,1}})); },"nonsymmetric definiteness accepted");

        // Random congruences T^T D T with unit upper-triangular T have the known inertia of D.
        // For nonsingular D the LDL^T factors are T^T and D, so the floating factorization is
        // exact and every proposal must certify. A singular D must never be certified.
        std::mt19937 generator(20261018);
        auto small=[&](int range) { return static_cast<long>(generator()%(2*range+1))-range; };
        for(std::size_t n=1;n<=10;++n)for(unsigned trial=0;trial<24;++trial) {
            auto t=identity<Z>(n);
            for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j)t(i,j)=small(3);
            ZMatrix d(n);
            Inertia expected;
            for(std::size_t i=0;i<n;++i) {
                d(i,i)=small(5);
                if(d(i,i)>0)++expected.positive;
                else if(d(i,i)<0)++expected.negative;
                else ++expected.zero;
            }
            const auto a=multiply(multiply(transpose(t),d),t);
            check(same(inertia(a),expected),"elimination agrees with the constructed inertia");
            for(unsigned precision: {64u,256u}) {
                const auto in=congruence_inertia(a,precision);
                if(expected.zero)check(!in,"congruence certified a singular matrix");
                else check(in&&same(*in,expected),"congruence missed or miscounted a nonsingular congruence");
            }
            const auto report=check_positive_definite(a);
            check(report.positive_definite==(expected.positive==n),"definiteness disagrees with the constructed inertia");
            if(expected.zero)check(report.method=="elimination"||report.method=="diagonal","a singular matrix was certified by congruence");
        }

        // Diagonal matrices are decided by their signs.
        auto diagonal_matrix=from_rows({{1,0,0},{0,2,0},{0,0,3}});
        auto report=check_positive_definite(diagonal_matrix);
        check(report.positive_definite&&report.method=="diagonal"&&report.precision_bits==0,"diagonal positive definite");
        diagonal_matrix(1,1)=0;
        check(!check_positive_definite(diagonal_matrix).positive_definite,"diagonal with a zero is not positive definite");
        diagonal_matrix(1,1)=-1;
        check(!check_positive_definite(diagonal_matrix).positive_definite,"diagonal with a negative is not positive definite");

        // The fast path proves both answers: a banded stiffness-like matrix with entries spread
        // over 30 decades is positive definite, and the coupled pair is not.
        {
            const std::size_t n=60,band=5;
            auto t=identity<Z>(n);
            ZMatrix d(n);
            for(std::size_t i=0;i<n;++i) {
                d(i,i)=pow10(static_cast<unsigned>(i%31))*(1+static_cast<long>(generator()%9));
                for(std::size_t j=i+1;j<n&&j<=i+band;++j)t(i,j)=small(4);
            }
            const auto banded=multiply(multiply(transpose(t),d),t);
            report=check_positive_definite(banded);
            check(report.positive_definite&&report.method=="congruence"&&report.precision_bits==64,
                  "banded positive definite matrix proved at the first precision");
            check(inertia(banded).positive==n,"banded matrix is positive definite by elimination");
        }
        report=check_positive_definite(coupled);
        check(!report.positive_definite&&report.method=="congruence"&&report.precision_bits==64,
              "an indefinite matrix is refuted by congruence");

        // Precision escalation: Hilbert matrices need more than 64 bits, and n=40 needs 256.
        report=check_positive_definite(hilbert(20));
        check(report.positive_definite&&report.method=="congruence"&&report.precision_bits>64,
              "Hilbert 20 escalates beyond 64 bits");
        check(!congruence_inertia(hilbert(20),64),"Hilbert 20 is not certified at 64 bits");
        report=check_positive_definite(hilbert(40));
        check(report.positive_definite&&report.method=="congruence"&&report.precision_bits==256,
              "Hilbert 40 is certified at 256 bits");

        // Fallback: when no proposal certifies, the fraction-free elimination decides.
        for(unsigned precision: {64u,128u,256u})
            check(!congruence_inertia(cancelling(1000),precision),"a lost floating pivot is not certified");
        report=check_positive_definite(cancelling(1000));
        check(report.positive_definite&&report.method=="elimination"&&report.precision_bits==0,
              "elimination proves the cancelling matrix positive definite");
        report=check_positive_definite(hilbert(70));
        check(report.positive_definite&&report.method=="elimination","elimination proves Hilbert 70 positive definite");
        report=check_positive_definite(from_rows({{1,1},{1,1}}));
        check(!report.positive_definite&&report.method=="elimination","a singular semidefinite matrix is refuted by elimination");
        report=check_positive_definite(from_rows({{0,1},{1,0}}));
        check(!report.positive_definite&&report.method=="elimination","a zero pivot falls back to elimination");

        // The oil-rig pencil (dense K, diagonal M): K is proved at 64 bits, and a shift inside
        // the spectrum has the elimination's inertia.
        if(argc>1) {
            const std::filesystem::path root(argv[1]);
            const auto pencil=integer_pencil(read_matrix_market(root/"matrices/oilrig/bcsstk02.mtx"),
                                             read_matrix_market(root/"matrices/oilrig/bcsstm02.mtx"));
            report=check_positive_definite(pencil.F);
            check(report.positive_definite&&report.method=="congruence"&&report.precision_bits==64,
                  "oil-rig K proved by congruence");
            report=check_positive_definite(pencil.E);
            check(report.positive_definite&&report.method=="diagonal","oil-rig M proved by its diagonal");
            ZMatrix shifted(pencil.F.n);
            Z trace=0,mass=0;
            for(std::size_t i=0;i<pencil.F.n;++i) {
                trace+=pencil.F(i,i);
                mass+=pencil.E(i,i);
            }
            // trace(F)/trace(E) lies strictly inside the spectrum of a non-scalar pencil.
            for(std::size_t i=0;i<shifted.data.size();++i)shifted.data[i]=mass*pencil.F.data[i]-trace*pencil.E.data[i];
            const auto expected=inertia(shifted);
            check(expected.positive>0&&expected.negative>0,"the oil-rig shift is indefinite");
            const auto in=congruence_inertia(shifted,256);
            check(in&&same(*in,expected),"oil-rig shifted inertia certified by congruence");
            check(!check_positive_definite(shifted).positive_definite,"the oil-rig shift is not positive definite");
        }
        std::cout<<"definiteness tests passed: "<<checks<<" checks\n";
        return 0;
    }
    catch(const std::exception& e) {
        std::cerr<<"FAILED: "<<e.what()<<"\n";
        return 1;
    }
}
