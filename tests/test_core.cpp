#include "secular/exact.hpp"
#include "secular/solver.hpp"
#include <iostream>
#include <functional>
#include <sstream>

using namespace secular;
static int count=0;
void require(bool b,const std::string& m) {
    ++count;
    if(!b)throw std::runtime_error(m);
}
template<class F> void rejected(F f) {
    bool yes=false;
    try {
        f();
    }
    catch(const std::exception&) {
        yes=true;
    }
    require(yes,"expected exception");
}
int main() {
    try {
        require(parse_rational("-1.25e-2")==Q(-1,80),"scientific exact parse");
        require(parse_rational(".75")==Q(3,4),"dyadic parse");
        require(parse_rational("1/-2")==Q(-1,2),"negative denominator");
        require(parse_rational("2D3")==2000,"Fortran exponent");
        for(auto s: {
            "1/0","nan","1.2.3","2e","--2","e9","1/2/3"
        })rejected([&] {
            parse_rational(s);
        });
        ZMatrix a(2);
        a(0,0)=2;
        a(0,1)=1;
        a(1,0)=1;
        a(1,1)=2;
        require(determinant(a)==3,"Bareiss determinant");
        auto adj=adjugate(a);
        require(adj(0,0)==2 && adj(0,1)==-1,"adjugate");
        ZMatrix swap(2);
        swap(0,1)=swap(1,0)=1;
        require(determinant(swap)==-1,"det row pivot");
        auto in=inertia(swap);
        require(in.positive==1 && in.negative==1 && in.zero==0,"2x2 inertia pivot");
        ZMatrix sing(3);
        sing(0,0)=1;
        sing(2,2)=-1;
        in=inertia(sing);
        require(in.positive==1 && in.negative==1 && in.zero==1,"singular inertia");
        ZMatrix e=identity<Z>(2);
        auto cp=characteristic(e,a);
        require(cp.coefficients==std::vector<Z>( {
            1,-4,3
        }),"Faddeev coefficients");
        require(cp.terminal_zero && verify_coefficients(e,a,cp.coefficients),"independent determinant oracle");
        auto bad=cp.coefficients;
        bad.back()+=1;
        require(!verify_coefficients(e,a,bad),"reject bad coefficients");
        auto af=affine_characteristic(e,a,identity<Z>(2));
        require(af[1]==std::vector<Z>( {
            -4,-2
        }) && af[2]==std::vector<Z>( {
            3,4,1
        }),"affine recurrence");
        require(legendre_value(3,Q(1,3))==Q(-11,27),"exact Legendre");
        require(format_decimal(round_scaled(Q(1,4),1),1)=="0.3","positive tie");
        require(format_decimal(round_scaled(Q(-1,4),1),1)=="-0.3","negative tie");
        require(has_repeated_roots( {
            1,-2,1
        }),"repeated detection");
        require(!has_repeated_roots( {
            1,-4,3
        }),"squarefree detection");
        std::istringstream mm("%%MatrixMarket matrix coordinate real symmetric\n% comment\n2 2 3\n1 1 2\n2 1 1\n2 2 2\n");
        auto q=read_matrix_market(mm);
        require(q(0,1)==1 && q(1,0)==1,"MM symmetry");
        std::istringstream invalid("%%MatrixMarket matrix coordinate real general\n2 3 0\n");
        rejected([&] {
            read_matrix_market(invalid);
        });
        Options o;
        o.digits=30;
        o.precision=128;
        o.max_precision=512;
        o.max_squarings=80;
        auto result=solve(q,identity<Q>(2),o);
        require(result.certified && result.roots.size()==2,"solve 2x2");
        require(result.roots[0].decimal=="3.000000000000000000000000000000","largest exact root");
        require(result.roots[1].decimal=="1.000000000000000000000000000000","smallest exact root");
        require(verify_root_certificates(q,identity<Q>(2),result.roots,o.digits),"replay roots");
        auto corrupt=result.roots;
        corrupt[0].decimal="4.000000000000000000000000000000";
        require(!verify_root_certificates(q,identity<Q>(2),corrupt,o.digits),"reject wrong decimal");
        QMatrix tied(2);
        tied(0,0)=Q(1,4);
        tied(1,1)=Q(9,4);
        o.digits=1;
        result=solve(tied,identity<Q>(2),o);
        require(result.certified && result.roots[1].decimal=="0.3","solve rational tie");
        result=solve(identity<Q>(2),identity<Q>(2),o);
        require(!result.certified && result.status=="REPEATED_SPECTRUM","reject repeated roots");
        q(0,0)=-4;
        rejected([&] {
            solve(q,identity<Q>(2),o);
        });
        std::cout<<"PASS "<<count<<" assertions\n";
        return 0;
    }
    catch(const std::exception& ex) {
        std::cerr<<"FAIL after "<<count<<" checks: "<<ex.what()<<"\n";
        return 1;
    }
}
