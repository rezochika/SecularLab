#include "secular/solver.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

using namespace secular;

namespace {
    void require(bool condition,const std::string& message) {
        if(!condition)throw std::runtime_error(message);
    }

    void check_diagonal(const std::string& label,const Q& small,const Q& large,
                        const std::string& small_decimal,const std::string& large_decimal,
                        unsigned precision,unsigned squarings=128) {
        QMatrix k(2);
        k(0,0)=small;
        k(1,1)=large;
        auto m=identity<Q>(2);
        Options o;
        o.digits=1;
        o.precision=o.max_precision=precision;
        o.max_squarings=squarings;
        auto result=solve(k,m,o);
        require(result.certified,label+": "+result.status+" after "+
                std::to_string(result.graeffe_steps)+" squarings: "+result.message);
        require(verify_root_certificates(k,m,result.roots,o.digits),label+": exact certificate replay");
        require(result.roots.size()==2,label+": complete spectrum");
        const auto& high=result.roots[0];
        const auto& low=result.roots[1];
        require(high.lower<=large&&large<=high.upper,label+": largest known root enclosed");
        require(low.lower<=small&&small<=low.upper,label+": smallest known root enclosed");
        require(high.decimal==large_decimal&&low.decimal==small_decimal,label+": expected decimals");
        require(low.upper<high.lower,label+": distinct certificates despite identical rounded output");
    }

    void cluster() {
        for(unsigned precision:{128U,512U}) {
            check_diagonal("cluster work="+std::to_string(precision),Q(1),
                           Q(1)+parse_rational("1e-20"),"1.0","1.0",precision);
        }
    }

    void scale() {
        for(unsigned precision:{128U,512U}) {
            check_diagonal("tiny work="+std::to_string(precision),parse_rational("1e-100"),
                           parse_rational("2e-100"),"0.0","0.0",precision);
            check_diagonal("tiny cluster work="+std::to_string(precision),parse_rational("1e-100"),
                           parse_rational("1e-100")+parse_rational("1e-120"),"0.0","0.0",precision);
            check_diagonal("large work="+std::to_string(precision),parse_rational("1e100"),
                           parse_rational("2e100"),"1"+std::string(100,'0')+".0",
                           "2"+std::string(100,'0')+".0",precision);
        }
    }

    void limits() {
        QMatrix k(2);
        k(0,0)=1;
        k(1,1)=Q(1)+parse_rational("1e-20");
        Options o;
        o.digits=1;
        o.precision=o.max_precision=128;
        o.max_squarings=2;
        auto result=solve(k,identity<Q>(2),o);
        require(!result.certified&&result.status=="NOT_CERTIFIED"&&result.roots.empty(),
                "insufficient squaring budget must remain NOT_CERTIFIED");
    }

    void exact_scalar() {
        // Output precision can exceed candidate working precision when the
        // recovered root is exact; interval proposals must still try it.
        Options o;
        o.digits=50;
        o.precision=o.max_precision=32;
        o.max_squarings=128;
        auto k=identity<Q>(1);
        auto result=solve(k,k,o);
        require(result.certified,"exact scalar must certify beyond working precision");
        require(verify_root_certificates(k,k,result.roots,o.digits),"exact scalar certificate replay");
        require(result.roots[0].decimal=="1."+std::string(50,'0'),"exact scalar expected decimal");
        require(result.roots[0].lower<=1&&result.roots[0].upper>=1,"exact scalar known root enclosed");
    }
}

int main(int argc,char** argv) {
    try {
        const std::string mode=argc>1?argv[1]:"all";
        require(mode=="all"||mode=="cluster"||mode=="scale"||mode=="limits"||mode=="scalar","unknown test mode");
        if(mode=="all"||mode=="cluster")cluster();
        if(mode=="all"||mode=="scale")scale();
        if(mode=="all"||mode=="limits")limits();
        if(mode=="all"||mode=="scalar")exact_scalar();
        std::cout<<"PASS recovery "<<mode<<"\n";
        return 0;
    }
    catch(const std::exception& ex) {
        std::cerr<<"FAIL recovery: "<<ex.what()<<"\n";
        return 1;
    }
}
