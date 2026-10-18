#include "secular/solver.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace secular;

namespace {
    void require(bool condition,const std::string& message) {
        if(!condition)throw std::runtime_error(message);
    }

    SolveResult check_diagonal(const std::string& label,const Q& small,const Q& large,
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
        return result;
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
            // Separating a relative gap of 1e-20 takes about 80 squarings, while the binary
            // exponents of roots near 1e-100 outgrow 64 bits after about 54.
            const auto cluster=check_diagonal("tiny cluster work="+std::to_string(precision),parse_rational("1e-100"),
                                              parse_rational("1e-100")+parse_rational("1e-120"),"0.0","0.0",precision);
            require(cluster.graeffe_unbounded_exponents&&cluster.graeffe_steps>64,
                    "tiny cluster did not continue with unbounded exponents");
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

    // Descending integer coefficients of prod_i (z - roots[i]) times a common denominator.
    std::vector<Z> monic_times(const std::vector<Q>& roots,const Z& denominator) {
        std::vector<Q> c {
            1
        };
        for(const auto& r:roots) {
            std::vector<Q> next(c.size()+1);
            for(std::size_t i=0;i<c.size();++i) {
                next[i]+=c[i];
                next[i+1]-=r*c[i];
            }
            c=std::move(next);
        }
        std::vector<Z> out;
        for(const auto& v:c) {
            Q scaled=v*denominator;
            require(scaled.get_den()==1,"test polynomial denominator does not clear");
            out.push_back(scaled.get_num());
        }
        return out;
    }

    void exponents() {
        // Roots 1e-100 and 1e-100+1e-120, cleared by 10^240.
        const std::vector<Q> cluster {
            parse_rational("1e-100")+parse_rational("1e-120"),parse_rational("1e-100")
        };
        const auto tiny=monic_times(cluster,pow10(240));
        for(unsigned squarings:{40U,80U,128U}) {
            const auto narrow=graeffe_candidates(tiny,128,squarings);
            require(narrow.unbounded_exponents==(squarings>=80),"int64 exponents switched at an unexpected step");
            // Forcing the switch to unbounded exponents at once, or after the first steps,
            // must give bit-identical candidates.
            for(unsigned bits:{1U,12U,40U}) {
                const auto forced=graeffe_candidates(tiny,128,squarings,bits);
                require(forced.unbounded_exponents,"small exponent limit did not switch");
                require(forced.roots==narrow.roots,"candidates depend on the exponent switch");
            }
            if(squarings>=80)for(std::size_t i=0;i<cluster.size();++i) {
                Q error=narrow.roots[i]/cluster[i]-1;
                if(error<0)error=-error;
                require(error<parse_rational("1e-100"),"separated tiny cluster candidate is inaccurate");
            }
        }
        // Separated roots over a wide range converge after a few squarings.
        const std::vector<Q> spread {
            parse_rational("3e40"),Q(7),Q(5),Q(1,3),parse_rational("2e-30")
        };
        const auto wide=monic_times(spread,Z(3)*pow10(30));
        for(unsigned bits:{63U,30U,5U}) {
            const auto c=graeffe_candidates(wide,256,24,bits);
            for(std::size_t i=0;i<spread.size();++i) {
                Q error=c.roots[i]/spread[i]-1;
                if(error<0)error=-error;
                require(error<parse_rational("1e-200"),"separated candidate is inaccurate");
            }
        }
        bool caught=false;
        try {
            (void)graeffe_candidates({1},128,4);
        }
        catch(const std::invalid_argument&) {
            caught=true;
        }
        require(caught,"constant polynomial accepted");
        caught=false;
        try {
            (void)graeffe_candidates({1,1},128,4);
        }
        catch(const std::invalid_argument&) {
            // Request validation must accept this shape; only the sign invariant rejects it.
        }
        catch(const std::logic_error&) {
            caught=true;
        }
        require(caught,"coefficients without the positive-root sign pattern accepted");
    }

    void stable_endpoints() {
        // Proposal widths depend on candidate magnitudes and gaps only through their nearest
        // powers of ten, so every working precision that certifies yields the same intervals.
        auto m=identity<Q>(3);
        QMatrix k(3);
        k(0,0)=parse_rational("1e-100");
        k(1,1)=parse_rational("2e-100");
        k(2,2)=parse_rational("1e-100")+parse_rational("3e-118");
        for(unsigned digits:{1U,8U}) {
            std::string reference;
            for(unsigned precision:{64U,128U,256U,512U,1024U}) {
                Options o;
                o.digits=digits;
                o.precision=o.max_precision=precision;
                auto result=solve(k,m,o);
                require(result.certified,"stable case did not certify at work="+std::to_string(precision));
                std::ostringstream text;
                write_certificates(text,result.roots,digits);
                if(reference.empty())reference=text.str();
                require(text.str()==reference,"certificate endpoints changed with working precision "+std::to_string(precision));
            }
        }
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
        require(mode=="all"||mode=="cluster"||mode=="scale"||mode=="limits"||mode=="scalar"||
                mode=="exponents"||mode=="stable","unknown test mode");
        if(mode=="all"||mode=="cluster")cluster();
        if(mode=="all"||mode=="scale")scale();
        if(mode=="all"||mode=="limits")limits();
        if(mode=="all"||mode=="scalar")exact_scalar();
        if(mode=="all"||mode=="exponents")exponents();
        if(mode=="all"||mode=="stable")stable_endpoints();
        std::cout<<"PASS recovery "<<mode<<"\n";
        return 0;
    }
    catch(const std::exception& ex) {
        std::cerr<<"FAIL recovery: "<<ex.what()<<"\n";
        return 1;
    }
}
