#include "secular/solver.hpp"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace secular;
static void require(bool b,const std::string& message) {
    if(!b)throw std::runtime_error(message);
}
static void compare_methods(const QMatrix& k,unsigned digits,const std::vector<std::string>& expected,bool early_failure=false) {
    auto m=identity<Q>(k.n);
    Options options;
    options.digits=digits;
    options.oracle_workers=2;
    options.precision=options.max_precision=128;
    auto fast=solve(k,m,options);
    require(fast.certified&&fast.certification_method=="polynomial","auto solve uses polynomial proof");
    require(fast.coefficient_identity_verified&&fast.polynomial_sign_evaluations>0,
            "polynomial proof requires verified identity and actual sign checks");
    // Earlier candidates of a cluster fail inside certify_one. Under the verified identity
    // exact inertia provably fails on them too, so it is not tried.
    require(fast.inertia_evaluations==2,"fast solve tried the futile inertia fallback");
    require(early_failure?fast.certification_attempts>1:fast.certification_attempts>=1,
            "cluster had no failed polynomial attempt");
    require(fast.coefficient_oracle=="modular"&&fast.oracle_moduli>=1&&
            fast.oracle_workers==std::min<std::size_t>(2,fast.oracle_moduli)&&fast.recurrence_workers==2&&
            fast.squarefree_modular_prime!=0,"oracle, recurrence workers and modular witness reported");
    for(std::size_t i=0;i<expected.size();++i)require(fast.roots[i].decimal==expected[i],"hand-known rounded eigenvalue");
    auto checked=verify_certificates(k,m,fast.roots,digits,options);
    require(checked.certified&&checked.status=="VERIFIED"&&checked.certification_method=="polynomial",
            "detailed verifier proves full polynomial rootset");
    require(checked.coefficient_identity_verified&&checked.coefficients==fast.coefficients,
            "verifier reconstructs original matrix coefficients");
    require(checked.inertia_evaluations==2&&checked.polynomial_sign_evaluations>0,"fast verify avoids shifted inertia");
    require(checked.roots.size()==k.n,"detailed verifier returns verified roots");
    options.polynomial_certification=false;
    auto slow=solve(k,m,options);
    require(slow.certified&&slow.certification_method=="exact_inertia"&&slow.inertia_evaluations>2,
            "explicit inertia solve uses exact shifts");
    std::ostringstream a,b;
    write_certificates(a,fast.roots,digits);
    write_certificates(b,slow.roots,digits);
    require(a.str()==b.str(),"proof choice preserves all certificate fields");
    auto replay=verify_certificates(k,m,fast.roots,digits,options);
    require(replay.certified&&replay.certification_method=="exact_inertia"&&replay.coefficients.empty(),
            "inertia verifier requires no polynomial reconstruction");
    options.polynomial_certification=true;
    options.full_coefficient_check=false;
    auto unchecked=solve(k,m,options);
    require(unchecked.certified&&unchecked.certification_method=="exact_inertia"&&!unchecked.coefficient_identity_verified,
            "disabled coefficient check forces inertia acceptance");
    checked=verify_certificates(k,m,fast.roots,digits,options);
    require(checked.certified&&checked.certification_method=="polynomial"&&checked.coefficient_identity_verified,
            "auto verifier cannot skip original-pencil coefficient check");
}
static void malformed_certificates() {
    QMatrix k(2);
    k(0,0)=1;
    k(1,1)=3;
    const auto m=identity<Q>(2);
    const std::vector<RootCertificate> good={{0,Q(2999,1000),Q(3001,1000),Q(1,1000),"3.00"},
                                             {1,Q(999,1000),Q(1001,1000),Q(1,1000),"1.00"}};
    for(bool polynomial:{true,false}) {
        Options o;
        o.polynomial_certification=polynomial;
        auto reject=[&](const std::vector<RootCertificate>& roots,const std::string& name) {
            auto result=verify_certificates(k,m,roots,2,o);
            require(!result.certified&&result.status=="FAILED"&&result.roots.empty(),name);
            // A failed exact sign under the verified identity also fails exact inertia, and
            // a failed shared row check stops both before any shifted inertia.
            require(!polynomial||result.inertia_evaluations==2,name+" after futile shifted inertia");
        };
        require(verify_certificates(k,m,good,2,o).certified,"hand-derived rational brackets accepted");
        auto bad=good; bad.pop_back(); reject(bad,"missing root rejected");
        bad=good; bad.push_back(good[1]); reject(bad,"extra root rejected");
        bad=good; std::swap(bad[0],bad[1]); reject(bad,"reordered rows rejected");
        bad=good; bad[1]=bad[0]; bad[1].index=1; reject(bad,"duplicate brackets rejected");
        bad=good; bad[1].upper=good[0].lower; reject(bad,"touching brackets rejected");
        bad=good; bad[0].index=1; reject(bad,"wrong index rejected");
        bad=good; bad[0].decimal="3.0"; reject(bad,"noncanonical decimal rejected");
        bad=good; bad[0].decimal="3.01"; reject(bad,"wrong rounded decimal rejected");
        bad=good; bad[0].radius=0; reject(bad,"insufficient radius rejected");
        bad=good; bad[0].radius=-1; reject(bad,"negative radius rejected");
        bad=good; bad[0].lower=4; reject(bad,"reversed interval rejected");
        bad=good; bad[0].lower=Q(30001,10000); reject(bad,"non-bracketing same-sign interval rejected");
        require(polynomial||verify_certificates(k,m,bad,2,o).inertia_evaluations>2,"inertia replay skipped shifted inertia");
        bad=good; bad[0].lower=bad[0].upper=Q(30001,10000); reject(bad,"false singleton rejected");
        bad=good; bad[0].lower=3; reject(bad,"zero nonsingleton endpoint rejected");
        auto changed=k; changed(1,1)=4;
        require(!verify_certificates(changed,m,good,2,o).certified,"matrix mismatch rejected");
        changed=k; changed(0,1)=1;
        require(!verify_certificates(changed,m,good,2,o).certified,"nonsymmetric matrix rejected");
        changed=m; changed(0,0)=-1;
        require(!verify_certificates(k,changed,good,2,o).certified,"indefinite mass rejected");
        changed=k; changed(0,0)=-1;
        require(!verify_certificates(changed,m,good,2,o).certified,"indefinite stiffness rejected");
        require(!verify_certificates(k,identity<Q>(1),good,2,o).certified,"matrix dimension mismatch rejected");
        require(!verify_certificates(k,m,good,0,o).certified,"invalid precision rejected");
        bad={{0,Q(3),Q(3),Q(0),"3.00"},{1,Q(1),Q(1),Q(0),"1.00"}};
        require(verify_certificates(k,m,bad,2,o).certified,"exact singleton full rootset accepted");
        // Denominators 300 are not decimal: the general evaluator decides those signs.
        const std::vector<RootCertificate> mixed={{0,Q(899,300),Q(901,300),Q(1,300),"3.00"},
                                                  {1,Q(999,1000),Q(1001,1000),Q(1,1000),"1.00"}};
        require(verify_certificates(k,m,mixed,2,o).certified,"nondecimal and decimal brackets accepted");
        bad=mixed; bad[0].lower=Q(9001,3000); bad[0].upper=Q(9002,3000); reject(bad,"nondecimal same-sign interval rejected");
    }
    QMatrix repeated=identity<Q>(2);
    auto result=solve(repeated,m,Options{});
    require(!result.certified&&result.status=="REPEATED_SPECTRUM","repeated solve preserves status");
    const std::vector<RootCertificate> duplicates={{0,Q(1),Q(1),Q(0),"1.00"},{1,Q(1),Q(1),Q(0),"1.00"}};
    require(!verify_certificates(repeated,m,duplicates,2).certified,"repeated verify cannot certify distinct full spectrum");
}
int main() {
    try {
        QMatrix k(2);
        k(0,0)=1;
        k(1,1)=3;
        Options options;
        options.digits=8;
        auto result=solve(k,identity<Q>(2),options);
        require(result.certified,"known simple roots must certify");
        std::ostringstream json;
        write_result_json(json,result,options.digits);
        require(json.str().find("\"certification_method\": \"polynomial\"")!=std::string::npos,
                "default solve must report a polynomial certificate proof");
        // Each root certifies at its first proposal: two endpoint signs, which the final
        // rootset gate reuses instead of evaluating again. A replay has no such memo.
        require(result.polynomial_sign_evaluations==4,"rootset gate re-evaluated certify_one endpoints");
        require(verify_certificates(k,identity<Q>(2),result.roots,options.digits).polynomial_sign_evaluations==4,
                "replay evaluates every endpoint once");
        compare_methods(k,8,{"3.00000000","1.00000000"});
        QMatrix irrational(2);
        irrational(0,0)=3; irrational(1,1)=1;
        irrational(0,1)=irrational(1,0)=1;
        compare_methods(irrational,8,{"3.41421356","0.58578644"});
        QMatrix clustered(2);
        clustered(0,0)=1; clustered(1,1)=Q(1)+parse_rational("1e-20");
        compare_methods(clustered,1,{"1.0","1.0"},true);
        QMatrix ties(2);
        ties(0,0)=Q(1,4); ties(1,1)=Q(5,4);
        compare_methods(ties,1,{"1.3","0.3"});
        Options tie_options;
        tie_options.digits=1;
        auto tied=solve(ties,identity<Q>(2),tie_options);
        require(tied.roots[0].lower==Q(5,4)&&tied.roots[0].upper==Q(5,4)&&
                tied.roots[1].lower==Q(1,4)&&tied.roots[1].upper==Q(1,4),"half ties use exact singleton certificates");
        malformed_certificates();
        std::cout<<"PASS certificate proof regression\n";
    }
    catch(const std::exception& ex) {
        std::cerr<<"FAIL certification: "<<ex.what()<<'\n';
        return 1;
    }
}
