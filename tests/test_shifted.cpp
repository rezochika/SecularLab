#include "secular/solver.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

using namespace secular;

namespace {
    unsigned checks=0;
    void check(bool value,const std::string& message) {
        ++checks;
        if(!value)throw std::runtime_error(message);
    }
    bool same(const Inertia& a,const Inertia& b) {
        return a.positive==b.positive&&a.negative==b.negative&&a.zero==b.zero;
    }
    Q power_of_two(unsigned long exponent) {
        Z p;
        mpz_setbit(p.get_mpz_t(),exponent);
        return Q(p);
    }
    ZMatrix shifted_matrix(const IntegerPencil& p,const Q& s) {
        ZMatrix a(p.F.n);
        for(std::size_t i=0;i<a.data.size();++i)a.data[i]=s.get_den()*p.F.data[i]-s.get_num()*p.E.data[i];
        return a;
    }
    struct Pencil {
        QMatrix k,m;
    };
    // K = T^T diag(d) T and M = T^T T for a random unit upper-triangular integer T. The generalized
    // eigenvalues are exactly d, because K-sM = T^T (diag(d)-s) T, yet neither matrix is diagonal.
    Pencil congruent_pencil(const std::vector<Q>& d,std::mt19937& generator) {
        const auto n=d.size();
        auto t=identity<Q>(n);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j)t(i,j)=static_cast<long>(generator()%5)-2;
        Pencil p {QMatrix(n),QMatrix(n)};
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)for(std::size_t l=0;l<=std::min(i,j);++l) {
            p.k(i,j)+=t(l,i)*d[l]*t(l,j);
            p.m(i,j)+=t(l,i)*t(l,j);
        }
        return p;
    }
    Inertia expected_inertia(const std::vector<Q>& d,const Q& s) {
        Inertia in;
        for(const auto& x:d) {
            if(x>s)++in.positive;
            else if(x<s)++in.negative;
            else ++in.zero;
        }
        return in;
    }
    std::string certificate_text(const std::vector<RootCertificate>& roots,unsigned digits) {
        std::ostringstream out;
        write_certificates(out,roots,digits);
        return out.str();
    }
    // Solves at several worker counts, where exact inertia must decide. Rows and every counter must
    // equal those of one worker, the serial order of evaluation.
    SolveResult inertia_solves_agree(const QMatrix& k,const QMatrix& m,Options options,const std::string& name) {
        options.oracle_workers=1;
        const auto serial=solve(k,m,options);
        check(serial.certified&&serial.certification_method=="exact_inertia",name+": serial inertia solve certifies");
        for(unsigned workers: {2u,3u,8u,16u}) {
            options.oracle_workers=workers;
            const auto parallel=solve(k,m,options);
            const auto label=name+" at "+std::to_string(workers)+" workers: ";
            check(parallel.certified&&parallel.certification_method=="exact_inertia",label+"certifies");
            check(certificate_text(parallel.roots,options.digits)==certificate_text(serial.roots,options.digits),label+"same rows");
            check(parallel.inertia_evaluations==serial.inertia_evaluations,label+"same inertia count");
            check(parallel.polynomial_sign_evaluations==serial.polynomial_sign_evaluations,label+"same sign count");
            check(parallel.certification_attempts==serial.certification_attempts&&parallel.attempts==serial.attempts&&
                  parallel.graeffe_steps==serial.graeffe_steps,label+"same attempts");
        }
        return serial;
    }
    // Moves row r off its root, keeping its decimal: the shared row checks pass, its inertia fails.
    std::vector<RootCertificate> displaced(std::vector<RootCertificate> roots,std::size_t r) {
        auto& row=roots[r];
        const Q width=row.upper>row.lower?row.upper-row.lower:Q(1,pow10(60));
        row.lower=row.upper+width;
        row.upper=row.lower+width;
        const Q printed=parse_rational(row.decimal);
        const Q below=abs(printed-row.lower),above=abs(printed-row.upper);
        row.radius=below>above?below:above;
        return roots;
    }
    void ordered_runs() {
        const std::size_t jobs=64;
        for(unsigned workers: {1u,2u,3u,8u,16u}) {
            const auto label="run_ordered at "+std::to_string(workers)+" workers: ";
            // fail and fault list the indices whose job returns false or throws "fault <i>".
            auto run=[&](std::vector<std::size_t> fail,std::vector<std::size_t> fault,std::size_t stop,int thrown) {
                std::vector<std::atomic<int>> runs(jobs);
                const auto outcome=run_ordered(jobs,workers,[&](std::size_t i) {
                    ++runs[i];
                    // The stopping jobs finish last, so later indices complete before them.
                    const bool stops=std::find(fail.begin(),fail.end(),i)!=fail.end()||
                                     std::find(fault.begin(),fault.end(),i)!=fault.end();
                    std::this_thread::sleep_for(std::chrono::microseconds(stops?3000:(i%5)*100));
                    if(std::find(fault.begin(),fault.end(),i)!=fault.end())throw std::runtime_error("fault "+std::to_string(i));
                    return std::find(fail.begin(),fail.end(),i)==fail.end();
                });
                check(outcome.stop==stop,label+"stop index");
                for(std::size_t i=0;i<jobs;++i) {
                    if(i<=stop&&i<jobs)check(runs[i]==1,label+"every job up to the stop runs once");
                    else check(runs[i]<=1,label+"no job runs twice");
                }
                if(thrown<0)check(!outcome.error,label+"no exception");
                else {
                    check(static_cast<bool>(outcome.error),label+"exception returned");
                    try {
                        std::rethrow_exception(outcome.error);
                    }
                    catch(const std::runtime_error& e) {
                        check(e.what()==std::string("fault ")+std::to_string(thrown),label+"lowest exception kept");
                    }
                }
            };
            run({},{},jobs,-1);
            run({40},{},40,-1);
            run({40,50},{},40,-1);
            run({40},{17},17,17);
            run({9},{17},9,-1);
            run({},{30,50},30,30);
            run({0},{1},0,-1);
            run({},{63},63,63);
        }
        check(run_ordered(0,8,[](std::size_t)->bool {throw std::logic_error("no job exists");}).stop==0,"empty run");
    }
}

int main(int argc,char** argv) {
    try {
        // The precision schedule: 64+bits(den), rounded up to whole limbs, then doubled twice,
        // never beyond the congruence limit.
        using List=std::vector<unsigned>;
        check(shifted_precisions(Q(3))==List({128,256,512}),"integer shift starts at two limbs");
        check(shifted_precisions(Q(1,pow10(58)))==List({320,640,1280}),"58-place endpoint starts at 320 bits");
        check(shifted_precisions(Q(1)/power_of_two(63))==List({128,256,512}),"64-bit denominator fits two limbs");
        check(shifted_precisions(Q(1)/power_of_two(64))==List({192,384,768}),"65-bit denominator needs three limbs");
        check(shifted_precisions(Q(1)/power_of_two(8127))==List({8192,16384}),"schedule capped at the congruence limit");
        check(shifted_precisions(Q(1)/power_of_two(16256))==List({16384}),"one proposal at the limit");
        check(shifted_precisions(Q(1)/power_of_two(16384)).empty(),"beyond the limit only elimination remains");

        std::mt19937 generator(20261018);
        const std::vector<Q> spectrum={Q(7),Q(5),Q(3),Q(2),Q(1),Q(1,2),Q(1,3),Q(1,7)};
        const auto pencil=congruent_pencil(spectrum,generator);
        const auto p=integer_pencil(pencil.k,pencil.m);
        check(!diagonal(p.F)&&!diagonal(p.E),"synthetic pencil is not diagonal");
        // Fast path: endpoints within 10^-50 and 10^-200 of each eigenvalue are proved by the first
        // proposal of the schedule, and agree with elimination.
        for(const auto& d:spectrum)for(unsigned places: {50u,200u})for(int side: {-1,1}) {
            const Q s=d+Q(side,pow10(places));
            const auto report=shifted_inertia_report(p.F,p.E,s);
            check(report.method=="congruence"&&report.precision_bits==shifted_precisions(s).front(),
                  "near-singular endpoint proved by the first proposal");
            check(same(report.inertia,expected_inertia(spectrum,s)),"congruence inertia counts the eigenvalues above the shift");
            check(same(report.inertia,inertia(shifted_matrix(p,s))),"congruence agrees with elimination");
            check(same(shifted_inertia(p.F,p.E,s),report.inertia),"shifted_inertia returns the report's inertia");
        }
        // Fallback: an exact eigenvalue makes the endpoint matrix singular, so no congruence can
        // pass and elimination finds the nullity.
        for(const auto& d:spectrum) {
            const auto report=shifted_inertia_report(p.F,p.E,d);
            check(report.method=="elimination"&&report.precision_bits==0,"singular endpoint falls back to elimination");
            check(same(report.inertia,expected_inertia(spectrum,d))&&report.inertia.zero==1,"elimination proves nullity one");
        }
        // Fallback for a nonsingular matrix: the pivot of [[1,N],[N,N^2+1]], N=2^1000, is lost at
        // every scheduled precision (128, 256, 512 bits for an integer shift).
        {
            IntegerPencil lost {ZMatrix(2),ZMatrix(2),1};
            Z big;
            mpz_setbit(big.get_mpz_t(),1000);
            lost.F(0,0)=1;
            lost.F(0,1)=lost.F(1,0)=big;
            lost.F(1,1)=big*big+1;
            lost.E(0,0)=lost.E(1,1)=1;
            const auto report=shifted_inertia_report(lost.F,lost.E,Q(0));
            check(report.method=="elimination"&&same(report.inertia,{2,0,0}),"lost pivot falls back to elimination");
            // A shift whose schedule is empty goes to elimination directly.
            IntegerPencil coupled {ZMatrix(2),ZMatrix(2),1};
            coupled.F(0,0)=coupled.F(1,1)=2;
            coupled.F(0,1)=coupled.F(1,0)=1;
            coupled.E(0,0)=coupled.E(1,1)=1;
            const auto tiny=shifted_inertia_report(coupled.F,coupled.E,Q(1)/power_of_two(16384));
            check(tiny.method=="elimination"&&tiny.precision_bits==0&&same(tiny.inertia,{2,0,0}),"empty schedule uses elimination");
            // At s=2 the endpoint matrix [[0,1],[1,0]] has a zero pivot; at s=3/2 it is [[1,2],[2,1]].
            auto between=shifted_inertia_report(coupled.F,coupled.E,Q(2));
            check(between.method=="elimination"&&same(between.inertia,{1,1,0}),"a zero floating pivot falls back to elimination");
            between=shifted_inertia_report(coupled.F,coupled.E,Q(3,2));
            check(between.method=="congruence"&&same(between.inertia,{1,1,0}),"coupled pencil between its eigenvalues 1 and 3");
            IntegerPencil diagonal_pencil {ZMatrix(2),ZMatrix(2),1};
            diagonal_pencil.F(0,0)=1;
            diagonal_pencil.F(1,1)=3;
            diagonal_pencil.E(0,0)=diagonal_pencil.E(1,1)=1;
            const auto signs=shifted_inertia_report(diagonal_pencil.F,diagonal_pencil.E,Q(2));
            check(signs.method=="diagonal"&&signs.precision_bits==0&&same(signs.inertia,{1,1,0}),"diagonal pencil counts signs");
        }

        ordered_runs();

        // Parallel inertia certification reproduces the serial rows and counters. Half ties at one
        // place are exact singletons: a midpoint sign, then a singular endpoint that falls back to
        // elimination. A cluster at the top fails its first attempts inside certify_one while the
        // roots below it would pass, so the counters must stop at the failing root.
        {
            Options options;
            options.digits=1;
            options.precision=options.max_precision=128;
            options.polynomial_certification=false;
            const std::vector<Q> ties={Q(13,4),Q(9,4),Q(2),Q(5,4),Q(1,2),Q(1,4)};
            const auto tied=congruent_pencil(ties,generator);
            const auto serial=inertia_solves_agree(tied.k,tied.m,options,"ties");
            check(serial.polynomial_sign_evaluations>0,"ties evaluate midpoint signs");
            check(serial.roots[0].lower==Q(13,4)&&serial.roots[0].upper==Q(13,4)&&serial.roots[0].decimal=="3.3",
                  "a tie is an exact singleton");
            options.full_coefficient_check=false;
            options.polynomial_certification=true;
            const auto unchecked=inertia_solves_agree(tied.k,tied.m,options,"unchecked ties");
            check(!unchecked.coefficient_identity_verified,"unchecked solve uses inertia");
            check(certificate_text(unchecked.roots,1)==certificate_text(serial.roots,1),"unchecked rows equal checked rows");
            options=Options {};
            options.digits=1;
            options.precision=options.max_precision=128;
            options.polynomial_certification=false;
            const std::vector<Q> cluster={Q(3)+Q(1,pow10(20)),Q(3),Q(2),Q(3,2),Q(1),Q(1,2),Q(1,3),Q(1,4)};
            const auto clustered=congruent_pencil(cluster,generator);
            const auto result=inertia_solves_agree(clustered.k,clustered.m,options,"cluster");
            check(result.certification_attempts>1,"cluster fails an attempt before certifying");
        }

        // The oil-rig pencil (n=66, dense K, diagonal M).
        if(argc>1) {
            const std::filesystem::path root(argv[1]);
            const auto k=read_matrix_market(root/"matrices/oilrig/bcsstk02.mtx");
            const auto m=read_matrix_market(root/"matrices/oilrig/bcsstm02.mtx");
            const auto oil=integer_pencil(k,m);
            Options options;
            options.oracle_workers=8;
            const auto fast=solve(k,m,options);
            check(fast.certified&&fast.certification_method=="polynomial"&&fast.roots.size()==66,"oil rig certifies");
            // Every certificate endpoint is proved by congruence, without elimination.
            for(const auto& r:fast.roots)for(const Q* x: {&r.lower,&r.upper}) {
                const auto report=shifted_inertia_report(oil.F,oil.E,*x);
                check(report.method=="congruence","oil-rig endpoint proved by congruence");
                check(report.inertia.zero==0&&report.inertia.positive==(x==&r.lower?r.index+1:r.index),
                      "oil-rig endpoint counts the roots above it");
            }
            const Q& first=fast.roots[0].lower;
            check(same(shifted_inertia(oil.F,oil.E,first),inertia(shifted_matrix(oil,first))),"oil-rig congruence agrees with elimination");
            // Exact inertia gives the same rows and counters for every worker count.
            options.polynomial_certification=false;
            const auto serial=inertia_solves_agree(k,m,options,"oil rig");
            check(certificate_text(serial.roots,50)==certificate_text(fast.roots,50),"inertia rows equal polynomial rows");
            check(serial.inertia_evaluations==2+2*66,"oil rig needs two inertias per root");
            // A replay stops counting at the first failing row, whatever the worker count.
            const std::size_t r=20;
            const auto bad=displaced(fast.roots,r);
            for(unsigned workers: {1u,3u,8u,16u}) {
                Options replay;
                replay.polynomial_certification=false;
                replay.oracle_workers=workers;
                const auto result=verify_certificates(k,m,bad,50,replay);
                check(!result.certified&&result.status=="FAILED","displaced row rejected");
                check(result.inertia_evaluations==2+2*(r+1),"replay counts the rows up to the failing one");
                replay.polynomial_certification=true;
                const auto signed_result=verify_certificates(k,m,bad,50,replay);
                check(!signed_result.certified&&signed_result.inertia_evaluations==2,"a failed sign skips inertia");
            }
            Options replay;
            replay.polynomial_certification=false;
            replay.oracle_workers=16;
            const auto good=verify_certificates(k,m,fast.roots,50,replay);
            check(good.certified&&good.certification_method=="exact_inertia"&&good.inertia_evaluations==2+2*66,
                  "parallel inertia replay verifies the oil rig");
        }
        std::cout<<"shifted inertia tests passed: "<<checks<<" checks\n";
        return 0;
    }
    catch(const std::exception& e) {
        std::cerr<<"FAILED: "<<e.what()<<"\n";
        return 1;
    }
}
