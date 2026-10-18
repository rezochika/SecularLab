#include "secular/solver.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <iomanip>
#include <optional>
#include <algorithm>

using namespace secular;
namespace fs=std::filesystem;
namespace {
    void help() {
        std::cout<<R"(SecularLab 1.0 - exact coefficients, Graeffe candidates, exact root certificates

  secular solve --matrix K.mtx [--mass M.mtx] [--digits 50]
                [--precision 256] [--max-precision 2048] [--max-squarings 128]
                [--affine K1.mtx --parameter 1/3]
                [--affine K1.mtx --legendre-degree 128 --argument 1/3]
                [--workers 0..64] [--certification auto|inertia]
                [--oracle auto|bareiss] [--no-coefficient-check] [--out results/new-run]
  secular coefficients --matrix K.mtx [--mass M.mtx] [--affine K1.mtx]
                       [--workers 0..64] [--oracle auto|bareiss] [--out results/new-coefficients]
  secular verify --matrix K.mtx [--mass M.mtx] --cert eigenvalues.tsv
                 [--workers 0..64] [--certification auto|inertia] [--oracle auto|bareiss]
                 [--out results/replay]

Matrix Market real/integer coordinate or array files; exact decimal/scientific
and rational tokens supported. Omitted mass is identity. Affine input means
K(t) = K + t*K1. Solve specializes the rational parameter BEFORE coefficient
construction. Coefficients mode without --parameter retains formal t.

Root solving currently requires symmetric positive definite K and M, and a simple
spectrum. Coefficients mode only requires nonsingular M. Candidate arithmetic is
NOT a certificate. Auto certification uses an exactly verified characteristic
polynomial and a complete disjoint root proof, with exact inertia fallback.
Inertia mode checks the original pencil directly. Decimal rounding is exact
half-away-from-zero. Auto verify derives its polynomial from the input matrices.
Oracle auto proves the coefficient identity modulo proven primes, using n+1 Bareiss
determinants when inconclusive or estimated cheaper; bareiss always uses them.
Workers run fixed-pencil recurrence rows and oracle jobs: 0 means automatic (at most 8);
explicit 1..64 is clamped to the jobs of each phase. Results never depend on it.
--no-coefficient-check forces solve to use inertia and is forbidden for verify.

Exit 0: successful computation/verification; exit 1: unresolved/repeated spectrum
or invalid certificate; exit 2: invalid input, resource or I/O error.
Output directory must be new or empty. Never overwrites existing run results.
)";
    }
    unsigned number(const std::string& s,unsigned max) {
        Q q=parse_rational(s);
        if(q.get_den()!=1||q<0||q>max)throw std::invalid_argument("invalid unsigned option: "+s);
        return static_cast<unsigned>(q.get_num().get_ui());
    }
    std::ofstream output(const fs::path& p) {
        std::ofstream f(p);
        f.exceptions(std::ios::badbit|std::ios::failbit);
        return f;
    }
    void init_output(const fs::path& p) {
        if(fs::exists(p)&&(!fs::is_directory(p)||!fs::is_empty(p)))throw std::runtime_error("output directory must be new or empty: "+p.string());
        fs::create_directories(p);
    }
    void coefficients_file(const fs::path& file,const std::vector<Z>& co,const Z& scale) {
        auto f=output(file);
        f<<"# det(z*(L*M)-L*K), L="<<scale<<"; divide by L^n for original determinant\n# descending powers\n";
        for(std::size_t i=0;i<co.size();++i)f<<co.size()-1-i<<"\t"<<co[i]<<"\n";
        f.close();
    }
}
int main(int argc,char** argv) {
    try {
        if(argc<2||std::string(argv[1])=="--help"||std::string(argv[1])=="help") {
            help();
            return 0;
        }
        std::string mode=argv[1];
        if(mode!="solve"&&mode!="coefficients"&&mode!="verify")throw std::invalid_argument("unknown command: "+mode);
        const std::set<std::string> valued= {
            "--matrix","--mass","--affine","--parameter","--legendre-degree","--argument","--digits","--precision","--max-precision","--max-squarings","--out","--cert","--workers","--certification","--oracle"
        };
        std::map<std::string,std::string> args;
        bool check=true;
        for(int i=2;i<argc;++i) {
            std::string a=argv[i];
            if(a=="--help") {
                help();
                return 0;
            }
            if(a=="--no-coefficient-check") {
                check=false;
                continue;
            }
            if(!valued.count(a)||i+1>=argc)throw std::invalid_argument("unknown option or missing value: "+a);
            if(args.count(a))throw std::invalid_argument("duplicate option: "+a);
            args[a]=argv[++i];
        }
        if(!args.count("--matrix"))throw std::invalid_argument("--matrix is required");
        Options o;
        o.full_coefficient_check=check;
        if(args.count("--workers"))o.oracle_workers=number(args.at("--workers"),64);
        if(args.count("--certification")) {
            const auto& method=args.at("--certification");
            if(method!="auto"&&method!="inertia")throw std::invalid_argument("--certification must be auto or inertia");
            if(mode=="coefficients")throw std::invalid_argument("--certification requires solve or verify");
            o.polynomial_certification=method=="auto";
        }
        if(args.count("--oracle")) {
            const auto& oracle=args.at("--oracle");
            if(oracle!="auto"&&oracle!="bareiss")throw std::invalid_argument("--oracle must be auto or bareiss");
            if(!check||(mode=="verify"&&!o.polynomial_certification))
                throw std::invalid_argument("--oracle requires a coefficient identity check");
            o.coefficient_oracle=oracle=="auto"?CoefficientOracle::automatic:CoefficientOracle::bareiss;
        }
        if(mode=="verify"&&!check)throw std::invalid_argument("--no-coefficient-check is not allowed for verify");
        if(args.count("--digits"))o.digits=number(args.at("--digits"),1000);
        if(args.count("--precision"))o.precision=number(args.at("--precision"),20000);
        if(args.count("--max-precision"))o.max_precision=number(args.at("--max-precision"),20000);
        if(args.count("--max-squarings"))o.max_squarings=number(args.at("--max-squarings"),512);
        if(mode=="solve")validate_options(o);
        auto k=read_matrix_market(fs::path(args.at("--matrix")));
        auto m=args.count("--mass")?read_matrix_market(fs::path(args.at("--mass"))):identity<Q>(k.n);
        if(m.n!=k.n)throw std::invalid_argument("mass dimension differs");
        std::optional<QMatrix> k1;
        if(args.count("--affine")) {
            k1=read_matrix_market(fs::path(args.at("--affine")));
            if(k1->n!=k.n)throw std::invalid_argument("affine dimension differs");
        }
        if(args.count("--parameter")&&args.count("--legendre-degree"))throw std::invalid_argument("choose --parameter or --legendre-degree, not both");
        std::optional<Q> parameter;
        if(args.count("--parameter"))parameter=parse_rational(args.at("--parameter"));
        if(args.count("--legendre-degree"))parameter=legendre_value(number(args.at("--legendre-degree"),10000),args.count("--argument")?parse_rational(args.at("--argument")):Q(1,3));
        if(args.count("--argument")&&!args.count("--legendre-degree"))throw std::invalid_argument("--argument requires --legendre-degree");
        if(parameter&&!k1)throw std::invalid_argument("parameter requires --affine");
        if(k1&&mode!="coefficients"&&!parameter)throw std::invalid_argument("solve/verify affine input requires parameter");
        if(parameter)for(std::size_t i=0;i<k.data.size();++i)k.data[i]+=*parameter*k1->data[i];
        if(mode=="verify") {
            if(!args.count("--cert"))throw std::invalid_argument("verify requires --cert");
            std::ifstream in(args.at("--cert"));
            if(!in)throw std::runtime_error("cannot open certificate file");
            unsigned digits;
            auto roots=read_certificates(in,digits);
            if(args.count("--out"))init_output(fs::path(args.at("--out")));
            auto result=verify_certificates(k,m,roots,digits,o);
            if(args.count("--out")) {
                auto f=output(fs::path(args.at("--out"))/"result.json");
                write_result_json(f,result,digits);
                f.close();
            }
            std::cout<<(result.certified?"VERIFIED":"FAILED")<<" "<<roots.size()<<" eigenvalue certificates, "<<digits<<" decimal places"
                     <<" certification_method="<<result.certification_method<<" total_seconds="<<std::setprecision(6)<<result.total_seconds<<'\n';
            return result.certified?0:1;
        }
        const fs::path dir=args.count("--out")?args.at("--out"):"results/new-run";
        init_output(dir);
        {
            auto f=output(dir/"input_K.mtx");
            write_matrix_market(f,k);
            f.close();
            auto g=output(dir/"input_M.mtx");
            write_matrix_market(g,m);
            g.close();
        }
        {
            auto f=output(dir/"command.txt");
            for(int i=0;i<argc;++i)f<<argv[i]<<'\n';
            if(parameter)f<<"exact_parameter="<<*parameter<<'\n';
            f.close();
        }
        if(mode=="coefficients") {
            if(k1&&!parameter) {
                Z l=1;
                for(const auto* a: {
                    &k,&m,&*k1
                })for(const auto& v:a->data)mpz_lcm(l.get_mpz_t(),l.get_mpz_t(),v.get_den_mpz_t());
                ZMatrix e(k.n),f0(k.n),f1(k.n);
                for(std::size_t i=0;i<k.data.size();++i) {
                    e.data[i]=m.data[i].get_num()*exact_divide(l,m.data[i].get_den());
                    f0.data[i]=k.data[i].get_num()*exact_divide(l,k.data[i].get_den());
                    f1.data[i]=k1->data[i].get_num()*exact_divide(l,k1->data[i].get_den());
                }
                auto co=affine_characteristic(e,f0,f1);
                bool verified=false;
                if(check) {
                    verified=true;
                    for(std::size_t t=0;t<=k.n;++t) {
                        Z value=static_cast<unsigned long>(t);
                        ZMatrix f=f0;
                        for(std::size_t i=0;i<f.data.size();++i)f.data[i]+=value*f1.data[i];
                        std::vector<Z> specialized;
                        for(auto row:co) {
                            std::reverse(row.begin(),row.end());
                            specialized.push_back(evaluate(row,value));
                        }
                        if(!check_coefficients(e,f,specialized,o.oracle_workers,o.coefficient_oracle).verified)
                            throw std::logic_error("affine determinant-grid identity failed");
                    }
                }
                auto f=output(dir/"coefficients_affine.tsv");
                f<<"# det(z*(L*M)-L*K0-t*L*K1), L="<<l<<"\nz_power\tt_power\tinteger_coefficient\n";
                for(std::size_t i=0;i<co.size();++i)for(std::size_t j=0;j<co[i].size();++j)f<<k.n-i<<'\t'<<j<<'\t'<<co[i][j]<<'\n';
                f.close();
                auto g=output(dir/"input_K1.mtx");
                write_matrix_market(g,*k1);
                g.close();
                std::cout<<"EXACT_AFFINE_COEFFICIENTS order="<<k.n<<" determinant_grid_verified="<<(verified?"yes":"no")<<'\n';
                return 0;
            }
            auto p=integer_pencil(k,m);
            // The recurrence runs on the reduced pencil; the oracle checks the original one.
            auto co=reduced_characteristic(p.E,p.F,o.oracle_workers);
            OracleReport oracle;
            if(check)oracle=check_coefficients(p.E,p.F,co.coefficients,o.oracle_workers,o.coefficient_oracle);
            if(check&&!oracle.verified)throw std::logic_error("coefficient identity failed");
            coefficients_file(dir/"coefficients.tsv",co.coefficients,p.scale);
            std::cout<<"EXACT_COEFFICIENTS order="<<k.n<<" determinant_identity_verified="<<(oracle.verified?"yes":"no")
                     <<" coefficient_oracle="<<oracle.method<<'\n';
            return 0;
        }
        o.progress=[](const std::string& s) {
            std::cerr<<s<<std::endl;
        };
        auto result=solve(k,m,o);
        {
            auto f=output(dir/"result.json");
            write_result_json(f,result,o.digits);
            f.close();
        }
        coefficients_file(dir/"coefficients.tsv",result.coefficients,result.determinant_scale);
        if(result.certified) {
            auto f=output(dir/"eigenvalues.tsv");
            write_certificates(f,result.roots,o.digits);
            f.close();
            for(const auto& r:result.roots)std::cout<<"lambda["<<r.index<<"] = "<<r.decimal<<'\n';
        }
        std::cout<<result.status<<" order="<<result.order<<" roots="<<result.roots.size()<<" digits="<<o.digits<<" squarings="<<result.graeffe_steps<<" work_digits="<<result.working_digits<<" certification_method="<<result.certification_method<<" total_seconds="<<std::setprecision(6)<<result.total_seconds<<'\n';
        return result.certified?0:1;
    }
    catch(const std::exception& ex) {
        std::cerr<<"ERROR: "<<ex.what()<<'\n';
        return 2;
    }
}
