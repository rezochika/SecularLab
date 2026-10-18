#include "secular/exact.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <sstream>
#include <limits>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>

namespace secular {
    namespace {
        std::string trim(std::string s) {
            auto a=s.find_first_not_of(" \t\r\n");
            if(a==std::string::npos)return {};
            return s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);
        }
        Z integer(const std::string& s) {
            if(s.empty())throw std::invalid_argument("empty integer");
            std::size_t i=(s[0]=='-' || s[0]=='+')?1:0;
            if(i==s.size())throw std::invalid_argument("missing integer digits");
            for(std::size_t j=i;j<s.size();++j)if(!std::isdigit(static_cast<unsigned char>(s[j])))throw std::invalid_argument("invalid integer: "+s);
            Z z;
            auto t=s[0]=='+'?s.substr(1):s;
            if(z.set_str(t,10)!=0)throw std::invalid_argument("invalid integer");
            return z;
        }
        std::string data_line(std::istream& s) {
            std::string line;
            while(std::getline(s,line)) {
                line=trim(line);
                if(!line.empty()&&line[0]!='%')return line;
            }
            throw std::invalid_argument("unexpected end of matrix file");
        }
        void swap_symmetric(ZMatrix& a,std::size_t i,std::size_t j) {
            if(i==j)return;
            for(std::size_t k=0;k<a.n;++k)std::swap(a(i,k),a(j,k));
            for(std::size_t k=0;k<a.n;++k)std::swap(a(k,i),a(k,j));
        }
        void swap_symmetric(QMatrix& a,std::size_t i,std::size_t j) {
            if(i==j)return;
            for(std::size_t k=0;k<a.n;++k)std::swap(a(i,k),a(j,k));
            for(std::size_t k=0;k<a.n;++k)std::swap(a(k,i),a(k,j));
        }
        Inertia rational_inertia(QMatrix a) {
            Inertia out;
            const auto n=a.n;
            for(std::size_t k=0;k<n;) {
                std::size_t p=k;
                while(p<n&&a(p,p)==0)++p;
                if(p<n) {
                    swap_symmetric(a,k,p);
                    Q pivot=a(k,k);
                    if(pivot>0)++out.positive;
                    else ++out.negative;
                    for(std::size_t i=k+1;i<n;++i)for(std::size_t j=i;j<n;++j) {
                        Q x=a(i,j)-a(i,k)*a(k,j)/pivot;
                        a(i,j)=x;
                        a(j,i)=x;
                    }
                    ++k;
                    continue;
                }
                // No nonzero diagonal: a nonzero off-diagonal gives a [[0,b],[b,0]]
                // pivot. Its inertia is (1,1,0); this also handles exact-zero pivots.
                std::size_t i=k,j=k;
                bool found=false;
                for(i=k;i<n&&!found;++i)for(j=i+1;j<n;++j)if(a(i,j)!=0) {
                    found=true;
                    break;
                }
                if(!found) {
                    out.zero+=n-k;
                    break;
                }
                --i;
                swap_symmetric(a,k,i);
                swap_symmetric(a,k+1,j);
                Q b=a(k,k+1);
                ++out.positive;
                ++out.negative;
                for(std::size_t r=k+2;r<n;++r)for(std::size_t c=r;c<n;++c) {
                    Q v=a(r,c)-(a(r,k)*a(k+1,c)+a(r,k+1)*a(k,c))/b;
                    a(r,c)=v;
                    a(c,r)=v;
                }
                k+=2;
            }
            return out;
        }
        ZMatrix apply_e_inverse(const ZMatrix& e,const ZMatrix& j,const Z& d,const ZMatrix& t,std::size_t* count) {
            ZMatrix out(e.n);
            if(diagonal(e)) {
                // Algebraically identical to adj(E)*t/d, but avoids multiplying by
                // det(E)/E_ii and cancelling it on every row of a diagonal mass matrix.
                for(std::size_t r=0;r<e.n;++r)for(std::size_t c=0;c<e.n;++c) {
                    out(r,c)=exact_divide(t(r,c),e(r,r));
                    if(count)++*count;
                }
            }
            else {
                out=multiply(j,t);
                for(auto& v:out.data) {
                    v=exact_divide(v,d);
                    if(count)++*count;
                }
            }
            return out;
        }
        void primitive(std::vector<Z>& p) {
            while(!p.empty()&&p.back()==0)p.pop_back();
            if(p.empty())return;
            Z g=0;
            for(const auto& v:p)mpz_gcd(g.get_mpz_t(),g.get_mpz_t(),v.get_mpz_t());
            if(p.back()<0)g=-g;
            if(g!=1)for(auto& v:p)v=exact_divide(v,g);
        }
        // All residues and products are uint64_t: unsigned long is only 32 bits
        // on Windows, and (65537-1)^2 already exceeds that range.
        using Residue=std::uint64_t;
        using ModularPolynomial=std::vector<Residue>;
        void trim_modular(ModularPolynomial& p) {
            while(!p.empty()&&p.back()==0)p.pop_back();
        }
        Residue modular_power(Residue base,Residue exponent,Residue prime) {
            Residue result=1;
            while(exponent) {
                if(exponent&1)result=result*base%prime;
                base=base*base%prime;
                exponent>>=1;
            }
            return result;
        }
        bool modular_squarefree(const std::vector<Z>& primitive_polynomial,unsigned long prime) {
            ModularPolynomial a;
            a.reserve(primitive_polynomial.size());
            for(const auto& v:primitive_polynomial)a.push_back(mpz_fdiv_ui(v.get_mpz_t(),prime));
            // A degree drop cannot prove anything about the original polynomial.
            if(a.back()==0)return false;
            ModularPolynomial b;
            b.reserve(a.size()-1);
            for(std::size_t i=1;i<a.size();++i)b.push_back((i%prime)*a[i]%prime);
            trim_modular(b);
            while(!b.empty()) {
                const auto inverse=modular_power(b.back(),prime-2,prime);
                while(a.size()>=b.size()) {
                    const auto offset=a.size()-b.size();
                    const auto factor=a.back()*inverse%prime;
                    for(std::size_t j=0;j<b.size();++j)
                        a[offset+j]=(a[offset+j]+prime-factor*b[j]%prime)%prime;
                    trim_modular(a);
                }
                a.swap(b);
            }
            return a.size()==1;
        }
        // A persistent team for indexed jobs: the calling thread plus size()-1 helpers.
        // Indices are claimed in increasing order and claiming stops soon after a failure.
        // Every claimed index runs, so every index below the lowest failing one has run, and
        // that exception is rethrown, as a serial loop would report it. Exact results never
        // depend on the schedule.
        class Team {
        public:
            explicit Team(unsigned workers) {
                try {
                    helpers.reserve(workers>1?workers-1:0);
                    for(unsigned slot=1;slot<workers;++slot)helpers.emplace_back([this,slot] {
                        serve(slot);
                    });
                }
                catch(...) {
                    stop();
                    throw;
                }
            }
            Team(const Team&)=delete;
            Team& operator=(const Team&)=delete;
            ~Team() {
                stop();
            }
            unsigned size()const {
                return static_cast<unsigned>(helpers.size()+1);
            }
            // Runs body(index,slot) for every index below count; slot < size() selects
            // per-thread scratch space. Returns after every claimed job has finished.
            void run(std::size_t count,const std::function<void(std::size_t,unsigned)>& body) {
                if(count==0)return;
                errors.assign(count,nullptr);
                job=&body;
                jobs=count;
                next.store(0);
                failed.store(false);
                {
                    std::lock_guard lock(mutex);
                    active=helpers.size();
                    ++generation;
                }
                start.notify_all();
                work(0);
                {
                    std::unique_lock lock(mutex);
                    done.wait(lock,[&] {
                        return active==0;
                    });
                }
                job=nullptr;
                for(const auto& error:errors)if(error)std::rethrow_exception(error);
            }
        private:
            void work(unsigned slot) {
                while(!failed.load()) {
                    const auto index=next.fetch_add(1);
                    if(index>=jobs)break;
                    try {
                        (*job)(index,slot);
                    }
                    catch(...) {
                        errors[index]=std::current_exception();
                        failed.store(true);
                    }
                }
            }
            void serve(unsigned slot) {
                std::size_t seen=0;
                for(;;) {
                    {
                        std::unique_lock lock(mutex);
                        start.wait(lock,[&] {
                            return stopping||generation!=seen;
                        });
                        if(stopping)return;
                        seen=generation;
                    }
                    work(slot);
                    std::lock_guard lock(mutex);
                    if(--active==0)done.notify_one();
                }
            }
            void stop() {
                {
                    std::lock_guard lock(mutex);
                    stopping=true;
                }
                start.notify_all();
                helpers.clear();
            }
            std::mutex mutex;
            std::condition_variable start,done;
            std::size_t generation=0,active=0,jobs=0;
            bool stopping=false;
            const std::function<void(std::size_t,unsigned)>* job=nullptr;
            std::atomic<std::size_t> next{0};
            std::atomic<bool> failed{false};
            std::vector<std::exception_ptr> errors;
            // Last member: helpers are joined before anything they use is destroyed.
            std::vector<std::jthread> helpers;
        };
        // Checked exact quotient written into an existing integer, keeping its allocation.
        void exact_divide_into(Z& q,const Z& a,const Z& b) {
            if(b==0||!mpz_divisible_p(a.get_mpz_t(),b.get_mpz_t()))throw std::logic_error("non-exact integer division");
            mpz_divexact(q.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
        }
        void add_product(Z& sum,const Z& a,const Z& b) {
            if(mpz_sgn(b.get_mpz_t()))mpz_addmul(sum.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
        }
        using SparseRows=std::vector<std::vector<std::pair<std::size_t,const Z*>>>;
        SparseRows sparse_rows(const ZMatrix& a) {
            SparseRows rows(a.n);
            for(std::size_t i=0;i<a.n;++i)for(std::size_t j=0;j<a.n;++j)if(a(i,j)!=0)rows[i].emplace_back(j,&a(i,j));
            return rows;
        }
        // |[z^m] det(zE-F)| <= [z^m] prod_i(ceil||E_i|| z+ceil||F_i||): expand det(zE-F)
        // multilinearly in its rows and bound each mixed determinant by Hadamard's
        // inequality. Returned in descending powers, like the coefficients.
        std::vector<Z> identity_bound(const ZMatrix& e,const ZMatrix& f) {
            auto norm=[](const ZMatrix& a,std::size_t i) {
                Z square=0,root,remainder;
                for(std::size_t j=0;j<a.n;++j)mpz_addmul(square.get_mpz_t(),a(i,j).get_mpz_t(),a(i,j).get_mpz_t());
                mpz_sqrtrem(root.get_mpz_t(),remainder.get_mpz_t(),square.get_mpz_t());
                if(remainder!=0)++root;
                return root;
            };
            std::vector<Z> ascending {
                1
            };
            for(std::size_t i=0;i<e.n;++i) {
                const Z row_e=norm(e,i),row_f=norm(f,i);
                std::vector<Z> product(ascending.size()+1);
                for(std::size_t k=0;k<ascending.size();++k) {
                    mpz_addmul(product[k].get_mpz_t(),row_f.get_mpz_t(),ascending[k].get_mpz_t());
                    mpz_addmul(product[k+1].get_mpz_t(),row_e.get_mpz_t(),ascending[k].get_mpz_t());
                }
                ascending=std::move(product);
            }
            return {ascending.rbegin(),ascending.rend()};
        }
        bool determinant_identity(const ZMatrix& e,const ZMatrix& f,const std::vector<Z>& c,unsigned workers) {
            auto check_point=[&](std::size_t k) {
                Z x=-Z(static_cast<unsigned long>(k+1));
                ZMatrix a(e.n);
                for(std::size_t i=0;i<a.data.size();++i)a.data[i]=x*e.data[i]-f.data[i];
                return determinant(std::move(a))==evaluate(c,x);
            };
            if(workers==1) {
                for(std::size_t k=0;k<=e.n;++k)if(!check_point(k))return false;
                return true;
            }
            std::atomic<std::size_t> next{0};
            std::atomic<bool> valid{true};
            std::mutex exception_mutex;
            std::exception_ptr first_exception;
            {
                // jthread destruction joins even if a later thread creation throws.
                // Inputs are read-only; determinant matrices and GMP temporaries are
                // local to each job. Success consumes all n+1 distinct points.
                std::vector<std::jthread> threads;
                threads.reserve(workers);
                for(unsigned i=0;i<workers;++i)threads.emplace_back([&] {
                    try {
                        for(;;) {
                            const auto k=next.fetch_add(1,std::memory_order_relaxed);
                            if(k>e.n)break;
                            if(!check_point(k))valid.store(false,std::memory_order_relaxed);
                        }
                    } catch(...) {
                        std::lock_guard lock(exception_mutex);
                        if(!first_exception)first_exception=std::current_exception();
                    }
                });
            }
            if(first_exception)std::rethrow_exception(first_exception);
            return valid.load(std::memory_order_relaxed);
        }
#if GMP_NUMB_BITS==64&&defined(__SIZEOF_INT128__)
#define SECULAR_MODULAR_ORACLE 1
        using Word=std::uint64_t;
        __extension__ typedef unsigned __int128 Wide;
        // Montgomery arithmetic modulo an odd p < 2^62. Values are residues below p
        // in Montgomery form, except the arguments of to() and results of from().
        struct Montgomery {
            Word p,negative_inverse,r2;
            explicit Montgomery(Word modulus):p(modulus) {
                Word inverse=modulus;
                // Newton iteration doubles the correct low bits: 3, 6, ..., 96 >= 64.
                for(int i=0;i<5;++i)inverse*=2-modulus*inverse;
                negative_inverse=0-inverse;
                const Wide r=(Wide(1)<<64)%modulus;
                r2=static_cast<Word>(r*r%modulus);
            }
            Word reduce(Wide t)const {
                const Word m=static_cast<Word>(t)*negative_inverse;
                const Word u=static_cast<Word>((t+Wide(m)*p)>>64);
                return u>=p?u-p:u;
            }
            Word mul(Word a,Word b)const {
                return reduce(Wide(a)*b);
            }
            Word to(Word a)const {
                return mul(a,r2);
            }
            Word from(Word a)const {
                return reduce(a);
            }
            Word add(Word a,Word b)const {
                const Word s=a+b;
                return s>=p?s-p:s;
            }
            Word sub(Word a,Word b)const {
                return a>=b?a-b:a+p-b;
            }
            Word power(Word a,Word exponent)const {
                Word result=to(1);
                for(;exponent;exponent>>=1) {
                    if(exponent&1)result=mul(result,a);
                    a=mul(a,a);
                }
                return result;
            }
            Word inverse(Word a)const {
                return power(a,p-2);
            }
        };
        // Miller-Rabin with the first twelve prime bases (2..37) is deterministic below
        // 318665857834031151167461 (about 3.18e23), far above every modulus below 2^62.
        bool proven_prime(Word n) {
            static constexpr Word bases[]= {
                2,3,5,7,11,13,17,19,23,29,31,37
            };
            for(Word q:bases)if(n%q==0)return n==q;
            const Montgomery m(n);
            Word d=n-1;
            unsigned s=0;
            for(;!(d&1);d>>=1)++s;
            const Word one=m.to(1),minus_one=m.to(n-1);
            for(Word a:bases) {
                Word x=m.power(m.to(a),d);
                if(x==one||x==minus_one)continue;
                bool composite=true;
                for(unsigned r=1;r<s&&composite;++r) {
                    x=m.mul(x,x);
                    if(x==minus_one)composite=false;
                }
                if(composite)return false;
            }
            return true;
        }
        Word prime_at(std::size_t index) {
            static std::mutex table_mutex;
            static std::vector<Word> table;
            std::lock_guard lock(table_mutex);
            Word candidate=table.empty()?(Word(1)<<62)-1:table.back()-2;
            for(;table.size()<=index;candidate-=2) {
                // Every modulus exceeds 2^61; the singular-E test below relies on it.
                if(candidate<(Word(1)<<61))throw std::length_error("oracle prime sequence exhausted");
                if(proven_prime(candidate))table.push_back(candidate);
            }
            return table[index];
        }
        // All residues stay in 64-bit words: unsigned long is only 32 bits on Windows.
        Word residue(const Z& z,Word p) {
            const auto size=mpz_size(z.get_mpz_t());
            if(size==0)return 0;
            const Word r=mpn_mod_1(mpz_limbs_read(z.get_mpz_t()),static_cast<mp_size_t>(size),p);
            return mpz_sgn(z.get_mpz_t())<0&&r?p-r:r;
        }
        // A=E^-1 F over F_p in Montgomery form, and det(E) mod p. False if p divides det(E).
        bool reduce_pencil(const ZMatrix& e,const ZMatrix& f,bool diagonal_e,const Montgomery& m,
                           std::vector<Word>& a,Word& det_e) {
            const auto n=e.n;
            a.resize(n*n);
            for(std::size_t i=0;i<a.size();++i)a[i]=m.to(residue(f.data[i],m.p));
            det_e=m.to(1);
            if(diagonal_e) {
                for(std::size_t i=0;i<n;++i) {
                    const Word pivot=m.to(residue(e(i,i),m.p));
                    if(pivot==0)return false;
                    det_e=m.mul(det_e,pivot);
                    const Word inverse=m.inverse(pivot);
                    for(std::size_t j=0;j<n;++j)a[i*n+j]=m.mul(inverse,a[i*n+j]);
                }
                return true;
            }
            // Gauss-Jordan on [E | F].
            std::vector<Word> left(n*n);
            for(std::size_t i=0;i<left.size();++i)left[i]=m.to(residue(e.data[i],m.p));
            for(std::size_t k=0;k<n;++k) {
                std::size_t row=k;
                while(row<n&&left[row*n+k]==0)++row;
                if(row==n)return false;
                if(row!=k) {
                    for(std::size_t j=0;j<n;++j) {
                        std::swap(left[k*n+j],left[row*n+j]);
                        std::swap(a[k*n+j],a[row*n+j]);
                    }
                    det_e=m.sub(0,det_e);
                }
                const Word pivot=left[k*n+k];
                det_e=m.mul(det_e,pivot);
                const Word inverse=m.inverse(pivot);
                for(std::size_t j=0;j<n;++j) {
                    left[k*n+j]=m.mul(left[k*n+j],inverse);
                    a[k*n+j]=m.mul(a[k*n+j],inverse);
                }
                for(std::size_t i=0;i<n;++i)if(i!=k&&left[i*n+k]) {
                    const Word factor=left[i*n+k];
                    for(std::size_t j=0;j<n;++j) {
                        left[i*n+j]=m.sub(left[i*n+j],m.mul(factor,left[k*n+j]));
                        a[i*n+j]=m.sub(a[i*n+j],m.mul(factor,a[k*n+j]));
                    }
                }
            }
            return true;
        }
        // det(zI-A) over F_p in descending powers: similarity reduction to upper Hessenberg
        // form, then the Hessenberg determinant recurrence. Exact over any field.
        std::vector<Word> hessenberg_charpoly(std::vector<Word>& h,std::size_t n,const Montgomery& m) {
            auto at=[&](std::size_t i,std::size_t j)->Word& {
                return h[i*n+j];
            };
            for(std::size_t j=0;j+2<n;++j) {
                std::size_t pivot=j+1;
                while(pivot<n&&at(pivot,j)==0)++pivot;
                if(pivot==n)continue;
                if(pivot!=j+1) {
                    for(std::size_t c=0;c<n;++c)std::swap(at(pivot,c),at(j+1,c));
                    for(std::size_t r=0;r<n;++r)std::swap(at(r,pivot),at(r,j+1));
                }
                const Word inverse=m.inverse(at(j+1,j));
                for(std::size_t r=j+2;r<n;++r) {
                    if(at(r,j)==0)continue;
                    // Row r -= u*row j+1, then column j+1 += u*column r: a similarity.
                    const Word u=m.mul(at(r,j),inverse);
                    for(std::size_t c=j;c<n;++c)at(r,c)=m.sub(at(r,c),m.mul(u,at(j+1,c)));
                    for(std::size_t q=0;q<n;++q)at(q,j+1)=m.add(at(q,j+1),m.mul(u,at(q,r)));
                }
            }
            // 1-based: p_k=(z-h_kk)p_(k-1)-sum_(i<k) h_ik*(h_(i+1,i)...h_(k,k-1))*p_(i-1).
            // Each p_k is stored in ascending powers.
            std::vector<std::vector<Word>> p(n+1);
            p[0]= {
                m.to(1)
            };
            for(std::size_t k=1;k<=n;++k) {
                auto& current=p[k];
                current.assign(k+1,0);
                const auto& previous=p[k-1];
                const Word diagonal_entry=at(k-1,k-1);
                for(std::size_t d=0;d<previous.size();++d) {
                    current[d+1]=m.add(current[d+1],previous[d]);
                    current[d]=m.sub(current[d],m.mul(diagonal_entry,previous[d]));
                }
                Word product=m.to(1);
                for(std::size_t i=k-1;i>=1;--i) {
                    product=m.mul(product,at(i,i-1));
                    if(product==0)break;
                    const Word coefficient=m.mul(at(i-1,k-1),product);
                    const auto& lower=p[i-1];
                    for(std::size_t d=0;d<lower.size();++d)current[d]=m.sub(current[d],m.mul(coefficient,lower[d]));
                }
            }
            return {p[n].rbegin(),p[n].rend()};
        }
        enum class ModularStatus {
            verified,refuted,inconclusive
        };
        struct ModularOutcome {
            ModularStatus status=ModularStatus::inconclusive;
            std::size_t moduli=0;
            unsigned workers=0;
        };
        // Size-only work estimates in limb operations. They only choose between two exact
        // proofs: residues cost every modulus a pass over every entry and coefficient, which for
        // tiny orders with very long entries exceeds the n+1 Bareiss determinants.
        bool bareiss_clearly_cheaper(const ZMatrix& e,const ZMatrix& f,const std::vector<Z>& c,std::size_t need_bits) {
            const double n=static_cast<double>(e.n);
            double residue_limbs=0,entry_limbs=0,coefficient_limbs=0;
            for(std::size_t i=0;i<e.data.size();++i) {
                const double le=static_cast<double>(mpz_size(e.data[i].get_mpz_t()));
                const double lf=static_cast<double>(mpz_size(f.data[i].get_mpz_t()));
                residue_limbs+=le+lf;
                entry_limbs=std::max(entry_limbs,std::max(le,lf)+1);
            }
            for(const auto& v:c) {
                const double l=static_cast<double>(mpz_size(v.get_mpz_t()));
                residue_limbs+=l;
                coefficient_limbs=std::max(coefficient_limbs,l);
            }
            const double modular=std::ceil(static_cast<double>(need_bits)/61)*(residue_limbs+2*n*n*n);
            if(modular<1e7)return false;
            // Per point: form xE-F and evaluate q, then (n-k)^2 updates of two products and an
            // exact division of k-minor-sized entries, with multiplication modelled as m^1.6.
            double elimination=0;
            for(std::size_t k=1;k<e.n;++k)elimination+=(n-k)*(n-k)*4*std::pow(static_cast<double>(k)*entry_limbs,1.6);
            const double bareiss=(n+1)*(n*n*entry_limbs+(n+1)*coefficient_limbs+elimination);
            return 2*bareiss<modular;
        }
        // For p not dividing det(E), det(zE-F) = det(E) det(zI-E^-1 F) holds over F_p.
        // Agreement modulo distinct primes whose product exceeds max_m(|c_m|+B_m) >= |c_m-P_m|
        // proves c = P, the coefficients of det(zE-F); one disagreement refutes c.
        ModularOutcome modular_identity(const ZMatrix& e,const ZMatrix& f,const std::vector<Z>& c,unsigned requested) {
            const auto n=e.n;
            const auto bound=identity_bound(e,f);
            Z need=0;
            for(std::size_t k=0;k<=n;++k) {
                Z v=abs(c[k])+bound[k];
                if(v>need)need=std::move(v);
            }
            ModularOutcome out;
            // Every modulus exceeds 2^61, so k matching moduli have a product above 2^(61k),
            // and 61k >= bits(need) gives a product above need.
            const std::size_t need_bits=bits(need);
            if(bareiss_clearly_cheaper(e,f,c,need_bits))return out;
            const bool diagonal_e=diagonal(e);
            constexpr unsigned char unusable=0,match=1,mismatch=2;
            std::size_t proven_bits=0,next_prime=0,unusable_count=0;
            std::optional<Team> team;
            for(;;) {
                std::vector<Word> batch;
                for(std::size_t planned=proven_bits;planned<need_bits;planned+=61)batch.push_back(prime_at(next_prime++));
                if(!team) {
                    out.workers=oracle_worker_count(requested,batch.size());
                    team.emplace(out.workers);
                }
                std::vector<unsigned char> outcome(batch.size(),unusable);
                team->run(batch.size(),[&](std::size_t index,unsigned) {
                    const Montgomery m(batch[index]);
                    std::vector<Word> a;
                    Word det_e=0;
                    if(!reduce_pencil(e,f,diagonal_e,m,a,det_e))return;
                    const auto chi=hessenberg_charpoly(a,n,m);
                    for(std::size_t k=0;k<=n;++k)if(residue(c[k],m.p)!=m.from(m.mul(det_e,chi[k]))) {
                        outcome[index]=mismatch;
                        return;
                    }
                    outcome[index]=match;
                });
                for(std::size_t index=0;index<batch.size();++index) {
                    if(outcome[index]==unusable) {
                        ++unusable_count;
                        continue;
                    }
                    ++out.moduli;
                    if(outcome[index]==mismatch) {
                        out.status=ModularStatus::refuted;
                        return out;
                    }
                    proven_bits+=61;
                }
                if(proven_bits>=need_bits) {
                    out.status=ModularStatus::verified;
                    return out;
                }
                // Distinct primes above 2^61 dividing a nonzero det(E) <= bound[0] number fewer
                // than bits(bound[0])/61. Reaching that count proves det(E)=0: inconclusive.
                if(unusable_count*61>=bits(bound[0]))return out;
            }
        }
#endif
    }
    Z pow10(unsigned d) {
        if(d>100000)throw std::invalid_argument("decimal scale exceeds limit");
        Z v;
        mpz_ui_pow_ui(v.get_mpz_t(),10,d);
        return v;
    }
    std::size_t bits(const Z& z) {
        return z==0?0:mpz_sizeinbase(z.get_mpz_t(),2);
    }
    Q parse_rational(const std::string& raw) {
        std::string s=trim(raw);
        if(s.empty()||s.size()>20000)throw std::invalid_argument("invalid numeric token length");
        auto slash=s.find('/');
        if(slash!=std::string::npos) {
            Z n=integer(s.substr(0,slash)),d=integer(s.substr(slash+1));
            if(d==0)throw std::invalid_argument("zero denominator");
            Q q(n,d);
            q.canonicalize();
            return q;
        }
        long exponent=0;
        auto pos=s.find_first_of("eEdD");
        if(pos!=std::string::npos) {
            Z ex=integer(s.substr(pos+1));
            if(ex>10000||ex< -10000)throw std::invalid_argument("decimal exponent exceeds 10000");
            exponent=ex.get_si();
            s=s.substr(0,pos);
        }
        int sign=1;
        if(!s.empty()&&(s[0]=='+'||s[0]=='-')) {
            if(s[0]=='-')sign=-1;
            s.erase(0,1);
        }
        bool has_digit=false,has_point=false;
        for(unsigned char ch:s) {
            if(std::isdigit(ch))has_digit=true;
            else if(ch=='.'&&!has_point)has_point=true;
            else throw std::invalid_argument("invalid numeric mantissa");
        }
        if(!has_digit)throw std::invalid_argument("missing mantissa digits");
        auto point=s.find('.');
        long places=0;
        if(point!=std::string::npos) {
            places=static_cast<long>(s.size()-point-1);
            s.erase(point,1);
        }
        Z num=integer(s);
        num*=sign;
        long power=exponent-places;
        Q q;
        if(power>=0)q=Q(num*pow10(static_cast<unsigned>(power)));
        else q=Q(num,pow10(static_cast<unsigned>(-power)));
        q.canonicalize();
        return q;
    }
    QMatrix read_matrix_market(std::istream& in) {
        std::string line;
        if(!std::getline(in,line))throw std::invalid_argument("empty matrix file");
        std::istringstream h(line);
        std::string magic,object,format,field,storage,extra;
        if(!(h>>magic>>object>>format>>field>>storage)||h>>extra||magic!="%%MatrixMarket"||object!="matrix")throw std::invalid_argument("Matrix Market header required");
        if(format!="coordinate"&&format!="array")throw std::invalid_argument("only coordinate/array Matrix Market supported");
        if(field!="real"&&field!="integer")throw std::invalid_argument("real or integer numeric values required (not pattern/complex)");
        if(storage!="symmetric"&&storage!="general")throw std::invalid_argument("only symmetric/general storage supported");
        std::istringstream dims(data_line(in));
        long rows=0,cols=0,entries=0;
        if(!(dims>>rows>>cols)||rows<=0||rows!=cols||rows>static_cast<long>(max_dimension))throw std::invalid_argument("matrix must be square, order 1..256");
        const auto n=static_cast<std::size_t>(rows);
        QMatrix a(n);
        if(format=="coordinate") {
            if(!(dims>>entries)||entries<0||entries>static_cast<long>(4*n*n)||dims>>extra)throw std::invalid_argument("invalid entry count");
            int triangle=0;
            for(long k=0;k<entries;++k) {
                std::istringstream row(data_line(in));
                long r=0,c=0;
                std::string value;
                if(!(row>>r>>c>>value)||row>>extra||r<1||c<1||r>rows||c>cols)throw std::invalid_argument("invalid matrix coordinate");
                Q q=parse_rational(value);
                auto i=static_cast<std::size_t>(r-1),j=static_cast<std::size_t>(c-1);
                if(storage=="symmetric"&&i!=j) {
                    int tr=i>j?1:-1;
                    if(triangle&&tr!=triangle)throw std::invalid_argument("symmetric file mixes upper and lower triangles");
                    triangle=tr;
                }
                a(i,j)+=q;
                if(storage=="symmetric"&&i!=j)a(j,i)+=q;
            }
        }
        else {
            if(dims>>extra)throw std::invalid_argument("extra array dimensions");
            // Matrix Market array entries are column-major, lower triangle if symmetric.
            for(std::size_t c=0;c<n;++c)for(std::size_t r=(storage=="symmetric"?c:0);r<n;++r) {
                std::istringstream item(data_line(in));
                std::string value;
                if(!(item>>value)||item>>extra)throw std::invalid_argument("array requires one value per line");
                a(r,c)=parse_rational(value);
                if(storage=="symmetric")a(c,r)=a(r,c);
            }
        }
        while(std::getline(in,line)) {
            line=trim(line);
            if(!line.empty()&&line[0]!='%')throw std::invalid_argument("extra data after matrix entries");
        }
        return a;
    }
    QMatrix read_matrix_market(const std::filesystem::path& p) {
        std::ifstream in(p);
        if(!in)throw std::runtime_error("cannot open matrix: "+p.string());
        return read_matrix_market(in);
    }
    void write_matrix_market(std::ostream& out,const QMatrix& a) {
        bool sym=symmetric(a);
        std::size_t count=0;
        for(std::size_t i=0;i<a.n;++i)for(std::size_t j=0;j<(sym?i+1:a.n);++j)if(a(i,j)!=0)++count;
        out<<"%%MatrixMarket matrix coordinate real "<<(sym?"symmetric":"general")<<"\n% SecularLab rational-token extension; values are exact.\n"<<a.n<<" "<<a.n<<" "<<count<<"\n";
        for(std::size_t i=0;i<a.n;++i)for(std::size_t j=0;j<(sym?i+1:a.n);++j)if(a(i,j)!=0)out<<i+1<<" "<<j+1<<" "<<a(i,j)<<"\n";
        if(!out)throw std::runtime_error("matrix write failed");
    }
    Z exact_divide(const Z& a,const Z& b) {
        if(b==0||!mpz_divisible_p(a.get_mpz_t(),b.get_mpz_t()))throw std::logic_error("non-exact integer division");
        Z q;
        mpz_divexact(q.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
        return q;
    }
    ZMatrix multiply(const ZMatrix& a,const ZMatrix& b) {
        if(a.n!=b.n)throw std::invalid_argument("matrix dimensions differ");
        ZMatrix c(a.n);
        for(std::size_t i=0;i<a.n;++i)for(std::size_t k=0;k<a.n;++k)if(a(i,k)!=0)for(std::size_t j=0;j<a.n;++j)if(b(k,j)!=0)mpz_addmul(c(i,j).get_mpz_t(),a(i,k).get_mpz_t(),b(k,j).get_mpz_t());
        return c;
    }
    Z determinant(ZMatrix a) {
        const auto n=a.n;
        Z previous=1;
        int sign=1;
        for(std::size_t k=0;k+1<n;++k) {
            std::size_t row=k;
            while(row<n&&a(row,k)==0)++row;
            if(row==n)return 0;
            if(row!=k) {
                for(std::size_t j=0;j<n;++j)std::swap(a(k,j),a(row,j));
                sign=-sign;
            }
            Z pivot=a(k,k);
            for(std::size_t i=k+1;i<n;++i) {
                for(std::size_t j=k+1;j<n;++j) {
                    Z v=pivot*a(i,j)-a(i,k)*a(k,j);
                    a(i,j)=exact_divide(v,previous);
                }
                a(i,k)=0;
            }
            previous=std::move(pivot);
        }
        return sign*a(n-1,n-1);
    }
    ZMatrix adjugate(const ZMatrix& a) {
        Z d=determinant(a);
        if(d==0)throw std::invalid_argument("singular leading matrix");
        ZMatrix result(a.n);
        if(diagonal(a)) {
            for(std::size_t i=0;i<a.n;++i)result(i,i)=exact_divide(d,a(i,i));
            return result;
        }
        QMatrix left(a.n),right=identity<Q>(a.n);
        for(std::size_t i=0;i<a.data.size();++i)left.data[i]=a.data[i];
        for(std::size_t k=0;k<a.n;++k) {
            std::size_t row=k;
            while(row<a.n&&left(row,k)==0)++row;
            if(row==a.n)throw std::logic_error("adjugate pivot missing");
            for(std::size_t j=0;j<a.n;++j) {
                std::swap(left(k,j),left(row,j));
                std::swap(right(k,j),right(row,j));
            }
            Q pivot=left(k,k);
            for(std::size_t j=0;j<a.n;++j) {
                left(k,j)/=pivot;
                right(k,j)/=pivot;
            }
            for(std::size_t i=0;i<a.n;++i)if(i!=k) {
                Q factor=left(i,k);
                for(std::size_t j=0;j<a.n;++j) {
                    left(i,j)-=factor*left(k,j);
                    right(i,j)-=factor*right(k,j);
                }
            }
        }
        for(std::size_t i=0;i<a.data.size();++i) {
            Q v=right.data[i]*d;
            if(v.get_den()!=1)throw std::logic_error("noninteger adjugate");
            result.data[i]=v.get_num();
        }
        return result;
    }
    Inertia inertia(const ZMatrix& original) {
        if(!symmetric(original))throw std::invalid_argument("inertia requires symmetric matrix");
        ZMatrix a=original;
        Inertia out;
        Z previous=1;
        for(std::size_t k=0;k<a.n;++k) {
            std::size_t p=k;
            while(p<a.n&&a(p,p)==0)++p;
            if(p==a.n) {
                bool zero=true;
                for(std::size_t i=k;i<a.n;++i)for(std::size_t j=k;j<a.n;++j)if(a(i,j)!=0)zero=false;
                if(zero) {
                    out.zero+=a.n-k;
                    return out;
                }
                QMatrix q(a.n);
                for(std::size_t i=0;i<original.data.size();++i)q.data[i]=original.data[i];
                return rational_inertia(std::move(q));
            }
            swap_symmetric(a,k,p);
            Z pivot=a(k,k);
            int sign=mpz_sgn(pivot.get_mpz_t())*mpz_sgn(previous.get_mpz_t());
            if(sign>0)++out.positive;
            else ++out.negative;
            for(std::size_t i=k+1;i<a.n;++i)for(std::size_t j=i;j<a.n;++j) {
                Z v=pivot*a(i,j)-a(i,k)*a(k,j);
                v=exact_divide(v,previous);
                a(i,j)=v;
                a(j,i)=v;
            }
            previous=std::move(pivot);
        }
        return out;
    }
    Inertia inertia(const QMatrix& a) {
        if(!symmetric(a))throw std::invalid_argument("inertia requires symmetric matrix");
        Z l=1;
        for(auto& v:a.data)mpz_lcm(l.get_mpz_t(),l.get_mpz_t(),v.get_den_mpz_t());
        ZMatrix b(a.n);
        for(std::size_t i=0;i<a.data.size();++i)b.data[i]=a.data[i].get_num()*exact_divide(l,a.data[i].get_den());
        return inertia(b);
    }
    IntegerPencil integer_pencil(const QMatrix& k,const QMatrix& m) {
        if(k.n!=m.n)throw std::invalid_argument("K and M dimensions differ");
        Z l=1;
        for(auto* a: {
            &k,&m
        })for(const auto& v:a->data)mpz_lcm(l.get_mpz_t(),l.get_mpz_t(),v.get_den_mpz_t());
        IntegerPencil p {
            ZMatrix(k.n),ZMatrix(k.n),l
        };
        for(std::size_t i=0;i<k.data.size();++i) {
            p.E.data[i]=m.data[i].get_num()*exact_divide(l,m.data[i].get_den());
            p.F.data[i]=k.data[i].get_num()*exact_divide(l,k.data[i].get_den());
        }
        return p;
    }
    Inertia shifted_inertia(const ZMatrix& f,const ZMatrix& e,const Q& s) {
        if(f.n!=e.n)throw std::invalid_argument("pencil dimensions differ");
        ZMatrix a(f.n);
        for(std::size_t i=0;i<f.data.size();++i)a.data[i]=s.get_den()*f.data[i]-s.get_num()*e.data[i];
        return inertia(a);
    }
    CoefficientResult characteristic(const ZMatrix& e,const ZMatrix& f,unsigned requested) {
        if(e.n!=f.n)throw std::invalid_argument("pencil dimensions differ");
        const auto n=e.n;
        if(n==0||e.data.size()!=n*n||f.data.size()!=n*n)throw std::invalid_argument("invalid pencil storage");
        const unsigned workers=oracle_worker_count(requested,n);
        const bool diagonal_e=diagonal(e);
        // The determinant and adjugate of a diagonal E are products and exact quotients.
        Z d=1;
        if(diagonal_e)for(std::size_t i=0;i<n;++i)d*=e(i,i);
        else d=determinant(e);
        if(d==0)throw std::invalid_argument("singular mass/leading matrix");
        ZMatrix j(n);
        if(diagonal_e)for(std::size_t i=0;i<n;++i)j(i,i)=exact_divide(d,e(i,i));
        else j=adjugate(e);
        // S_k is the z^(n-1-k) coefficient of adj(zE-F), hence symmetric for a symmetric
        // pencil. Then only entries (i,c) with c>=i of S_k are computed and stored.
        const bool symmetric_pencil=symmetric(e)&&symmetric(f);
        const bool fused=symmetric_pencil&&diagonal_e;
        const auto rows=sparse_rows(f);
        ZMatrix s(n),next(n),t;
        for(std::size_t i=0;i<n;++i)for(std::size_t c=symmetric_pencil?i:0;c<n;++c)s(i,c)=j(i,c);
        if(!fused)t=ZMatrix(n);
        std::vector<std::vector<Z>> accumulators(workers,std::vector<Z>(n));
        std::vector<Z> products(workers),diagonal_t(fused?n:0);
        std::vector<std::size_t> row_bits(n),row_divisions(n);
        auto stored=[&](const ZMatrix& a,std::size_t row,std::size_t column)->const Z& {
            return symmetric_pencil&&column<row?a(column,row):a(row,column);
        };
        // Rows of a stage are independent jobs: every worker count gives identical integers.
        Team team(workers);
        CoefficientResult r;
        r.workers=workers;
        r.coefficients.push_back(d);
        for(std::size_t k=1;k<=n;++k) {
            if(fused)team.run(n,[&](std::size_t i,unsigned slot) {
                // Row i of T_k=F*S_(k-1) on and right of the diagonal. Off-diagonal entries
                // are divided straight into S_k. The mirrored entry T_k(c,i)=E_cc*S_k(i,c)
                // enters the workspace statistic, formed only when it could raise the maximum.
                auto& acc=accumulators[slot];
                for(std::size_t c=i;c<n;++c)mpz_set_ui(acc[c].get_mpz_t(),0);
                for(const auto& [l,value]:rows[i])for(std::size_t c=i;c<n;++c)add_product(acc[c],*value,stored(s,l,c));
                std::size_t width=0;
                for(std::size_t c=i+1;c<n;++c) {
                    width=std::max(width,bits(acc[c]));
                    exact_divide_into(next(i,c),acc[c],e(i,i));
                    const Z& entry=next(i,c);
                    if(entry!=0&&bits(e(c,c))+bits(entry)>width) {
                        mpz_mul(products[slot].get_mpz_t(),e(c,c).get_mpz_t(),entry.get_mpz_t());
                        width=std::max(width,bits(products[slot]));
                    }
                }
                mpz_swap(diagonal_t[i].get_mpz_t(),acc[i].get_mpz_t());
                row_bits[i]=width;
                row_divisions[i]=n-1-i;
            });
            else team.run(n,[&](std::size_t i,unsigned) {
                for(std::size_t c=0;c<n;++c)mpz_set_ui(t(i,c).get_mpz_t(),0);
                for(const auto& [l,value]:rows[i])for(std::size_t c=0;c<n;++c)add_product(t(i,c),*value,stored(s,l,c));
                std::size_t width=0;
                for(std::size_t c=0;c<n;++c)if(c!=i)width=std::max(width,bits(t(i,c)));
                row_bits[i]=width;
            });
            Z trace=0;
            for(std::size_t i=0;i<n;++i)trace+=fused?diagonal_t[i]:t(i,i);
            const Z c=exact_divide(trace,Z(static_cast<unsigned long>(k)));
            ++r.exact_divisions;
            r.coefficients.push_back(-c);
            for(std::size_t i=0;i<n;++i) {
                Z& entry=fused?diagonal_t[i]:t(i,i);
                entry-=c;
                r.max_workspace_bits=std::max({r.max_workspace_bits,row_bits[i],bits(entry)});
            }
            if(fused)for(std::size_t i=0;i<n;++i) {
                exact_divide_into(next(i,i),diagonal_t[i],e(i,i));
                ++r.exact_divisions;
            }
            else team.run(n,[&](std::size_t i,unsigned slot) {
                // S_k=E^-1(T_k-c_k I): row division for diagonal E, otherwise J(T_k-c_k I)/d.
                const std::size_t first=symmetric_pencil?i:0;
                if(diagonal_e)for(std::size_t c=first;c<n;++c)exact_divide_into(next(i,c),t(i,c),e(i,i));
                else {
                    auto& acc=accumulators[slot];
                    for(std::size_t c=first;c<n;++c)mpz_set_ui(acc[c].get_mpz_t(),0);
                    for(std::size_t l=0;l<n;++l)if(j(i,l)!=0)for(std::size_t c=first;c<n;++c)add_product(acc[c],j(i,l),t(l,c));
                    for(std::size_t c=first;c<n;++c)exact_divide_into(next(i,c),acc[c],d);
                }
                row_divisions[i]=n-first;
            });
            for(const auto count:row_divisions)r.exact_divisions+=count;
            std::swap(s,next);
        }
        r.terminal_zero=std::all_of(s.data.begin(),s.data.end(),[](const Z& z) {
            return z==0;
        });
        if(!r.terminal_zero)throw std::logic_error("coefficient recurrence did not terminate at S_n=0");
        for(auto& c:r.coefficients)r.max_coefficient_bits=std::max(r.max_coefficient_bits,bits(c));
        return r;
    }
    std::vector<std::vector<Z>> affine_characteristic(const ZMatrix& e,const ZMatrix& f0,const ZMatrix& f1) {
        if(e.n!=f0.n||e.n!=f1.n)throw std::invalid_argument("affine dimensions differ");
        if(e.n>24)throw std::invalid_argument("symbolic affine mode limited to order 24");
        Z d=determinant(e);
        if(d==0)throw std::invalid_argument("singular affine leading matrix");
        ZMatrix j=adjugate(e);
        std::vector<ZMatrix> old {
            j
        };
        std::vector<std::vector<Z>> co {
            {
                d
            }
        };
        for(std::size_t k=1;k<=e.n;++k) {
            std::vector<ZMatrix> next;
            std::vector<Z> row;
            for(std::size_t t=0;t<=k;++t) {
                ZMatrix v(e.n);
                if(t<old.size())v=multiply(f0,old[t]);
                if(t>0) {
                    auto w=multiply(f1,old[t-1]);
                    for(std::size_t u=0;u<v.data.size();++u)v.data[u]+=w.data[u];
                }
                Z tr=0;
                for(std::size_t i=0;i<e.n;++i)tr+=v(i,i);
                Z c=exact_divide(tr,Z(static_cast<unsigned long>(k)));
                row.push_back(-c);
                for(std::size_t i=0;i<e.n;++i)v(i,i)-=c;
                next.push_back(apply_e_inverse(e,j,d,v,nullptr));
            }
            co.push_back(std::move(row));
            old=std::move(next);
        }
        for(auto& a:old)for(auto& z:a.data)if(z!=0)throw std::logic_error("affine S_n is not zero");
        return co;
    }
    Z evaluate(const std::vector<Z>& c,const Z& x) {
        Z v=0;
        for(const auto& a:c)v=v*x+a;
        return v;
    }
    Q evaluate(const std::vector<Z>& c,const Q& x) {
        Q v=0;
        for(const auto& a:c)v=v*x+a;
        return v;
    }
    int polynomial_sign(const std::vector<Z>& c,const Q& x) {
        if(c.empty())throw std::invalid_argument("empty polynomial");
        // At stage k, v = denominator^k times the partial Horner value.
        // A canonical rational has positive denominator, so no sign is lost.
        const Z numerator=x.get_num(),denominator=x.get_den();
        Z value=c.front(),denominator_power=1;
        for(std::size_t i=1;i<c.size();++i) {
            denominator_power*=denominator;
            value*=numerator;
            mpz_addmul(value.get_mpz_t(),c[i].get_mpz_t(),denominator_power.get_mpz_t());
        }
        int sign=mpz_sgn(value.get_mpz_t());
        if(denominator<0&&(c.size()-1)%2)sign=-sign;
        return sign;
    }
    unsigned oracle_worker_count(unsigned requested,std::size_t jobs) {
        if(requested>64)throw std::invalid_argument("oracle workers must be 0..64");
        if(jobs==0)return 0;
        unsigned workers=requested;
        if(workers==0)workers=std::min(8u,std::max(1u,std::thread::hardware_concurrency()));
        return static_cast<unsigned>(std::min<std::size_t>(workers,jobs));
    }
    OracleReport check_coefficients(const ZMatrix& e,const ZMatrix& f,const std::vector<Z>& c,unsigned requested,CoefficientOracle oracle) {
        // Validate before scheduling; public Matrix fields may have been edited.
        const auto workers=oracle_worker_count(requested,e.n<max_dimension+1?e.n+1:0);
        OracleReport report;
        if(e.n==0||e.n>max_dimension||e.n!=f.n||c.size()!=e.n+1||
           e.data.size()!=e.n*e.n||f.data.size()!=e.n*e.n)return report;
#ifdef SECULAR_MODULAR_ORACLE
        if(oracle==CoefficientOracle::automatic) {
            const auto modular=modular_identity(e,f,c,requested);
            if(modular.status!=ModularStatus::inconclusive) {
                report.verified=modular.status==ModularStatus::verified;
                report.method="modular";
                report.moduli=modular.moduli;
                report.workers=modular.workers;
                return report;
            }
        }
#else
        (void)oracle;
#endif
        // Explicit request, det(E)=0 (every modulus divides det E), or a tiny order with very
        // long entries where the size estimate favours them: n+1 Bareiss points.
        report.verified=determinant_identity(e,f,c,workers);
        report.method="bareiss";
        report.workers=workers;
        return report;
    }
    bool verify_coefficients(const ZMatrix& e,const ZMatrix& f,const std::vector<Z>& c,unsigned requested) {
        return check_coefficients(e,f,c,requested).verified;
    }
    std::uint64_t oracle_prime(std::size_t index) {
#ifdef SECULAR_MODULAR_ORACLE
        return prime_at(index);
#else
        (void)index;
        throw std::logic_error("modular oracle requires 64-bit GMP limbs and 128-bit integers");
#endif
    }
    SquarefreeResult check_squarefree(const std::vector<Z>& c) {
        if(c.size()<2||c.front()==0)throw std::invalid_argument("invalid polynomial");
        std::vector<Z> a(c.rbegin(),c.rend()),b;
        primitive(a);
        // Fixed proven primes; a modular gcd of one is sufficient, never a
        // nontrivial modular gcd (which may merely reflect colliding roots).
        for(unsigned long prime:{2UL,3UL,5UL,7UL,11UL,65537UL})
            if(modular_squarefree(a,prime))return {false,prime};
        for(std::size_t i=1;i<a.size();++i)b.push_back(a[i]*static_cast<unsigned long>(i));
        primitive(b);
        while(!b.empty()) {
            auto r=a;
            while(r.size()>=b.size()) {
                const auto offset=r.size()-b.size();
                Z lc=r.back(),blc=b.back();
                for(auto& x:r)x*=blc;
                for(std::size_t j=0;j<b.size();++j)r[offset+j]-=lc*b[j];
                primitive(r);
            }
            a=std::move(b);
            b=std::move(r);
        }
        return {a.size()>1,0};
    }
    bool has_repeated_roots(const std::vector<Z>& c) {
        return check_squarefree(c).repeated;
    }
    Z round_scaled(const Q& x,unsigned d) {
        Q q=x*pow10(d);
        Z num=q.get_num(),den=q.get_den();
        bool neg=num<0;
        if(neg)num=-num;
        Z a=2*num+den,b=2*den,v;
        mpz_fdiv_q(v.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
        return neg?-v:v;
    }
    std::string format_decimal(const Z& z,unsigned d) {
        std::string v=(z<0?-z:z).get_str();
        if(d) {
            if(v.size()<=d)v=std::string(d+1-v.size(),'0')+v;
            v.insert(v.size()-d,1,'.');
        }
        return z<0?"-"+v:v;
    }
    Q legendre_value(unsigned degree,const Q& x) {
        if(degree>10000)throw std::invalid_argument("Legendre degree limit 10000");
        Q p0=1,p1=x;
        if(degree==0)return p0;
        for(unsigned i=1;i<degree;++i) {
            Q q=((2*i+1)*x*p1-i*p0)/(i+1);
            p0=p1;
            p1=q;
        }
        return p1;
    }
}
