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
#include <system_error>
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
        void poll(const CancelCheck& cancelled) {
            if(cancelled&&cancelled())throw Cancelled {};
        }
        Inertia rational_inertia(QMatrix a,const CancelCheck& cancelled) {
            Inertia out;
            const auto n=a.n;
            for(std::size_t k=0;k<n;) {
                poll(cancelled);
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
        // One truncating division and a zero-remainder test cost less than a divisibility
        // test followed by divexact. q and the scratch remainder r must be distinct objects.
        void exact_divide_into(Z& q,const Z& a,const Z& b,Z& r) {
            if(b==0)throw std::logic_error("non-exact integer division");
            mpz_tdiv_qr(q.get_mpz_t(),r.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
            if(r!=0)throw std::logic_error("non-exact integer division");
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
        // The work of one modulus, shared by the oracle and Chinese remaindering: det(zE-F) modulo
        // p as det(E) det(zI-E^-1 F), in descending powers, as plain residues below p. False when
        // p divides det(E).
        bool modular_characteristic(const ZMatrix& e,const ZMatrix& f,bool diagonal_e,Word p,std::vector<Word>& q) {
            const Montgomery m(p);
            std::vector<Word> a;
            Word det_e=0;
            if(!reduce_pencil(e,f,diagonal_e,m,a,det_e))return false;
            const auto chi=hessenberg_charpoly(a,e.n,m);
            q.resize(chi.size());
            for(std::size_t k=0;k<chi.size();++k)q[k]=m.from(m.mul(det_e,chi[k]));
            return true;
        }
        // The batch that would bring `proven_bits` to `need_bits` if every modulus in it is usable.
        std::vector<Word> prime_batch(std::size_t proven_bits,std::size_t need_bits,std::size_t& next_prime) {
            std::vector<Word> batch;
            for(std::size_t planned=proven_bits;planned<need_bits;planned+=61)batch.push_back(prime_at(next_prime++));
            return batch;
        }
        Z word(Word w) {
            Z z;
            mpz_import(z.get_mpz_t(),1,-1,sizeof w,0,0,&w);
            return z;
        }
        // Chinese remaindering for fixed distinct odd primes p_i with product M. A residue vector
        // r_i gives x = sum_i ((r_i w_i) mod p_i) M/p_i, with w_i = (M/p_i)^-1 mod p_i, so x = r_i
        // modulo every p_i and 0 <= x < count*M. The sum is formed up a product tree, where a node
        // with children (x_L,M_L) and (x_R,M_R) holds x_L M_R + x_R M_L, and then reduced to the
        // symmetric range |x| <= (M-1)/2.
        class Remainderer {
        public:
            Remainderer(const std::vector<Word>& moduli,Team& team):primes(moduli),weights(moduli.size()) {
                arithmetic.reserve(primes.size());
                tree.emplace_back();
                for(const auto p:primes) {
                    arithmetic.emplace_back(p);
                    tree.back().push_back(word(p));
                }
                while(tree.back().size()>1) {
                    const auto& below=tree.back();
                    std::vector<Z> level((below.size()+1)/2);
                    for(std::size_t j=0;j<level.size();++j)
                        level[j]=2*j+1<below.size()?Z(below[2*j]*below[2*j+1]):below[2*j];
                    tree.push_back(std::move(level));
                }
                half=tree.back().front()/2;
                team.run(primes.size(),[&](std::size_t i,unsigned) {
                    // Every prime is below 2^62 < 2*p_i, so p_j mod p_i is one conditional subtraction.
                    const Montgomery& m=arithmetic[i];
                    Word others=m.to(1);
                    for(std::size_t j=0;j<primes.size();++j)if(j!=i)
                        others=m.mul(others,m.to(primes[j]>=primes[i]?primes[j]-primes[i]:primes[j]));
                    weights[i]=m.inverse(others);
                });
            }
            // residues[i*stride] is the residue modulo p_i.
            Z combine(const Word* residues,std::size_t stride)const {
                std::vector<Z> x(primes.size());
                for(std::size_t i=0;i<primes.size();++i) {
                    // A plain residue times a Montgomery-form weight gives their plain product.
                    x[i]=word(arithmetic[i].mul(residues[i*stride],weights[i]));
                }
                Z sum;
                for(std::size_t h=0;h+1<tree.size();++h) {
                    const auto& moduli=tree[h];
                    const auto size=x.size();
                    // Node j reads only children 2j and 2j+1 >= j, so the level is formed in place.
                    for(std::size_t j=0;2*j+1<size;++j) {
                        mpz_mul(sum.get_mpz_t(),x[2*j].get_mpz_t(),moduli[2*j+1].get_mpz_t());
                        mpz_addmul(sum.get_mpz_t(),x[2*j+1].get_mpz_t(),moduli[2*j].get_mpz_t());
                        mpz_swap(x[j].get_mpz_t(),sum.get_mpz_t());
                    }
                    if(size%2)mpz_swap(x[size/2].get_mpz_t(),x[size-1].get_mpz_t());
                    x.resize((size+1)/2);
                }
                Z& value=x.front();
                mpz_mod(value.get_mpz_t(),value.get_mpz_t(),tree.back().front().get_mpz_t());
                // M is odd, so half=(M-1)/2 and the representatives are -(M-1)/2..(M-1)/2.
                if(value>half)value-=tree.back().front();
                return std::move(value);
            }
        private:
            std::vector<Word> primes;
            std::vector<Montgomery> arithmetic;
            // w_i in Montgomery form modulo p_i.
            std::vector<Word> weights;
            // tree[0] holds the primes, tree[h+1][j] the product of tree[h][2j] and tree[h][2j+1]
            // (a lone last node is carried up); the root is M.
            std::vector<std::vector<Z>> tree;
            Z half;
        };
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
                const auto batch=prime_batch(proven_bits,need_bits,next_prime);
                if(!team) {
                    out.workers=oracle_worker_count(requested,batch.size());
                    team.emplace(out.workers);
                }
                std::vector<unsigned char> outcome(batch.size(),unusable);
                team->run(batch.size(),[&](std::size_t index,unsigned) {
                    std::vector<Word> q;
                    if(!modular_characteristic(e,f,diagonal_e,batch[index],q))return;
                    for(std::size_t k=0;k<=n;++k)if(residue(c[k],batch[index])!=q[k]) {
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
        Z q,r;
        exact_divide_into(q,a,b,r);
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
        Z previous=1,update,remainder;
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
                    mpz_mul(update.get_mpz_t(),pivot.get_mpz_t(),a(i,j).get_mpz_t());
                    mpz_submul(update.get_mpz_t(),a(i,k).get_mpz_t(),a(k,j).get_mpz_t());
                    exact_divide_into(a(i,j),update,previous,remainder);
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
    namespace {
        // The fraction-free inertia, polling `cancelled` once per pivot.
        Inertia eliminate(const ZMatrix& original,const CancelCheck& cancelled) {
            if(!symmetric(original))throw std::invalid_argument("inertia requires symmetric matrix");
            Inertia out;
            if(diagonal(original)) {
                // The eigenvalues of a diagonal matrix are its diagonal entries.
                for(std::size_t i=0;i<original.n;++i) {
                    const int sign=mpz_sgn(original(i,i).get_mpz_t());
                    if(sign>0)++out.positive;
                    else if(sign<0)++out.negative;
                    else ++out.zero;
                }
                return out;
            }
            ZMatrix a=original;
            Z previous=1,update,remainder;
            for(std::size_t k=0;k<a.n;++k) {
                poll(cancelled);
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
                    return rational_inertia(std::move(q),cancelled);
                }
                swap_symmetric(a,k,p);
                Z pivot=a(k,k);
                int sign=mpz_sgn(pivot.get_mpz_t())*mpz_sgn(previous.get_mpz_t());
                if(sign>0)++out.positive;
                else ++out.negative;
                for(std::size_t i=k+1;i<a.n;++i)for(std::size_t j=i;j<a.n;++j) {
                    mpz_mul(update.get_mpz_t(),pivot.get_mpz_t(),a(i,j).get_mpz_t());
                    mpz_submul(update.get_mpz_t(),a(i,k).get_mpz_t(),a(k,j).get_mpz_t());
                    exact_divide_into(a(i,j),update,previous,remainder);
                    a(j,i)=a(i,j);
                }
                previous=std::move(pivot);
            }
            return out;
        }
    }
    Inertia inertia(const ZMatrix& a) {
        return eliminate(a,{});
    }
    Inertia inertia(const QMatrix& a) {
        if(!symmetric(a))throw std::invalid_argument("inertia requires symmetric matrix");
        Z l=1;
        for(auto& v:a.data)mpz_lcm(l.get_mpz_t(),l.get_mpz_t(),v.get_den_mpz_t());
        ZMatrix b(a.n);
        for(std::size_t i=0;i<a.data.size();++i)b.data[i]=a.data[i].get_num()*exact_divide(l,a.data[i].get_den());
        return inertia(b);
    }
    namespace {
        std::optional<Inertia> dominant_congruence(const ZMatrix& a,const ZMatrix& x,const CancelCheck& cancelled) {
            if(a.n!=x.n)throw std::invalid_argument("congruence dimensions differ");
            if(!symmetric(a))throw std::invalid_argument("inertia requires symmetric matrix");
            const auto n=a.n;
            // Upper triangular with a nonzero diagonal, so det(x)!=0 and x^T*a*x is a congruence.
            for(std::size_t i=0;i<n;++i) {
                if(x(i,i)==0)return std::nullopt;
                for(std::size_t j=0;j<i;++j)if(x(i,j)!=0)return std::nullopt;
            }
            // B = x^T*c with c = a*x. Both are needed only on and above the diagonal: b(i,j), i<=j,
            // reads c(m,j) for m<=i, and c(m,j) reads x(l,j) for l<=j.
            ZMatrix c(n);
            for(std::size_t m=0;m<n;++m) {
                poll(cancelled);
                for(std::size_t l=0;l<n;++l) {
                    const auto aml=a(m,l).get_mpz_t();
                    if(mpz_sgn(aml)==0)continue;
                    for(std::size_t j=std::max(l,m);j<n;++j) {
                        const auto xlj=x(l,j).get_mpz_t();
                        if(mpz_sgn(xlj)!=0)mpz_addmul(c(m,j).get_mpz_t(),aml,xlj);
                    }
                }
            }
            // B is formed one row at a time and not kept. Row i of the upper triangle completes the
            // off-diagonal sum of row i, because b(k,i), k<i, came with the earlier rows, and adds
            // its entries to the sums of the later rows. A failing row ends the check at once.
            // A strictly diagonally dominant symmetric B is nonsingular, and so is every
            // diag(B)+t(B-diag(B)), 0<=t<=1. No eigenvalue crosses zero along that path, so B,
            // and by Sylvester's law a, has the inertia of diag(B), with nullity zero.
            Inertia out;
            std::vector<Z> row(n),off(n);
            for(std::size_t i=0;i<n;++i) {
                poll(cancelled);
                for(std::size_t j=i;j<n;++j)row[j]=0;
                for(std::size_t m=0;m<=i;++m) {
                    const auto xmi=x(m,i).get_mpz_t();
                    if(mpz_sgn(xmi)==0)continue;
                    for(std::size_t j=i;j<n;++j) {
                        const auto cmj=c(m,j).get_mpz_t();
                        if(mpz_sgn(cmj)!=0)mpz_addmul(row[j].get_mpz_t(),xmi,cmj);
                    }
                }
                for(std::size_t j=i+1;j<n;++j) {
                    const auto v=row[j].get_mpz_t();
                    for(Z* sum: {&off[i],&off[j]}) {
                        if(mpz_sgn(v)<0)mpz_sub(sum->get_mpz_t(),sum->get_mpz_t(),v);
                        else mpz_add(sum->get_mpz_t(),sum->get_mpz_t(),v);
                    }
                }
                if(mpz_cmpabs(row[i].get_mpz_t(),off[i].get_mpz_t())<=0)return std::nullopt;
                if(mpz_sgn(row[i].get_mpz_t())>0)++out.positive;
                else ++out.negative;
            }
            return out;
        }
    }
    std::optional<Inertia> certify_congruence(const ZMatrix& a,const ZMatrix& x) {
        return dominant_congruence(a,x,{});
    }
    namespace {
        std::optional<Inertia> proposed_congruence(const ZMatrix& a,unsigned precision_bits,const CancelCheck& cancelled) {
            if(!symmetric(a))throw std::invalid_argument("inertia requires symmetric matrix");
            if(precision_bits<16||precision_bits>max_congruence_precision)
                throw std::invalid_argument("congruence precision must be 16.."+std::to_string(max_congruence_precision)+" bits");
            const auto n=a.n;
            // Untrusted proposal: an LDL^T of a without pivoting, in floating point whose precision
            // is set on every object; the process-global default precision is never read or changed.
            // The upper triangle holds the active matrix, the strict lower triangle L's multipliers.
            const mpf_class zero(0,precision_bits),one(1,precision_bits);
            std::vector<mpf_class> f(n*n,zero);
            mpf_class product=zero,sum=zero;
            auto at=[n](std::vector<mpf_class>& v,std::size_t i,std::size_t j) {
                return v[i*n+j].get_mpf_t();
            };
            std::size_t entry_bits=0;
            for(std::size_t i=0;i<n;++i)for(std::size_t j=i;j<n;++j) {
                mpf_set_z(at(f,i,j),a(i,j).get_mpz_t());
                entry_bits=std::max(entry_bits,bits(a(i,j)));
            }
            for(std::size_t k=0;k<n;++k) {
                poll(cancelled);
                const auto pivot=at(f,k,k);
                if(mpf_sgn(pivot)==0)return std::nullopt;
                for(std::size_t i=k+1;i<n;++i) {
                    if(mpf_sgn(at(f,k,i))==0)continue;
                    const auto l=at(f,i,k);
                    mpf_div(l,at(f,k,i),pivot);
                    for(std::size_t j=i;j<n;++j)if(mpf_sgn(at(f,k,j))!=0) {
                        mpf_mul(product.get_mpf_t(),l,at(f,k,j));
                        mpf_sub(at(f,i,j),at(f,i,j),product.get_mpf_t());
                    }
                }
            }
            // W = L^-1, unit lower triangular, one column at a time. Only L and the pivots are needed
            // from now on, so W(i,j), i>j, is stored at f(j,i) in the upper triangle; that is W^T, the
            // layout of x. The unit diagonal of W is implicit.
            for(std::size_t j=0;j<n;++j) {
                poll(cancelled);
                for(std::size_t i=j+1;i<n;++i) {
                    mpf_set_ui(sum.get_mpf_t(),0);
                    for(std::size_t m=j;m<i;++m) {
                        const auto wmj=m==j?one.get_mpf_t():at(f,j,m);
                        if(mpf_sgn(at(f,i,m))!=0&&mpf_sgn(wmj)!=0) {
                            mpf_mul(product.get_mpf_t(),at(f,i,m),wmj);
                            mpf_sub(sum.get_mpf_t(),sum.get_mpf_t(),product.get_mpf_t());
                        }
                    }
                    mpf_set(at(f,j,i),sum.get_mpf_t());
                }
            }
            // x = 2^s*W^T, truncated to integers: 2^s on the diagonal and zero below it by
            // construction. Any integers would do for soundness; only the dominance test decides.
            // Entries far beyond the input's scale mean a failed factorization, and are refused
            // before they are converted, which bounds the exact work.
            const mp_bitcnt_t scale=precision_bits-8;
            const long limit=static_cast<long>(precision_bits+entry_bits);
            ZMatrix x(n);
            for(std::size_t r=0;r<n;++r) {
                mpz_setbit(x(r,r).get_mpz_t(),scale);
                for(std::size_t c=r+1;c<n;++c) {
                    const auto v=at(f,r,c);
                    if(mpf_sgn(v)==0)continue;
                    long exponent=0;
                    mpf_get_d_2exp(&exponent,v);
                    if(exponent>limit)return std::nullopt;
                    mpf_mul_2exp(product.get_mpf_t(),v,scale);
                    mpz_set_f(x(r,c).get_mpz_t(),product.get_mpf_t());
                }
            }
            return dominant_congruence(a,x,cancelled);
        }
    }
    std::optional<Inertia> congruence_inertia(const ZMatrix& a,unsigned precision_bits) {
        return proposed_congruence(a,precision_bits,{});
    }
    InertiaReport scheduled_inertia(const ZMatrix& a,const std::vector<unsigned>& precisions,const CancelCheck& cancelled) {
        if(!symmetric(a))throw std::invalid_argument("inertia requires symmetric matrix");
        InertiaReport report;
        if(diagonal(a)) {
            report.method="diagonal";
            report.inertia=inertia(a);
            return report;
        }
        for(unsigned precision:precisions)if(const auto in=proposed_congruence(a,precision,cancelled)) {
            report.method="congruence";
            report.precision_bits=precision;
            report.inertia=*in;
            return report;
        }
        report.method="elimination";
        report.inertia=eliminate(a,cancelled);
        return report;
    }
    DefinitenessReport check_positive_definite(const ZMatrix& a) {
        const auto decided=scheduled_inertia(a,{64u,128u,256u});
        DefinitenessReport report;
        report.positive_definite=decided.inertia.positive==a.n;
        report.method=decided.method;
        report.precision_bits=decided.precision_bits;
        return report;
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
    std::vector<unsigned> shifted_precisions(const Q& s) {
        std::vector<unsigned> out;
        // A canonical denominator is positive, so bits() is at least one.
        const std::size_t first=(64+bits(s.get_den())+63)/64*64;
        for(std::size_t precision=first;precision<=max_congruence_precision&&out.size()<3;precision*=2)
            out.push_back(static_cast<unsigned>(precision));
        return out;
    }
    InertiaReport shifted_inertia_report(const ZMatrix& f,const ZMatrix& e,const Q& s,const CancelCheck& cancelled) {
        if(f.n!=e.n||f.data.size()!=e.data.size())throw std::invalid_argument("pencil dimensions differ");
        ZMatrix a(f.n);
        for(std::size_t i=0;i<f.data.size();++i)a.data[i]=s.get_den()*f.data[i]-s.get_num()*e.data[i];
        return scheduled_inertia(a,shifted_precisions(s),cancelled);
    }
    Inertia shifted_inertia(const ZMatrix& f,const ZMatrix& e,const Q& s,const CancelCheck& cancelled) {
        return shifted_inertia_report(f,e,s,cancelled).inertia;
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
        std::vector<Z> products(workers),remainders(workers),diagonal_t(fused?n:0);
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
                    exact_divide_into(next(i,c),acc[c],e(i,i),remainders[slot]);
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
                exact_divide_into(next(i,i),diagonal_t[i],e(i,i),remainders[0]);
                ++r.exact_divisions;
            }
            else team.run(n,[&](std::size_t i,unsigned slot) {
                // S_k=E^-1(T_k-c_k I): row division for diagonal E, otherwise J(T_k-c_k I)/d.
                const std::size_t first=symmetric_pencil?i:0;
                if(diagonal_e)for(std::size_t c=first;c<n;++c)exact_divide_into(next(i,c),t(i,c),e(i,i),remainders[slot]);
                else {
                    auto& acc=accumulators[slot];
                    for(std::size_t c=first;c<n;++c)mpz_set_ui(acc[c].get_mpz_t(),0);
                    for(std::size_t l=0;l<n;++l)if(j(i,l)!=0)for(std::size_t c=first;c<n;++c)add_product(acc[c],j(i,l),t(l,c));
                    for(std::size_t c=first;c<n;++c)exact_divide_into(next(i,c),acc[c],d,remainders[slot]);
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
    namespace {
        // Exponents a_i >= 0 with a_i+a_j <= v_q(A_ij) for every nonzero entry of E and F. The
        // start a_i = floor(min_j v_q / 2) over row and column i is feasible; each coordinate is
        // then raised to its largest feasible value, in a fixed order. Every feasible choice
        // gives an exact reduction; the ascent only aims at a large sum.
        std::vector<unsigned long> congruence_exponents(const ZMatrix& e,const ZMatrix& f,unsigned long q) {
            const auto n=e.n;
            constexpr auto unbounded=std::numeric_limits<unsigned long>::max();
            std::vector<std::vector<std::pair<std::size_t,unsigned long>>> limits(n);
            std::vector<unsigned long> cap(n,unbounded);
            const Z prime=q;
            Z rest;
            for(const ZMatrix* a: {
                &e,&f
            })for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
                const Z& entry=(*a)(i,j);
                if(entry==0)continue;
                const unsigned long v=q==2?mpz_scan1(entry.get_mpz_t(),0):
                    mpz_remove(rest.get_mpz_t(),entry.get_mpz_t(),prime.get_mpz_t());
                limits[i].emplace_back(j,v);
                if(i!=j)limits[j].emplace_back(i,v);
                cap[i]=std::min(cap[i],v);
                cap[j]=std::min(cap[j],v);
            }
            std::vector<unsigned long> exponent(n);
            for(std::size_t i=0;i<n;++i)exponent[i]=cap[i]==unbounded?0:cap[i]/2;
            for(bool changed=true;changed;) {
                changed=false;
                for(std::size_t i=0;i<n;++i) {
                    unsigned long best=unbounded;
                    // Feasibility gives v >= exponent[i]+exponent[j], so v-exponent[j] cannot wrap.
                    for(const auto& [j,v]:limits[i])best=std::min(best,j==i?v/2:v-exponent[j]);
                    if(best!=unbounded&&best>exponent[i]) {
                        exponent[i]=best;
                        changed=true;
                    }
                }
            }
            return exponent;
        }
        Z content(const ZMatrix& a) {
            Z g=0;
            for(const auto& v:a.data)mpz_gcd(g.get_mpz_t(),g.get_mpz_t(),v.get_mpz_t());
            return g==0?Z(1):g;
        }
    }
    PencilReduction reduce_pencil(const ZMatrix& e,const ZMatrix& f) {
        if(e.n!=f.n)throw std::invalid_argument("pencil dimensions differ");
        const auto n=e.n;
        if(n==0||e.data.size()!=n*n||f.data.size()!=n*n)throw std::invalid_argument("invalid pencil storage");
        PencilReduction r {
            e,f,std::vector<Z>(n,1),1,1,true
        };
        for(unsigned long q: {
            2UL,5UL
        }) {
            const auto exponent=congruence_exponents(e,f,q);
            for(std::size_t i=0;i<n;++i)if(exponent[i]) {
                Z power;
                mpz_ui_pow_ui(power.get_mpz_t(),q,exponent[i]);
                r.scale[i]*=power;
                r.trivial=false;
            }
        }
        // E' = P^-1 E P^-1 and F' = P^-1 F P^-1 by checked exact divisions, then their contents.
        if(!r.trivial)for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
            const Z d=r.scale[i]*r.scale[j];
            if(d==1)continue;
            r.E(i,j)=exact_divide(e(i,j),d);
            r.F(i,j)=exact_divide(f(i,j),d);
        }
        r.e_content=content(r.E);
        r.f_content=content(r.F);
        for(auto [a,g]: {
            std::pair<ZMatrix*,const Z*>{&r.E,&r.e_content},std::pair<ZMatrix*,const Z*>{&r.F,&r.f_content}
        })if(*g!=1) {
            for(auto& v:a->data)v=exact_divide(v,*g);
            r.trivial=false;
        }
        return r;
    }
    CoefficientResult reduced_characteristic(const ZMatrix& e,const ZMatrix& f,unsigned workers) {
        const auto reduction=reduce_pencil(e,f);
        if(reduction.trivial)return characteristic(e,f,workers);
        auto r=characteristic(reduction.E,reduction.F,workers);
        // c_k = c'_k alpha^(n-k) beta^k prod(scale)^2, by exact integer multiplication only.
        const auto n=e.n;
        Z square=1;
        for(const auto& s:reduction.scale)square*=s;
        square*=square;
        std::vector<Z> alpha(n+1,Z(1)),beta(n+1,Z(1));
        for(std::size_t k=1;k<=n;++k) {
            alpha[k]=alpha[k-1]*reduction.e_content;
            beta[k]=beta[k-1]*reduction.f_content;
        }
        r.max_coefficient_bits=0;
        for(std::size_t k=0;k<=n;++k) {
            r.coefficients[k]*=alpha[n-k]*beta[k]*square;
            r.max_coefficient_bits=std::max(r.max_coefficient_bits,bits(r.coefficients[k]));
        }
        return r;
    }
    CoefficientResult crt_characteristic(const ZMatrix& e,const ZMatrix& f,unsigned requested) {
        if(e.n!=f.n)throw std::invalid_argument("pencil dimensions differ");
        const auto n=e.n;
        if(n==0||e.data.size()!=n*n||f.data.size()!=n*n)throw std::invalid_argument("invalid pencil storage");
        oracle_worker_count(requested,1);
#ifdef SECULAR_MODULAR_ORACLE
        // |q_m| <= B_m < M/2 once M > 2*max_m B_m; 61k >= bits(max B)+1 ensures it for k moduli.
        const auto bound=identity_bound(e,f);
        const std::size_t need_bits=bits(*std::max_element(bound.begin(),bound.end()))+1;
        const bool diagonal_e=diagonal(e);
        CoefficientResult r;
        std::vector<Word> primes,residues;
        std::size_t proven_bits=0,next_prime=0;
        std::optional<Team> team;
        while(proven_bits<need_bits) {
            const auto batch=prime_batch(proven_bits,need_bits,next_prime);
            if(!team) {
                r.workers=oracle_worker_count(requested,batch.size());
                team.emplace(r.workers);
            }
            std::vector<std::vector<Word>> found(batch.size());
            std::vector<unsigned char> usable(batch.size(),0);
            team->run(batch.size(),[&](std::size_t index,unsigned) {
                usable[index]=modular_characteristic(e,f,diagonal_e,batch[index],found[index]);
            });
            for(std::size_t index=0;index<batch.size();++index) {
                if(!usable[index]) {
                    ++r.skipped_moduli;
                    continue;
                }
                primes.push_back(batch[index]);
                residues.insert(residues.end(),found[index].begin(),found[index].end());
                proven_bits+=61;
            }
            // As in the oracle: that many skipped moduli above 2^61 prove det(E)=0.
            if(proven_bits<need_bits&&r.skipped_moduli*61>=bits(bound[0]))
                throw std::invalid_argument("singular mass/leading matrix");
        }
        r.moduli=primes.size();
        const Remainderer crt(primes,*team);
        r.coefficients.resize(n+1);
        team->run(n+1,[&](std::size_t k,unsigned) {
            r.coefficients[k]=crt.combine(residues.data()+k,n+1);
        });
        // q_0=det(E) is nonzero modulo every usable prime, so it cannot be zero; checked anyway.
        if(r.coefficients.front()==0)throw std::logic_error("Chinese remaindering gave det(E)=0");
        for(const auto& c:r.coefficients)r.max_coefficient_bits=std::max(r.max_coefficient_bits,bits(c));
        return r;
#else
        throw std::logic_error("Chinese remaindering requires the modular oracle: 64-bit GMP limbs and 128-bit integers");
#endif
    }
    Z polynomial_content(const std::vector<Z>& coefficients) {
        Z g=0;
        for(const auto& c:coefficients)mpz_gcd(g.get_mpz_t(),g.get_mpz_t(),c.get_mpz_t());
        return g==0?Z(1):g;
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
    ScaledPolynomial::ScaledPolynomial(const std::vector<Z>& c,const Z& denominator):q(denominator) {
        if(c.empty())throw std::invalid_argument("empty polynomial");
        if(q<=0)throw std::invalid_argument("scaled polynomial denominator must be positive");
        table.resize(c.size());
        Z power=1;
        for(std::size_t i=0;i<c.size();++i) {
            mpz_mul(table[i].get_mpz_t(),c[i].get_mpz_t(),power.get_mpz_t());
            if(i+1<c.size())power*=q;
        }
    }
    bool ScaledPolynomial::covers(const Q& x)const {
        return mpz_sgn(x.get_den_mpz_t())>0&&mpz_divisible_p(q.get_mpz_t(),x.get_den_mpz_t());
    }
    int ScaledPolynomial::sign(const Q& x)const {
        if(!covers(x))throw std::invalid_argument("denominator does not divide the scaled polynomial denominator");
        // p=x*q is an integer, and v ends as sum_i c_i q^i p^(n-i) = q^n P(x) with q>0.
        Z p,value=table.front();
        mpz_divexact(p.get_mpz_t(),q.get_mpz_t(),x.get_den_mpz_t());
        p*=x.get_num();
        for(std::size_t i=1;i<table.size();++i) {
            mpz_mul(value.get_mpz_t(),value.get_mpz_t(),p.get_mpz_t());
            mpz_add(value.get_mpz_t(),value.get_mpz_t(),table[i].get_mpz_t());
        }
        return mpz_sgn(value.get_mpz_t());
    }
    unsigned oracle_worker_count(unsigned requested,std::size_t jobs) {
        if(requested>64)throw std::invalid_argument("oracle workers must be 0..64");
        if(jobs==0)return 0;
        unsigned workers=requested;
        if(workers==0)workers=std::min(64u,std::max(1u,std::thread::hardware_concurrency()));
        return static_cast<unsigned>(std::min<std::size_t>(workers,jobs));
    }
    OrderedRun run_ordered(std::size_t jobs,unsigned workers,const std::function<bool(std::size_t,const CancelCheck&)>& job) {
        std::vector<std::exception_ptr> errors(jobs);
        std::atomic<std::size_t> next{0},stop{jobs};
        // Indices are handed out in increasing order. Once one exceeds the lowest stop seen so far,
        // so does every later one, and every index below the final stop has been handed out and
        // run; the lowest stop is thus the serial one, whatever the timing.
        auto drain=[&] {
            for(;;) {
                const auto i=next.fetch_add(1);
                if(i>=jobs||i>stop.load())return;
                bool passed=false;
                // The stop only decreases, so once it is below i it stays there: job i is past
                // the final stop and its outcome, a Cancelled exception included, is discarded.
                const CancelCheck superseded=[&stop,i] {
                    return stop.load()<i;
                };
                try {
                    passed=job(i,superseded);
                }
                catch(...) {
                    errors[i]=std::current_exception();
                }
                if(!passed)for(auto lowest=stop.load();i<lowest&&!stop.compare_exchange_weak(lowest,i);) {}
            }
        };
        {
            // jthread destruction joins every helper before the shared state goes away. A helper
            // that cannot be created only leaves more of the work to the others and the caller.
            std::vector<std::jthread> helpers;
            const auto threads=std::min<std::size_t>(workers,jobs);
            helpers.reserve(threads);
            for(std::size_t w=1;w<threads;++w) {
                try {
                    helpers.emplace_back(drain);
                }
                catch(...) {
                    break;
                }
            }
            drain();
        }
        OrderedRun out;
        out.stop=stop.load();
        if(out.stop<jobs)out.error=errors[out.stop];
        return out;
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
        // The three primes below 2^32 (2^31-1, 2^32-5, 2^32-17) follow, in this order,
        // so an earlier witness is still reported; they make the exact fallback rare
        // when small primes divide the leading coefficient or discriminant.
        for(unsigned long prime:{2UL,3UL,5UL,7UL,11UL,65537UL,2147483647UL,4294967291UL,4294967279UL})
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
