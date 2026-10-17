#include "secular/exact.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <limits>

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
    CoefficientResult characteristic(const ZMatrix& e,const ZMatrix& f) {
        if(e.n!=f.n)throw std::invalid_argument("pencil dimensions differ");
        Z d=determinant(e);
        if(d==0)throw std::invalid_argument("singular mass/leading matrix");
        ZMatrix j=adjugate(e),s=j;
        CoefficientResult r;
        r.coefficients.push_back(d);
        for(std::size_t k=1;k<=e.n;++k) {
            ZMatrix t=multiply(f,s);
            Z trace=0;
            for(std::size_t i=0;i<e.n;++i)trace+=t(i,i);
            Z c=exact_divide(trace,Z(static_cast<unsigned long>(k)));
            ++r.exact_divisions;
            r.coefficients.push_back(-c);
            for(std::size_t i=0;i<e.n;++i)t(i,i)-=c;
            for(auto& v:t.data)r.max_workspace_bits=std::max(r.max_workspace_bits,bits(v));
            s=apply_e_inverse(e,j,d,t,&r.exact_divisions);
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
    bool verify_coefficients(const ZMatrix& e,const ZMatrix& f,const std::vector<Z>& c) {
        if(e.n!=f.n||c.size()!=e.n+1)return false;
        for(std::size_t k=0;k<=e.n;++k) {
            Z x=-Z(static_cast<unsigned long>(k+1));
            ZMatrix a(e.n);
            for(std::size_t i=0;i<a.data.size();++i)a.data[i]=x*e.data[i]-f.data[i];
            if(determinant(a)!=evaluate(c,x))return false;
        }
        return true;
    }
    bool has_repeated_roots(const std::vector<Z>& c) {
        if(c.size()<2||c.front()==0)throw std::invalid_argument("invalid polynomial");
        std::vector<Z> a(c.rbegin(),c.rend()),b;
        for(std::size_t i=1;i<a.size();++i)b.push_back(a[i]*static_cast<unsigned long>(i));
        primitive(a);
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
        return a.size()>1;
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
