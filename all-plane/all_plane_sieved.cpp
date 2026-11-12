// Accelerated exact whole-plane search with CRT generation and local obstructions.
// See all_plane_sieved.md for the generator proof and verification details.
// The exact reconstruction core below is shared in mathematical form with all_plane.cpp.
// Original reference core:
// Build: g++ -O3 -std=c++20 all_plane_sieved.cpp -o all_plane_sieved
// Run:   ./all_plane_sieved 100000 --validate
// This is an independent reference search, not the fast interior search.
//
// Symmetry: X <= Y <= N/2. Let A=distance to (0,0),
// p=distance to (N,0)-A, q=distance to (0,N)-A.
// Then 0 <= q <= p <= N and
//   Q=p^2+q^2-N^2, T=p(N^2-p^2)+q(N^2-q^2),
//   R=(N^2-p^2)^2+(N^2-q^2)^2,
//   4 Q A^2 - 4 T A + R = 0.
// K=T^2-Q R=(N^2-p^2)(N^2-q^2)(2N^2-(p-q)^2).
// For Q!=0, A=(T +/- sqrt(K))/(2Q).
// For Q==0,T!=0, A=R/(4T). Q==T==0 has R>0.
// Coordinates follow from X=(N^2-p^2-2Ap)/(2N), likewise Y.
// All roots, including boundary and median cases, are handled.
//
// Only in FOUR-distance search mode, proven global necessary conditions are
// used: 24|N; all four distances odd; X,Y have opposite parity, the even
// coordinate is divisible by 4; precisely one of X,Y is divisible by 3.
// Thus p,q are even. The validation instead enumerates ALL integer
// three-distance points, without these four-distance-only restrictions.
//
// Bound: K<=2N^6 and every nonnegative A root is <=2N^3.
// N<=1,000,000 keeps every intermediate inside signed __int128.
// Integer square tests use a floating estimate followed by exact correction.

#include <algorithm>
#include <array>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

using i128 = __int128_t;
using u128 = __uint128_t;
using i64 = int64_t;
using u64 = uint64_t;
using Points = std::set<std::pair<i64,i64>>;

struct Stats {
    u64 sides=0, parameter_pairs=0, linear_pairs=0;
    u64 square_discriminants=0, integral_A=0, integer_coordinates=0;
    u64 primitive_three_candidates=0, solutions=0;
};

static u64 isqrt(u128 z) {
    u64 r=static_cast<u64>(std::sqrt(static_cast<long double>(z)));
    while (u128(r)*r>z) --r;
    while (u128(r+1)*(r+1)<=z) ++r;
    return r;
}

static bool square(u128 z, u64 &root) {
    // These filters are exact. No floating-point equality is used.
    static constexpr unsigned q5=(1u<<0)|(1u<<1)|(1u<<4);
    static constexpr unsigned q7=(1u<<0)|(1u<<1)|(1u<<2)|(1u<<4);
    if (!((q5 >> unsigned(z%5))&1) || !((q7 >> unsigned(z%7))&1))
        return false;
    root=isqrt(z);
    return u128(root)*root==z;
}

static i64 magnitude(i64 x) { return x<0 ? -x : x; }

// Callback receives exact integer THREE-distance points.
template<class Callback>
static void evaluate(i64 N, i64 p, i64 q, bool four_filters,
                     Stats &st, Callback accept) {
    ++st.parameter_pairs;
    const i128 nn=i128(N)*N, pp=i128(p)*p, qq=i128(q)*q;
    const i128 U=nn-pp, V=nn-qq;
    const i128 Q=pp+qq-nn;
    const i128 T=i128(p)*U+i128(q)*V;
    const i128 R=U*U+V*V;
    i128 roots[2]; int nr=0;
    if (Q==0) {
        ++st.linear_pairs;
        if (T==0) {
            if (R==0) throw std::logic_error("identically zero quadratic");
            return;
        }
        if (R%(4*T)!=0) return;
        roots[nr++]=R/(4*T);
    } else {
        const i128 delta=i128(p)-q;
        const u128 K=u128(U*V*(2*nn-delta*delta));
        u64 sr;
        if (!square(K,sr)) return;
        ++st.square_discriminants;
        const i128 den=2*Q;
        const i128 plus=T+i128(sr), minus=T-i128(sr);
        if (plus%den==0) roots[nr++]=plus/den;
        if (sr!=0 && minus%den==0) roots[nr++]=minus/den;
    }
    for (int k=0;k<nr;++k) {
        const i128 A=roots[k];
        if (A<0) continue;
        ++st.integral_A;
        if (four_filters && A%2==0) continue;
        const i128 xn=U-2*A*p, yn=V-2*A*q, den=2*i128(N);
        if (xn%den || yn%den) continue;
        const i128 X=xn/den,Y=yn/den;
        if (X>Y || 2*Y>N) throw std::logic_error("symmetry failure");
        if (X*X+Y*Y!=A*A ||
            (N-X)*(N-X)+Y*Y!=(A+p)*(A+p) ||
            X*X+(N-Y)*(N-Y)!=(A+q)*(A+q))
            throw std::logic_error("three-distance identity failure");
        if (A>2*i128(N)*N*N) throw std::logic_error("radius bound failure");
        ++st.integer_coordinates;
        const i64 x=i64(X), y=i64(Y), a=i64(A);
        if (four_filters) {
            if (std::gcd(std::gcd(magnitude(x),magnitude(y)),N)!=1) continue;
            const bool xodd=(x%2!=0), yodd=(y%2!=0);
            if (xodd==yodd) continue;
            if ((xodd ? y : x)%4!=0) continue;
            if ((x%3==0)==(y%3==0)) continue;
            ++st.primitive_three_candidates;
        }
        accept(x,y,a);
    }
}

static void validate() {
    u64 total=0;
    // Independent coordinate enumeration over the PROVED complete bound.
    // It includes exterior points and does not use four-distance parity rules.
    for (i64 N=1;N<=12;++N) {
        Points algebraic,brute;
        Stats ignored;
        for (i64 p=0;p<=N;++p) for (i64 q=0;q<=p;++q)
            evaluate(N,p,q,false,ignored,[&](i64 x,i64 y,i64){
                algebraic.emplace(x,y);
            });
        const i64 radius=2*N*N*N;
        for (i64 y=-radius;y<=N/2;++y) for (i64 x=-radius;x<=y;++x) {
            u64 a,b,c;
            if (!square(u128(i128(x)*x+i128(y)*y),a)) continue;
            if (!square(u128(i128(N-x)*(N-x)+i128(y)*y),b)) continue;
            if (!square(u128(i128(x)*x+i128(N-y)*(N-y)),c)) continue;
            brute.emplace(x,y);
        }
        if (algebraic!=brute)
            throw std::runtime_error("independent validation mismatch at N="+std::to_string(N));
        total+=algebraic.size();
    }
    // Q=0 with an actual integral three-distance point, not merely no roots.
    Stats ignored;
    Points linear;
    evaluate(8400,6720,5040,false,ignored,[&](i64 x,i64 y,i64 a){
        if (a!=1685) throw std::logic_error("linear-case radius failure");
        linear.emplace(x,y);
    });
    if (linear!=Points{{164,1677}})
        throw std::runtime_error("nontrivial Q=0 validation failure");
    Points impossible;
    evaluate(24,24,0,false,ignored,[&](i64 x,i64 y,i64){impossible.emplace(x,y);});
    if (!impossible.empty()) throw std::runtime_error("Q=T=0 validation failure");
    std::cerr << "Validation passed: all integer three-distance points in the entire plane "
              << "for N=1..12 (" << total << " canonical points); "
              << "Q=0 integral example and Q=T=0 exception.\n";
}

using u32=std::uint32_t;

struct CRTPlan {
    std::array<u32,7> prime{}, power{}, idempotent{};
    unsigned count=0;
};

static u32 inverse_mod(u32 a,u32 modulus) {
    i64 r=modulus,s=a,t=0,w=1;
    while(s) {
        i64 q=r/s;
        i64 next_r=r-q*s; r=s; s=next_r;
        i64 next_t=t-q*w; t=w; w=next_t;
    }
    if(r!=1) throw std::logic_error("inverse of nonunit");
    t%=modulus; if(t<0) t+=modulus;
    return u32(t);
}

struct ArithmeticTables {
    std::vector<u32> spf,bad_part;
    std::vector<CRTPlan> plans;
    explicit ArithmeticTables(u32 bound):spf(bound+1),bad_part(bound+1,1),plans(bound+1) {
        for(u32 p=2;p<=bound;++p) if(!spf[p]) {
            spf[p]=p;
            if(u64(p)*p<=bound)
                for(u32 j=p*p;j<=bound;j+=p) if(!spf[j]) spf[j]=p;
        }
        bad_part[0]=0;
        for(u32 m=2;m<=bound;++m) {
            const u32 p=spf[m];
            bad_part[m]=bad_part[m/p]*((p%8==3 || p%8==5)?p:1);
            u32 rem=m;
            auto &plan=plans[m];
            while(rem>1) {
                const u32 prime=spf[rem];
                u32 power=1;
                do {power*=prime;rem/=prime;} while(rem%prime==0);
                if(plan.count>=plan.prime.size()) throw std::logic_error("CRT factor capacity");
                const unsigned i=plan.count++;
                plan.prime[i]=prime; plan.power[i]=power;
                const u32 other=m/power;
                plan.idempotent[i]=u32(u64(other)*inverse_mod(other%power,power)%m);
            }
        }
    }
};

// Integer coordinates require n | u(A+u), n | v(A+v), where
// n=N/2,u=p/2,v=q/2. Eliminating A gives the congruence compatibility
// gcd(n/gcd(n,u),n/gcd(n,v)) | (u-v).
// Set m=n/gcd(n,u). For every prime power t || m this is equivalent
// to v==0 or v==u (mod t). This generator enumerates precisely that set.
template<class Callback>
static u64 generate_crt(u32 n,u32 u,const ArithmeticTables &tables,Callback emit,u32 known_gcd=0) {
    const u32 m=n/(known_gcd?known_gcd:std::gcd(n,u));
    const auto &plan=tables.plans[m];
    std::array<u32,128> residues{};
    unsigned count=1;
    for(unsigned i=0;i<plan.count;++i) {
        if(u%plan.power[i]==0) continue; // The two local choices coincide.
        const u32 shift=u32(u64(u)*plan.idempotent[i]%m);
        const unsigned old=count;
        if(2*old>residues.size()) throw std::logic_error("CRT residue capacity");
        for(unsigned j=0;j<old;++j) {
            u32 r=residues[j]+shift;
            if(r>=m) r-=m;
            residues[count++]=r;
        }
    }
    for(unsigned i=0;i<count;++i)
        for(u32 v=residues[i];v<=u;v+=m) emit(v);
    return count;
}

// For primes 2 or 3 mod 4 dividing n, A must be a unit in a primitive
// four-distance configuration. Thus u (and v) is either a unit or a
// multiple of the entire prime power in n; partial valuations fail.
static std::vector<std::uint8_t> permitted_valuations(u32 n,const CRTPlan &plan) {
    std::vector<std::uint8_t> permitted(n+1,1);
    for(unsigned i=0;i<plan.count;++i) {
        const u32 p=plan.prime[i],power=plan.power[i];
        if(!(p==2 || p%4==3) || power==p) continue;
        for(u32 base=0;base<=n;base+=power)
            for(u32 j=p;j<power && base+j<=n;j+=p) permitted[base+j]=0;
    }
    // If 2^e||n, e>=2, valuation exactly e+1 would make X==2 mod4.
    if(n%4==0) {
        const u32 two_part=n&(~n+1);
        for(u32 u=2*two_part;u<=n;u+=4*two_part) permitted[u]=0;
    }
    return permitted;
}

// X-Y=-(u-v)(A+u+v)/n. X-Y is odd and has no bad (3,5 mod8)
// prime factor. Therefore bad_part(u-v) divides n; if u-v is even,
// its 2-adic valuation must equal that of n. Delta zero is a diagonal.
static bool permitted_difference(u32 n,u32 u,u32 v,const ArithmeticTables &tables) {
    const u32 delta=u-v;
    if(!delta) return false;
    if(!(delta&1)) {
        const u32 two_part=n&(~n+1);
        if((delta & (2*two_part-1))!=two_part) return false;
    }
    const u32 bad=tables.bad_part[delta];
    return bad==1 || n%bad==0;
}

// Let g=gcd(n,u,v), n0=n/g,u0=u/g,v0=v/g and Q0=u0^2+v0^2-n0^2.
// With h=gcd(A,g), coordinate primitivity implies h|n0. Dividing the norm
// equation by h^2 and reducing mod g/h gives g/h|Q0. Thus
// g/gcd(g,n0)|Q0 is necessary. A is odd, so if 2^e||g then 2^(e+1)|Q0.
static bool permitted_scale(u32 n,u32 u,u32 v,u32 gcd_nu) {
    if(gcd_nu==1) return true;
    const u32 g=std::gcd(gcd_nu,v);
    if(g==1) return true;
    const i64 n0=n/g,u0=u/g,v0=v/g;
    const i64 Q0=u0*u0+v0*v0-n0*n0;
    const u32 required=g/std::gcd(g,u32(n0));
    if(required>1 && Q0%required!=0) return false;
    if(!(g&1)) {
        const u32 two_part=g&(~g+1);
        if(Q0%(2*two_part)!=0) return false;
    }
    return true;
}

// Full finite-field feasibility, not just a discriminant square test.
// For N=1 enumerate all x,y and all signed roots A,B,C; D needs to exist.
// Record (B-A,C-A). Nonzero N follows by scaling. For N=0 the necessary
// and sufficient field condition is p*q*(p-q)==0.
template<unsigned P>
struct PrimeProjection {
    std::array<std::uint8_t,P*P*P> cube{};
    std::array<std::uint8_t,P*P> any{};
    std::array<std::uint8_t,P> qr{};
    PrimeProjection() {
        std::array<std::vector<unsigned>,P> roots;
        for(unsigned a=0;a<P;++a) {roots[a*a%P].push_back(a);qr[a*a%P]=1;}
        std::array<std::uint8_t,P*P> normal{};
        for(unsigned x=0;x<P;++x) for(unsigned y=0;y<P;++y) {
            const unsigned xx=(1+P-x)%P,yy=(1+P-y)%P;
            const unsigned aa=(x*x+y*y)%P,bb=(xx*xx+y*y)%P;
            const unsigned cc=(x*x+yy*yy)%P,dd=(xx*xx+yy*yy)%P;
            if(!qr[dd]) continue;
            for(unsigned a:roots[aa]) for(unsigned b:roots[bb]) for(unsigned c:roots[cc])
                normal[((b+P-a)%P)*P+(c+P-a)%P]=1;
        }
        for(unsigned n=0;n<P;++n) {
            const unsigned inv=n?inverse_mod(n,P):0;
            for(unsigned u=0;u<P;++u) for(unsigned v=0;v<P;++v) {
                const bool okay=n?normal[((u*inv)%P)*P+(v*inv)%P]
                    :(u*v*(u+P-v))%P==0;
                cube[(n*P+u)*P+v]=okay;
                any[n*P+u]|=okay;
            }
        }
    }
    const std::uint8_t *row(u32 n,u32 u) const {return cube.data()+((n%P)*P+u%P)*P;}
    bool row_exists(u32 n,u32 u) const {return any[(n%P)*P+u%P];}
    bool accepts(u32 n,u32 u,u32 v) const {return row(n,u)[v%P];}
    void verify() const {
        // Independent reconstruction from (n,u,v,A), rather than x,y enumeration.
        auto mod=[](int z){z%=int(P);return unsigned(z<0?z+P:z);};
        for(unsigned n=0;n<P;++n) for(unsigned u=0;u<P;++u) for(unsigned v=0;v<P;++v) {
            bool expected=false;
            const unsigned inv=n?inverse_mod((2*n)%P,P):0;
            for(unsigned a=0;a<P && !expected;++a) {
                if(!n) {expected=(u*(2*a+u))%P==0 && (v*(2*a+v))%P==0;continue;}
                const unsigned x=mod(int(n*n)-int(u*u)-int(2*a*u))*inv%P;
                const unsigned y=mod(int(n*n)-int(v*v)-int(2*a*v))*inv%P;
                const unsigned xx=(n+P-x)%P,yy=(n+P-y)%P;
                expected=(x*x+y*y)%P==a*a%P && qr[(xx*xx+yy*yy)%P];
            }
            if(expected!=bool(cube[(n*P+u)*P+v]))
                throw std::runtime_error("finite-field projection mismatch mod "+std::to_string(P));
        }
    }
};

struct LocalSieve {
    PrimeProjection<5> p5;
    PrimeProjection<7> p7;
    PrimeProjection<11> p11;
    PrimeProjection<13> p13;
    PrimeProjection<17> p17;
    PrimeProjection<19> p19;
    PrimeProjection<23> p23;
    PrimeProjection<29> p29;
    PrimeProjection<31> p31;
    PrimeProjection<37> p37;
    PrimeProjection<41> p41;
    PrimeProjection<43> p43;
    bool rows_possible(u32 n,u32 u) const {
        return p11.row_exists(n,u) && p19.row_exists(n,u);
    }
    bool accepts(u32 n,u32 u,u32 v) const {
        return p5.accepts(n,u,v) && p11.accepts(n,u,v) && p19.accepts(n,u,v)
            && p7.accepts(n,u,v) && p29.accepts(n,u,v) && p37.accepts(n,u,v)
            && p13.accepts(n,u,v) && p17.accepts(n,u,v) && p23.accepts(n,u,v)
            && p31.accepts(n,u,v) && p41.accepts(n,u,v) && p43.accepts(n,u,v);
    }
    void verify() const {
        p5.verify();p7.verify();p11.verify();p13.verify();p17.verify();p19.verify();
        p23.verify();p29.verify();p31.verify();p37.verify();p41.verify();p43.verify();
    }
};

struct SearchOptions {
    u32 bound=10000,min_denominator=1;
    bool local=true,difference=true,valuation=true,scale=true,dense=false;
};

struct SearchStats {
    u64 dense_even_pairs=0,rows=0,valuation_rows=0,local_rows=0;
    u64 crt_residues=0,generated_pairs=0,valuation_pairs=0,difference_pairs=0,scale_pairs=0,local_pairs=0;
    Stats exact;
};

template<class Callback>
static SearchStats run_search(const SearchOptions &opt,const ArithmeticTables &arithmetic,
                              const LocalSieve &local,Callback accept) {
    SearchStats out;
    const u32 first_n=12*((opt.min_denominator+23)/24);
    for(u32 n=first_n;n<=opt.bound/2;n+=12) {
        const u32 N=2*n;
        ++out.exact.sides;
        out.dense_even_pairs+=u64(n+1)*(n+2)/2;
        auto permitted=permitted_valuations(n,arithmetic.plans[n]);
        std::vector<std::uint8_t> permitted_delta(n+1,1);
        if(opt.difference) for(u32 delta=0;delta<=n;++delta)
            permitted_delta[delta]=permitted_difference(n,delta,0,arithmetic);
        for(u32 u=0;u<=n;++u) {
            ++out.rows;
            if(opt.valuation && !permitted[u]) continue;
            ++out.valuation_rows;
            if(opt.local && !local.rows_possible(n,u)) continue;
            ++out.local_rows;
            const u32 gcd_nu=std::gcd(n,u);
            auto visit=[&](u32 v) {
                ++out.generated_pairs;
                if(opt.valuation && !permitted[v]) return;
                ++out.valuation_pairs;
                if(opt.difference && !permitted_delta[u-v]) return;
                ++out.difference_pairs;
                if(opt.local && !local.accepts(n,u,v)) return;
                ++out.local_pairs;
                if(opt.scale && !permitted_scale(n,u,v,gcd_nu)) return;
                ++out.scale_pairs;
                evaluate(N,2*u,2*v,true,out.exact,[&](i64 x,i64 y,i64 a){
                    accept(N,2*u,2*v,x,y,a,out.exact);
                });
            };
            if(opt.dense) for(u32 v=0;v<=u;++v) visit(v);
            else out.crt_residues+=generate_crt(n,u,arithmetic,visit,gcd_nu);
        }
    }
    return out;
}

static void validate_accelerated(const ArithmeticTables &arithmetic,const LocalSieve &local) {
    // Regression: g|Q0 is too strong; the correct requirement is
    // g/gcd(g,n0)|Q0. This real primitive THREE-distance point must survive.
    if(!permitted_scale(10140,3900,1235,std::gcd(10140u,3900u)))
        throw std::runtime_error("common-scale regression failure");
    // Check the actual generated sets, including duplicates, for every small n,u.
    u64 crt_pairs=0;
    for(u32 n=1;n<=256;++n) for(u32 u=0;u<=n;++u) {
        std::vector<u32> actual,expected;
        generate_crt(n,u,arithmetic,[&](u32 v){actual.push_back(v);});
        std::sort(actual.begin(),actual.end());
        for(u32 v=0;v<=u;++v) {
            const u32 d=std::gcd(n/std::gcd(n,u),n/std::gcd(n,v));
            if((u-v)%d==0) expected.push_back(v);
        }
        if(actual!=expected) throw std::runtime_error("CRT candidate set mismatch");
        crt_pairs+=actual.size();
    }
    local.verify();
    // Retain the baseline's nonzero three-distance examples: the extra GOOD
    // and finite-field four-distance restrictions are disabled for this check.
    using Key=std::array<i64,6>;
    std::set<Key> reference,accelerated;
    Stats ignored;
    for(i64 N=24;N<=1000;N+=24)
        for(i64 p=0;p<=N;p+=2) for(i64 q=0;q<=p;q+=2)
            evaluate(N,p,q,true,ignored,[&](i64 x,i64 y,i64 a){
                reference.insert({N,p,q,x,y,a});
            });
    SearchOptions opt;
    opt.bound=1000;opt.local=false;opt.difference=false;
    run_search(opt,arithmetic,local,[&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats &){
        accelerated.insert({N,p,q,x,y,a});
    });
    if(reference!=accelerated || reference.empty())
        throw std::runtime_error("baseline three-distance candidate mismatch");
    // Extra filters must yield the identical parameter set whether generation
    // is dense with compatibility checked or CRT based.
    u64 filtered_pairs=0;
    for(u32 n=12;n<=1000;n+=12) {
        auto permitted=permitted_valuations(n,arithmetic.plans[n]);
        for(u32 u=0;u<=n;++u) {
            if(!permitted[u] || !local.rows_possible(n,u)) continue;
            std::vector<u32> actual,expected;
            auto allowed=[&](u32 v){return permitted[v] &&
                permitted_difference(n,u,v,arithmetic) &&
                permitted_scale(n,u,v,std::gcd(n,u)) && local.accepts(n,u,v);};
            generate_crt(n,u,arithmetic,[&](u32 v){if(allowed(v))actual.push_back(v);});
            for(u32 v=0;v<=u;++v) {
                const u32 d=std::gcd(n/std::gcd(n,u),n/std::gcd(n,v));
                if((u-v)%d==0 && allowed(v)) expected.push_back(v);
            }
            std::sort(actual.begin(),actual.end());
            if(actual!=expected) throw std::runtime_error("filtered generator mismatch");
            filtered_pairs+=actual.size();
        }
    }
    std::cerr << "Accelerated validation passed: CRT exact sets for every n=1..256 ("
              << crt_pairs << " pairs); finite-field tables independently checked; "
              << reference.size() << " baseline primitive three-distance examples retained; "
              << filtered_pairs << " small-bound filtered parameter pairs matched.\n";
}

int main(int argc,char **argv) try {
    SearchOptions opt;
    bool check=false;
    std::string audit_path;
    for(int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if(arg=="--validate") check=true;
        else if(arg=="--no-local") opt.local=false;
        else if(arg=="--no-difference") opt.difference=false;
        else if(arg=="--no-valuation") opt.valuation=false;
        else if(arg=="--no-scale") opt.scale=false;
        else if(arg=="--dense-reference") opt.dense=true;
        else if(arg=="--min-denominator") {
            if(++i>=argc) throw std::runtime_error("missing minimum denominator");
            const auto value=std::stoull(argv[i]);
            if(value<1 || value>1000000) throw std::runtime_error("minimum denominator must be in [1,1000000]");
            opt.min_denominator=u32(value);
        }
        else if(arg=="--audit-three") {
            if(++i>=argc) throw std::runtime_error("missing audit filename");
            audit_path=argv[i];
        } else if(arg=="--help") {
            std::cout << "all_plane_sieved [B] [--min-denominator L] [--validate] [--audit-three CSV] "
                         "[--no-local] [--no-difference] [--no-valuation] [--no-scale] [--dense-reference]\n";
            return 0;
        } else {
            const auto value=std::stoull(arg);
            if(value<1 || value>1000000) throw std::runtime_error("B must be in [1,1000000]");
            opt.bound=u32(value);
        }
    }
    if(opt.min_denominator>opt.bound) throw std::runtime_error("minimum denominator exceeds maximum");
    const auto start=std::chrono::steady_clock::now();
    ArithmeticTables arithmetic(std::max<u32>(1000,opt.bound/2));
    LocalSieve local;
    const auto setup_done=std::chrono::steady_clock::now();
    if(check) {validate();validate_accelerated(arithmetic,local);}
    std::ofstream audit;
    if(!audit_path.empty()) {
        audit.open(audit_path);
        if(!audit) throw std::runtime_error("cannot open audit CSV");
        audit << "N,p,q,X,Y,A\n";
    }
    const auto search_started=std::chrono::steady_clock::now();
    const auto stats=run_search(opt,arithmetic,local,
        [&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats &exact) {
            if(audit.is_open()) audit << N << ',' << p << ',' << q << ',' << x << ',' << y << ',' << a << '\n';
            u64 d;
            if(!square(u128(i128(N-x)*(N-x)+i128(N-y)*(N-y)),d)) return;
            ++exact.solutions;
            std::cout << "{\"solution\":{\"N\":" << N << ",\"X\":" << x
                      << ",\"Y\":" << y << ",\"A\":" << a << ",\"D\":" << d << "}}\n";
        });
    const auto finished=std::chrono::steady_clock::now();
    auto elapsed=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
    const auto &exact=stats.exact;
    std::cout << "{\"scope\":\"whole_plane\",\"max_minimal_coordinate_denominator\":" << opt.bound
              << ",\"min_denominator\":" << opt.min_denominator
              << ",\"generator\":\"" << (opt.dense?"dense_reference":"CRT_coordinate_integrality") << "\""
              << ",\"local_filters\":" << (opt.local?"true":"false")
              << ",\"difference_filters\":" << (opt.difference?"true":"false")
              << ",\"valuation_filters\":" << (opt.valuation?"true":"false")
              << ",\"scale_filter\":" << (opt.scale?"true":"false")
              << ",\"dense_even_pairs\":" << stats.dense_even_pairs
              << ",\"u_rows\":" << stats.rows
              << ",\"valuation_rows\":" << stats.valuation_rows
              << ",\"local_rows\":" << stats.local_rows
              << ",\"crt_residue_choices\":" << stats.crt_residues
              << ",\"generated_pairs\":" << stats.generated_pairs
              << ",\"valuation_pairs\":" << stats.valuation_pairs
              << ",\"difference_pairs\":" << stats.difference_pairs
              << ",\"local_pairs\":" << stats.local_pairs
              << ",\"scale_pairs\":" << stats.scale_pairs
              << ",\"exact_evaluations\":" << exact.parameter_pairs
              << ",\"square_discriminants\":" << exact.square_discriminants
              << ",\"linear_pairs\":" << exact.linear_pairs
              << ",\"primitive_three_candidates\":" << exact.primitive_three_candidates
              << ",\"solutions\":" << exact.solutions
              << ",\"seconds_setup\":" << elapsed(start,setup_done)
              << ",\"seconds_search\":" << elapsed(search_started,finished)
              << ",\"seconds_total_excluding_validation\":" << elapsed(start,setup_done)+elapsed(search_started,finished)
              << "}\n";
    return 0;
} catch(const std::exception &e) {std::cerr << e.what() << '\n';return 1;}

