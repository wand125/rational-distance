// Version 3: critical 5/13-content gap rows and complementary
// denominator gamma conditions. Derived from the frozen v2 generator; the
// reconstruction core and all pre-existing mathematical predicates are shared.
// Do not mistake surviving THREE-distance diagnostics for all three-distance points.
// New FOUR-distance conditions can legitimately remove old diagnostic points.
// The mod-19 progression wheel is the independently checked implementation
// from research_wheel/wheel_suffix.cpp, copied here to keep one source dependency.
// Arithmetic refinements are consequences of the prior denominator and critical
// content results; see prior_arithmetic_overlap.md and research_family_refinement/.
// Build: g++ -O3 -std=c++17 all_plane_gap_v3.cpp -o all_plane_gap_v3
// The reconstruction and necessary-condition predicates are shared with the
// independently verified u-first implementation; this file owns its generator.
// This entry point allows B<=2,000,000. The included reference's old 1,000,000
// comment/entry-point limit is conservative: universally A<=N^3. For Q>0,
// A<=T/Q<=N^3 because t(N-t)(N+t)<=N^3/2 for 0<=t<=N; for Q<0,
// A<=sqrt(K)/2<=N^3/sqrt(2); for Q=0, A<=N^2/4. Thus coordinates
// fit int64, and K<=2N^6 plus all norm sums fit signed int128 at this cap.
// The first power-of-two failure for signed K is N=2^21; the chosen cap is lower.
#define main included_sieved_main
#include "all_plane_sieved.cpp"
#undef main

template<unsigned P>
struct GapProjection {
    std::array<std::uint8_t,P*P> possible{};
    explicit GapProjection(const PrimeProjection<P> &projection) {
        for(unsigned n=0;n<P;++n) for(unsigned h=0;h<P;++h)
            for(unsigned u=0;u<P;++u)
                possible[n*P+h]|=projection.accepts(n,u,(u+P-h)%P);
    }
    bool accepts(u32 n,u32 h) const {return possible[(n%P)*P+h%P];}
};

struct GapWheel {
    static constexpr u32 modulus=5*11*19;
    std::array<std::vector<u32>,modulus> residues;
    GapProjection<5> p5;
    GapProjection<11> p11;
    GapProjection<19> p19;
    explicit GapWheel(const LocalSieve &local):p5(local.p5),p11(local.p11),p19(local.p19) {
        for(u32 n=0;n<modulus;++n) for(u32 h=0;h<modulus;++h)
            if(p5.accepts(n,h) && p11.accepts(n,h) && p19.accepts(n,h))
                residues[n].push_back(h);
    }
    template<class Callback>
    void generate(u32 n,Callback emit) const {
        const auto &row=residues[n%modulus];
        for(u32 block=0;block<=n;block+=modulus)
            for(u32 r:row) {
                const u32 h=block+r;
                if(h>n) break;
                emit(h);
            }
    }
};

// Fix h=u-v. At each prime ell^e || n, compatibility of the two A
// congruences says e-max(val(u),val(v)) <= val(h). Equivalently, with
// t=ell^max(e-val(h),0), u==0 or u==h (mod t). Combine those choices
// by CRT, and generate u in [h,n]. This is an exact change of variables.
template<class Callback>
static u64 generate_gap_crt(u32 n,u32 h,const ArithmeticTables &tables,
                            Callback emit,u32 known_gcd=0) {
    const u32 m=n/(known_gcd?known_gcd:std::gcd(n,h));
    const auto &plan=tables.plans[m];
    std::array<u32,128> residues{};
    unsigned count=1;
    for(unsigned i=0;i<plan.count;++i) {
        if(h%plan.power[i]==0) continue;
        const u32 shift=u32(u64(h)*plan.idempotent[i]%m);
        const unsigned old=count;
        if(2*old>residues.size()) throw std::logic_error("gap CRT capacity");
        for(unsigned j=0;j<old;++j) {
            u32 r=residues[j]+shift;
            if(r>=m) r-=m;
            residues[count++]=r;
        }
    }
    for(unsigned i=0;i<count;++i) {
        u32 first=residues[i];
        if(first<h) first+=((h-first+m-1)/m)*m;
        for(u32 u=first;u<=n;u+=m) emit(u);
    }
    return count;
}

// Append after the frozen v2 definitions, excluding v2's main function.
template<unsigned W> struct ProgressionWheel {
    std::vector<std::uint8_t> next;
    explicit ProgressionWheel(const LocalSieve &local):next(std::size_t(W)*W*W*W,255) {
        for(unsigned n=0;n<W;++n) for(unsigned h=0;h<W;++h) {
            std::array<bool,W> allowed{};
            for(unsigned u=0;u<W;++u) {
                const unsigned v=(u+W-h)%W;
                bool yes=true;
                if constexpr(W%5==0) yes=yes && local.p5.accepts(n,u,v);
                if constexpr(W%11==0) yes=yes && local.p11.accepts(n,u,v);
                if constexpr(W%19==0) yes=yes && local.p19.accepts(n,u,v);
                allowed[u]=yes;
            }
            for(unsigned step=0;step<W;++step) {
                auto *row=next.data()+(((std::size_t(n)*W+h)*W+step)*W);
                const unsigned cycles=std::gcd(step,W),len=W/cycles;
                for(unsigned first=0;first<cycles;++first) {
                    int nearest=-1;
                    for(int j=int(2*len)-1;j>=0;--j) {
                        const unsigned u=(first+unsigned(j)*step)%W;
                        if(allowed[u]) nearest=j;
                        if(j<int(len)) row[u]=nearest<0?255:std::uint8_t(nearest-j);
                    }
                }
            }
        }
    }
    const std::uint8_t *row(u32 n,u32 h,u32 step) const {
        return next.data()+(((std::size_t(n%W)*W+h%W)*W+step%W)*W);
    }
};

template<unsigned W,class Callback>
static u64 generate_field_gap_crt(u32 n,u32 h,const ArithmeticTables &tables,
                                 const ProgressionWheel<W> &wheel,Callback emit,u32 known_gcd=0) {
    const u32 m=n/(known_gcd?known_gcd:std::gcd(n,h));
    const auto &plan=tables.plans[m];
    std::array<u32,128> residues{};
    unsigned count=1;
    for(unsigned i=0;i<plan.count;++i) {
        if(h%plan.power[i]==0) continue;
        const u32 shift=u32(u64(h)*plan.idempotent[i]%m);
        const unsigned old=count;
        if(2*old>residues.size()) throw std::logic_error("CRT capacity");
        for(unsigned j=0;j<old;++j) {
            u32 r=residues[j]+shift;
            if(r>=m) r-=m;
            residues[count++]=r;
        }
    }
    const auto *row=wheel.row(n,h,m);
    for(unsigned i=0;i<count;++i) {
        u32 first=residues[i];
        if(first<h) first+=((h-first+m-1)/m)*m;
        u64 u=first;
        while(u<=n) {
            const unsigned skip=row[u%W];
            if(skip==255) break;
            u+=u64(skip)*m;
            if(u>n) break;
            emit(u32(u));
            u+=m;
        }
    }
    return count;
}


template<unsigned W>
static void validate_field_wheel(const ArithmeticTables &arithmetic,const LocalSieve &local,
                                 const ProgressionWheel<W> &wheel) {
    u64 total=0;
    for(u32 n=1;n<=128;++n) for(u32 h=0;h<=n;++h) {
        std::vector<u32> actual,expected;
        generate_field_gap_crt(n,h,arithmetic,wheel,[&](u32 u){actual.push_back(u);});
        generate_gap_crt(n,h,arithmetic,[&](u32 u){
            bool yes=true;
            if constexpr(W%5==0) yes=yes && local.p5.accepts(n,u,u-h);
            if constexpr(W%11==0) yes=yes && local.p11.accepts(n,u,u-h);
            if constexpr(W%19==0) yes=yes && local.p19.accepts(n,u,u-h);
            if(yes) expected.push_back(u);
        });
        std::sort(actual.begin(),actual.end());std::sort(expected.begin(),expected.end());
        if(actual!=expected) throw std::logic_error("progression wheel mismatch");
        total+=actual.size();
    }
    std::cerr<<"Wheel validation passed: "<<total<<" emitted CRT members.\n";
}

// For g=gcd(n,t), the common divisor of the complementary denominators is
// gamma_n(t)=gcd(g,n/g). Any primitive four-distance point has the same gamma
// for u, v and h=u-v, and that common content has primes only 1 (mod 4).
static u32 gamma_value(u32 n,u32 t) {
    const u32 g=std::gcd(n,t);
    return std::gcd(g,n/g);
}
static bool permitted_content(u32 g,const ArithmeticTables &tables) {
    const auto &plan=tables.plans[g];
    for(unsigned i=0;i<plan.count;++i)
        if(plan.prime[i]%4!=1) return false;
    return true;
}

// At p=5 or 13, v_p(N)=2*v_p(G)>0 is impossible. This is an algebraic
// finite-field obstruction and does not use the point being inside the square.
// gamma_n(h) predicts v_p(G)=min(k,e-k). For even e, k=e/2 is forbidden.
struct CriticalContentRows {
    std::array<u32,2> midpoint{},next{};
    unsigned count=0;
    explicit CriticalContentRows(const CRTPlan &plan) {
        for(unsigned i=0;i<plan.count;++i) {
            const u32 prime=plan.prime[i];
            if(prime!=5 && prime!=13) continue;
            u32 power=plan.power[i];unsigned exponent=0;
            for(u32 t=power;t>1;t/=prime) ++exponent;
            if(exponent%2 || exponent==0) continue;
            u32 mid=1;
            for(unsigned j=0;j<exponent/2;++j) mid*=prime;
            midpoint[count]=mid;next[count]=mid*prime;++count;
        }
    }
    void mark(std::vector<std::uint8_t> &permitted,u32 n) const {
        for(unsigned i=0;i<count;++i)
            for(u32 t=midpoint[i];t<=n;t+=midpoint[i])
                if(t%next[i]) permitted[t]=0;
    }
    bool accepts(u32 h) const {
        for(unsigned i=0;i<count;++i)
            if(h%midpoint[i]==0 && h%next[i]!=0) return false;
        return true;
    }
};

// Only squared primes 1 mod4 can contribute nontrivial gamma after the
// existing CRT and valuation filters. At every other prime the gamma is 1.
static u32 gamma_relevant_modulus(const CRTPlan &plan) {
    u32 modulus=1;
    for(unsigned i=0;i<plan.count;++i)
        if(plan.prime[i]%4==1 && plan.power[i]>plan.prime[i]) modulus*=plan.power[i];
    return modulus;
}

// At a critical prime p with e=v_p(n)=2s>0 and v_p(h)=s, four normalized
// squared radii force z=(u/p^s)/(h/p^s) into the masks below. These masks
// are independently obtained from the complete finite-field critical model;
// see research_family_refinement/critical_projection.py and its certificate.
// Unlike the existing mod-p sieve, this restriction observes digits after
// removing the common p^s factor.
static u64 critical_mask(u32 prime) {
    if(prime==17) return 0x15554ULL; // {2,4,6,8,10,12,14,16}
    if(prime==29) return ((1ULL<<29)-1) & ~((1ULL<<0)|(1ULL<<1)|(1ULL<<4)|(1ULL<<8)|(1ULL<<22)|(1ULL<<26));
    if(prime==37) return ((1ULL<<37)-1) & ~((1ULL<<0)|(1ULL<<1)|(1ULL<<2)|(1ULL<<4)|(1ULL<<8)|(1ULL<<12)|(1ULL<<17)|(1ULL<<21)|(1ULL<<26)|(1ULL<<30)|(1ULL<<34)|(1ULL<<36));
    if(prime==41) return ((1ULL<<41)-1) & ~((1ULL<<0)|(1ULL<<1)|(1ULL<<2)|(1ULL<<3)|(1ULL<<4)|(1ULL<<5)|(1ULL<<6)|(1ULL<<9)|(1ULL<<10)|(1ULL<<16)|(1ULL<<21)|(1ULL<<26)|(1ULL<<32)|(1ULL<<33)|(1ULL<<36)|(1ULL<<37)|(1ULL<<38)|(1ULL<<39)|(1ULL<<40));
    return 0;
}
struct CriticalProjectionPlan {
    std::array<u32,4> prime{},midpoint{};
    std::array<u64,4> mask{};
    unsigned count=0;
    explicit CriticalProjectionPlan(const CRTPlan &plan) {
        for(unsigned i=0;i<plan.count;++i) {
            const u32 p=plan.prime[i];const u64 bits=critical_mask(p);
            if(!bits || plan.power[i]==p) continue;
            unsigned e=0;for(u32 t=plan.power[i];t>1;t/=p) ++e;
            if(e%2) continue;
            u32 mid=1;for(unsigned j=0;j<e/2;++j) mid*=p;
            prime[count]=p;midpoint[count]=mid;mask[count]=bits;++count;
        }
    }
    bool accepts(u32 u,u32 h,u32 *bad_prime=nullptr,u32 *bad_z=nullptr) const {
        for(unsigned i=0;i<count;++i) {
            const u32 p=prime[i],mid=midpoint[i];
            if(h%mid || h%(mid*p)==0) continue;
            if(u%mid) {if(bad_prime)*bad_prime=p;if(bad_z)*bad_z=p;return false;}
            const u32 z=u64((u/mid)%p)*inverse_mod((h/mid)%p,p)%p;
            if(!((mask[i]>>z)&1)) {if(bad_prime)*bad_prime=p;if(bad_z)*bad_z=z;return false;}
        }
        return true;
    }
};

struct GapStats {
    u64 dense_even_pairs=0,gap_rows=0,local_gap_rows=0,valuation_gap_rows=0,difference_gap_rows=0;
    u64 crt_residues=0,generated_pairs=0,valuation_pairs=0,scale_pairs=0,local_pairs=0;
    u64 critical_gap_rows=0,gamma_pairs=0,gamma_rejected=0,critical_projection_pairs=0,critical_projection_rejected=0;
    std::vector<std::array<u32,6>> gamma_rejection_examples;
    std::vector<std::array<u32,5>> critical_projection_rejection_examples;
    Stats exact;
};

template<class ParameterCallback,class PointCallback>
static GapStats run_gap_search(const SearchOptions &opt,const ArithmeticTables &arithmetic,
                               const LocalSieve &local,const GapWheel &wheel,
                               ParameterCallback parameter,PointCallback accept,
                               bool gap_valuation=false,bool gamma=true,bool critical=true,
                               const ProgressionWheel<19>* fieldwheel=nullptr,bool critical_projections=true) {
    GapStats out;
    const u32 first_n=12*((opt.min_denominator+23)/24);
    for(u32 n=first_n;n<=opt.bound/2;n+=12) {
        const u32 N=2*n;
        ++out.exact.sides;
        out.dense_even_pairs+=u64(n+1)*(n+2)/2;
        out.gap_rows+=n+1;
        auto permitted=permitted_valuations(n,arithmetic.plans[n]);
        const CriticalContentRows critical_rows(arithmetic.plans[n]);
        if(!opt.valuation) std::fill(permitted.begin(),permitted.end(),1);
        if(critical) critical_rows.mark(permitted,n);
        const u32 relevant=opt.valuation?gamma_relevant_modulus(arithmetic.plans[n]):n;
        const CriticalProjectionPlan critical_projection(arithmetic.plans[n]);
        auto visit_gap=[&](u32 h) {
            ++out.local_gap_rows;
            // At primes 2 or 3 mod4 dividing n, A is a unit and each of
            // u,v is 0 or -A modulo the full prime power. Thus h=u-v
            // is either a unit or a multiple of that full prime power.
            // The table's extra 2-adic rejection is also safe for h:
            // when h is even, n(X-Y)=-h(A+u+v) and X-Y,A+u+v are odd,
            // so val_2(h)=val_2(n), never val_2(n)+1.
            // With all default filters this precheck is implied by the
            // existing later predicates, but it avoids constructing CRT rows.
            if(opt.valuation && gap_valuation && !permitted[h]) return;
            ++out.valuation_gap_rows;
            if(opt.difference && !permitted_difference(n,h,0,arithmetic)) return;
            ++out.difference_gap_rows;
            if(critical && !critical_rows.accepts(h)) return;
            ++out.critical_gap_rows;
            const u32 gcd_nh=std::gcd(n,h);
            auto visit_pair=[&](u32 u) {
                const u32 v=u-h;
                ++out.generated_pairs;
                if((opt.valuation || critical) && (!permitted[u] || !permitted[v])) return;
                ++out.valuation_pairs;
                if(opt.local && !local.accepts(n,u,v)) return;
                ++out.local_pairs;
                // gcd(gcd(n,h),v)=gcd(n,u,v), so the cached gcd_nh can
                // safely replace the predicate's usual gcd(n,u) argument.
                // Only field-compatible candidates pay for these gcds.
                if(opt.scale && !permitted_scale(n,u,v,gcd_nh)) return;
                ++out.scale_pairs;
                if(gamma && relevant>1) {
                    const u32 gh=gamma_value(relevant,h);
                    const u32 gu=gamma_value(relevant,u),gv=gamma_value(relevant,v);
                    if(gh!=gu || gh!=gv || !permitted_content(gh,arithmetic)) {
                        ++out.gamma_rejected;
                        if(out.gamma_rejection_examples.size()<20)
                            out.gamma_rejection_examples.push_back({n,u,v,gh,gu,gv});
                        return;
                    }
                }
                ++out.gamma_pairs;
                if(critical_projections) {
                    u32 bad_prime=0,bad_z=0;
                    if(!critical_projection.accepts(u,h,&bad_prime,&bad_z)) {
                        ++out.critical_projection_rejected;
                        if(out.critical_projection_rejection_examples.size()<20)
                            out.critical_projection_rejection_examples.push_back({n,u,v,bad_prime,bad_z});
                        return;
                    }
                }
                ++out.critical_projection_pairs;
                parameter(N,2*u,2*v);
                evaluate(N,2*u,2*v,true,out.exact,[&](i64 x,i64 y,i64 a) {
                    accept(N,2*u,2*v,x,y,a,out.exact);
                });
            };
            if(opt.local && fieldwheel)
                out.crt_residues+=generate_field_gap_crt(n,h,arithmetic,*fieldwheel,visit_pair,gcd_nh);
            else out.crt_residues+=generate_gap_crt(n,h,arithmetic,visit_pair,gcd_nh);
        };
        if(opt.local) wheel.generate(n,visit_gap);
        else for(u32 h=0;h<=n;++h) visit_gap(h);
    }
    return out;
}

static void validate_gap(const ArithmeticTables &arithmetic,const LocalSieve &local,
                         const GapWheel &wheel) {
    u64 all_crt_pairs=0;
    for(u32 n=1;n<=256;++n) {
        std::vector<std::array<u32,2>> actual,expected;
        for(u32 h=0;h<=n;++h)
            generate_gap_crt(n,h,arithmetic,[&](u32 u){actual.push_back({u,u-h});});
        std::sort(actual.begin(),actual.end());
        for(u32 u=0;u<=n;++u) for(u32 v=0;v<=u;++v) {
            const u32 d=std::gcd(n/std::gcd(n,u),n/std::gcd(n,v));
            if((u-v)%d==0) expected.push_back({u,v});
        }
        if(actual!=expected) throw std::runtime_error("gap CRT set/uniqueness mismatch");
        all_crt_pairs+=actual.size();
    }
    using Parameter=std::array<u32,3>;
    std::set<Parameter> reference_parameters,gap_parameters;
    for(u32 n=12;n<=1000;n+=12) {
        const auto permitted=permitted_valuations(n,arithmetic.plans[n]);
        for(u32 u=0;u<=n;++u) {
            if(!permitted[u] || !local.rows_possible(n,u)) continue;
            const u32 g=std::gcd(n,u);
            generate_crt(n,u,arithmetic,[&](u32 v) {
                if(permitted[v] && permitted_difference(n,u,v,arithmetic) &&
                   permitted_scale(n,u,v,g) && local.accepts(n,u,v))
                    reference_parameters.insert({2*n,2*u,2*v});
            },g);
        }
    }
    SearchOptions filtered;filtered.bound=2000;
    run_gap_search(filtered,arithmetic,local,wheel,
        [&](u32 N,u32 p,u32 q){gap_parameters.insert({N,p,q});},
        [](i64,i64,i64,i64,i64,i64,Stats&){},false,false,false,nullptr,false);
    if(reference_parameters!=gap_parameters)
        throw std::runtime_error("gap filtered parameter set mismatch with u-first");

    using Key=std::array<i64,6>;
    std::set<Key> reference_points,gap_points;
    SearchOptions positive;positive.bound=1000;positive.local=false;positive.difference=false;
    run_search(positive,arithmetic,local,[&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats&){
        reference_points.insert({N,p,q,x,y,a});
    });
    run_gap_search(positive,arithmetic,local,wheel,[](u32,u32,u32){},
        [&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats&){gap_points.insert({N,p,q,x,y,a});},false,false,false,nullptr,false);
    if(reference_points.empty() || reference_points!=gap_points)
        throw std::runtime_error("gap primitive three-distance positive examples mismatch");
    std::cerr << "Gap validation passed: " << all_crt_pairs
              << " compatible pairs for every n=1..256; " << gap_parameters.size()
              << " fully filtered parameter tuples at N<=2000 exactly match u-first; "
              << gap_points.size() << " positive primitive three-distance examples retained.\n";
}

#ifndef RATIONAL_SQUARE_V3_NO_MAIN
int main(int argc,char**argv) try {
    SearchOptions opt;
    bool check=false;
    bool gap_valuation=false;
    bool gamma=true,critical=true,field_wheel19=true,critical_projections=true;
    std::string audit_path;
    for(int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if(arg=="--validate") check=true;
        else if(arg=="--no-gamma") gamma=false;
        else if(arg=="--no-critical-content") critical=false;
        else if(arg=="--field-wheel19") field_wheel19=true;
        else if(arg=="--no-field-wheel19") field_wheel19=false;
        else if(arg=="--no-critical-projections") critical_projections=false;
        else if(arg=="--no-local") opt.local=false;
        else if(arg=="--no-difference") opt.difference=false;
        else if(arg=="--no-valuation") opt.valuation=false;
        else if(arg=="--gap-valuation") gap_valuation=true;
        else if(arg=="--no-gap-valuation") gap_valuation=false;
        else if(arg=="--no-scale") opt.scale=false;
        else if(arg=="--min-denominator") {
            if(++i>=argc) throw std::runtime_error("missing minimum denominator");
            const u64 value=std::stoull(argv[i]);
            if(value<1 || value>2000000) throw std::runtime_error("minimum must be1..2000000");
            opt.min_denominator=u32(value);
        } else if(arg=="--audit-three") {
            if(++i>=argc) throw std::runtime_error("missing audit path");
            audit_path=argv[i];
        } else if(arg=="--help") {
            std::cout<<"all_plane_gap_v3 [B] [--min-denominator L] [--validate] "
                       "[--audit-three CSV] [--no-local] [--no-difference] "
                       "[--no-valuation] [--gap-valuation] [--no-scale] [--no-gamma] [--no-critical-content] [--no-field-wheel19] [--no-critical-projections]\n";
            return 0;
        } else {
            const u64 value=std::stoull(arg);
            if(value<1 || value>2000000) throw std::runtime_error("B must be1..2000000");
            opt.bound=u32(value);
        }
    }
    if(opt.min_denominator>opt.bound) throw std::runtime_error("minimum exceeds maximum");
    const auto started=std::chrono::steady_clock::now();
    ArithmeticTables arithmetic(std::max<u32>(1000,opt.bound/2));
    LocalSieve local;
    GapWheel wheel(local);
    ProgressionWheel<19> progression_wheel(local);
    const auto setup_done=std::chrono::steady_clock::now();
    if(check) {validate_gap(arithmetic,local,wheel);validate_field_wheel(arithmetic,local,progression_wheel);}
    std::ofstream audit;
    if(!audit_path.empty()) {
        audit.open(audit_path);
        if(!audit) throw std::runtime_error("cannot open audit file");
        audit<<"N,p,q,X,Y,A\n";
    }
    const auto search_started=std::chrono::steady_clock::now();
    const auto stats=run_gap_search(opt,arithmetic,local,wheel,[](u32,u32,u32){},
        [&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats &exact) {
            if(audit.is_open()) audit<<N<<','<<p<<','<<q<<','<<x<<','<<y<<','<<a<<'\n';
            u64 d;
            if(!square(u128(i128(N-x)*(N-x)+i128(N-y)*(N-y)),d)) return;
            ++exact.solutions;
            std::cout<<"{\"solution\":{\"N\":"<<N<<",\"X\":"<<x<<",\"Y\":"<<y
                     <<",\"A\":"<<a<<",\"D\":"<<d<<"}}\n";
        },gap_valuation,gamma,critical,field_wheel19?&progression_wheel:nullptr,critical_projections);
    const auto finished=std::chrono::steady_clock::now();
    auto seconds=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
    const auto &exact=stats.exact;
    std::cout<<"{\"scope\":\"whole_plane\",\"generator\":\"gap_CRT_mod1045_v3\""
             <<",\"max_minimal_coordinate_denominator\":"<<opt.bound
             <<",\"min_denominator\":"<<opt.min_denominator
             <<",\"local_filters\":"<<(opt.local?"true":"false")
             <<",\"difference_filters\":"<<(opt.difference?"true":"false")
             <<",\"valuation_filters\":"<<(opt.valuation?"true":"false")
             <<",\"gap_valuation_precheck\":"<<(gap_valuation?"true":"false")
             <<",\"scale_filter\":"<<(opt.scale?"true":"false")
             <<",\"gamma_filter\":"<<(gamma?"true":"false")
             <<",\"critical_5_13_content_filter\":"<<(critical?"true":"false")
             <<",\"field_wheel19\":"<<(field_wheel19 && opt.local?"true":"false")
             <<",\"critical_projections\":"<<(critical_projections?"true":"false")
             <<",\"sides\":"<<exact.sides
             <<",\"dense_even_pairs\":"<<stats.dense_even_pairs
             <<",\"gap_rows\":"<<stats.gap_rows
             <<",\"local_gap_rows\":"<<stats.local_gap_rows
             <<",\"valuation_gap_rows\":"<<stats.valuation_gap_rows
             <<",\"difference_gap_rows\":"<<stats.difference_gap_rows
             <<",\"critical_gap_rows\":"<<stats.critical_gap_rows
             <<",\"crt_residue_choices\":"<<stats.crt_residues
             <<",\"generated_pairs\":"<<stats.generated_pairs
             <<",\"valuation_pairs\":"<<stats.valuation_pairs
             <<",\"scale_pairs\":"<<stats.scale_pairs
             <<",\"local_pairs\":"<<stats.local_pairs
             <<",\"gamma_pairs\":"<<stats.gamma_pairs
             <<",\"gamma_rejected\":"<<stats.gamma_rejected
             <<",\"critical_projection_pairs\":"<<stats.critical_projection_pairs
             <<",\"critical_projection_rejected\":"<<stats.critical_projection_rejected
             <<",\"exact_parameter_pairs\":"<<exact.parameter_pairs
             <<",\"square_discriminants\":"<<exact.square_discriminants
             <<",\"linear_pairs\":"<<exact.linear_pairs
             <<",\"primitive_three_candidates\":"<<exact.primitive_three_candidates
             <<",\"solutions\":"<<exact.solutions
             <<",\"seconds_setup\":"<<seconds(started,setup_done)
             <<",\"seconds_search\":"<<seconds(search_started,finished)
             <<",\"seconds_total_excluding_validation\":"<<seconds(started,setup_done)+seconds(search_started,finished)
             <<",\"gamma_rejection_examples\":[";
    for(std::size_t i=0;i<stats.gamma_rejection_examples.size();++i) {
        if(i)std::cout<<',';
        const auto &e=stats.gamma_rejection_examples[i];
        std::cout<<"{\"n\":"<<e[0]<<",\"u\":"<<e[1]<<",\"v\":"<<e[2]
                 <<",\"gamma_h\":"<<e[3]<<",\"gamma_u\":"<<e[4]<<",\"gamma_v\":"<<e[5]<<'}';
    }
    std::cout<<"],\"critical_projection_rejection_examples\":[";
    for(std::size_t i=0;i<stats.critical_projection_rejection_examples.size();++i) {
        if(i)std::cout<<',';
        const auto &e=stats.critical_projection_rejection_examples[i];
        std::cout<<"{\"n\":"<<e[0]<<",\"u\":"<<e[1]<<",\"v\":"<<e[2]
                 <<",\"prime\":"<<e[3]<<",\"z\":"<<e[4]<<'}';
    }
    std::cout<<"]}\n";
    return 0;
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}

#endif // RATIONAL_SQUARE_V3_NO_MAIN
