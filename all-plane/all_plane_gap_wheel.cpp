// Experimental progression wheel; v2 prefix is copied unchanged.
// Gap-first CRT search, version 2: inexpensive field tests precede gcd work.
// Build: g++ -O3 -std=c++17 all_plane_gap_v2.cpp -o all_plane_gap_v2
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

struct GapStats {
    u64 dense_even_pairs=0,gap_rows=0,local_gap_rows=0,valuation_gap_rows=0,difference_gap_rows=0;
    u64 crt_residues=0,generated_pairs=0,valuation_pairs=0,scale_pairs=0,local_pairs=0;
    Stats exact;
};

template<class ParameterCallback,class PointCallback>
static GapStats run_gap_search(const SearchOptions &opt,const ArithmeticTables &arithmetic,
                               const LocalSieve &local,const GapWheel &wheel,
                               ParameterCallback parameter,PointCallback accept,
                               bool gap_valuation=false) {
    GapStats out;
    const u32 first_n=12*((opt.min_denominator+23)/24);
    for(u32 n=first_n;n<=opt.bound/2;n+=12) {
        const u32 N=2*n;
        ++out.exact.sides;
        out.dense_even_pairs+=u64(n+1)*(n+2)/2;
        out.gap_rows+=n+1;
        auto permitted=permitted_valuations(n,arithmetic.plans[n]);
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
            const u32 gcd_nh=std::gcd(n,h);
            out.crt_residues+=generate_gap_crt(n,h,arithmetic,[&](u32 u) {
                const u32 v=u-h;
                ++out.generated_pairs;
                if(opt.valuation && (!permitted[u] || !permitted[v])) return;
                ++out.valuation_pairs;
                if(opt.local && !local.accepts(n,u,v)) return;
                ++out.local_pairs;
                // gcd(gcd(n,h),v)=gcd(n,u,v), so the cached gcd_nh can
                // safely replace the predicate's usual gcd(n,u) argument.
                // Only field-compatible candidates pay for these gcds.
                if(opt.scale && !permitted_scale(n,u,v,gcd_nh)) return;
                ++out.scale_pairs;
                parameter(N,2*u,2*v);
                evaluate(N,2*u,2*v,true,out.exact,[&](i64 x,i64 y,i64 a) {
                    accept(N,2*u,2*v,x,y,a,out.exact);
                });
            },gcd_nh);
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
        [](i64,i64,i64,i64,i64,i64,Stats&){});
    if(reference_parameters!=gap_parameters)
        throw std::runtime_error("gap filtered parameter set mismatch with u-first");

    using Key=std::array<i64,6>;
    std::set<Key> reference_points,gap_points;
    SearchOptions positive;positive.bound=1000;positive.local=false;positive.difference=false;
    run_search(positive,arithmetic,local,[&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats&){
        reference_points.insert({N,p,q,x,y,a});
    });
    run_gap_search(positive,arithmetic,local,wheel,[](u32,u32,u32){},
        [&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats&){gap_points.insert({N,p,q,x,y,a});});
    if(reference_points.empty() || reference_points!=gap_points)
        throw std::runtime_error("gap primitive three-distance positive examples mismatch");
    std::cerr << "Gap validation passed: " << all_crt_pairs
              << " compatible pairs for every n=1..256; " << gap_parameters.size()
              << " fully filtered parameter tuples at N<=2000 exactly match u-first; "
              << gap_points.size() << " positive primitive three-distance examples retained.\n";
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

template<unsigned W,class PointCallback>
static GapStats run_field_gap(const SearchOptions &opt,const ArithmeticTables &arithmetic,
                             const LocalSieve &local,const GapWheel &gapwheel,
                             const ProgressionWheel<W> &wheel,PointCallback accept) {
    GapStats out;
    const u32 first_n=12*((opt.min_denominator+23)/24);
    for(u32 n=first_n;n<=opt.bound/2;n+=12) {
        const u32 N=2*n;
        ++out.exact.sides;
        out.dense_even_pairs+=u64(n+1)*(n+2)/2;
        out.gap_rows+=n+1;
        auto permitted=permitted_valuations(n,arithmetic.plans[n]);
        gapwheel.generate(n,[&](u32 h){
            ++out.local_gap_rows;
            if(!permitted_difference(n,h,0,arithmetic)) return;
            ++out.difference_gap_rows;
            const u32 gcd_nh=std::gcd(n,h);
            out.crt_residues+=generate_field_gap_crt(n,h,arithmetic,wheel,[&](u32 u) {
                const u32 v=u-h;
                ++out.generated_pairs;
                if(!permitted[u] || !permitted[v]) return;
                ++out.valuation_pairs;
                if(!local.accepts(n,u,v)) return;
                ++out.local_pairs;
                if(!permitted_scale(n,u,v,gcd_nh)) return;
                ++out.scale_pairs;
                evaluate(N,2*u,2*v,true,out.exact,[&](i64 x,i64 y,i64 a){
                    accept(N,2*u,2*v,x,y,a,out.exact);
                });
            },gcd_nh);
        });
    }
    return out;
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

int main(int argc,char**argv) try {
    SearchOptions opt;unsigned mode=55;bool check=false;std::string auditpath;
    for(int i=1;i<argc;++i) {
        std::string arg=argv[i];
        if(arg=="--wheel") {if(++i==argc) throw std::runtime_error("missing wheel");mode=std::stoul(argv[i]);}
        else if(arg=="--validate") check=true;
        else if(arg=="--audit-three") {if(++i==argc) throw std::runtime_error("missing audit");auditpath=argv[i];}
        else {opt.bound=std::stoul(arg);}
    }
    if(opt.bound<1 || opt.bound>2000000) throw std::runtime_error("bound1..2000000");
    const auto start=std::chrono::steady_clock::now();
    ArithmeticTables arithmetic(std::max<u32>(1000,opt.bound/2));LocalSieve local;GapWheel gapwheel(local);
    std::ofstream audit;
    if(!auditpath.empty()){audit.open(auditpath);if(!audit) throw std::runtime_error("audit open");audit<<"N,p,q,X,Y,A\n";}
    auto point=[&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats &stats){
        if(audit.is_open()) audit<<N<<','<<p<<','<<q<<','<<x<<','<<y<<','<<a<<'\n';
        u64 d;if(square(u128(i128(N-x)*(N-x)+i128(N-y)*(N-y)),d)) {
            ++stats.solutions;
            std::cout<<"{\"solution\":{\"N\":"<<N<<",\"X\":"<<x<<",\"Y\":"<<y<<",\"A\":"<<a<<",\"D\":"<<d<<"}}\n";
        }
    };
    GapStats stats;auto search_start=start;
    auto run=[&](auto &wheel){
        if(check) validate_field_wheel(arithmetic,local,wheel);
        search_start=std::chrono::steady_clock::now();
        stats=run_field_gap(opt,arithmetic,local,gapwheel,wheel,point);
    };
    if(mode==0){search_start=std::chrono::steady_clock::now();stats=run_gap_search(opt,arithmetic,local,gapwheel,[](u32,u32,u32){},point);}
    else if(mode==11){ProgressionWheel<11> wheel(local);run(wheel);}
    else if(mode==19){ProgressionWheel<19> wheel(local);run(wheel);}
    else if(mode==55){ProgressionWheel<55> wheel(local);run(wheel);}
    else throw std::runtime_error("wheel must be0,11,19,55");
    const auto end=std::chrono::steady_clock::now();
    std::cout<<"{\"bound\":"<<opt.bound<<",\"wheel_modulus\":"<<mode
             <<",\"generated_pairs\":"<<stats.generated_pairs<<",\"valuation_pairs\":"<<stats.valuation_pairs
             <<",\"local_pairs\":"<<stats.local_pairs<<",\"exact_parameter_pairs\":"<<stats.exact.parameter_pairs
             <<",\"primitive_three_candidates\":"<<stats.exact.primitive_three_candidates
             <<",\"solutions\":"<<stats.exact.solutions
             <<",\"seconds_search\":"<<std::chrono::duration<double>(end-search_start).count()
             <<",\"seconds_total\":"<<std::chrono::duration<double>(end-start).count()<<"}\n";
} catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
