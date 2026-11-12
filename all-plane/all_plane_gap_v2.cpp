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

int main(int argc,char**argv) try {
    SearchOptions opt;
    bool check=false;
    bool gap_valuation=false;
    std::string audit_path;
    for(int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if(arg=="--validate") check=true;
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
            std::cout<<"all_plane_gap_v2 [B] [--min-denominator L] [--validate] "
                       "[--audit-three CSV] [--no-local] [--no-difference] "
                       "[--no-valuation] [--gap-valuation] [--no-scale]\n";
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
    const auto setup_done=std::chrono::steady_clock::now();
    if(check) validate_gap(arithmetic,local,wheel);
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
        },gap_valuation);
    const auto finished=std::chrono::steady_clock::now();
    auto seconds=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
    const auto &exact=stats.exact;
    std::cout<<"{\"scope\":\"whole_plane\",\"generator\":\"gap_CRT_mod1045_v2\""
             <<",\"max_minimal_coordinate_denominator\":"<<opt.bound
             <<",\"min_denominator\":"<<opt.min_denominator
             <<",\"local_filters\":"<<(opt.local?"true":"false")
             <<",\"difference_filters\":"<<(opt.difference?"true":"false")
             <<",\"valuation_filters\":"<<(opt.valuation?"true":"false")
             <<",\"gap_valuation_precheck\":"<<(gap_valuation?"true":"false")
             <<",\"scale_filter\":"<<(opt.scale?"true":"false")
             <<",\"sides\":"<<exact.sides
             <<",\"dense_even_pairs\":"<<stats.dense_even_pairs
             <<",\"gap_rows\":"<<stats.gap_rows
             <<",\"local_gap_rows\":"<<stats.local_gap_rows
             <<",\"valuation_gap_rows\":"<<stats.valuation_gap_rows
             <<",\"difference_gap_rows\":"<<stats.difference_gap_rows
             <<",\"crt_residue_choices\":"<<stats.crt_residues
             <<",\"generated_pairs\":"<<stats.generated_pairs
             <<",\"valuation_pairs\":"<<stats.valuation_pairs
             <<",\"scale_pairs\":"<<stats.scale_pairs
             <<",\"local_pairs\":"<<stats.local_pairs
             <<",\"exact_parameter_pairs\":"<<exact.parameter_pairs
             <<",\"square_discriminants\":"<<exact.square_discriminants
             <<",\"linear_pairs\":"<<exact.linear_pairs
             <<",\"primitive_three_candidates\":"<<exact.primitive_three_candidates
             <<",\"solutions\":"<<exact.solutions
             <<",\"seconds_setup\":"<<seconds(started,setup_done)
             <<",\"seconds_search\":"<<seconds(search_started,finished)
             <<",\"seconds_total_excluding_validation\":"<<seconds(started,setup_done)+seconds(search_started,finished)
             <<"}\n";
    return 0;
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
