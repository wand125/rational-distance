// Whole-plane search v5: the v4 generator and every v4 filter, unchanged,
// with the exact reconstruction widened for minimal coordinate denominators
// N <= 10^7 (see arithmetic_bounds_v5.md and wide_arith_v5.hpp).
// Frozen v4, v3 and the shared core are included, not edited.
#define RATIONAL_SQUARE_V4_NO_MAIN
#include "all_plane_gap_v4.cpp"
#undef RATIONAL_SQUARE_V4_NO_MAIN
#include <random>
#include "wide_arith_v5.hpp"

static constexpr u32 V5_BOUND=10000000;
static constexpr u32 FROZEN_BOUND=2000000; // safe range of the v2..v4 reconstruction

// Identical to run_v4 except that evaluate_v5 replaces evaluate and the
// point callback receives i128 coordinates and radius.
template<GammaStage stage,bool skip19,class AdditionalSieve,class ParameterCallback,class PointCallback>
static V4Stats run_v5(const SearchOptions &opt,const ArithmeticTables &arithmetic,
                     const LocalSieve &local,const GapWheel &gapwheel,
                     const ProgressionWheel<19> &wheel,AdditionalSieve &additional,
                     ParameterCallback parameter,PointCallback accept) {
    V4Stats out;
    const bool trim_endpoints=additional.trim_endpoints();
    const u32 first_n=12*((opt.min_denominator+23)/24);
    for(u32 n=first_n;n<=opt.bound/2;n+=12) {
        const u32 N=2*n;
        ++out.exact.sides;
        out.dense_even_pairs+=u64(n+1)*(n+2)/2;
        out.gap_rows+=n+1;
        auto permitted=permitted_valuations(n,arithmetic.plans[n]);
        const CriticalContentRows critical_rows(arithmetic.plans[n]);
        critical_rows.mark(permitted,n);
        additional.mark(n,permitted);
        const GammaPowerPlan gamma_plan(arithmetic.plans[n]);
        const GammaLabels labels(stage==GammaStage::table?n:0,gamma_plan);
        const CriticalProjectionPlan critical_projection(arithmetic.plans[n]);
        gapwheel.generate(n,[&](u32 h) {
            ++out.local_gap_rows;++out.valuation_gap_rows;
            if(trim_endpoints && h>=n-1) {++out.endpoint_empty_rows;return;}
            if(!permitted_difference(n,h,0,arithmetic)) return;
            ++out.difference_gap_rows;
            if(!critical_rows.accepts(h)) return;
            ++out.critical_gap_rows;
            const u32 gcd_nh=std::gcd(n,h);
            const GammaPowerRow gamma_row(h,gamma_plan,stage!=GammaStage::late && stage!=GammaStage::table);
            auto visit_pair=[&](u32 u) {
                const u32 v=u-h;
                ++out.generated_pairs;
                if(!permitted[u] || !permitted[v]) return;
                ++out.valuation_pairs;
                if((stage==GammaStage::table && !labels.accepts(u,v,h)) ||
                   (stage==GammaStage::rules && !gamma_row.accepts(u))) {
                    ++out.gamma_early_rejected;return;
                }
                ++out.gamma_early_passed;
                if(skip19 ? !local_without19(local,n,u,v) : !local.accepts(n,u,v)) return;
                ++out.local_pairs;
                if(!additional.accepts(n,u,v)) {++out.additional_rejected;return;}
                ++out.additional_pairs;
                if(!permitted_scale(n,u,v,gcd_nh)) return;
                ++out.scale_pairs;
                if(stage==GammaStage::late && gamma_plan.relevant>1) {
                    const u32 gh=gamma_value(gamma_plan.relevant,h);
                    if(gh!=gamma_value(gamma_plan.relevant,u) || gh!=gamma_value(gamma_plan.relevant,v)) {
                        ++out.gamma_rejected;return;
                    }
                }
                ++out.gamma_pairs;
                if(!critical_projection.accepts(u,h)) {++out.critical_projection_rejected;return;}
                ++out.critical_projection_pairs;
                parameter(N,2*u,2*v);
                evaluate_v5(N,2*u,2*v,true,out.exact,[&](i128 x,i128 y,i128 a) {
                    accept(N,2*u,2*v,x,y,a,out.exact);
                });
            };
            if((stage==GammaStage::generation || stage==GammaStage::lift) && gamma_row.count)
                out.crt_residues+=generate_gamma_field_crt(n,h,arithmetic,wheel,gamma_row,
                    stage==GammaStage::lift,visit_pair,out,gcd_nh,trim_endpoints);
            else if(trim_endpoints)
                out.crt_residues+=generate_field_gap_crt_range(n,h,arithmetic,wheel,visit_pair,gcd_nh,h+1,n-1);
            else out.crt_residues+=generate_field_gap_crt(n,h,arithmetic,wheel,visit_pair,gcd_nh);
        });
    }
    return out;
}

template<bool skip19,class AdditionalSieve,class ParameterCallback,class PointCallback>
static V4Stats dispatch_gamma_v5(const SearchOptions &opt,const ArithmeticTables &arithmetic,
                                const LocalSieve &local,const GapWheel &gapwheel,
                                const ProgressionWheel<19> &wheel,GammaStage stage,
                                AdditionalSieve &additional,ParameterCallback parameter,PointCallback accept) {
    switch(stage) {
        case GammaStage::late:return run_v5<GammaStage::late,skip19>(opt,arithmetic,local,gapwheel,wheel,additional,parameter,accept);
        case GammaStage::table:return run_v5<GammaStage::table,skip19>(opt,arithmetic,local,gapwheel,wheel,additional,parameter,accept);
        case GammaStage::rules:return run_v5<GammaStage::rules,skip19>(opt,arithmetic,local,gapwheel,wheel,additional,parameter,accept);
        case GammaStage::generation:return run_v5<GammaStage::generation,skip19>(opt,arithmetic,local,gapwheel,wheel,additional,parameter,accept);
        case GammaStage::lift:return run_v5<GammaStage::lift,skip19>(opt,arithmetic,local,gapwheel,wheel,additional,parameter,accept);
    }
    throw std::logic_error("unknown gamma stage");
}

template<class AdditionalSieve,class ParameterCallback,class PointCallback>
static V4Stats run_v5_dispatch(const SearchOptions &opt,const ArithmeticTables &arithmetic,
                              const LocalSieve &local,const GapWheel &gapwheel,
                              const ProgressionWheel<19> &wheel,GammaStage stage,bool skip19,
                              AdditionalSieve &additional,ParameterCallback parameter,PointCallback accept) {
    if(skip19) return dispatch_gamma_v5<true>(opt,arithmetic,local,gapwheel,wheel,stage,additional,parameter,accept);
    return dispatch_gamma_v5<false>(opt,arithmetic,local,gapwheel,wheel,stage,additional,parameter,accept);
}

// ---------------------------------------------------------------- validation
static bool stats_equal(const Stats &a,const Stats &b) {
    return a.sides==b.sides && a.parameter_pairs==b.parameter_pairs && a.linear_pairs==b.linear_pairs
        && a.square_discriminants==b.square_discriminants && a.integral_A==b.integral_A
        && a.integer_coordinates==b.integer_coordinates
        && a.primitive_three_candidates==b.primitive_three_candidates && a.solutions==b.solutions;
}

static void validate_wide_arithmetic() {
    std::mt19937_64 rng(20261105);
    u64 checks=0;
    auto check_root=[&](u128 r) {
        const U256 s=mul128(r,r);u128 root=0;
        if(!square256(s,root) || root!=r || isqrt256(s)!=r) throw std::runtime_error("wide square of root failed");
        if(r>0) {
            const U256 below=sub256(s,U256::from(1));
            if(isqrt256(below)!=r-1 || (r>1 && square256(below,root))) throw std::runtime_error("wide r^2-1 failed");
        }
        const U256 above=add256(s,U256::from(u128(2)*r)); // r^2+2r < (r+1)^2
        if(isqrt256(above)!=r) throw std::runtime_error("wide r^2+2r root failed");
        if(r>0 && square256(above,root)) throw std::runtime_error("wide r^2+2r square failed");
        checks+=4;
    };
    for(u128 r=0;r<5000;++r) check_root(r);
    for(unsigned bits=1;bits<=127;++bits) for(int t=0;t<200;++t) {
        u128 r=(u128(rng())<<64)|rng();
        r&=(bits==128?~u128(0):((u128(1)<<bits)-1));
        check_root(r);
        check_root((u128(1)<<(bits-1)));
        check_root((u128(1)<<bits)-1);
    }
    // Product check against the u128 product where it fits, and against
    // the identity (a*b)*(c) computed two ways for the K factorisation.
    for(int t=0;t<200000;++t) {
        const u128 a=rng(),b=rng();
        if(mul128(a,b)!=U256::from(a*b)) throw std::runtime_error("mul128 small product");
        const u128 c=(u128(rng())<<40)|(rng()>>24),d=rng()>>16;
        const U256 left=mul128(c,d);
        U256 right=add256(mul128(u128(u64(c)),d),U256{});
        const U256 high=mul128(c>>64,d); // shifted by 64 bits
        U256 shifted;shifted.w[1]=high.w[0];shifted.w[2]=high.w[1];shifted.w[3]=high.w[2];
        if(high.w[3]) throw std::runtime_error("unexpected high limb");
        right=add256(right,shifted);
        if(left!=right) throw std::runtime_error("mul128 limb decomposition");
        checks+=2;
    }
    // i128 printing round trip.
    for(int t=0;t<20000;++t) {
        const i128 x=i128((u128(rng())<<60)^rng())*((t&1)?-1:1);
        if(std::stoll(i128_string(i128(i64(x>>70))))!=i64(x>>70)) throw std::runtime_error("i128 print");
    }
    std::cerr<<"V5 wide arithmetic validation passed: "<<checks<<" exact checks.\n";
}

static void validate_v5_reconstruction() {
    // (a) Independent brute force over the proved complete radius, N=1..12,
    // including exterior points and without four-distance rules.
    u64 total=0;
    for(i64 N=1;N<=12;++N) {
        Points algebraic,brute;Stats ignored;
        for(i64 p=0;p<=N;++p) for(i64 q=0;q<=p;++q)
            evaluate_v5(N,p,q,false,ignored,[&](i128 x,i128 y,i128){algebraic.emplace(i64(x),i64(y));});
        const i64 radius=2*N*N*N;
        for(i64 y=-radius;y<=N/2;++y) for(i64 x=-radius;x<=y;++x) {
            u64 a,b,c;
            if(!square(u128(i128(x)*x+i128(y)*y),a)) continue;
            if(!square(u128(i128(N-x)*(N-x)+i128(y)*y),b)) continue;
            if(!square(u128(i128(x)*x+i128(N-y)*(N-y)),c)) continue;
            brute.emplace(x,y);
        }
        if(algebraic!=brute) throw std::runtime_error("v5 brute validation mismatch at N="+std::to_string(N));
        total+=algebraic.size();
    }
    // (b) Equality with the frozen evaluate() inside its proved range,
    // in both modes: every pair for N<=400, random pairs up to N=2,000,000.
    using Key=std::array<i64,3>;
    u64 compared=0,points=0;
    auto compare=[&](i64 N,i64 p,i64 q) {
        for(bool four:{false,true}) {
            Stats s_old,s_new;std::vector<Key> old_points,new_points;
            evaluate(N,p,q,four,s_old,[&](i64 x,i64 y,i64 a){old_points.push_back({x,y,a});});
            evaluate_v5(N,p,q,four,s_new,[&](i128 x,i128 y,i128 a){new_points.push_back({i64(x),i64(y),i64(a)});});
            if(old_points!=new_points || !stats_equal(s_old,s_new))
                throw std::runtime_error("v5 differs from frozen evaluate at N="+std::to_string(N)
                                         +" p="+std::to_string(p)+" q="+std::to_string(q));
            points+=new_points.size();
        }
        ++compared;
    };
    for(i64 N=1;N<=400;++N) for(i64 p=0;p<=N;++p) for(i64 q=0;q<=p;++q) compare(N,p,q);
    std::mt19937_64 rng(5);
    for(int t=0;t<3000000;++t) {
        const i64 N=1+i64(rng()%FROZEN_BOUND),p=i64(rng()%u64(N+1)),q=i64(rng()%u64(p+1));
        compare(N,p,q);
    }
    // Known nontrivial three-distance parameters, then the same points scaled
    // close to N=10^7, where K, sqrt(K), A and the norms leave 64/128 bits.
    std::vector<std::array<i64,3>> seeds;
    for(i64 N=1;N<=12;++N) for(i64 p=0;p<=N;++p) for(i64 q=0;q<=p;++q) {
        Stats ignored;
        evaluate(N,p,q,false,ignored,[&](i64,i64,i64){seeds.push_back({N,p,q});});
    }
    seeds.push_back({8400,6720,5040});
    u64 scaled=0;
    for(const auto &s:seeds) {
        const i64 N=s[0],p=s[1],q=s[2];
        Stats st;std::vector<std::array<i128,3>> base;
        evaluate_v5(N,p,q,false,st,[&](i128 x,i128 y,i128 a){base.push_back({x,y,a});});
        for(i64 k:{i64(1),i64(V5_BOUND/N),i64(V5_BOUND/N)-1,i64(V5_BOUND/N/2)+1}) {
            if(k<1 || k*N>V5_BOUND) continue;
            std::vector<std::array<i128,3>> got;
            evaluate_v5(k*N,k*p,k*q,false,st,[&](i128 x,i128 y,i128 a){got.push_back({x,y,a});});
            std::vector<std::array<i128,3>> expected;
            for(auto b:base) expected.push_back({k*b[0],k*b[1],k*b[2]});
            std::sort(got.begin(),got.end());std::sort(expected.begin(),expected.end());
            if(got!=expected) throw std::runtime_error("scaled three-distance point lost at N="+std::to_string(k*N));
            ++scaled;
        }
    }
    // Extreme parameters at the top of the range must not overflow.
    {
        Stats st;const i64 N=V5_BOUND;
        for(i64 p:{i64(0),i64(1),N/2,N-1,N}) for(i64 q:{i64(0),i64(1),N/3,N-1,N}) if(q<=p)
            evaluate_v5(N,p,q,false,st,[&](i128,i128,i128){});
    }
    std::cerr<<"V5 reconstruction validation passed: brute N=1..12 ("<<total<<" points), "
             <<compared<<" parameter triples equal to the frozen core ("<<points<<" points), "
             <<scaled<<" scaled seed families up to N=10^7.\n";
}

static void write_wide_selftest(const std::string &path) {
    // Random cases for an independent Python big-integer check.
    std::ofstream out(path);if(!out) throw std::runtime_error("selftest open");
    std::mt19937_64 rng(77);
    out<<"a,b,product,isqrt_product_plus_c,c\n";
    for(int t=0;t<20000;++t) {
        const unsigned bits_a=1+rng()%127,bits_b=1+rng()%127;
        u128 a=((u128(rng())<<64)|rng())&((u128(1)<<bits_a)-1);
        u128 b=((u128(rng())<<64)|rng())&((u128(1)<<bits_b)-1);
        const U256 product=mul128(a,b);
        const u128 c=rng()%3==0?0:rng();
        const U256 z=add256(product,U256::from(c));
        auto u256_string=[](U256 v) {
            std::string s;
            while(v.w[0]||v.w[1]||v.w[2]||v.w[3]) {
                u128 rem=0;
                for(int i=3;i>=0;--i) {const u128 cur=(rem<<64)|v.w[i];v.w[i]=u64(cur/10);rem=cur%10;}
                s.push_back(char('0'+unsigned(rem)));
            }
            if(s.empty()) s="0";
            return std::string(s.rbegin(),s.rend());
        };
        out<<i128_string(i128(a>>1))<<"*2+"<<unsigned(a&1)<<','<<i128_string(i128(b>>1))<<"*2+"<<unsigned(b&1)<<','
           <<u256_string(product)<<','<<i128_string(i128(isqrt256(z)>>1))<<"*2+"<<unsigned(isqrt256(z)&1)<<','
           <<i128_string(i128(c))<<'\n';
    }
}

#ifndef RATIONAL_SQUARE_V5_NO_MAIN
int main(int argc,char**argv) try {
    SearchOptions opt;
    GammaStage stage=GammaStage::late;
    bool skip19=true,check=false,new_filters=true;
    AdditionalOptions additional_options;
    std::string audit_path,parameter_path,selftest_path;
    for(int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        auto value=[&]()->std::string {if(++i>=argc) throw std::runtime_error("missing option value");return argv[i];};
        if(arg=="--validate") check=true;
        else if(arg=="--no-new-filters") new_filters=false;
        else if(arg=="--no-five-adic") additional_options.fifth=FifthStage::off;
        else if(arg=="--five-adic-late") additional_options.fifth=FifthStage::late;
        else if(arg=="--slow-five-adic-mark") additional_options.fast_fifth=false;
        else if(arg=="--no-prime-power-projections") additional_options.lifted=false;
        else if(arg=="--no-endpoint-pruning") additional_options.endpoints=false;
        else if(arg=="--no-critical-extension") additional_options.critical=false;
        else if(arg=="--keep-redundant19") skip19=false;
        else if(arg=="--gamma-stage") {
            const auto name=value();
            if(name=="late") stage=GammaStage::late;
            else if(name=="table") stage=GammaStage::table;
            else if(name=="rules") stage=GammaStage::rules;
            else if(name=="generation") stage=GammaStage::generation;
            else if(name=="lift") stage=GammaStage::lift;
            else throw std::runtime_error("gamma stage: late|table|rules|generation|lift");
        } else if(arg=="--audit-three") audit_path=value();
        else if(arg=="--audit-parameters") parameter_path=value();
        else if(arg=="--wide-selftest") selftest_path=value();
        else if(arg=="--min-denominator") {
            const u64 lower=std::stoull(value());
            if(lower<1 || lower>V5_BOUND) throw std::runtime_error("minimum1..10000000");
            opt.min_denominator=u32(lower);
        } else if(arg=="--help") {
            std::cout<<"all_plane_gap_v5 [B<=10000000] [--min-denominator L] [--validate] [--audit-three CSV] "
                       "[--audit-parameters CSV] [--wide-selftest CSV] [--keep-redundant19] [--no-new-filters] "
                       "[--no-five-adic] [--five-adic-late] [--slow-five-adic-mark] "
                       "[--no-prime-power-projections] [--no-endpoint-pruning] [--no-critical-extension] "
                       "[--gamma-stage late|table|rules|generation|lift]\n";return 0;
        } else {
            const u64 bound=std::stoull(arg);
            if(bound<1 || bound>V5_BOUND) throw std::runtime_error("bound1..10000000");
            opt.bound=u32(bound);
        }
    }
    if(!selftest_path.empty()) {write_wide_selftest(selftest_path);return 0;}
    if(opt.min_denominator>opt.bound) throw std::runtime_error("minimum exceeds bound");
    const auto started=std::chrono::steady_clock::now();
    ArithmeticTables arithmetic(std::max<u32>(1000,opt.bound/2));LocalSieve local;
    GapWheel gapwheel(local);ProgressionWheel<19> wheel(local);
    std::unique_ptr<AdditionalSieve> additional;
    if(new_filters) additional=std::make_unique<AdditionalSieve>(additional_options);
    const auto setup_done=std::chrono::steady_clock::now();
    if(check) {
        validate_gap(arithmetic,local,gapwheel);
        validate_field_wheel(arithmetic,local,wheel);
        validate_v4_generators(local,wheel);
        if(additional) additional->verify();
        validate_wide_arithmetic();
        validate_v5_reconstruction();
    }
    std::ofstream audit,parameters;
    if(!audit_path.empty()) {audit.open(audit_path);if(!audit)throw std::runtime_error("audit open");audit<<"N,p,q,X,Y,A\n";}
    if(!parameter_path.empty()) {parameters.open(parameter_path);if(!parameters)throw std::runtime_error("parameter audit open");parameters<<"N,p,q\n";}
    auto accept=[&](i64 N,i64 p,i64 q,i128 x,i128 y,i128 a,Stats &stats) {
        if(audit.is_open())audit<<N<<','<<p<<','<<q<<','<<i128_string(x)<<','<<i128_string(y)<<','<<i128_string(a)<<'\n';
        const U256 d2=add256(square_i128(N-x),square_i128(N-y));
        u128 d;
        if(square256(d2,d)) {
            ++stats.solutions;
            std::cout<<"{\"solution\":{\"N\":"<<N<<",\"X\":"<<i128_string(x)<<",\"Y\":"<<i128_string(y)
                     <<",\"A\":"<<i128_string(a)<<",\"D\":"<<i128_string(i128(d))<<"}}\n"<<std::flush;
        }
    };
    auto parameter=[&](u32 N,u32 p,u32 q) {if(parameters.is_open())parameters<<N<<','<<p<<','<<q<<'\n';};
    const auto search_started=std::chrono::steady_clock::now();
    V4Stats stats;
    if(additional) stats=run_v5_dispatch(opt,arithmetic,local,gapwheel,wheel,stage,skip19,*additional,parameter,accept);
    else {NoAdditionalSieve none;stats=run_v5_dispatch(opt,arithmetic,local,gapwheel,wheel,stage,skip19,none,parameter,accept);}
    const auto finished=std::chrono::steady_clock::now();
    if(audit.is_open()) {audit.flush();if(!audit) throw std::runtime_error("audit write failed");}
    if(parameters.is_open()) {parameters.flush();if(!parameters) throw std::runtime_error("parameter write failed");}
    auto seconds=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
    std::cout<<"{\"scope\":\"whole_plane\",\"generator\":\"gap_CRT_v5\""
             <<",\"reconstruction\":\"exact_u256_v5\""
             <<",\"max_minimal_coordinate_denominator\":"<<opt.bound
             <<",\"min_denominator\":"<<opt.min_denominator
             <<",\"reference_v3\":false"
             <<",\"gamma_stage\":\""<<gamma_stage_name(stage)<<"\""
             <<",\"skip_redundant19\":"<<(skip19?"true":"false")
             <<",\"additional_filters\":"<<(additional?"true":"false")
             <<",\"endpoint_pruning\":"<<(additional && additional_options.endpoints?"true":"false")
             <<",\"fifth_power_condition\":\""<<(!additional || additional_options.fifth==FifthStage::off?"off":
                                                      additional_options.fifth==FifthStage::early?"early":"late")<<"\""
             <<",\"prime_power_projections\":"<<(additional && additional_options.lifted?"true":"false")
             <<",\"critical_extension\":"<<(additional && additional_options.critical?"true":"false")
             <<",\"fast_fifth_mark\":"<<(additional && additional_options.fifth==FifthStage::early && additional_options.fast_fifth?"true":"false")
             <<",\"sides\":"<<stats.exact.sides
             <<",\"dense_even_pairs\":"<<stats.dense_even_pairs
             <<",\"gap_rows\":"<<stats.gap_rows
             <<",\"local_gap_rows\":"<<stats.local_gap_rows
             <<",\"difference_gap_rows\":"<<stats.difference_gap_rows
             <<",\"critical_gap_rows\":"<<stats.critical_gap_rows
             <<",\"endpoint_empty_rows\":"<<stats.endpoint_empty_rows
             <<",\"crt_residue_choices\":"<<stats.crt_residues
             <<",\"generated_pairs\":"<<stats.generated_pairs
             <<",\"valuation_pairs\":"<<stats.valuation_pairs
             <<",\"gamma_early_rejected\":"<<stats.gamma_early_rejected
             <<",\"gamma_lifted_progressions\":"<<stats.gamma_lifted_progressions
             <<",\"gamma_lift_residues\":"<<stats.gamma_lift_residues
             <<",\"examined_generation_rejections\":"<<stats.examined_generation_rejections
             <<",\"local_pairs\":"<<stats.local_pairs
             <<",\"additional_pairs\":"<<stats.additional_pairs
             <<",\"additional_rejected\":"<<stats.additional_rejected
             <<",\"scale_pairs\":"<<stats.scale_pairs
             <<",\"gamma_pairs\":"<<stats.gamma_pairs
             <<",\"gamma_rejected\":"<<stats.gamma_rejected
             <<",\"critical_projection_pairs\":"<<stats.critical_projection_pairs
             <<",\"critical_projection_rejected\":"<<stats.critical_projection_rejected
             <<",\"exact_parameter_pairs\":"<<stats.exact.parameter_pairs
             <<",\"square_discriminants\":"<<stats.exact.square_discriminants
             <<",\"linear_pairs\":"<<stats.exact.linear_pairs
             <<",\"primitive_three_candidates\":"<<stats.exact.primitive_three_candidates
             <<",\"solutions\":"<<stats.exact.solutions
             <<",\"seconds_setup\":"<<seconds(started,setup_done)
             <<",\"seconds_search\":"<<seconds(search_started,finished)
             <<",\"seconds_total_excluding_validation\":"<<seconds(started,setup_done)+seconds(search_started,finished)<<"}\n";
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
#endif
