#define RATIONAL_SQUARE_V4_NO_MAIN
#include "../all_plane_gap_v4.cpp"

// Independent 5-adic oracle: inspect the complete exact discriminant product,
// not either one-variable marking routine. The n-multiple-of-5 case is outside
// the additional condition's proved scope and is deliberately retained.
static bool full_five_square(u32 n,u32 u,u32 v) {
    if(n%5==0) return true;
    const u64 nn=u64(n)*n,uu=u64(u)*u,vv=u64(v)*v,hh=u64(u-v)*(u-v);
    u128 value=u128(nn-uu)*(nn-vv)*(2*nn-hh);
    if(!value) return true;
    unsigned exponent=0;
    while(value%5==0) {value/=5;++exponent;}
    return exponent%2==0 && (value%5==1 || value%5==4);
}

int main(int argc,char **argv) try {
    const u32 bound=argc>1?u32(std::stoul(argv[1])):100000;
    if(bound<24 || bound>200000) throw std::runtime_error("validation bound24..200000");
    ArithmeticTables arithmetic(std::max<u32>(1000,bound/2));LocalSieve local;
    GapWheel gaps(local);ProgressionWheel<19> wheel(local);
    AdditionalOptions options;AdditionalSieve extra(options);
    extra.verify();
    using Key=std::array<u32,3>;
    using Point=std::array<i64,6>;
    std::set<Key> baseline,expected,actual,late_actual,optimization_only;
    std::set<Point> expected_points,actual_points;
    u64 endpoint_rejected=0,five_rejected=0,lift_rejected=0,critical_rejected=0,five_oracle_comparisons=0;
    u64 sequential_endpoint=0,sequential_five=0,sequential_lift=0,sequential_critical=0,old_points=0;
    CriticalExtension critical;u32 last_n=0;
    struct Witness {u32 n,u,v;bool endpoints,five,lifted;};
    std::vector<Witness> critical_witnesses;
    SearchOptions search;search.bound=bound;
    run_gap_search(search,arithmetic,local,gaps,[&](u32 N,u32 p,u32 q) {
        if(!baseline.insert({N,p,q}).second) throw std::runtime_error("duplicate frozen-v3 parameter");
        const u32 n=N/2,u=p/2,v=q/2;
        const bool endpoints=(v>0 && u<n),five=full_five_square(n,u,v);
        if(five!=FifthPowerCondition::accepts(n,u,v))
            throw std::runtime_error("one-variable fifth marking differs from full discriminant oracle");
        ++five_oracle_comparisons;
        const bool lifted=extra.p125.accepts(n,u,v) && extra.p169.accepts(n,u,v) &&
                          extra.p121.accepts(n,u,v) && extra.p49.accepts(n,u,v);
        if(n!=last_n) {critical.prepare(n);last_n=n;}
        const bool critical_ok=critical.accepts(u,u-v);
        endpoint_rejected+=!endpoints;five_rejected+=!five;lift_rejected+=!lifted;critical_rejected+=!critical_ok;
        if(!critical_ok)critical_witnesses.push_back({n,u,v,endpoints,five,lifted});
        if(!endpoints) {++sequential_endpoint;return;}
        if(!five) {++sequential_five;return;}
        if(!lifted) {++sequential_lift;return;}
        if(!critical_ok) {++sequential_critical;return;}
        expected.insert({N,p,q});
    },[&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats &) {
        ++old_points;
        if(expected.count({u32(N),u32(p),u32(q)})) expected_points.insert({N,p,q,x,y,a});
    },false,true,true,&wheel,true);

    auto no_point=[](i64,i64,i64,i64,i64,i64,Stats&){};
    NoAdditionalSieve none;
    run_v4_dispatch(search,arithmetic,local,gaps,wheel,GammaStage::late,true,none,
        [&](u32 N,u32 p,u32 q) {
            if(!optimization_only.insert({N,p,q}).second) throw std::runtime_error("duplicate unextended-v4 parameter");
        },no_point);
    if(optimization_only!=baseline) throw std::runtime_error("unextended v4 differs from frozen v3");
    const auto refined=run_v4_dispatch(search,arithmetic,local,gaps,wheel,GammaStage::late,true,extra,
        [&](u32 N,u32 p,u32 q) {
            if(!actual.insert({N,p,q}).second) throw std::runtime_error("duplicate extended-v4 parameter");
        },[&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats&){actual_points.insert({N,p,q,x,y,a});});
    options.fifth=FifthStage::late;AdditionalSieve late(options);
    run_v4_dispatch(search,arithmetic,local,gaps,wheel,GammaStage::late,true,late,
        [&](u32 N,u32 p,u32 q) {
            if(!late_actual.insert({N,p,q}).second) throw std::runtime_error("duplicate late-fifth parameter");
        },no_point);
    if(expected.empty() || expected!=actual || expected!=late_actual)
        throw std::runtime_error("optimized additional-filter set differs from frozen-v3 late oracle");
    if(expected_points!=actual_points) throw std::runtime_error("expected three-distance diagnostic set mismatch");
    std::cout<<"{\"bound\":"<<bound<<",\"projection_tables_independently_verified\":true"
             <<",\"unextended_v4_equals_frozen_v3\":true,\"nonempty_exact_parameter_sets_equal\":true"
             <<",\"early_and_late_fifth_sets_equal\":true,\"duplicates\":0"
             <<",\"fast_fifth_mark_equals_exact_definition\":true"
             <<",\"baseline_parameters\":"<<baseline.size()<<",\"extended_parameters\":"<<actual.size()
             <<",\"full_discriminant_fifth_oracle_comparisons\":"<<five_oracle_comparisons
             <<",\"independent_rejected\":{\"endpoints\":"<<endpoint_rejected<<",\"fifth\":"<<five_rejected
             <<",\"prime_power_projections\":"<<lift_rejected<<",\"critical_extension\":"<<critical_rejected<<"}"
             <<",\"sequential_rejected\":{\"endpoints\":"<<sequential_endpoint<<",\"fifth\":"<<sequential_five
             <<",\"prime_power_projections\":"<<sequential_lift<<",\"critical_extension\":"<<sequential_critical<<"}"
             <<",\"old_three_distance_points\":"<<old_points<<",\"retained_three_distance_points\":"<<actual_points.size()
             <<",\"generated_pairs\":"<<refined.generated_pairs<<",\"valuation_pairs\":"<<refined.valuation_pairs
             <<",\"critical_witnesses\":[";
    for(std::size_t j=0;j<critical_witnesses.size();++j) {
        if(j)std::cout<<',';const auto &w=critical_witnesses[j];
        std::cout<<"{\"n\":"<<w.n<<",\"u\":"<<w.u<<",\"v\":"<<w.v
                 <<",\"passes_endpoints\":"<<(w.endpoints?"true":"false")
                 <<",\"passes_fifth\":"<<(w.five?"true":"false")
                 <<",\"passes_lifted\":"<<(w.lifted?"true":"false")<<'}';
    }
    std::cout<<"]}\n";
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
