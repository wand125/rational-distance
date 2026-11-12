#define RATIONAL_SQUARE_V3_NO_MAIN
#include "../all_plane_gap_v3.cpp"

// Independent definition of the critical residue mask, by all signed radii.
// b^2-a^2 and c^2-a^2 are nonzero and have ratio a square root of -1.
// The fourth squared radius is b^2+c^2-a^2. The 2 factors in u and h cancel.
static u64 enumerate_critical_mask(u32 p) {
    std::vector<bool> squares(p,false);
    for(u32 a=0;a<p;++a) squares[a*a%p]=true;
    u64 mask=0;
    for(u32 a=0;a<p;++a) for(u32 b=0;b<p;++b) for(u32 c=0;c<p;++c) {
        const u32 first=(b*b+p-a*a%p)%p,second=(c*c+p-a*a%p)%p;
        if(!first || !second || (first*first+second*second)%p) continue;
        const u32 fourth=(b*b+c*c+p-a*a%p)%p;
        if(!squares[fourth] || b==c) continue;
        const u32 z=u64((b+p-a)%p)*inverse_mod((b+p-c)%p,p)%p;
        mask|=1ULL<<z;
    }
    return mask;
}

int main() try {
    std::vector<std::pair<u32,u64>> masks;
    for(u32 p:{5u,13u,17u,29u,37u,41u}) {
        const u64 actual=enumerate_critical_mask(p);
        const u64 expected=(p==5 || p==13)?0:critical_mask(p);
        if(actual!=expected) throw std::runtime_error("critical mask mismatch p="+std::to_string(p));
        masks.emplace_back(p,actual);
    }

    constexpr u32 bound=50000;
    ArithmeticTables arithmetic(bound/2);
    LocalSieve local;
    GapWheel gapwheel(local);
    ProgressionWheel<19> progression(local);
    validate_field_wheel(arithmetic,local,progression);

    using Key=std::array<u32,3>;
    std::set<Key> baseline,expected,actual;
    u64 rejected=0,gamma_reduction_comparisons=0;
    SearchOptions options;options.bound=bound;
    auto no_point=[](i64,i64,i64,i64,i64,i64,Stats&){};
    const auto old=run_gap_search(options,arithmetic,local,gapwheel,
        [&](u32 N,u32 p,u32 q) {
            const u32 n=N/2,u=p/2,v=q/2,h=u-v;
            if(!baseline.insert({N,p,q}).second) throw std::runtime_error("duplicate baseline parameter");
            const u32 gh=gamma_value(n,h),gu=gamma_value(n,u),gv=gamma_value(n,v);
            const bool full_gamma=(gh==gu && gh==gv && permitted_content(gh,arithmetic));
            const u32 relevant=gamma_relevant_modulus(arithmetic.plans[n]);
            const u32 reduced_h=gamma_value(relevant,h);
            const bool reduced_gamma=(reduced_h==gamma_value(relevant,u) &&
                                      reduced_h==gamma_value(relevant,v));
            if(full_gamma!=reduced_gamma) throw std::runtime_error("reduced gamma mismatch");
            ++gamma_reduction_comparisons;
            const CriticalContentRows critical(arithmetic.plans[n]);
            const CriticalProjectionPlan projection(arithmetic.plans[n]);
            const bool allowed=full_gamma && critical.accepts(h) && critical.accepts(u) &&
                               critical.accepts(v) && projection.accepts(u,h);
            if(allowed) expected.insert({N,p,q});
            else ++rejected;
        },no_point,false,false,false,nullptr,false);
    const auto refined=run_gap_search(options,arithmetic,local,gapwheel,
        [&](u32 N,u32 p,u32 q) {
            if(!actual.insert({N,p,q}).second) throw std::runtime_error("duplicate refined parameter");
        },no_point,false,true,true,&progression,true);
    if(expected.empty() || rejected==0 || expected!=actual)
        throw std::runtime_error("optimized set differs from independent late postfilter");
    if(old.exact.parameter_pairs!=baseline.size() || refined.exact.parameter_pairs!=actual.size())
        throw std::runtime_error("parameter counters differ from set cardinalities");
    std::cout<<"{\"bound\":"<<bound<<",\"critical_masks_independently_enumerated\":true"
             <<",\"progression_wheel_generator_verified\":true"
             <<",\"nonempty_parameter_sets_equal\":true"
             <<",\"duplicates\":0,\"baseline_parameters\":"<<baseline.size()
             <<",\"refined_parameters\":"<<actual.size()<<",\"rejected_parameters\":"<<rejected
             <<",\"full_vs_reduced_gamma_comparisons\":"<<gamma_reduction_comparisons
             <<",\"masks\":[";
    for(std::size_t i=0;i<masks.size();++i) {
        if(i)std::cout<<',';
        std::cout<<"{\"prime\":"<<masks[i].first<<",\"mask_decimal\":"<<masks[i].second<<'}';
    }
    std::cout<<"]}\n";
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
