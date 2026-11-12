// Bounded diagnostic: the frozen full finite-field sieve is reused unchanged.
// This does not search for or assert a four-distance solution.
#define main frozen_reference_main
#include "../all_plane_sieved.cpp"
#undef main

static u32 gamma(u32 n,u32 t) {
    const u32 d=std::gcd(n,t);
    return std::gcd(d,n/d);
}

static unsigned vp(u32 n,u32 p) {
    unsigned e=0;
    while(n && n%p==0) {++e;n/=p;}
    return e;
}

static bool passes_critical_5_13(u32 n,u32 h) {
    for(u32 p:{5u,13u}) {
        unsigned e=vp(n,p),k=vp(h,p);
        if(e && e%2==0 && 2*k==e) return false;
    }
    return true;
}

int main(int argc,char**argv) {
    const u32 limit=argc>1?std::stoul(argv[1]):12000;
    if(limit>50000) throw std::runtime_error("bounded diagnostic cap50000");
    const ArithmeticTables arithmetic(std::max<u32>(1000,limit));
    const LocalSieve local;
    u64 old_pass=0, gamma_fail=0,gamma_noncritical=0,gamma_no_critical_coordinate=0;
    std::array<u32,3> noncritical_example{};
    bool got_noncritical=false;
    bool first=true;
    std::cout << "{\"max_n\":" << limit << ",\"examples\":[";
    for(u32 n=12;n<=limit;n+=12) {
        // A gamma mismatch beyond current2/3mod4 valuations needs a squared
        // prime1mod4 dividing n. Skipping other n does not claim full coverage.
        bool interesting=false;
        for(unsigned i=0;i<arithmetic.plans[n].count;++i) {
            u32 p=arithmetic.plans[n].prime[i];
            if(p%4==1 && arithmetic.plans[n].power[i]>p) interesting=true;
        }
        if(!interesting) continue;
        const auto permitted=permitted_valuations(n,arithmetic.plans[n]);
        for(u32 u=0;u<=n;++u) {
            if(!permitted[u]) continue;
            const u32 d=std::gcd(n,u);
            generate_crt(n,u,arithmetic,[&](u32 v){
                if(!permitted[v] || !permitted_difference(n,u,v,arithmetic) ||
                   !local.accepts(n,u,v) || !permitted_scale(n,u,v,d)) return;
                ++old_pass;
                const u32 a=gamma(n,u),b=gamma(n,v),c=gamma(n,u-v);
                if(a==b && b==c) return;
                ++gamma_fail;
                if(passes_critical_5_13(n,u-v)) {
                    ++gamma_noncritical;
                    if(passes_critical_5_13(n,u) && passes_critical_5_13(n,v)) {
                        ++gamma_no_critical_coordinate;
                        if(!got_noncritical && 0<v && v<u && u<n) {
                            noncritical_example={n,u,v};got_noncritical=true;
                        }
                    }
                }
                if(gamma_fail<=12) {
                    if(!first) std::cout << ',';
                    first=false;
                    std::cout << "{\"n\":" << n << ",\"u\":" << u << ",\"v\":" << v
                              << ",\"h\":" << u-v << ",\"gamma\":[" << a << ',' << b << ',' << c
                              << "],\"all_existing_filters_pass\":true}";
                }
            },d);
        }
    }
    std::cout << "],\"old_pass_on_selected_denominators\":" << old_pass
              << ",\"gamma_rejected\":" << gamma_fail
              << ",\"gamma_rejected_after_critical5_13_h_only\":" << gamma_noncritical
              << ",\"gamma_rejected_after_critical5_13_all_three\":" << gamma_no_critical_coordinate
              << ",\"noncritical_nondegenerate_example\":";
    if(!got_noncritical) std::cout << "null";
    else {
        const auto [n,u,v]=noncritical_example;
        std::cout << "{\"n\":" << n << ",\"u\":" << u << ",\"v\":" << v
                  << ",\"h\":" << u-v << ",\"gamma\":[" << gamma(n,u)
                  << ',' << gamma(n,v) << ',' << gamma(n,u-v) << "]}";
    }
    std::cout << "}\n";
}
