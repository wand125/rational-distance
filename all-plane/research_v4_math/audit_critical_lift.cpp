#define RATIONAL_SQUARE_V3_NO_MAIN
#include "../all_plane_gap_v3.cpp"

// Read-only v3 callback: the frozen generator and reconstruction are unchanged.
static bool new41(u32 n,u32 u,u32 v,bool &eligible,u32 &z,u32 &r) {
    eligible=false; z=r=0;
    u32 e=0,temp=n,power=1;
    while(temp%41==0) {temp/=41;++e;}
    if(!e || (e&1)) return true;
    for(u32 k=0;k<e/2;++k) power*=41;
    const u32 h=u-v;
    if(h%power || (h/power)%41==0) return true;
    if(u%power) throw std::runtime_error("v3 gamma did not force content");
    eligible=true;
    const u32 U=u/power,H=h/power;
    z=u64(U%41)*inverse_mod(H%41,41)%41;
    r=u64(U%1681)*inverse_mod(H%1681,1681)%1681;
    switch(z) {
        case 7:return r==48;
        case 12:return r==750;
        case 15:return r==1327;
        case 27:return r==355;
        case 30:return r==932;
        case 35:return r==1634;
        default:return true;
    }
}

int main(int argc,char **argv) try {
    const u32 bound=argc>1?std::stoul(argv[1]):100000;
    if(bound>100000) throw std::runtime_error("bounded verification only: N<=100000");
    ArithmeticTables arithmetic(bound/2);
    LocalSieve local;
    GapWheel gapwheel(local);
    ProgressionWheel<19> progression(local);
    SearchOptions options;options.bound=bound;
    u64 tested=0,eligible_count=0,rejected=0;
    std::vector<std::array<u32,5>> examples;
    const auto stats=run_gap_search(options,arithmetic,local,gapwheel,
        [&](u32 N,u32 twice_u,u32 twice_v) {
            ++tested;
            const u32 n=N/2,u=twice_u/2,v=twice_v/2;
            bool eligible;u32 z,r;
            const bool allowed=new41(n,u,v,eligible,z,r);
            eligible_count+=eligible;
            if(!allowed) {
                ++rejected;
                if(examples.size()<20) examples.push_back({n,u,v,z,r});
            }
        },[](i64,i64,i64,i64,i64,i64,Stats&){},false,true,true,&progression,true);
    std::cout<<"{\"bound\":"<<bound<<",\"v3_parameters\":"<<tested
             <<",\"critical_41_eligible\":"<<eligible_count<<",\"new_rejected\":"<<rejected
             <<",\"primitive_three_diagnostics\":"<<stats.exact.primitive_three_candidates
             <<",\"examples\":[";
    for(std::size_t j=0;j<examples.size();++j) {
        if(j)std::cout<<',';
        const auto &x=examples[j];
        std::cout<<"{\"n\":"<<x[0]<<",\"u\":"<<x[1]<<",\"v\":"<<x[2]
                 <<",\"r_mod_41\":"<<x[3]<<",\"r_mod_1681\":"<<x[4]<<'}';
    }
    std::cout<<"]}\n";
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
