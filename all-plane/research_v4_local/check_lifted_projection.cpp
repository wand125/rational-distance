#define RATIONAL_SQUARE_V3_NO_MAIN
#include "../all_plane_gap_v3.cpp"
#include "lifted_projection.hpp"

int main(int argc,char**argv) {
    const u32 bound=argc>1?std::stoul(argv[1]):100000;
    if(bound<24 || bound>200000) throw std::runtime_error("diagnostic bound 24..200000");
    LiftedProjection<25,5> p25;
    LiftedProjection<49,7> p49;
    LiftedProjection<121,11> p121;
    LiftedProjection<169,13> p169;
    LiftedProjection<125,5> p125;
    p25.verify();p49.verify();p121.verify();p169.verify();p125.verify();
    ArithmeticTables arithmetic(std::max(bound/2,u32(1045)));
    LocalSieve local;GapWheel wheel(local);ProgressionWheel<19> fieldwheel(local);
    SearchOptions opt;opt.bound=bound;
    std::array<u64,5> independent_rejected{},sequential_rejected{};
    std::array<std::vector<std::array<u32,5>>,5> examples;
    u64 input=0,survivors=0;
    std::vector<std::array<i64,6>> three_survivors;
    const auto stats=run_gap_search(opt,arithmetic,local,wheel,
        [&](u32 N,u32 p,u32 q){
            const u32 n=N/2,u=p/2,v=q/2;
            const std::array<bool,5> okay={p25.accepts(n,u,v),p49.accepts(n,u,v),p121.accepts(n,u,v),p169.accepts(n,u,v),p125.accepts(n,u,v)};
            const std::array<unsigned,5> moduli={25,49,121,169,125};
            ++input;
            bool passed=true;
            for(unsigned i=0;i<okay.size();++i) {
                if(!okay[i]) {
                    ++independent_rejected[i];
                    if(examples[i].size()<4) {
                        const unsigned m=moduli[i],inv=inverse_mod(n%m,m);
                        examples[i].push_back({n,u,v,unsigned(u64(u)*inv%m),unsigned(u64(v)*inv%m)});
                    }
                    if(passed) {++sequential_rejected[i];passed=false;}
                }
            }
            survivors+=passed;
        },
        [&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats&){
            if(p25.accepts(N/2,p/2,q/2) && p49.accepts(N/2,p/2,q/2) && p121.accepts(N/2,p/2,q/2) && p169.accepts(N/2,p/2,q/2) && p125.accepts(N/2,p/2,q/2))
                three_survivors.push_back({N,p,q,x,y,a});
        },false,true,true,&fieldwheel,true);
    std::cout<<"{\"bound\":"<<bound<<",\"projection_verification\":true,\"moduli\":[25,49,121,169,125],\"normal_pair_counts\":["<<p25.size()<<","<<p49.size()<<","<<p121.size()<<","<<p169.size()<<","<<p125.size()<<"],\"input_parameters\":"<<input<<",\"independent_rejected\":[";
    for(unsigned i=0;i<5;++i) std::cout<<(i?",":"")<<independent_rejected[i];
    std::cout<<"],\"sequential_rejected\":[";
    for(unsigned i=0;i<5;++i) std::cout<<(i?",":"")<<sequential_rejected[i];
    std::cout<<"],\"surviving_parameters\":"<<survivors<<",\"examples\":[";
    for(unsigned i=0;i<5;++i) {
        std::cout<<(i?",":"")<<"[";
        for(unsigned j=0;j<examples[i].size();++j) {
            std::cout<<(j?",":"")<<"[";
            for(unsigned k=0;k<5;++k) std::cout<<(k?",":"")<<examples[i][j][k];
            std::cout<<"]";
        }
        std::cout<<"]";
    }
    std::cout<<"],\"input_three_points\":"<<stats.exact.primitive_three_candidates<<",\"surviving_three_points\":[";
    for(unsigned j=0;j<three_survivors.size();++j) {
        std::cout<<(j?",":"")<<"[";
        for(unsigned k=0;k<6;++k) std::cout<<(k?",":"")<<three_survivors[j][k];
        std::cout<<"]";
    }
    std::cout<<"]}\n";
}
