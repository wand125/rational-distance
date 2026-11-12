#define RATIONAL_SQUARE_V3_NO_MAIN
#include "../all_plane_gap_v3.cpp"
#include "lifted_projection.hpp"
int main() {
    LiftedProjection<25,5> p25;LiftedProjection<125,5> p125;
    LiftedProjection<49,7> p49;LiftedProjection<121,11> p121;LiftedProjection<169,13> p169;
    // The depth-one rejection equals the entire mod25 refinement of mod5.
    LocalSieve local;
    u64 checked=0;
    for(u32 n=1;n<=25;++n) if(n%5) for(u32 u=0;u<25;++u) for(u32 v=0;v<25;++v) {
        auto first_digit=[&](u32 t){return (t%5!=n%5 && t%5!=(5-n%5)%5) || t%25==n%25 || t%25==(25-n%25)%25;};
        if(p25.accepts(n,u,v)!=(local.p5.accepts(n,u,v) && first_digit(u) && first_digit(v)))
            throw std::logic_error("mod25 one-variable equivalence failed");
        ++checked;
    }
    ArithmeticTables arithmetic(50000);GapWheel wheel(local);ProgressionWheel<19> fw(local);
    SearchOptions opt;opt.bound=100000;
    u64 total=0,fifth_rejected=0,all_before=0,all_after=0,actual_three_rejected=0;
    std::vector<std::array<u32,3>> examples;
    run_gap_search(opt,arithmetic,local,wheel,[&](u32 N,u32 p,u32 q){
        u32 n=N/2,u=p/2,v=q/2;++total;
        bool fifth=FifthPowerCondition::accepts(n,u,v);
        if(!fifth){++fifth_rejected;if(examples.size()<10)examples.push_back({n,u,v});}
        bool lift=p25.accepts(n,u,v)&&p125.accepts(n,u,v)&&p49.accepts(n,u,v)&&p121.accepts(n,u,v)&&p169.accepts(n,u,v);
        all_before+=lift;all_after+=lift&&fifth;
    },[&](i64 N,i64 p,i64 q,i64,i64,i64,Stats&){actual_three_rejected+=!FifthPowerCondition::accepts(N/2,p/2,q/2);},false,true,true,&fw,true);
    std::cout<<"{\"mod25_full_equivalence_checks\":"<<checked<<",\"v3_parameters\":"<<total<<",\"fifth_rejected\":"<<fifth_rejected<<",\"lifted_survivors_before\":"<<all_before<<",\"lifted_survivors_after\":"<<all_after<<",\"actual_three_points_rejected\":"<<actual_three_rejected<<",\"examples\":[";
    for(unsigned j=0;j<examples.size();++j){std::cout<<(j?",":"")<<"[";for(unsigned k=0;k<3;++k)std::cout<<(k?",":"")<<examples[j][k];std::cout<<"]";}
    std::cout<<"]}\n";
}
