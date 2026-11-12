// Extra range audit; never edits or changes the original search programs.
// Build with -fsanitize=undefined,float-cast-overflow -fno-sanitize-recover=all.
#define main included_sieved_main
#include "all_plane_sieved.cpp"
#undef main

int main() {
    Stats stats;
    u64 positive=0;
    const std::array<i64,3> bounds={1999992,2000000,2097151};
    for(i64 N:bounds) {
        std::array<std::pair<i64,i64>,9> pairs={{{0,0},{N,0},{N,N},{N,1},
            {N,2},{N-1,N-1},{N-2,0},{N-2,N-4},{N/2,N/3}}};
        for(auto [p,q]:pairs) evaluate(N,p,q,false,stats,[&](i64 x,i64 y,i64 a) {
            if(a<0 || i128(a)>i128(N)*N*N) throw std::logic_error("A bound failed");
            if(i128(magnitude(x))>i128(N)*N*N || i128(magnitude(y))>i128(N)*N*N)
                throw std::logic_error("coordinate bound failed");
            u64 d;
            square(u128(i128(N-x)*(N-x)+i128(N-y)*(N-y)),d);
            if(a>0) ++positive;
        });
    }
    // Near-zero Q with comparatively large positive radicands.
    evaluate(1331714,941664,941664,false,stats,[](i64,i64,i64){});
    evaluate(551614,390050,390050,false,stats,[](i64,i64,i64){});
    // An actual exterior primitive three-distance point near the new round cap.
    bool exterior_seen=false;
    const i64 N=1999992, expected=N*N/4-1;
    evaluate(N,N,2,true,stats,[&](i64 x,i64 y,i64 a) {
        if(x!=-expected || y!=0 || a!=expected) throw std::logic_error("exterior regression");
        exterior_seen=true;
    });
    if(!exterior_seen || !positive) throw std::logic_error("empty positive verification");
    const u128 max_signed=(u128(1)<<127)-1;
    for(u128 value:std::array<u128,4>{max_signed,max_signed-1,u128(1)<<126,(u128(1)<<126)-1}) {
        const u64 r=isqrt(value);
        if(u128(r)*r>value || u128(r+1)*(r+1)<=value) throw std::logic_error("wide square root");
    }
    ArithmeticTables arithmetic(1048575);
    unsigned count=0;
    const u64 residue_count=generate_crt(510510,1,arithmetic,[&](u32){++count;});
    if(residue_count!=128 || arithmetic.plans[510510].count!=7 || !count)
        throw std::logic_error("maximum CRT capacity");
    std::cout << "Range audit passed: N=1999992,2000000,2097151; extreme K; "
                 "wide integer square roots; positive exterior point; 7 factors/128 CRT residues. "
              << positive << " additional positive three-distance reconstructions.\n";
}
