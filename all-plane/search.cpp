// Exact bounded search for the rational-distance problem, square INTERIOR.
// Build: g++ -O3 -DNDEBUG -std=c++20 search.cpp -o square_search
// Run:   ./square_search --limit 1000000 --near near.csv
// Every final square decision is verified with integer arithmetic.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <tuple>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;
using Clock = std::chrono::steady_clock;

struct Edge { u32 leg, scale; };
struct Near { u32 n, a, c; u64 left_low, left_high, cross_square_root;
              char missing; u64 missing_squared, nearest_root, residual; };

u64 integer_root(u64 n) {
    // Floating point is only a seed, never a decision. Both correction loops
    // and the eventual equality check use 128-bit integer multiplication.
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((u128)(r + 1) * (r + 1) <= n) ++r;
    while ((u128)r * r > n) --r;
    return r;
}
bool square(u64 n, u64 &r) { r = integer_root(n); return (u128)r*r == n; }
u64 sumsq(u32 a, u32 b) { return (u64)a*a + (u64)b*b; }

std::vector<std::uint8_t> good_numbers(u32 limit) {
    // GOOD = positive odd integers with no prime factor congruent to 3 or 5 mod 8.
    std::vector<std::uint8_t> good((std::size_t)limit + 1, 1);
    std::vector<bool> composite((std::size_t)limit + 1, false);
    for (u64 j=0; j<=limit; j+=2) good[j]=0;
    for (u64 p=3; p<=limit; p+=2) {
        if (composite[p]) continue;
        if (p*p<=limit) for (u64 j=p*p; j<=limit; j+=2*p) composite[j]=true;
        if (p%8==3 || p%8==5)
            for (u64 j=p; j<=limit; j+=2*p) good[j]=0;
    }
    return good;
}
template<u32 M> std::array<std::uint8_t,M> qr_table() {
    std::array<std::uint8_t,M> q{};
    for (u64 i=0; i<M; ++i) q[i*i%M]=1;
    return q;
}

int main(int argc, char **argv) try {
    u32 bound=100000, side_multiple=24;
    bool use_good=true, use_sieve=true;
    std::string audit_path, near_path;
    for (int i=1; i<argc; ++i) {
        std::string s=argv[i];
        auto arg = [&]() -> std::string { if (++i>=argc) throw std::runtime_error("missing argument"); return argv[i]; };
        if (s=="--limit") { auto v=std::stoull(arg()); if(v<4 || v>100000000) throw std::runtime_error("limit must be 4..100000000"); bound=(u32)v; }
        else if(s=="--side-multiple") side_multiple=(u32)std::stoul(arg());
        else if(s=="--no-good") use_good=false;
        else if(s=="--no-sieve") use_sieve=false;
        else if(s=="--audit-wedges") audit_path=arg();
        else if(s=="--near") near_path=arg();
        else if(s=="--help") {
            std::cout << "--limit B [--side-multiple 4|12|24] [--no-good] [--no-sieve] [--audit-wedges CSV] [--near CSV]\n";
            return 0;
        } else throw std::runtime_error("unknown argument: "+s);
    }
    if(side_multiple!=4 && side_multiple!=12 && side_multiple!=24)
        throw std::runtime_error("side multiple must be 4, 12, or 24");
    auto t0=Clock::now();
    auto good=use_good ? good_numbers(2*bound) : std::vector<std::uint8_t>();
    auto t1=Clock::now();
    const u32 max_odd=bound/2, nodes=(max_odd+1)/2, groups=side_multiple/4;
    std::vector<u32> counts(nodes,0);
    u64 primitive_triples=0;
    auto generate = [&](auto emit, bool count_triples) {
        // m^2 <= 2B is a safe bound whenever both primitive legs are <= B.
        for (u64 m=2; m*m<=2ull*bound; ++m) {
            for (u64 n=(m%2 ? 2 : 1); n<m; n+=2) {
                const u64 odd=m*m-n*n, even=2*m*n;
                if(even>bound) break;
                if(odd>max_odd || std::gcd(m,n)!=1) continue;
                if(count_triples) ++primitive_triples;
                u64 last=std::min((u64)max_odd/odd,(u64)bound/even);
                for(u64 k=1; k<=last; k+=2) {
                    if(use_good && !good[k]) continue;
                    emit((u32)(k*odd),(u32)(k*even),(u32)k);
                }
            }
        }
    };
    generate([&](u32 a,u32,u32){ ++counts[a/2]; },true);
    std::vector<u64> start((std::size_t)nodes+1,0);
    for(u32 i=0;i<nodes;++i) start[i+1]=start[i]+counts[i];
    const u64 edge_count=start.back();
    std::vector<Edge> edges(edge_count);
    std::fill(counts.begin(),counts.end(),0);
    generate([&](u32 a,u32 c,u32 k){ edges[start[a/2]+counts[a/2]++]={c,k}; },false);
    for(u32 i=0;i<nodes;++i) {
        auto lo=edges.begin()+start[i], hi=edges.begin()+start[i+1];
        std::sort(lo,hi,[groups](const Edge &a,const Edge &b){
            const u32 ra=(a.leg/4)%groups, rb=(b.leg/4)%groups;
            return ra==rb ? a.leg<b.leg : ra<rb;
        });
    }
    auto t2=Clock::now();
    const auto qr1=qr_table<5005>();             // 5*7*11*13
    const auto qr2=qr_table<215441>();           // 17*19*23*29
    std::ofstream audit,near_file;
    if(!audit_path.empty()) { audit.open(audit_path); if(!audit) throw std::runtime_error("cannot open audit file"); audit<<"N,X,Y\n"; }
    u64 raw_wedges=0,residue_wedges=0,bounded_wedges=0,primitive_wedges=0,good_wedges=0;
    u64 pass_qr1=0,pass_qr2=0,exact_square_s=0,exact_square_t=0,three_after_sieves=0,solutions=0;
    u64 tuple_xor=0,tuple_sum=0;
    std::vector<std::array<u64,7>> found;
    std::vector<Near> near;
    auto mix=[](u64 z) { z+=0x9e3779b97f4a7c15ull; z=(z^(z>>30))*0xbf58476d1ce4e5b9ull; z=(z^(z>>27))*0x94d049bb133111ebull; return z^(z>>31); };
    for(u32 node=0;node<nodes;++node) {
        const u32 a=2*node+1;
        const u64 begin=start[node],end=start[node+1],degree=end-begin;
        raw_wedges+=degree*(degree+1)/2;
        if(!degree) continue;
        std::array<u64,7> cut{};
        u64 z=begin;
        for(u32 r=0;r<groups;++r) {
            cut[r]=z;
            while(z<end && (edges[z].leg/4)%groups==r) ++z;
        }
        cut[groups]=end;
        for(u32 r=0;r<groups;++r) {
            u32 s=(groups-r)%groups;
            if(r>s) continue;
            u64 nr=cut[r+1]-cut[r],ns=cut[s+1]-cut[s];
            residue_wedges+=r==s ? nr*(nr+1)/2 : nr*ns;
            for(u64 i=cut[r];i<cut[r+1];++i) {
                const u32 e=edges[i].leg;
                if(e>=bound) continue;
                const u32 low=2ull*a>e ? 2*a-e : 0;
                const u32 high=bound-e;
                if(low>high) continue;
                auto sb=edges.begin()+cut[s],se=edges.begin()+cut[s+1];
                auto jl=std::lower_bound(sb,se,low,[](const Edge &e,u32 x){return e.leg<x;});
                auto jh=std::upper_bound(jl,se,high,[](u32 x,const Edge &e){return x<e.leg;});
                if(r==s && jl<edges.begin()+i) jl=edges.begin()+i;
                for(auto it=jl;it<jh;++it) {
                    ++bounded_wedges;
                    const u32 f=it->leg,N=e+f,b=N-a,c=std::min(e,f),d=std::max(e,f);
                    if(edges[i].scale!=1 && it->scale!=1 && std::gcd(edges[i].scale,it->scale)!=1) continue;
                    ++primitive_wedges;
                    if(use_good && (!good[b+c] || !good[b+d])) continue;
                    ++good_wedges;
                    // Diagnostic fingerprint only, not a proof certificate.
                    // Mix the full value of each coordinate without truncation.
                    u64 hash=mix(mix(mix(N)^a)^c);
                    tuple_xor^=hash; tuple_sum+=hash;
                    if(audit.is_open()) audit<<N<<','<<a<<','<<c<<'\n';
                    const u64 S=sumsq(b,c),T=sumsq(b,d);
                    if(use_sieve && (!qr1[S%5005] || !qr1[T%5005])) continue;
                    ++pass_qr1;
                    if(use_sieve && (!qr2[S%215441] || !qr2[T%215441])) continue;
                    ++pass_qr2;
                    u64 rs,rt;
                    bool ss=square(S,rs),tt=square(T,rt);
                    exact_square_s+=ss; exact_square_t+=tt;
                    if(ss || tt) ++three_after_sieves;
                    if(ss && tt) {
                        ++solutions;
                        found.push_back({N,a,c,integer_root(sumsq(a,c)),rs,integer_root(sumsq(a,d)),rt});
                    } else if(!near_path.empty() && (ss || tt)) {
                        u64 value=ss?T:S, k=ss?rt:rs;
                        u64 below=value-k*k,above=(k+1)*(k+1)-value;
                        if(above<below) ++k;
                        u64 residual=std::min(below,above);
                        near.push_back({N,a,c,integer_root(sumsq(a,c)),integer_root(sumsq(a,d)),ss?rs:rt,ss?'T':'S',value,k,residual});
                    }
                }
            }
        }
    }
    auto t3=Clock::now();
    if(!near_path.empty()) {
        std::sort(near.begin(),near.end(),[](const Near &a,const Near &b){ return std::tie(a.residual,a.n,a.a,a.c)<std::tie(b.residual,b.n,b.a,b.c); });
        near_file.open(near_path); if(!near_file) throw std::runtime_error("cannot open near file");
        near_file<<"N,X,Y,A,C,other_integer_distance,missing,S,nearest_root,integer_residual\n";
        for(const auto &v:near) near_file<<v.n<<','<<v.a<<','<<v.c<<','<<v.left_low<<','<<v.left_high<<','<<v.cross_square_root<<','<<v.missing<<','<<v.missing_squared<<','<<v.nearest_root<<','<<v.residual<<'\n';
    }
    struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
    auto seconds=[](auto a,auto b){ return std::chrono::duration<double>(b-a).count(); };
    std::cout<<std::fixed<<std::setprecision(6);
    std::cout<<"{\n  \"scope\":\"strict_interior_primitive_common_denominator\",\n"
      <<"  \"limit\":"<<bound<<",\n  \"side_multiple\":"<<side_multiple<<",\n"
      <<"  \"good_filter\":"<<(use_good?"true":"false")<<",\n  \"modular_sieve\":"<<(use_sieve?"true":"false")<<",\n"
      <<"  \"primitive_triples\":"<<primitive_triples<<",\n  \"edges\":"<<edge_count<<",\n"
      <<"  \"raw_wedges\":"<<raw_wedges<<",\n  \"residue_wedges\":"<<residue_wedges<<",\n"
      <<"  \"bounded_wedges\":"<<bounded_wedges<<",\n  \"primitive_wedges\":"<<primitive_wedges<<",\n"
      <<"  \"good_wedges\":"<<good_wedges<<",\n  \"pass_mod_5005\":"<<pass_qr1<<",\n"
      <<"  \"pass_mod_215441\":"<<pass_qr2<<",\n  \"exact_square_S\":"<<exact_square_s<<",\n"
      <<"  \"exact_square_T\":"<<exact_square_t<<",\n  \"three_distances_after_sieves\":"<<three_after_sieves<<",\n"
      <<"  \"solutions\":"<<solutions<<",\n  \"tuple_xor\":"<<tuple_xor<<",\n  \"tuple_sum\":"<<tuple_sum<<",\n"
      <<"  \"seconds_good_sieve\":"<<seconds(t0,t1)<<",\n  \"seconds_graph\":"<<seconds(t1,t2)<<",\n"
      <<"  \"seconds_search\":"<<seconds(t2,t3)<<",\n  \"seconds_total\":"<<seconds(t0,t3)<<",\n"
      <<"  \"peak_rss_kib\":"<<usage.ru_maxrss<<",\n  \"witnesses\":[";
    for(std::size_t i=0;i<found.size();++i) {
        if(i) std::cout<<',';
        std::cout<<'[';
        for(int j=0;j<7;++j) {if(j)std::cout<<',';std::cout<<found[i][j];}
        std::cout<<']';
    }
    std::cout<<"]\n}\n";
    return 0;
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
