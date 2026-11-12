// Exact whole-plane search for primitive integer square configurations.
// Build: g++ -O3 -std=c++17 all_plane.cpp -o all_plane
// Run:   ./all_plane 1000 --validate
// This is an independent reference search, not the fast interior search.
//
// Symmetry: X <= Y <= N/2. Let A=distance to (0,0),
// p=distance to (N,0)-A, q=distance to (0,N)-A.
// Then 0 <= q <= p <= N and
//   Q=p^2+q^2-N^2, T=p(N^2-p^2)+q(N^2-q^2),
//   R=(N^2-p^2)^2+(N^2-q^2)^2,
//   4 Q A^2 - 4 T A + R = 0.
// K=T^2-Q R=(N^2-p^2)(N^2-q^2)(2N^2-(p-q)^2).
// For Q!=0, A=(T +/- sqrt(K))/(2Q).
// For Q==0,T!=0, A=R/(4T). Q==T==0 has R>0.
// Coordinates follow from X=(N^2-p^2-2Ap)/(2N), likewise Y.
// All roots, including boundary and median cases, are handled.
//
// Only in FOUR-distance search mode, proven global necessary conditions are
// used: 24|N; all four distances odd; X,Y have opposite parity, the even
// coordinate is divisible by 4; precisely one of X,Y is divisible by 3.
// Thus p,q are even. The validation instead enumerates ALL integer
// three-distance points, without these four-distance-only restrictions.
//
// Bound: K<=2N^6 and every nonnegative A root is <=2N^3.
// N<=1,000,000 keeps every intermediate inside signed __int128.
// Integer square tests use a floating estimate followed by exact correction.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

using i128 = __int128_t;
using u128 = __uint128_t;
using i64 = int64_t;
using u64 = uint64_t;
using Points = std::set<std::pair<i64,i64>>;

struct Stats {
    u64 sides=0, parameter_pairs=0, linear_pairs=0;
    u64 square_discriminants=0, integral_A=0, integer_coordinates=0;
    u64 primitive_three_candidates=0, solutions=0;
};

static u64 isqrt(u128 z) {
    u64 r=static_cast<u64>(std::sqrt(static_cast<long double>(z)));
    while (u128(r)*r>z) --r;
    while (u128(r+1)*(r+1)<=z) ++r;
    return r;
}

static bool square(u128 z, u64 &root) {
    // These filters are exact. No floating-point equality is used.
    static constexpr unsigned q5=(1u<<0)|(1u<<1)|(1u<<4);
    static constexpr unsigned q7=(1u<<0)|(1u<<1)|(1u<<2)|(1u<<4);
    if (!((q5 >> unsigned(z%5))&1) || !((q7 >> unsigned(z%7))&1))
        return false;
    root=isqrt(z);
    return u128(root)*root==z;
}

static i64 magnitude(i64 x) { return x<0 ? -x : x; }

// Callback receives exact integer THREE-distance points.
template<class Callback>
static void evaluate(i64 N, i64 p, i64 q, bool four_filters,
                     Stats &st, Callback accept) {
    ++st.parameter_pairs;
    const i128 nn=i128(N)*N, pp=i128(p)*p, qq=i128(q)*q;
    const i128 U=nn-pp, V=nn-qq;
    const i128 Q=pp+qq-nn;
    const i128 T=i128(p)*U+i128(q)*V;
    const i128 R=U*U+V*V;
    i128 roots[2]; int nr=0;
    if (Q==0) {
        ++st.linear_pairs;
        if (T==0) {
            if (R==0) throw std::logic_error("identically zero quadratic");
            return;
        }
        if (R%(4*T)!=0) return;
        roots[nr++]=R/(4*T);
    } else {
        const i128 delta=i128(p)-q;
        const u128 K=u128(U*V*(2*nn-delta*delta));
        u64 sr;
        if (!square(K,sr)) return;
        ++st.square_discriminants;
        const i128 den=2*Q;
        const i128 plus=T+i128(sr), minus=T-i128(sr);
        if (plus%den==0) roots[nr++]=plus/den;
        if (sr!=0 && minus%den==0) roots[nr++]=minus/den;
    }
    for (int k=0;k<nr;++k) {
        const i128 A=roots[k];
        if (A<0) continue;
        ++st.integral_A;
        if (four_filters && A%2==0) continue;
        const i128 xn=U-2*A*p, yn=V-2*A*q, den=2*i128(N);
        if (xn%den || yn%den) continue;
        const i128 X=xn/den,Y=yn/den;
        if (X>Y || 2*Y>N) throw std::logic_error("symmetry failure");
        if (X*X+Y*Y!=A*A ||
            (N-X)*(N-X)+Y*Y!=(A+p)*(A+p) ||
            X*X+(N-Y)*(N-Y)!=(A+q)*(A+q))
            throw std::logic_error("three-distance identity failure");
        if (A>2*i128(N)*N*N) throw std::logic_error("radius bound failure");
        ++st.integer_coordinates;
        const i64 x=i64(X), y=i64(Y), a=i64(A);
        if (four_filters) {
            if (std::gcd(std::gcd(magnitude(x),magnitude(y)),N)!=1) continue;
            const bool xodd=(x%2!=0), yodd=(y%2!=0);
            if (xodd==yodd) continue;
            if ((xodd ? y : x)%4!=0) continue;
            if ((x%3==0)==(y%3==0)) continue;
            ++st.primitive_three_candidates;
        }
        accept(x,y,a);
    }
}

static void validate() {
    u64 total=0;
    // Independent coordinate enumeration over the PROVED complete bound.
    // It includes exterior points and does not use four-distance parity rules.
    for (i64 N=1;N<=12;++N) {
        Points algebraic,brute;
        Stats ignored;
        for (i64 p=0;p<=N;++p) for (i64 q=0;q<=p;++q)
            evaluate(N,p,q,false,ignored,[&](i64 x,i64 y,i64){
                algebraic.emplace(x,y);
            });
        const i64 radius=2*N*N*N;
        for (i64 y=-radius;y<=N/2;++y) for (i64 x=-radius;x<=y;++x) {
            u64 a,b,c;
            if (!square(u128(i128(x)*x+i128(y)*y),a)) continue;
            if (!square(u128(i128(N-x)*(N-x)+i128(y)*y),b)) continue;
            if (!square(u128(i128(x)*x+i128(N-y)*(N-y)),c)) continue;
            brute.emplace(x,y);
        }
        if (algebraic!=brute)
            throw std::runtime_error("independent validation mismatch at N="+std::to_string(N));
        total+=algebraic.size();
    }
    // Q=0 with an actual integral three-distance point, not merely no roots.
    Stats ignored;
    Points linear;
    evaluate(8400,6720,5040,false,ignored,[&](i64 x,i64 y,i64 a){
        if (a!=1685) throw std::logic_error("linear-case radius failure");
        linear.emplace(x,y);
    });
    if (linear!=Points{{164,1677}})
        throw std::runtime_error("nontrivial Q=0 validation failure");
    Points impossible;
    evaluate(24,24,0,false,ignored,[&](i64 x,i64 y,i64){impossible.emplace(x,y);});
    if (!impossible.empty()) throw std::runtime_error("Q=T=0 validation failure");
    std::cerr << "Validation passed: all integer three-distance points in the entire plane "
              << "for N=1..12 (" << total << " canonical points); "
              << "Q=0 integral example and Q=T=0 exception.\n";
}

int main(int argc,char**argv) {
    try {
        i64 B=1000; bool check=false;
        for (int i=1;i<argc;++i) {
            if (std::string(argv[i])=="--validate") check=true;
            else B=std::stoll(argv[i]);
        }
        if (B<1 || B>1000000)
            throw std::runtime_error("B must be in [1,1000000] for audited arithmetic bounds");
        if (check) validate();
        Stats st;
        const auto started=std::chrono::steady_clock::now();
        for (i64 N=24;N<=B;N+=24) {
            ++st.sides;
            for (i64 p=0;p<=N;p+=2) for (i64 q=0;q<=p;q+=2)
                evaluate(N,p,q,true,st,[&](i64 x,i64 y,i64 a){
                    u64 d;
                    if (!square(u128(i128(N-x)*(N-x)+i128(N-y)*(N-y)),d)) return;
                    ++st.solutions;
                    std::cout << "{\"solution\":{\"N\":" << N << ",\"X\":" << x
                              << ",\"Y\":" << y << ",\"A\":" << a
                              << ",\"D\":" << d << "}}\n";
                });
        }
        const double elapsed=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();
        std::cout << "{\"scope\":\"whole_plane\",\"max_minimal_coordinate_denominator\":" << B
                  << ",\"sides\":" << st.sides
                  << ",\"parameter_pairs\":" << st.parameter_pairs
                  << ",\"linear_pairs\":" << st.linear_pairs
                  << ",\"square_discriminants\":" << st.square_discriminants
                  << ",\"nonnegative_integral_A\":" << st.integral_A
                  << ",\"integer_coordinates_after_odd_A_filter\":" << st.integer_coordinates
                  << ",\"primitive_three_candidates_after_four_distance_filters\":" << st.primitive_three_candidates
                  << ",\"solutions\":" << st.solutions
                  << ",\"elapsed_seconds\":" << elapsed << "}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
