#ifndef RATIONAL_SQUARE_WIDE_ARITH_V5_HPP
#define RATIONAL_SQUARE_WIDE_ARITH_V5_HPP
// v5 exact reconstruction for minimal coordinate denominators N <= 10^7.
// Include after all_plane_sieved.cpp (needs i64,u64,i128,u128,Stats).
//
// Bounds (arithmetic_bounds_v5.md), N <= 10^7, 0 <= q <= p <= N:
//   N^2 < 2^47, T <= N^3 < 2^70, R <= 2N^4 < 2^95, U*V <= N^4 < 2^94,
//   W = 2N^2-(p-q)^2 <= 2N^2 < 2^48, K = U*V*W <= 2N^6 < 2^141,
//   sqrt(K) < 2^71, every nonnegative root A <= N^3 < 2^70,
//   |X|,|Y| <= N^3, all squared norms <= 2(N^3+N)^2 < 2^141.
// Everything below 2^126 stays in signed __int128. K and the squared norms
// are formed exactly in a 256-bit unsigned type. No floating point is used:
// the square test is an exact residue prefilter followed by an exact
// digit-by-digit integer square root. No search cut such as A <= N^2/2 is
// applied, so the three-distance diagnostic set is unchanged.

#include <array>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>

struct U256 {
    std::array<u64,4> w{}; // little endian limbs
    static U256 from(u128 x) {U256 r;r.w[0]=u64(x);r.w[1]=u64(x>>64);return r;}
    bool operator==(const U256 &o) const {return w==o.w;}
    bool operator!=(const U256 &o) const {return w!=o.w;}
    bool operator<(const U256 &o) const {
        for(int i=3;i>=0;--i) if(w[i]!=o.w[i]) return w[i]<o.w[i];
        return false;
    }
    bool operator<=(const U256 &o) const {return !(o<*this);}
    bool high_zero() const {return w[2]==0 && w[3]==0;}
    u128 low128() const {return (u128(w[1])<<64)|w[0];}
};

static U256 add256(const U256 &a,const U256 &b) {
    U256 r;u128 carry=0;
    for(int i=0;i<4;++i) {const u128 s=u128(a.w[i])+b.w[i]+carry;r.w[i]=u64(s);carry=s>>64;}
    if(carry) throw std::logic_error("U256 add overflow");
    return r;
}
static U256 sub256(const U256 &a,const U256 &b) { // requires b<=a
    U256 r;u64 borrow=0;
    for(int i=0;i<4;++i) {
        const u128 d=u128(a.w[i])-b.w[i]-borrow;
        r.w[i]=u64(d);borrow=u64(d>>64)&1;
    }
    if(borrow) throw std::logic_error("U256 subtract underflow");
    return r;
}
static U256 shr256(const U256 &a,unsigned s) { // s in {1,2}
    U256 r;
    for(int i=0;i<4;++i) {
        r.w[i]=a.w[i]>>s;
        if(i<3) r.w[i]|=a.w[i+1]<<(64-s);
    }
    return r;
}
static U256 mul128(u128 a,u128 b) {
    const u64 a0=u64(a),a1=u64(a>>64),b0=u64(b),b1=u64(b>>64);
    const u128 p00=u128(a0)*b0,p01=u128(a0)*b1,p10=u128(a1)*b0,p11=u128(a1)*b1;
    U256 r;
    r.w[0]=u64(p00);
    const u128 mid=(p00>>64)+u64(p01)+u64(p10);
    r.w[1]=u64(mid);
    const u128 high=(mid>>64)+(p01>>64)+(p10>>64)+u64(p11);
    r.w[2]=u64(high);
    r.w[3]=u64((high>>64)+(p11>>64));
    return r;
}
static u128 magnitude128(i128 x) {return x<0 ? u128(0)-u128(x) : u128(x);}
static U256 square_i128(i128 x) {const u128 m=magnitude128(x);return mul128(m,m);}

static unsigned mod256(const U256 &a,unsigned m) {
    // sum of limb*2^(64i) mod m, all exact in u128.
    const u128 base=(u128(1)<<64)%m;
    u128 r=0;
    for(int i=3;i>=0;--i) r=(r*base+a.w[i]%m)%m;
    return unsigned(r);
}

// Exact floor square root by the restoring binary method.
static u128 isqrt256(const U256 &x) {
    if(!x.w[3] && !x.w[2] && !x.w[1] && !x.w[0]) return 0;
    int top=255;
    while(!((x.w[top/64]>>(top%64))&1)) --top;
    U256 bit;const int b=top&~1; bit.w[b/64]=u64(1)<<(b%64);
    U256 rem=x,res;
    while(true) {
        const U256 trial=add256(res,bit);
        if(trial<=rem) {rem=sub256(rem,trial);res=add256(shr256(res,1),bit);}
        else res=shr256(res,1);
        if(bit.w[0]==1 && !bit.w[1] && !bit.w[2] && !bit.w[3]) break;
        bit=shr256(bit,2);
    }
    if(!res.high_zero()) throw std::logic_error("isqrt256 root exceeds 128 bits");
    return res.low128();
}

struct SquareResidues {
    std::array<std::uint8_t,64> q64{};std::array<std::uint8_t,63> q63{};
    std::array<std::uint8_t,65> q65{};std::array<std::uint8_t,11> q11{};
    SquareResidues() {
        for(unsigned a=0;a<64;++a) q64[a*a%64]=1;
        for(unsigned a=0;a<63;++a) q63[a*a%63]=1;
        for(unsigned a=0;a<65;++a) q65[a*a%65]=1;
        for(unsigned a=0;a<11;++a) q11[a*a%11]=1;
    }
};
static const SquareResidues &square_residues() {static const SquareResidues s;return s;}

// Exact: true iff z is a perfect square; then root=floor(sqrt(z)).
static bool square256(const U256 &z,u128 &root) {
    const auto &s=square_residues();
    if(!s.q64[z.w[0]&63]) return false;
    // One reduction mod 45045=63*65*11 serves three filters.
    const unsigned r=mod256(z,45045u);
    if(!s.q63[r%63] || !s.q65[r%65] || !s.q11[r%11]) return false;
    root=isqrt256(z);
    return mul128(root,root)==z;
}

static std::string i128_string(i128 x) {
    if(x==0) return "0";
    const bool negative=x<0;u128 m=magnitude128(x);
    std::string s;
    while(m) {s.push_back(char('0'+unsigned(m%10)));m/=10;}
    if(negative) s.push_back('-');
    return std::string(s.rbegin(),s.rend());
}

static constexpr i64 V5_MAX_DENOMINATOR=10000000;

// Same mathematics, root order, counters and filters as evaluate() in
// all_plane_sieved.cpp; only the integer widths differ. Callback receives
// (x,y,a) as i128.
template<class Callback>
static void evaluate_v5(i64 N,i64 p,i64 q,bool four_filters,Stats &st,Callback accept) {
    if(N<1 || N>V5_MAX_DENOMINATOR || q<0 || q>p || p>N)
        throw std::logic_error("evaluate_v5 input outside proved range");
    ++st.parameter_pairs;
    const i128 nn=i128(N)*N,pp=i128(p)*p,qq=i128(q)*q;
    const i128 U=nn-pp,V=nn-qq;
    const i128 Q=pp+qq-nn;
    const i128 T=i128(p)*U+i128(q)*V;
    const i128 R=U*U+V*V;
    i128 roots[2];int nr=0;
    if(Q==0) {
        ++st.linear_pairs;
        if(T==0) {
            if(R==0) throw std::logic_error("identically zero quadratic");
            return;
        }
        if(R%(4*T)!=0) return;
        roots[nr++]=R/(4*T);
    } else {
        const i128 delta=i128(p)-q;
        const i128 W=2*nn-delta*delta;
        if(U<0 || V<0 || W<=0) throw std::logic_error("discriminant factor sign");
        const U256 K=mul128(u128(U)*u128(V),u128(W));
        u128 sr;
        if(!square256(K,sr)) return;
        ++st.square_discriminants;
        if(sr>>100) throw std::logic_error("sqrt(K) beyond proved bound");
        const i128 den=2*Q;
        const i128 plus=T+i128(sr),minus=T-i128(sr);
        if(plus%den==0) roots[nr++]=plus/den;
        if(sr!=0 && minus%den==0) roots[nr++]=minus/den;
    }
    for(int k=0;k<nr;++k) {
        const i128 A=roots[k];
        if(A<0) continue;
        ++st.integral_A;
        if(four_filters && A%2==0) continue;
        const i128 xn=U-2*A*p,yn=V-2*A*q,den=2*i128(N);
        if(xn%den || yn%den) continue;
        const i128 X=xn/den,Y=yn/den;
        if(X>Y || 2*Y>N) throw std::logic_error("symmetry failure");
        if(add256(square_i128(X),square_i128(Y))!=square_i128(A) ||
           add256(square_i128(N-X),square_i128(Y))!=square_i128(A+p) ||
           add256(square_i128(X),square_i128(N-Y))!=square_i128(A+q))
            throw std::logic_error("three-distance identity failure");
        if(A>2*i128(N)*N*N) throw std::logic_error("radius bound failure");
        ++st.integer_coordinates;
        if(four_filters) {
            const i64 gx=i64(magnitude128(X)%u128(N)),gy=i64(magnitude128(Y)%u128(N));
            if(std::gcd(std::gcd(N,gx),gy)!=1) continue;
            const bool xodd=(X%2!=0),yodd=(Y%2!=0);
            if(xodd==yodd) continue;
            if((xodd ? Y : X)%4!=0) continue;
            if((X%3==0)==(Y%3==0)) continue;
            ++st.primitive_three_candidates;
        }
        accept(X,Y,A);
    }
}

#endif
