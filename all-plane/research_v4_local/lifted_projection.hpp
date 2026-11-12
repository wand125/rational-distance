#pragma once
// Given the existing mod-5 four-distance projection and 5 not dividing n,
// u/n=+/-1 mod5 forces v/n=+/-2 in the matching class. Consequently the
// other two factors of K=(n^2-u^2)(n^2-v^2)(2*n^2-(u-v)^2) are units
// whose product is a square mod5. The first factor must be a Q_5-square.
// This is a one-variable necessary condition; zero is retained.
struct FifthPowerCondition {
    static bool accepts_one(u32 n,u32 u) {
        if(n%5==0 || (u%5!=n%5 && u%5!=(5-n%5)%5)) return true;
        u64 value=u64(n-u)*(n+u);
        if(!value) return true;
        unsigned exponent=0;
        while(value%5==0) {value/=5;++exponent;}
        return exponent%2==0 && (value%5==1 || value%5==4);
    }
    static bool accepts(u32 n,u32 u,u32 v) {
        return accepts_one(n,u) && accepts_one(n,v);
    }
    static void mark(u32 n,std::vector<std::uint8_t> &permitted) {
        if(n%5==0) return;
        for(u32 start:{n%5,(5-n%5)%5})
            for(u32 u=start;u<=n;u+=5)
                if(!accepts_one(n,u)) permitted[u]=0;
    }
};

// Exact necessary condition over Z/MZ for odd prime powers M. It is used
// only when n is a unit. Scale the whole square by (2*n)^-1: the side
// becomes 1 and adjacent radius differences become u/n and v/n.
// All signed roots, including zero and nonunit roots, are retained.
// No assumptions about an interior point are used.
template<unsigned M, unsigned P>
struct LiftedProjection {
    std::array<std::uint8_t,M*M> normal{};
    std::array<std::uint16_t,M> inverse{};
    std::array<std::uint8_t,M> qr{};
    LiftedProjection() {
        static_assert(M % P == 0 && M % 2 == 1, "odd prime power required");
        std::array<std::vector<unsigned>,M> roots;
        for(unsigned a=0;a<M;++a) {
            roots[a*a%M].push_back(a);
            qr[a*a%M]=1;
            if(a%P) inverse[a]=std::uint16_t(inverse_mod(a,M));
        }
        for(unsigned x=0;x<M;++x) for(unsigned y=0;y<M;++y) {
            const unsigned xx=(1+M-x)%M, yy=(1+M-y)%M;
            const unsigned aa=(x*x+y*y)%M, bb=(xx*xx+y*y)%M;
            const unsigned cc=(x*x+yy*yy)%M, dd=(xx*xx+yy*yy)%M;
            if(!qr[dd]) continue;
            for(unsigned a:roots[aa]) for(unsigned b:roots[bb]) for(unsigned c:roots[cc])
                normal[((b+M-a)%M)*M+(c+M-a)%M]=1;
        }
    }
    bool accepts(u32 n,u32 u,u32 v) const {
        const unsigned inv=inverse[n%M];
        if(!inv) return true; // Nonunits need a different primitive local model.
        return normal[((u%M)*inv%M)*M+((v%M)*inv%M)];
    }
    unsigned size() const {
        return std::count(normal.begin(),normal.end(),std::uint8_t(1));
    }
    void verify() const {
        // Independent algebraic elimination: choose the first signed radius,
        // reconstruct x,y from the two differences, then check both norms.
        const unsigned inv2=(M+1)/2;
        auto mod=[](int value){value%=int(M);return unsigned(value<0?value+M:value);};
        for(unsigned u=0;u<M;++u) for(unsigned v=0;v<M;++v) {
            bool exists=false;
            for(unsigned a=0;a<M && !exists;++a) {
                const unsigned x=mod(1-int(u*u)-int(2*a*u))*inv2%M;
                const unsigned y=mod(1-int(v*v)-int(2*a*v))*inv2%M;
                const unsigned xx=(1+M-x)%M,yy=(1+M-y)%M;
                exists=(x*x+y*y)%M==a*a%M && qr[(xx*xx+yy*yy)%M];
            }
            if(exists!=bool(normal[u*M+v]))
                throw std::logic_error("lifted projection differs from radius reconstruction");
        }
    }
};
