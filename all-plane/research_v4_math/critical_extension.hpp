#ifndef RATIONAL_SQUARE_CRITICAL_EXTENSION_HPP
#define RATIONAL_SQUARE_CRITICAL_EXTENSION_HPP

// Whole-plane necessary conditions, proved in critical_lift.md and
// critical_sparse.md.  Include after u32, u64, and inverse_mod are defined.
// No reconstruction or global search bounds are changed by this header.
#include <array>
#include <cstdint>
#include <vector>

struct CriticalExtension {
    struct Lift {u32 residue,root;};
    struct PrimeData {
        u32 prime=0;
        std::array<std::uint8_t,282> allowed{};
        std::array<Lift,8> lifts{};
        u32 lift_count=0;
    };
    struct Active {const PrimeData* table;u32 mid;};
    std::vector<Active> active;

    static u32 quartic(u32 p,u32 iota,u32 r) {
        const u32 t=u64(r)*(r+p-1)%p;
        const u32 left=u64((2*r+p-1)%p)*((2*t+p-1)%p)%p;
        const u32 right=u64(2)*iota%p*t%p*(t+2)%p;
        return (left+p-right)%p;
    }

    static const std::vector<PrimeData>& tables() {
        static const std::vector<PrimeData> data=[] {
            std::vector<PrimeData> answer;
            for(u32 p=17;p<=281;p+=4) {
                bool prime=true;
                for(u32 d=3;d*d<=p;d+=2)
                    if(p%d==0) {prime=false;break;}
                if(!prime) continue;
                PrimeData table;table.prime=p;
                std::array<std::uint8_t,282> squares{};
                u32 iota=0;
                for(u32 a=0;a<p;++a) {
                    const u32 square=a*a%p;
                    squares[square]=1;
                    if(square==p-1 && !iota) iota=a;
                }
                for(u32 r=2;r<p;++r) {
                    const u32 denominator=(r+p-1+p-u64(iota)*r%p)%p;
                    if(denominator && squares[quartic(p,iota,r)]) {
                        table.allowed[r]=1;
                        // F_i(1-r)=-F_{-i}(r), and -1 is a square.
                        table.allowed[(1+p-r)%p]=1;
                    }
                }
                if(p==41) {
                    table.lifts={{{7,48},{12,750},{15,1327},{27,355},
                                  {30,932},{35,1634},{0,0},{0,0}}};
                    table.lift_count=6;
                } else if(p==113) {
                    table.lifts={{{9,10518},{50,3892},{64,8878},{105,2252},
                                  {0,0},{0,0},{0,0},{0,0}}};
                    table.lift_count=4;
                } else if(p==137) {
                    table.lifts={{{43,15113},{55,9234},{83,9536},{95,3657},
                                  {0,0},{0,0},{0,0},{0,0}}};
                    table.lift_count=4;
                } else if(p==257) {
                    table.lifts={{{23,36260},{62,44523},{117,3715},{141,62335},
                                  {196,21527},{235,29790},{0,0},{0,0}}};
                    table.lift_count=6;
                }
                answer.push_back(table);
            }
            return answer;
        }();
        return data;
    }

    void prepare(u32 n) {
        active.clear();
        if(!n) return;
        for(const auto &table:tables()) {
            const u32 p=table.prime;
            if(n%(p*p)) continue;
            u32 remaining=n,e=0,mid=1;
            do {remaining/=p;++e;} while(remaining%p==0);
            if(e&1) continue;
            for(u32 j=0;j<e/2;++j) mid*=p;
            active.push_back({&table,mid});
        }
    }

    bool accepts(u32 u,u32 h) const {
        for(const auto &item:active) {
            const auto &table=*item.table;
            const u32 p=table.prime,mid=item.mid;
            if(h%mid || (h/mid)%p==0) continue; // Only v_p(h)=s.
            if(u%mid) return false; // Actual critical radii force this.
            const u32 U=u/mid,H=h/mid;
            const u32 residue=u64(U%p)*inverse_mod(H%p,p)%p;
            if(!table.allowed[residue]) return false;
            for(u32 j=0;j<table.lift_count;++j) {
                const Lift &lift=table.lifts[j];
                if(residue==lift.residue) {
                    const u32 modulus=p*p;
                    if(U%modulus!=u64(lift.root)*(H%modulus)%modulus) return false;
                    break;
                }
            }
        }
        return true;
    }
};

#endif
