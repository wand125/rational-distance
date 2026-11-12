#define main baseline_all_plane_main
#include "../all_plane.cpp"
#undef main
#include <vector>

static std::vector<std::pair<i64,i64>> factors(i64 n) {
    std::vector<std::pair<i64,i64>> f;
    for(i64 p=2;p*p<=n;++p) if(n%p==0) {
        i64 power=1;
        do { n/=p; power*=p; } while(n%p==0);
        f.push_back({p,power});
    }
    if(n>1) f.push_back({n,n});
    return f;
}

static bool unit_filter(i64 n,i64 u,i64 v) {
    for(auto [ell,power]:factors(n)) {
        if(ell!=2 && ell%4!=3) continue;
        if((u%ell==0 && u%power!=0) || (v%ell==0 && v%power!=0)) return false;
    }
    return true;
}
static bool scale_filter(i64 n,i64 u,i64 v) {
    i64 g=std::gcd(std::gcd(n,u),v);
    i64 n0=n/g,u0=u/g,v0=v/g;
    i64 Q0=u0*u0+v0*v0-n0*n0;
    if((n0*Q0)%g!=0) return false;
    if(g%2==0 && n0%2==1) {
        i64 two_part=1;
        while(g%2==0) {g/=2;two_part*=2;}
        if(Q0%(2*two_part)!=0) return false;
    }
    return true;
}

int main(int argc,char**argv) {
    i64 bound=argc>1 ? std::stoll(argv[1]) : 1000;
    Stats st;
    u64 primitive=0,large_g=0,unit_fail=0,scale_fail=0;
    for(i64 N=2;N<=bound;N+=2) {
        for(i64 p=0;p<=N;p+=2) for(i64 q=0;q<=p;q+=2) {
            evaluate(N,p,q,false,st,[&](i64 x,i64 y,i64 a) {
                if(std::gcd(std::gcd(magnitude(x),magnitude(y)),N)!=1) return;
                ++primitive;
                if(a%2==0) throw std::logic_error("primitive even-N A must be odd");
                auto n=N/2,u=p/2,v=q/2;
                large_g+=std::gcd(std::gcd(n,u),v)>1;
                if(!unit_filter(n,u,v)) ++unit_fail;
                if(!scale_filter(n,u,v)) ++scale_fail;
            });
        }
    }
    std::cout<<"{\"bound\":"<<bound<<",\"even_N_primitive_three_distance_points\":"<<primitive
             <<",\"gcd_nuv_greater_than_one\":"<<large_g
             <<",\"unit_filter_false_rejections\":"<<unit_fail
             <<",\"scale_filter_false_rejections\":"<<scale_fail<<"}\n";
    return unit_fail || scale_fail;
}
