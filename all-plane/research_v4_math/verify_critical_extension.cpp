#define RATIONAL_SQUARE_V3_NO_MAIN
#include "../all_plane_gap_v3.cpp"
#include "critical_extension.hpp"

// Independent mod-p model: no quartic or square root of -1 is used.
static std::vector<bool> direct_base(u32 p) {
    std::vector<bool> squares(p,false),allowed(p,false);
    for(u32 a=0;a<p;++a) squares[a*a%p]=true;
    for(u32 r=2;r<p;++r) for(u32 a=0;a<p;++a) {
        const u32 first=u64(4)*r%p*(a+r)%p;
        const u32 second=u64(4)*(r-1)%p*(a+r-1)%p;
        if(!first || !second || (first*first+second*second)%p) continue;
        const u32 fourth=(a*a+first+second)%p;
        if(squares[fourth]) allowed[r]=true;
    }
    return allowed;
}

// Independent eliminated-radius model modulo p^2; does not use F_i.
static std::vector<bool> radius_reference(u32 p) {
    const u32 m=p*p;
    std::vector<bool> squares(m,false),allowed(m,false);
    std::vector<u32> iotas;
    for(u32 a=0;a<m;++a) {
        const u32 square=u64(a)*a%m;
        squares[square]=true;
        if(square==m-1) iotas.push_back(a);
    }
    if(iotas.size()!=2) throw std::runtime_error("bad roots of -1");
    for(u32 r=0;r<m;++r) {
        if(r%p==0 || r%p==1) continue;
        const u32 u=r,v=(r+m-1)%m;
        for(u32 i:iotas) {
            const u32 den=(v+m-u64(i)*u%m)%m;
            if(den%p==0) continue;
            const u32 num=(u64(i)*u%m*u%m+m-u64(v)*v%m)%m;
            const u32 a=u64(num)*inverse_mod(den,m)%m;
            const u32 b=(a+2*u)%m,c=(a+2*v)%m;
            const u32 fourth=(u64(b)*b%m+u64(c)*c%m+m-u64(a)*a%m)%m;
            if(squares[fourth]) allowed[r]=true;
        }
    }
    return allowed;
}

int main() try {
    u64 basic_checked=0,lift_checked=0;
    std::vector<std::array<u32,3>> counts;
    for(const auto &table:CriticalExtension::tables()) {
        const u32 p=table.prime;
        const auto direct=direct_base(p);
        for(u32 r=0;r<p;++r) {
            if(bool(table.allowed[r])!=direct[r])
                throw std::runtime_error("base mismatch p="+std::to_string(p));
            ++basic_checked;
        }
        CriticalExtension extension;extension.prepare(12*p*p);
        if(extension.active.size()!=1) throw std::runtime_error("bad active plan");
        if(!extension.accepts(1,1) || !extension.accepts(1,p*p))
            throw std::runtime_error("noncritical h was filtered");
        if(table.lift_count) {
            const auto reference=radius_reference(p);
            u32 old_count=0,new_count=0;
            for(u32 r=0;r<p*p;++r) {
                const bool accepted=extension.accepts(p*r,p);
                if(accepted!=reference[r])
                    throw std::runtime_error("p^2 mismatch p="+std::to_string(p)+",r="+std::to_string(r));
                old_count+=table.allowed[r%p];new_count+=accepted;++lift_checked;
            }
            counts.push_back({p,old_count,new_count});
        }
        extension.prepare(12*p*p*p);
        if(!extension.active.empty() || !extension.accepts(1,p))
            throw std::runtime_error("odd denominator exponent was filtered");
    }
    constexpr u32 bound=100000;
    ArithmeticTables arithmetic(bound/2);LocalSieve local;GapWheel wheel(local);
    ProgressionWheel<19> progression(local);
    SearchOptions options;options.bound=bound;
    CriticalExtension extension;u32 last_n=0;u64 parameters=0,rejected=0;
    std::vector<std::array<u32,3>> examples;
    run_gap_search(options,arithmetic,local,wheel,[&](u32 N,u32 twice_u,u32 twice_v) {
        const u32 n=N/2,u=twice_u/2,v=twice_v/2;
        if(n!=last_n) {extension.prepare(n);last_n=n;}
        ++parameters;
        if(!extension.accepts(u,u-v)) {
            ++rejected;
            if(examples.size()<20) examples.push_back({n,u,v});
        }
    },[](i64,i64,i64,i64,i64,i64,Stats&){},false,true,true,&progression,true);
    std::cout<<"{\"basic_prime_count\":"<<CriticalExtension::tables().size()
             <<",\"basic_residues_checked\":"<<basic_checked
             <<",\"lift_residues_checked\":"<<lift_checked
             <<",\"independent_models_match\":true,\"bound\":"<<bound
             <<",\"v3_parameters\":"<<parameters<<",\"new_rejected\":"<<rejected
             <<",\"lift_counts\":[";
    for(std::size_t j=0;j<counts.size();++j) {
        if(j) std::cout<<',';
        std::cout<<"{\"p\":"<<counts[j][0]<<",\"old\":"<<counts[j][1]
                 <<",\"new\":"<<counts[j][2]<<'}';
    }
    std::cout<<"],\"examples\":[";
    for(std::size_t j=0;j<examples.size();++j) {
        if(j) std::cout<<',';
        std::cout<<"{\"n\":"<<examples[j][0]<<",\"u\":"<<examples[j][1]
                 <<",\"v\":"<<examples[j][2]<<'}';
    }
    std::cout<<"]}\n";
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
