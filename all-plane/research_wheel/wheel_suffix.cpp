// Append after the frozen v2 definitions, excluding v2's main function.
template<unsigned W> struct ProgressionWheel {
    std::vector<std::uint8_t> next;
    explicit ProgressionWheel(const LocalSieve &local):next(std::size_t(W)*W*W*W,255) {
        for(unsigned n=0;n<W;++n) for(unsigned h=0;h<W;++h) {
            std::array<bool,W> allowed{};
            for(unsigned u=0;u<W;++u) {
                const unsigned v=(u+W-h)%W;
                bool yes=true;
                if constexpr(W%5==0) yes=yes && local.p5.accepts(n,u,v);
                if constexpr(W%11==0) yes=yes && local.p11.accepts(n,u,v);
                if constexpr(W%19==0) yes=yes && local.p19.accepts(n,u,v);
                allowed[u]=yes;
            }
            for(unsigned step=0;step<W;++step) {
                auto *row=next.data()+(((std::size_t(n)*W+h)*W+step)*W);
                const unsigned cycles=std::gcd(step,W),len=W/cycles;
                for(unsigned first=0;first<cycles;++first) {
                    int nearest=-1;
                    for(int j=int(2*len)-1;j>=0;--j) {
                        const unsigned u=(first+unsigned(j)*step)%W;
                        if(allowed[u]) nearest=j;
                        if(j<int(len)) row[u]=nearest<0?255:std::uint8_t(nearest-j);
                    }
                }
            }
        }
    }
    const std::uint8_t *row(u32 n,u32 h,u32 step) const {
        return next.data()+(((std::size_t(n%W)*W+h%W)*W+step%W)*W);
    }
};

template<unsigned W,class Callback>
static u64 generate_field_gap_crt(u32 n,u32 h,const ArithmeticTables &tables,
                                 const ProgressionWheel<W> &wheel,Callback emit,u32 known_gcd=0) {
    const u32 m=n/(known_gcd?known_gcd:std::gcd(n,h));
    const auto &plan=tables.plans[m];
    std::array<u32,128> residues{};
    unsigned count=1;
    for(unsigned i=0;i<plan.count;++i) {
        if(h%plan.power[i]==0) continue;
        const u32 shift=u32(u64(h)*plan.idempotent[i]%m);
        const unsigned old=count;
        if(2*old>residues.size()) throw std::logic_error("CRT capacity");
        for(unsigned j=0;j<old;++j) {
            u32 r=residues[j]+shift;
            if(r>=m) r-=m;
            residues[count++]=r;
        }
    }
    const auto *row=wheel.row(n,h,m);
    for(unsigned i=0;i<count;++i) {
        u32 first=residues[i];
        if(first<h) first+=((h-first+m-1)/m)*m;
        u64 u=first;
        while(u<=n) {
            const unsigned skip=row[u%W];
            if(skip==255) break;
            u+=u64(skip)*m;
            if(u>n) break;
            emit(u32(u));
            u+=m;
        }
    }
    return count;
}

template<unsigned W,class PointCallback>
static GapStats run_field_gap(const SearchOptions &opt,const ArithmeticTables &arithmetic,
                             const LocalSieve &local,const GapWheel &gapwheel,
                             const ProgressionWheel<W> &wheel,PointCallback accept) {
    GapStats out;
    const u32 first_n=12*((opt.min_denominator+23)/24);
    for(u32 n=first_n;n<=opt.bound/2;n+=12) {
        const u32 N=2*n;
        ++out.exact.sides;
        out.dense_even_pairs+=u64(n+1)*(n+2)/2;
        out.gap_rows+=n+1;
        auto permitted=permitted_valuations(n,arithmetic.plans[n]);
        gapwheel.generate(n,[&](u32 h){
            ++out.local_gap_rows;
            if(!permitted_difference(n,h,0,arithmetic)) return;
            ++out.difference_gap_rows;
            const u32 gcd_nh=std::gcd(n,h);
            out.crt_residues+=generate_field_gap_crt(n,h,arithmetic,wheel,[&](u32 u) {
                const u32 v=u-h;
                ++out.generated_pairs;
                if(!permitted[u] || !permitted[v]) return;
                ++out.valuation_pairs;
                if(!local.accepts(n,u,v)) return;
                ++out.local_pairs;
                if(!permitted_scale(n,u,v,gcd_nh)) return;
                ++out.scale_pairs;
                evaluate(N,2*u,2*v,true,out.exact,[&](i64 x,i64 y,i64 a){
                    accept(N,2*u,2*v,x,y,a,out.exact);
                });
            },gcd_nh);
        });
    }
    return out;
}

template<unsigned W>
static void validate_field_wheel(const ArithmeticTables &arithmetic,const LocalSieve &local,
                                 const ProgressionWheel<W> &wheel) {
    u64 total=0;
    for(u32 n=1;n<=128;++n) for(u32 h=0;h<=n;++h) {
        std::vector<u32> actual,expected;
        generate_field_gap_crt(n,h,arithmetic,wheel,[&](u32 u){actual.push_back(u);});
        generate_gap_crt(n,h,arithmetic,[&](u32 u){
            bool yes=true;
            if constexpr(W%5==0) yes=yes && local.p5.accepts(n,u,u-h);
            if constexpr(W%11==0) yes=yes && local.p11.accepts(n,u,u-h);
            if constexpr(W%19==0) yes=yes && local.p19.accepts(n,u,u-h);
            if(yes) expected.push_back(u);
        });
        std::sort(actual.begin(),actual.end());std::sort(expected.begin(),expected.end());
        if(actual!=expected) throw std::logic_error("progression wheel mismatch");
        total+=actual.size();
    }
    std::cerr<<"Wheel validation passed: "<<total<<" emitted CRT members.\n";
}

int main(int argc,char**argv) try {
    SearchOptions opt;unsigned mode=55;bool check=false;std::string auditpath;
    for(int i=1;i<argc;++i) {
        std::string arg=argv[i];
        if(arg=="--wheel") {if(++i==argc) throw std::runtime_error("missing wheel");mode=std::stoul(argv[i]);}
        else if(arg=="--validate") check=true;
        else if(arg=="--audit-three") {if(++i==argc) throw std::runtime_error("missing audit");auditpath=argv[i];}
        else {opt.bound=std::stoul(arg);}
    }
    if(opt.bound<1 || opt.bound>2000000) throw std::runtime_error("bound1..2000000");
    const auto start=std::chrono::steady_clock::now();
    ArithmeticTables arithmetic(std::max<u32>(1000,opt.bound/2));LocalSieve local;GapWheel gapwheel(local);
    std::ofstream audit;
    if(!auditpath.empty()){audit.open(auditpath);if(!audit) throw std::runtime_error("audit open");audit<<"N,p,q,X,Y,A\n";}
    auto point=[&](i64 N,i64 p,i64 q,i64 x,i64 y,i64 a,Stats &stats){
        if(audit.is_open()) audit<<N<<','<<p<<','<<q<<','<<x<<','<<y<<','<<a<<'\n';
        u64 d;if(square(u128(i128(N-x)*(N-x)+i128(N-y)*(N-y)),d)) {
            ++stats.solutions;
            std::cout<<"{\"solution\":{\"N\":"<<N<<",\"X\":"<<x<<",\"Y\":"<<y<<",\"A\":"<<a<<",\"D\":"<<d<<"}}\n";
        }
    };
    GapStats stats;auto search_start=start;
    auto run=[&](auto &wheel){
        if(check) validate_field_wheel(arithmetic,local,wheel);
        search_start=std::chrono::steady_clock::now();
        stats=run_field_gap(opt,arithmetic,local,gapwheel,wheel,point);
    };
    if(mode==0){search_start=std::chrono::steady_clock::now();stats=run_gap_search(opt,arithmetic,local,gapwheel,[](u32,u32,u32){},point);}
    else if(mode==11){ProgressionWheel<11> wheel(local);run(wheel);}
    else if(mode==19){ProgressionWheel<19> wheel(local);run(wheel);}
    else if(mode==55){ProgressionWheel<55> wheel(local);run(wheel);}
    else throw std::runtime_error("wheel must be0,11,19,55");
    const auto end=std::chrono::steady_clock::now();
    std::cout<<"{\"bound\":"<<opt.bound<<",\"wheel_modulus\":"<<mode
             <<",\"generated_pairs\":"<<stats.generated_pairs<<",\"valuation_pairs\":"<<stats.valuation_pairs
             <<",\"local_pairs\":"<<stats.local_pairs<<",\"exact_parameter_pairs\":"<<stats.exact.parameter_pairs
             <<",\"primitive_three_candidates\":"<<stats.exact.primitive_three_candidates
             <<",\"solutions\":"<<stats.exact.solutions
             <<",\"seconds_search\":"<<std::chrono::duration<double>(end-search_start).count()
             <<",\"seconds_total\":"<<std::chrono::duration<double>(end-start).count()<<"}\n";
} catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
