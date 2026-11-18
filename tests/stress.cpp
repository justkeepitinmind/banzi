#include "../lichao.cpp"
#include "../gp_hash_table.cpp"
#include "../scc.cpp"
#include "../bcc.cpp"
#include "../flow.cpp"
#include "../hlpp.cpp"
#include "../ntt.cpp"

mt19937_64 rng(20261116);
int rnd(int l,int r){return l+rng()%(r-l+1);}

void test_lines(){
    for(int tc=0;tc<1200;++tc){
        vector<ll> x;
        for(int i=0;i<20;++i)x.push_back(rnd(-40,40));
        LiChao lc(x);
        LineContainer hull;
        struct Seg{ll l,r,k,b;};
        vector<Seg> a;
        vector<pair<ll,ll>> lines;
        for(int it=0;it<60;++it){
            ll k = rnd(-100,100),b = rnd(-100,100);
            hull.add(k,b),lines.push_back({k,b});
            if(rnd(0,1))lc.add(k,b),a.push_back({LLONG_MIN,LLONG_MAX,k,b});
            else{
                ll l = rnd(-60,60),r = rnd(-60,60);
                if(l > r)swap(l,r);
                lc.add_segment(l,r,k,b),a.push_back({l,r,k,b});
            }
            for(ll q:x){
                i128 mn = LiChao::inf,mx = -LiChao::inf;
                for(auto s:a)if(s.l <= q && q <= s.r)mn = min(mn,(i128)s.k*q+s.b);
                for(auto [kk,bb]:lines)mx = max(mx,(i128)kk*q+bb);
                assert(lc.query(q) == mn && hull.query(q) == mx);
            }
        }
    }
    // ll 安全范围：大整数、相同斜率、负交点向下取整。
    static_assert(is_same_v<decltype(LineContainer().query(0)),ll>);
    vector<ll> v = {-1000000000,-1,0,1,1000000000};
    LineContainer h;
    vector<pair<ll,ll>> lines;
    for(ll k:v)for(ll b:v){
        h.add(k,b),lines.push_back({k,b});
        for(ll x:v){
            i128 ans = -LiChao::inf;
            for(auto [kk,bb]:lines)ans = max(ans,(i128)kk*x+bb);
            assert(h.query(x) == ans);
        }
    }
    LineContainer wide;
    wide.add(-1,-4000000000000000000LL);
    wide.add(0,4000000000000000000LL);
    wide.add(1,-4000000000000000000LL);
    for(ll x:v)assert(wide.query(x) == 4000000000000000000LL);
    LineContainer near_limit;
    near_limit.add(1,LLONG_MAX-2);
    assert(near_limit.query(1) == LLONG_MAX-1);
    LiChao one({-1000000000000000000LL,1000000000000000000LL});
    assert(one.query(1000000000000000000LL) == LiChao::inf);
    one.add(1000000000000000000LL,-1000000000000000000LL);
    assert(one.query(1000000000000000000LL) == (i128)1000000000000000000LL*1000000000000000000LL-1000000000000000000LL);
    LiChao single({7,7});single.add_segment(8,9,1,2);
    assert(single.query(7) == LiChao::inf);
    single.add_segment(7,7,3,4);assert(single.query(7) == 25);
    cout << "lines: 1200 random cases + integer boundaries OK\n";
}

void test_hash(){
    hashmap h({},{},{},{},{1<<16});
    map<ll,int> a;
    for(int i=0;i<100000;++i){
        ll x = rnd(-10000,10000);
        if(rnd(0,2))++h[x],++a[x];
        else h.erase(x),a.erase(x);
        auto p = h.find(x);
        assert((p == h.end()) == (a.find(x) == a.end()));
        if(p != h.end())assert(p->second == a[x]);
        assert(h.size() == a.size());
    }
    for(auto [x,y]:a)assert(h.find(x)->second == y);
    gp_hash_table<ull,int,chash> pairs;
    for(int x=-5;x<=5;++x)for(int y=-5;y<=5;++y)pairs[pack(x,y)] = x*20+y;
    assert(pairs.size() == 121);
    h.clear();assert(h.empty());
    cout << "hash: 100000 operations OK\n";
}

int components(int n,const vector<pair<int,int>>& edges,int removed_v=0,int removed_e=-1,int mask=-1){
    vector<vector<int>> g(n+1);
    for(int i=0;i<(int)edges.size();++i)if(i != removed_e){
        auto [u,v] = edges[i];
        g[u].push_back(v),g[v].push_back(u);
    }
    int cnt = 0;
    vector<char> vis(n+1);
    for(int u=1;u<=n;++u)if(u != removed_v && (mask == -1 || (mask>>(u-1)&1)) && !vis[u]){
        ++cnt;
        vector<int> st = {u};vis[u] = 1;
        while(!st.empty()){
            int x = st.back();st.pop_back();
            for(int v:g[x])if(v != removed_v && (mask == -1 || (mask>>(v-1)&1)) && !vis[v])
                vis[v] = 1,st.push_back(v);
        }
    }
    return cnt;
}

void test_graphs(){
    for(int tc=0;tc<1800;++tc){
        int n = rnd(1,8),m = rnd(0,20);
        SCC s(n);VBCC v(n);EBCC e(n);
        vector<pair<int,int>> edges;
        vector<vector<int>> reach(n+1,vector<int>(n+1));
        for(int u=1;u<=n;++u)reach[u][u] = 1;
        for(int i=0;i<m;++i){
            int x = rnd(1,n),y = rnd(1,n);
            edges.push_back({x,y});
            s.add(x,y),v.add(x,y),e.add(x,y),reach[x][y] = 1;
        }
        s.run(),v.run(),e.run();
        for(int k=1;k<=n;++k)for(int i=1;i<=n;++i)for(int j=1;j<=n;++j)
            reach[i][j] |= reach[i][k]&reach[k][j];
        for(int i=1;i<=n;++i)for(int j=1;j<=n;++j)
            assert((s.bel[i] == s.bel[j]) == (reach[i][j] && reach[j][i]));
        for(auto [x,y]:edges)if(s.bel[x] != s.bel[y])assert(s.bel[x] > s.bel[y]);
        int base = components(n,edges);
        for(int i=1;i<=n;++i)assert(bool(v.cut[i]) == (components(n,edges,i) > base));
        for(int i=0;i<m;++i)assert(bool(e.bridge[i]) == (components(n,edges,0,i) > base));
        vector<pair<int,int>> rest;
        for(int i=0;i<m;++i)if(!e.bridge[i])rest.push_back(edges[i]);
        assert(components(n,rest) == e.cnt);
        for(auto [x,y]:rest)assert(e.bel[x] == e.bel[y]);
        // 枚举所有诱导子图，选连通且无割点的极大顶点集合，独立验证点双。
        vector<int> good,expected,actual;
        for(int mask=1;mask<(1<<n);++mask){
            if(components(n,edges,0,-1,mask) != 1)continue;
            bool ok = true;
            for(int u=1;u<=n;++u)if(mask>>(u-1)&1)
                if(components(n,edges,u,-1,mask) > 1)ok = false;
            if(ok)good.push_back(mask);
        }
        for(int mask:good){
            bool maximal = true;
            for(int other:good)if(mask != other && (mask&other) == mask)maximal = false;
            if(maximal)expected.push_back(mask);
        }
        for(auto block:v.bcc){
            int mask = 0;
            for(int u:block){assert(!(mask>>(u-1)&1));mask |= 1<<(u-1);}
            actual.push_back(mask);
        }
        sort(expected.begin(),expected.end()),sort(actual.begin(),actual.end());
        assert(expected == actual);
        auto tr = v.block_tree();
        int N = tr.size()-1,cc = 0,twice_edges = 0;
        vector<int> seen(N+1);
        for(int u=1;u<=N;++u){
            twice_edges += tr[u].size();
            for(int w:tr[u])assert((u <= n) != (w <= n));
            if(seen[u])continue;
            ++cc;
            vector<int> st = {u};seen[u] = cc;
            while(!st.empty()){
                int x = st.back();st.pop_back();
                for(int w:tr[x])if(!seen[w])seen[w] = cc,st.push_back(w);
            }
        }
        assert(cc == base && twice_edges == 2*(N-cc));
        int bridge_edges = 0;
        for(auto list:e.bridge_tree())bridge_edges += list.size();
        assert(bridge_edges == 2*(e.cnt-base));
        auto old_s = s.bel,old_e = e.bel;
        s.run(),e.run(),v.run(); // 多次 run 必须重置状态。
        assert(s.bel == old_s && e.bel == old_e && v.bcc.size() == actual.size());
    }
    cout << "SCC/BCC/block forest: 1800 random multigraphs OK\n";
}

void test_dinic(){
    for(int tc=0;tc<2200;++tc){
        int n = rnd(2,8),m = rnd(0,30);
        Dinic d(n);
        vector<tuple<int,int,ll>> edges;
        for(int i=0;i<m;++i){
            int u = rnd(1,n),v = rnd(1,n);ll c = rnd(0,12);
            d.add(u,v,c),edges.push_back({u,v,c});
        }
        ll best = LLONG_MAX;
        for(int mask=0;mask<(1<<n);++mask)if((mask&1) && !(mask>>(n-1)&1)){
            ll cut = 0;
            for(auto [u,v,c]:edges)if((mask>>(u-1)&1) && !(mask>>(v-1)&1))cut += c;
            best = min(best,cut);
        }
        ll limit = rnd(0,20),first = d.run(1,n,limit),second = d.run(1,n);
        assert(first == min(limit,best) && first+second == best);
        ll cut = 0;
        for(auto [u,v,c]:edges)if(d.dep[u] != -1 && d.dep[v] == -1)cut += c;
        assert(cut == best);
    }
    Dinic big(3);
    big.add(1,2,3000000000000000000LL),big.add(2,3,3000000000000000000LL);
    assert(big.run(1,3) == 3000000000000000000LL);
    cout << "Dinic: 2200 exhaustive min-cut comparisons + 64-bit capacity OK\n";
}

void test_hlpp(){
    for(int tc=0;tc<4500;++tc){
        int n = rnd(2,9),m = rnd(0,50),s = rnd(1,n),t = rnd(1,n);
        while(s == t)t = rnd(1,n);
        HLPP f(n);
        vector<tuple<int,int,ll>> edges;
        vector<int> ids;
        for(int i=0;i<m;++i){
            int u = rnd(1,n),v = rnd(1,n);ll c = rnd(0,30);
            ids.push_back(f.add(u,v,c)),edges.push_back({u,v,c});
        }
        ll want = LLONG_MAX;
        for(int mask=0;mask<(1<<n);++mask)if((mask>>(s-1)&1) && !(mask>>(t-1)&1)){
            ll cut = 0;
            for(auto [u,v,c]:edges)if((mask>>(u-1)&1) && !(mask>>(v-1)&1))cut += c;
            want = min(want,cut);
        }
        assert(f.run(s,t) == want);
        ll cut = 0;
        vector<ll> balance(n+1);
        for(int i=0;i<m;++i){
            auto [u,v,c] = edges[i];
            if(f.source_side(u) && !f.source_side(v))cut += c;
            int id = ids[i];
            if(u == v){assert(id == -1);continue;}
            ll flow = f.e[id^1].cap;
            assert(0 <= flow && flow <= c && f.e[id].cap+flow == c);
            balance[u] -= flow,balance[v] += flow;
        }
        assert(cut == want && f.source_side(s) && !f.source_side(t));
        for(int u=1;u<=n;++u){
            ll expected = u == s ? -want : (u == t ? want : 0);
            assert(balance[u] == expected && f.ex[u] == expected);
        }
    }
    for(int tc=0;tc<500;++tc){
        int n = rnd(10,150),m = rnd(0,4000),s = rnd(1,n),t = rnd(1,n);
        while(s == t)t = rnd(1,n);
        HLPP f(n);Dinic d(n);
        for(int i=0;i<m;++i){
            int u = rnd(1,n),v = rnd(1,n);ll cap = rnd(0,1000000000);
            f.add(u,v,cap),d.add(u,v,cap);
        }
        assert(f.run(s,t) == d.run(s,t));
    }
    HLPP big(2);
    big.add(1,2,LLONG_MAX);
    assert(big.run(1,2) == LLONG_MAX);
    HLPP disconnected(4);
    disconnected.add(1,2,LLONG_MAX),disconnected.add(2,3,LLONG_MAX);
    assert(disconnected.run(1,4) == 0);
    assert(disconnected.ex[1] == 0 && disconnected.ex[2] == 0 && disconnected.ex[3] == 0);
    cout << "HLPP: 4500 exhaustive min-cut + 500 Dinic comparisons + 64-bit boundaries OK\n";
}

struct CostEdge{int u,v,cap,cost;};
void enumerate_flow(int i,int n,const vector<CostEdge>& e,vector<int>& bal,ll cost,ll limit,pair<ll,ll>& best){
    if(i == (int)e.size()){
        if(bal[1] < 0 || bal[1] > limit || bal[n] != -bal[1])return;
        for(int v=2;v<n;++v)if(bal[v])return;
        if(bal[1] > best.first || (bal[1] == best.first && cost < best.second))best = {bal[1],cost};
        return;
    }
    auto a = e[i];
    for(int f=0;f<=a.cap;++f){
        bal[a.u] += f,bal[a.v] -= f;
        enumerate_flow(i+1,n,e,bal,cost+(ll)f*a.cost,limit,best);
        bal[a.u] -= f,bal[a.v] += f;
    }
}
void test_mcmf(){
    for(int tc=0;tc<2400;++tc){
        int n = rnd(2,6),m = rnd(0,9);
        vector<int> p(n+1),bal(n+1);
        for(int& x:p)x = rnd(-5,5);
        vector<CostEdge> edges;
        MCMF flow(n);
        for(int i=0;i<m;++i){
            int u = rnd(1,n),v = rnd(1,n),cap = rnd(0,2);
            int cost = rnd(0,5)+p[v]-p[u]; // 保证有负边但无负环。
            if(tc%2 == 0)cost = rnd(0,8); // 也测完全非负、不需初始化的分支。
            edges.push_back({u,v,cap,cost}),flow.add(u,v,cap,cost);
        }
        ll limit = rnd(0,6);
        pair<ll,ll> best = {0,0};
        enumerate_flow(0,n,edges,bal,0,limit,best);
        auto copy = flow;
        auto ans = flow.run(1,n,limit);
        assert(ans == best);
        auto a = copy.run(1,n,min(1LL,limit));
        auto b = copy.run(1,n,limit-a.first);
        assert(a.first+b.first == best.first && a.second+b.second == best.second);
    }
    cout << "MCMF: 2400 exhaustive integral-flow comparisons OK\n";
}

void test_ntt(){
    using namespace NTT;
    for(int tc=0;tc<2500;++tc){
        int n = rnd(0,160),m = rnd(0,160);
        vector<int> a(n),b(m),want(n && m ? n+m-1 : 0);
        for(int& x:a)x = (int)rng();
        for(int& x:b)x = (int)rng();
        for(int i=0;i<n;++i)for(int j=0;j<m;++j){
            ll x = (a[i]%mod+mod)%mod,y = (b[j]%mod+mod)%mod;
            want[i+j] = (want[i+j]+x*y)%mod;
        }
        assert(multiply(a,b) == want);
    }
    // 变长顺序调用检查根缓存，覆盖长度 1 和 inverse。
    for(int k:{0,1,8,3,15,2,17,0}){
        vector<int> a(1<<k);
        for(int& x:a)x = rnd(0,mod-1);
        auto copy = a;
        ntt(a),ntt(a,true);assert(a == copy);
    }
    vector<int> a(1<<17,1),b(1<<17,1);
    auto c = multiply(a,b);
    for(int i=0;i<(int)c.size();++i)assert(c[i] == min(i+1,(int)c.size()-i));
    cout << "NTT: 2500 naive comparisons + round trips + large convolution OK\n";
}

int main(){
    test_lines();test_hash();test_graphs();test_dinic();test_hlpp();test_mcmf();test_ntt();
    cout << "ALL OK\n";
}
