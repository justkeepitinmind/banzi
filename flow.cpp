#include<bits/stdc++.h>
using namespace std;
using ll = long long;

// 点 1..n；有向边 add 一次。反向边自动添加，正反边下标 i 和 i^1。
// 容量非负，允许重边、自环；s != t。容量和、总流量必须装进 ll。
// add 返回正向边下标；e[id^1].cap 是该条边当前流量。
// run 返回本次新增流量，修改残量网络。新题/新测试重新构造对象。
struct Dinic{
    struct Edge{int to;ll cap;};
    int n;
    vector<Edge> e;
    vector<vector<int>> g;
    vector<int> dep,cur,q;
    Dinic(int n):n(n),g(n+1),dep(n+1),cur(n+1),q(n+1){}
    int add(int u,int v,ll cap){
        assert(cap >= 0);
        int id = e.size();
        g[u].push_back(id),e.push_back({v,cap});
        g[v].push_back(id^1),e.push_back({u,0});
        return id;
    }
    bool bfs(int s,int t){
        fill(dep.begin(),dep.end(),-1);
        int l = 0,r = 0;
        q[r++] = s,dep[s] = 0;
        while(l < r){
            int u = q[l++];
            for(int id:g[u])if(e[id].cap && dep[e[id].to] == -1){
                int v = e[id].to;
                dep[v] = dep[u]+1,q[r++] = v;
            }
        }
        return dep[t] != -1;
    }
    ll dfs(int u,int t,ll f){
        if(u == t)return f;
        for(int& i=cur[u];i<(int)g[u].size();++i){
            int id = g[u][i],v = e[id].to;
            if(!e[id].cap || dep[v] != dep[u]+1)continue;
            ll z = dfs(v,t,min(f,e[id].cap));
            if(z){e[id].cap -= z,e[id^1].cap += z;return z;}
        }
        dep[u] = -1; // 本轮层次图中已经无法增广。
        return 0;
    }
    ll run(int s,int t,ll limit = LLONG_MAX){
        assert(s != t && limit >= 0);
        ll flow = 0;
        while(flow < limit && bfs(s,t)){
            fill(cur.begin(),cur.end(),0);
            while(flow < limit){
                ll z = dfs(s,t,limit-flow);
                if(!z)break;
                flow += z;
            }
        }
        return flow;
    }
}; // 一般图 O(n^2 m)；DFS 最大递归深度 O(n)。

// 最小费用最大流：一次初始 SPFA 求势能，之后每次 Dijkstra 增广。
// 不存在负费用残量边时跳过 SPFA。支持负边，不支持源点可达的负费用环。
// run 返回 {本次新增流量,本次新增费用}；limit 可指定需求流量。
// 容量/总流量装进 ll；费用、距离、势能及其运算中间值绝对值须小于 inf。
struct MCMF{
    static constexpr ll inf = 4000000000000000000LL;
    struct Edge{int to;ll cap,cost;};
    int n;
    vector<Edge> e;
    vector<vector<int>> g;
    vector<ll> h,dis;
    vector<int> pre;
    MCMF(int n):n(n),g(n+1),h(n+1),dis(n+1),pre(n+1){}
    int add(int u,int v,ll cap,ll cost){
        assert(cap >= 0 && -inf < cost && cost < inf);
        int id = e.size();
        g[u].push_back(id),e.push_back({v,cap,cost});
        g[v].push_back(id^1),e.push_back({u,0,-cost});
        return id;
    }
    void init(int s){
        fill(h.begin(),h.end(),0);
        bool neg = false;
        for(auto a:e)if(a.cap && a.cost < 0){neg = true;break;}
        if(!neg)return;
        fill(h.begin(),h.end(),inf);
        vector<char> in(n+1);
        vector<int> len(n+1);
        queue<int> q;
        h[s] = 0,q.push(s),in[s] = 1;
        while(!q.empty()){
            int u = q.front();q.pop();in[u] = 0;
            for(int id:g[u]){
                auto a = e[id];
                if(!a.cap || h[a.to] <= h[u]+a.cost)continue;
                h[a.to] = h[u]+a.cost,len[a.to] = len[u]+1;
                assert(len[a.to] < n); // 可达负环不在本板子的适用范围。
                if(!in[a.to])q.push(a.to),in[a.to] = 1;
            }
        }
        for(int u=1;u<=n;++u)if(h[u] == inf)h[u] = 0;
    }
    bool dijkstra(int s,int t){
        fill(dis.begin(),dis.end(),inf);
        priority_queue<pair<ll,int>,vector<pair<ll,int>>,greater<pair<ll,int>>> q;
        dis[s] = 0,q.push({0,s});
        while(!q.empty()){
            auto [d,u] = q.top();q.pop();
            if(d != dis[u])continue;
            for(int id:g[u]){
                auto a = e[id];
                if(!a.cap)continue;
                ll nd = d+a.cost+h[u]-h[a.to];
                if(nd < dis[a.to]){
                    dis[a.to] = nd,pre[a.to] = id;
                    q.push({nd,a.to});
                }
            }
        }
        if(dis[t] == inf)return false;
        for(int u=1;u<=n;++u)if(dis[u] != inf)h[u] += dis[u];
        return true;
    }
    pair<ll,ll> run(int s,int t,ll limit = LLONG_MAX){
        assert(s != t && limit >= 0);
        if(!limit)return {0,0};
        init(s);
        ll flow = 0,cost = 0;
        while(flow < limit && dijkstra(s,t)){
            ll z = limit-flow;
            for(int v=t;v!=s;v=e[pre[v]^1].to)z = min(z,e[pre[v]].cap);
            for(int v=t;v!=s;v=e[pre[v]^1].to){
                int id = pre[v];
                cost += z*e[id].cost;
                e[id].cap -= z,e[id^1].cap += z;
            }
            flow += z;
        }
        return {flow,cost};
    }
}; // O(nm + A m log m)，A 为增广次数；首项仅负费用初始化需要。
