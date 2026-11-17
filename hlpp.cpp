#include<bits/stdc++.h>
using namespace std;
using ll = long long;

// HLPP 备用：最高标号预流推进 + gap 优化，无递归，不使用泛型模板。
// 改写自 KACTL PushRelabel，Simon Lindholm，CC0。
// https://github.com/kth-competitive-programming/kactl/blob/main/content/graph/PushRelabel.h
// 点 1..n；容量非负；有向边 add 一次，重边允许，自环忽略（返回 -1）。
// 源点出边容量之和须 <= LLONG_MAX：初始化会先把这些边全部推满。
// 同一对象只 run 一次，run 后不再加边；换题/换源汇重新构造。
// O(n^2 sqrt(m))，实际速度依赖图结构，不保证总比 Dinic 快。
struct HLPP{
    struct Edge{int to;ll cap;};
    int n,s = 0,t = 0;
    vector<Edge> e;
    vector<vector<int>> g,bucket;
    vector<ll> ex; // 点的盈余流量；结束时只有源汇可非零。
    vector<int> h,cur,cnt; // 高度、当前弧、各高度的点数。
    HLPP(int n):n(n),g(n+1),bucket(2*n+2),ex(n+1),h(n+1),cur(n+1),cnt(2*n+2){}
    int add(int u,int v,ll cap){
        assert(!s && cap >= 0);
        if(u == v)return -1;
        int id = e.size();
        g[u].push_back(id),e.push_back({v,cap});
        g[v].push_back(id^1),e.push_back({u,0});
        return id;
    }
    void push(int u,int id,ll f){
        if(!f)return;
        int v = e[id].to;
        if(!ex[v] && v != s && v != t)bucket[h[v]].push_back(v);
        e[id].cap -= f,e[id^1].cap += f;
        ex[u] -= f,ex[v] += f;
    }
    ll run(int S,int T){
        assert(!s && 1 <= S && S <= n && 1 <= T && T <= n && S != T);
        s = S,t = T,h[s] = n;
        cnt[0] = n-1,cnt[n] = 1;
        for(int id:g[s])push(s,id,e[id].cap); // 建立初始预流。
        for(int hi=0;;){
            while(bucket[hi].empty()){
                if(!hi)return ex[t];
                --hi;
            }
            int u = bucket[hi].back();bucket[hi].pop_back();
            while(ex[u] > 0){
                if(cur[u] == (int)g[u].size()){
                    int old = h[u];
                    h[u] = 2*n;
                    // 没有可推的当前弧，抬高到残量邻点最低高度 + 1。
                    for(int i=0;i<(int)g[u].size();++i){
                        int id = g[u][i];
                        if(e[id].cap && h[u] > h[e[id].to]+1)
                            h[u] = h[e[id].to]+1,cur[u] = i;
                    }
                    --cnt[old],++cnt[h[u]];
                    // old 层空了：更高且低于 n 的点都无法再到汇点。
                    if(!cnt[old] && old < n)for(int v=1;v<=n;++v)
                        if(old < h[v] && h[v] < n)
                            --cnt[h[v]],h[v] = n+1,++cnt[h[v]];
                    hi = h[u];
                }else{
                    int id = g[u][cur[u]],v = e[id].to;
                    if(e[id].cap && h[u] == h[v]+1)push(u,id,min(ex[u],e[id].cap));
                    else ++cur[u];
                }
            }
        }
    }
    bool source_side(int u)const{return h[u] >= n;} // run 后：最小割源侧。
};

// 用法：HLPP f(n); int id = f.add(u,v,cap); ll ans = f.run(s,t);
// id != -1 时，f.e[id^1].cap 为该条输入边的最终流量。
// 无向容量边：分别 add(u,v,c)、add(v,u,c)。本板子不提供 limit 参数。
