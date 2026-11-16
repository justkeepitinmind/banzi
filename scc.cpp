#include<bits/stdc++.h>
using namespace std;

// Tarjan 强连通分量，有向图，点 1..n，允许重边、自环、不连通。
// O(n+m)。bel[u] 为 1..cnt；跨分量边 u->v 满足 bel[u] > bel[v]。
// 递归版方便手写；长链深度 O(n)，大图注意运行环境的栈大小。
struct SCC{
    int n,tim = 0,cnt = 0;
    vector<vector<int>> g;
    vector<int> dfn,low,bel,sz,st;
    SCC(int n):n(n),g(n+1),dfn(n+1),low(n+1),bel(n+1),sz(n+1){}
    void add(int u,int v){g[u].push_back(v);}
    void dfs(int u){
        dfn[u] = low[u] = ++tim;
        st.push_back(u);
        for(int v:g[u]){
            if(!dfn[v])dfs(v),low[u] = min(low[u],low[v]);
            else if(!bel[v])low[u] = min(low[u],dfn[v]); // 只看仍在栈中的点。
        }
        if(low[u] != dfn[u])return;
        ++cnt;
        int v;
        do{
            v = st.back(),st.pop_back();
            bel[v] = cnt,++sz[cnt];
        }while(v != u);
    }
    void run(){
        tim = cnt = 0,st.clear();
        fill(dfn.begin(),dfn.end(),0);
        fill(low.begin(),low.end(),0);
        fill(bel.begin(),bel.end(),0);
        fill(sz.begin(),sz.end(),0);
        for(int u=1;u<=n;++u)if(!dfn[u])dfs(u);
    }
    vector<vector<int>> dag()const{ // run 后调用；保留重边，按需自行去重。
        vector<vector<int>> a(cnt+1);
        for(int u=1;u<=n;++u)for(int v:g[u])
            if(bel[u] != bel[v])a[bel[u]].push_back(bel[v]);
        return a;
    }
};
