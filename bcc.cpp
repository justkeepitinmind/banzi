#include<bits/stdc++.h>
using namespace std;

// 无向图，点 1..n；每条输入边只调用一次 add，边编号从 0 开始。
// 支持重边、不连通；自环不影响连通性，VBCC 直接忽略自环。
// 均为 O(n+m)。递归 DFS 长链深度 O(n)，注意栈大小。

// 点双 + 割点 + 圆方树。点双按顶点集合输出，桥两端也算一个点双。
// 可选：需要每条边所属点双时，取消下面各处「边归属」实现代码块的注释。
// 默认代码与原版一致；启用后 edge_bel[id] 是 bcc 下标（从 0 开始），自环为 -1。
struct VBCC{
    int n,m = 0,tim = 0;
    vector<vector<pair<int,int>>> g; // {邻点,无向边编号}
    vector<vector<int>> bcc; // bcc[0..size-1]，割点可在多个集合中。
    vector<int> dfn,low,st;
    /* 边归属：额外记录边的归属，以及尚未分块的边栈。
    vector<int> edge_bel,est;
    */
    vector<char> cut;
    VBCC(int n):n(n),g(n+1),dfn(n+1),low(n+1),cut(n+1){}
    int add(int u,int v){
        int id = m++;
        /* 边归属：自环不进图，归属保持 -1。
        edge_bel.push_back(-1);
        */
        if(u != v)g[u].push_back({v,id}),g[v].push_back({u,id});
        return id;
    }
    void dfs(int u,int pe){
        dfn[u] = low[u] = ++tim;
        st.push_back(u);
        int child = 0;
        for(auto [v,id]:g[u]){
            if(id == pe)continue; // 跳过父边，不能跳过所有通向父点的边。
            if(!dfn[v]){
                ++child;
                /* 边归属：树边在递归前入栈。
                est.push_back(id);
                */
                dfs(v,id);
                low[u] = min(low[u],low[v]);
                if(low[v] >= dfn[u]){
                    if(pe != -1 || child > 1)cut[u] = 1;
                    bcc.push_back({u});
                    /* 边归属：弹到当前树边为止，这些边恰属于刚建出的点双。
                    int eid;
                    do{
                        eid = est.back(),est.pop_back();
                        edge_bel[eid] = (int)bcc.size()-1;
                    }while(eid != id);
                    */
                    int x;
                    do{
                        x = st.back(),st.pop_back();
                        bcc.back().push_back(x);
                    }while(x != v);
                }
            }else low[u] = min(low[u],dfn[v]);
            /* 边归属：返祖边只在较深端入栈一次，重边分别记录。
            if(dfn[v] < dfn[u])est.push_back(id);
            */
        }
        if(pe == -1 && child == 0)bcc.push_back({u}); // 孤立点单独成块。
    }
    void run(){
        tim = 0,st.clear(),bcc.clear();
        /* 边归属：每次 run 重新计算。
        est.clear();
        fill(edge_bel.begin(),edge_bel.end(),-1);
        */
        fill(dfn.begin(),dfn.end(),0);
        fill(low.begin(),low.end(),0);
        fill(cut.begin(),cut.end(),0);
        for(int u=1;u<=n;++u)if(!dfn[u])dfs(u,-1),st.pop_back();
    }
    vector<vector<int>> block_tree()const{ // run 后调用；不连通时是森林。
        // 圆点 1..n；bcc[i] 对应方点 n+1+i；含孤立点的单点块。
        vector<vector<int>> tr(n+bcc.size()+1);
        for(int i=0;i<(int)bcc.size();++i)for(int u:bcc[i]){
            int v = n+1+i;
            tr[u].push_back(v),tr[v].push_back(u);
        }
        return tr;
    }
};

/* 输出示例：启用边归属后，将此段放在主函数中，v.run() 之后。
for(int id=0;id<v.m;++id)cout << v.edge_bel[id] << '\n';
// 非自环边的点双为 v.bcc[v.edge_bel[id]]，对应方点为 v.n+1+v.edge_bel[id]。
*/

// 边双 + 桥 + 桥树。删掉所有桥后，剩余每个连通块是一个边双。
struct EBCC{
    int n,tim = 0,cnt = 0;
    vector<vector<pair<int,int>>> g;
    vector<pair<int,int>> edges;
    vector<int> dfn,low,bel,sz;
    vector<char> bridge; // bridge[id] 表示输入的第 id 条边是否为桥。
    EBCC(int n):n(n),g(n+1),dfn(n+1),low(n+1),bel(n+1),sz(n+1){}
    int add(int u,int v){
        int id = edges.size();
        edges.push_back({u,v}),bridge.push_back(0);
        g[u].push_back({v,id}),g[v].push_back({u,id});
        return id;
    }
    void dfs(int u,int pe){
        dfn[u] = low[u] = ++tim;
        for(auto [v,id]:g[u]){
            if(id == pe)continue;
            if(!dfn[v]){
                dfs(v,id);
                low[u] = min(low[u],low[v]);
                if(low[v] > dfn[u])bridge[id] = 1;
            }else low[u] = min(low[u],dfn[v]);
        }
    }
    void run(){
        tim = cnt = 0;
        fill(dfn.begin(),dfn.end(),0);
        fill(low.begin(),low.end(),0);
        fill(bel.begin(),bel.end(),0);
        fill(sz.begin(),sz.end(),0);
        fill(bridge.begin(),bridge.end(),0);
        for(int u=1;u<=n;++u)if(!dfn[u])dfs(u,-1);
        vector<int> st;
        for(int u=1;u<=n;++u)if(!bel[u]){
            bel[u] = ++cnt,st.push_back(u);
            while(!st.empty()){
                int x = st.back();st.pop_back();
                ++sz[cnt];
                for(auto [v,id]:g[x])if(!bridge[id] && !bel[v])
                    bel[v] = cnt,st.push_back(v);
            }
        }
    }
    vector<vector<int>> bridge_tree()const{ // run 后调用，点编号 1..cnt。
        vector<vector<int>> tr(cnt+1);
        for(int i=0;i<(int)edges.size();++i)if(bridge[i]){
            auto [u,v] = edges[i];
            tr[bel[u]].push_back(bel[v]),tr[bel[v]].push_back(bel[u]);
        }
        return tr;
    }
};
