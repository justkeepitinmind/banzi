#include<bits/stdc++.h>
using namespace std;
using ll = long long;
using i128 = __int128_t;

// 两种写法按需抄：LineContainer 全用 ll；LiChao 的乘法/答案用 i128。
// 1. 动态凸包：任意顺序加整条直线，查询最大值；无需预知横坐标。
// KACTL LineContainer，Simon Lindholm，CC0。
// https://github.com/kth-competitive-programming/kactl/blob/main/content/data-structures/LineContainer.h
// add 均摊 O(log n)，query O(log n)，空间 O(n)。
// 差值、交点除法及 k*x+b 的中间计算须装进 ll，查询满足 |x| < LLONG_MAX。
// 例如 |k|、|b|、|x| <= 1e9 时可以安全使用。
struct HullLine{
    ll k,b;
    mutable ll p; // 本条直线胜过下一条的最右整数横坐标。
    bool operator<(const HullLine& a)const{return k < a.k;}
    bool operator<(ll x)const{return p < x;}
};
struct LineContainer:multiset<HullLine,less<>>{
    static constexpr ll inf = LLONG_MAX;
    ll div(ll a,ll b){ // 向下取整，不能直接用 C++ 的 a/b。
        return a/b-((a < 0) != (b < 0) && a%b != 0);
    }
    bool isect(iterator x,iterator y){
        if(y == end())return x->p = inf,false;
        if(x->k == y->k)x->p = x->b > y->b ? inf : -inf;
        else x->p = div(y->b-x->b,x->k-y->k);
        return x->p >= y->p;
    }
    void add(ll k,ll b){
        auto z = insert({k,b,0}),y = z++,x = y;
        while(isect(y,z))z = erase(z);
        if(x != begin() && isect(--x,y))isect(x,y = erase(y));
        while((y = x) != begin() && (--x)->p >= y->p)isect(x,erase(y));
    }
    ll query(ll x)const{
        assert(!empty());
        auto l = *lower_bound(x);
        return l.k*x+l.b;
    }
};

// 2. 离散李超树：查询最小值，所有查询横坐标预先给定；值仍在线插入。
// |k|、|b|、|x| <= 1e18。空位置返回 inf；不支持删除。
// 整线 add O(log q)，区间加线 add_segment O(log^2 q)，查询 O(log q)，空间 O(q)。
// 区间是原横坐标闭区间 [L,R]；只对预先登记的横坐标生效。
struct LiChao{
    static constexpr i128 inf = (i128)1<<126;
    struct Line{
        ll k = 0,b = 0;
        i128 at(ll x)const{return (i128)k*x+b;}
    };
    int n;
    vector<ll> x;
    vector<Line> tr;
    vector<char> has;
    LiChao(vector<ll> a):x(move(a)){
        sort(x.begin(),x.end());
        x.erase(unique(x.begin(),x.end()),x.end());
        n = x.size();
        assert(n > 0);
        tr.resize(4*n),has.resize(4*n);
    }
    void insert(int p,int l,int r,Line a){
        if(!has[p]){tr[p] = a,has[p] = 1;return;}
        int mid = (l+r)>>1;
        if(a.at(x[mid]) < tr[p].at(x[mid]))swap(a,tr[p]);
        if(l == r)return;
        if(a.at(x[l]) < tr[p].at(x[l]))insert(p*2,l,mid,a);
        else if(a.at(x[r]) < tr[p].at(x[r]))insert(p*2+1,mid+1,r,a);
    }
    void modify(int p,int l,int r,int L,int R,Line a){
        if(L <= l && r <= R){insert(p,l,r,a);return;}
        int mid = (l+r)>>1;
        if(L <= mid)modify(p*2,l,mid,L,R,a);
        if(R > mid)modify(p*2+1,mid+1,r,L,R,a);
    }
    void add(ll k,ll b){insert(1,0,n-1,{k,b});}
    void add_segment(ll L,ll R,ll k,ll b){
        int l = lower_bound(x.begin(),x.end(),L)-x.begin();
        int r = upper_bound(x.begin(),x.end(),R)-x.begin()-1;
        if(l <= r)modify(1,0,n-1,l,r,{k,b});
    }
    i128 query_index(int p,int l,int r,int q)const{
        i128 ans = has[p] ? tr[p].at(x[q]) : inf;
        if(l == r)return ans;
        int mid = (l+r)>>1;
        if(q <= mid)return min(ans,query_index(p*2,l,mid,q));
        return min(ans,query_index(p*2+1,mid+1,r,q));
    }
    i128 query(ll v)const{
        int q = lower_bound(x.begin(),x.end(),v)-x.begin();
        assert(q < n && x[q] == v);
        return query_index(1,0,n-1,q);
    }
};

// LineContainer 求最小值：add(-k,-b)，答案 -query(x)；k、b 不可为 LLONG_MIN。
// LiChao 求最大值同理取负；先判断返回值是否为 inf，再取负。
