#include<bits/stdc++.h>
using namespace std;
using ll = long long;

// 动态凸包：任意顺序加整条直线，查询最大值；无需预知横坐标。
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

// 求最小值：add(-k,-b)，答案 -query(x)；k、b 不可为 LLONG_MIN。
// 只支持整条直线插入，不支持区间加线或删除。
