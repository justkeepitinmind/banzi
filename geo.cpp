#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<random>
#include<cassert>
using namespace std;
using ll = long long;
using i128 = __int128_t;
using db = long double;
using T = ll; // 整数判定用 ll；浮点构造用 db，按题切换。

const db eps = 1e-10,pi = acosl(-1);
T sq(T x){return x*x;}
int sgn(int x){return (x > 0) - (x < 0);}
int sgn(ll x){return (x > 0) - (x < 0);}
int sgn(i128 x){return (x > 0) - (x < 0);}
int sgn(db x){return (x > eps) - (x < -eps);}

struct pt{
    T x,y;
    pt operator+(pt a)const{return {x+a.x,y+a.y};}
    pt operator-(pt a)const{return {x-a.x,y-a.y};}
    pt operator-()const{return {-x,-y};}
    pt operator*(T k)const{return {x*k,y*k};}
    T operator*(pt a)const{return x*a.x+y*a.y;}
    T operator^(pt a)const{return x*a.y-y*a.x;}
    int operator==(pt a)const{return x == a.x && y == a.y;}
    int operator<(pt a)const{
        if(x != a.x)return x < a.x;
        return y < a.y;
    }
    int toleft(pt a)const{return sgn(*this^a);} // a 在当前向量左侧为 1
    int quad()const{
        if(x == 0 && y == 0)return -1;
        if(y > 0 || (y == 0 && x > 0))return 0;
        return 1;
    }
    pt rot()const{return {-y,x};}
    pt rot(db a)const{return {(T)(x*cosl(a)-y*sinl(a)),(T)(x*sinl(a)+y*cosl(a))};}
    db len()const{return hypotl(x,y);}
    db arg()const{return atan2l(y,x);}
};

struct argcmp{
    int operator()(pt a,pt b)const{
        if(a.quad() != b.quad())return a.quad() < b.quad();
        return a.toleft(b) > 0;
    }
}cmp; // +x 轴起逆时针，同向等价，零向量最前。

db arg(pt a,pt b){return atan2l(a^b,a*b);} // 有向夹角；旋转、角度均为弧度。

struct line{pt p,v;}; // p+t*v；v 非零，半平面取左侧。

struct Cir{pt c;T r;};

// 按需抄函数；下方“依赖”默认已抄类型别名、pt、sq、sgn 和所用标准库头文件。
// 判定部分：T=ll，输入坐标/半径范围 [-1e9,1e9]，半径非负。
// 依赖：line。
int relation(line a,pt p){return a.v.toleft(p-a.p);}

// 无额外依赖。
int onseg(pt a,pt b,pt p){return !(a-p).toleft(b-p) && (a-p)*(b-p) <= 0;}

// 依赖：onseg。
int seg_inter(pt a,pt b,pt c,pt d){ // 含端点、重合、退化。
    if(onseg(a,b,c) || onseg(a,b,d) || onseg(c,d,a) || onseg(c,d,b))return 1;
    return (b-a).toleft(c-a)*(b-a).toleft(d-a) < 0
        && (d-c).toleft(a-c)*(d-c).toleft(b-c) < 0;
}
// 严格逆时针凸包，不重复首点，最小点在首位；全共线只留两端。
// 无额外依赖。
vector<pt> convexhull(vector<pt> a){
    sort(a.begin(),a.end());
    a.erase(unique(a.begin(),a.end()),a.end());
    int n = a.size(),k = 0;
    if(n <= 1)return a;
    vector<pt> b(2*n);
    for(auto p:a){
        while(k > 1 && (b[k-1]-b[k-2]).toleft(p-b[k-1]) <= 0)--k;
        b[k++] = p;
    }
    int t = k;
    for(int i=n-2;i>=0;--i){
        while(k > t && (b[k-1]-b[k-2]).toleft(a[i]-b[k-1]) <= 0)--k;
        b[k++] = a[i];
    }
    b.resize(k-1);
    return b;
}
// 有向面积，逆时针为正，累加用 db。
// 无额外依赖。
db area(const vector<pt>& a){
    db ans = 0;
    for(int i=1;i+1<(int)a.size();++i){
        auto p = a[i]-a[0],q = a[i+1]-a[0];
        ans += (db)p.x*q.y-(db)p.y*q.x;
    }
    return ans/2;
}
// O(n log n)，附加 [-inf,inf]^2 边界；返回逆时针边界直线。
// 空集/零面积返回 {}，相邻直线交点为结果顶点。同向线在入队时保留最左侧。
// 依赖：line、argcmp 和 cmp；不依赖 inter。
vector<line> halfinter(vector<line> a,const T inf=1000000000){
    auto side = [](line a,line b,line c){
        ll d = a.v^b.v,t = (b.p-a.p)^b.v;
        return sgn((i128)(c.v^(a.p-c.p))*d+(i128)(c.v^a.v)*t)*sgn(d);
    };
    
    a.push_back({{-inf,-inf},{1,0}});
    a.push_back({{inf,-inf},{0,1}});
    a.push_back({{inf,inf},{-1,0}});
    a.push_back({{-inf,inf},{0,-1}});
    sort(a.begin(),a.end(),[](line a,line b){return cmp(a.v,b.v);});
    vector<line> q(a.size());
    int l = 0,r = -1; // 有效队列 q[l..r]。
    for(auto x:a){
        if(l <= r && !cmp(q[r].v,x.v)){
            if(q[r].v.toleft(x.p-q[r].p) <= 0)continue;
            --r;
        }
        while(l < r && side(q[r-1],q[r],x) <= 0)--r;
        while(l < r && side(q[l],q[l+1],x) <= 0)++l;
        if(l <= r && !q[r].v.toleft(x.v))return {};
        q[++r] = x;
    }
    while(r-l > 1 && side(q[r-1],q[r],q[l]) <= 0)--r;
    while(r-l > 1 && side(q[l],q[l+1],q[r]) <= 0)++l;
    if(r-l < 2 || !q[r].v.toleft(q[l].v))return {};
    return {q.begin()+l,q.begin()+r+1};
}
// O(n+m)，输入严格逆时针凸包，起点任意；允许空集、点、线段。
// 无额外依赖；给任意点集时，先各自调用 convexhull。
vector<pt> convex_add(vector<pt> a,vector<pt> b){
    if(a.empty() || b.empty())return {};
    auto prep = [](vector<pt>& a){
        rotate(a.begin(),min_element(a.begin(),a.end(),[](pt a,pt b){
            if(a.y != b.y)return a.y < b.y;
            return a.x < b.x;
        }),a.end());
    };
    prep(a),prep(b);
    int n = a.size(),m = b.size(),i = 0,j = 0;
    pt p = a[0]+b[0];
    vector<pt> c = {p};
    while(i < n || j < m){
        pt u = a[(i+1)%n]-a[i%n],v = b[(j+1)%m]-b[j%m];
        int z = u.toleft(v),x = i < n && (j == m || z >= 0),y = j < m && (i == n || z <= 0);
        if(x)p = p+u,++i;
        if(y)p = p+v,++j;
        c.push_back(p);
    }
    c.pop_back();
    c.erase(unique(c.begin(),c.end()),c.end());
    if(c.size() > 1 && c.front() == c.back())c.pop_back();
    return c;
}

// O(log n)，严格逆时针凸包，起点任意；dir 返回每个顶点处的方向。
// 支撑点 extreme(a,[&](pt){return v;}) 最大化 v^a[i]；v 非零。
// 外点切点 extreme(a,[&](pt p){return p-u;})，凸包在 u->a[i] 右侧。
// 参考 Al.Cash: https://codeforces.com/blog/entry/48868 （extreme vertex）。
// 无额外依赖；下面三个封装按需抄一个即可。

template<class F>
int extreme(const vector<pt>& a,F dir){
    int n = a.size(),s;
    if(n <= 1)return n-1;
    auto compare = [&](int i,int j){return dir(a[j]).toleft(a[j]-a[i]);};
    auto check = [&](int i,int& s) -> int {
        s = compare((i+1)%n,i);
        return s >= 0 && compare(i,(i+n-1)%n) < 0;
    };
    if(check(0,s))return 0;
    int l = 1,r = n-1; // 候选闭区间 [l,r]，s 维护 l-1 处的符号。
    while(l < r){
        int mid = (l+r)>>1,t;
        if(check(mid,t))return mid;
        if(s < t || (s == t && s == compare(l-1,mid)))r = mid-1;
        else l = mid+1,s = t;
    }
    return l;
}
// 依赖：extreme。平行于 v 的支撑线接触点下标 {左侧最远点,右侧最远点}，v 非零。
pair<int,int> targent_line(const vector<pt>& a,pt v){
    return {extreme(a,[&](pt){return v;}),extreme(a,[&](pt){return -v;})};
}
// 依赖：extreme。外点 u 的切点下标 {右切点,左切点}，u 必须严格在凸包外。
pair<int,int> targent_pt(const vector<pt>& a,pt u){
    return {extreme(a,[&](pt p){return p-u;}),extreme(a,[&](pt p){return u-p;})};
}
// 依赖：extreme。投影最小/最大点的下标 {min,max}，即最小/最大化 a[i]*v；v 非零。
pair<int,int> min_max_proj(const vector<pt>& a,pt v){
    return {extreme(a,[&](pt){return v.rot();}),extreme(a,[&](pt){return -v.rot();})};
}
// O(log n)，严格逆时针凸包；内部 1，边界 0，外部 -1。
// 依赖：onseg；不需要另抄普通多边形的 relation。
int relation_convex(const vector<pt>& a,pt p){
    int n = a.size();
    if(n < 3){
        if(n && onseg(a[0],a.back(),p))return 0;
        return -1;
    }
    int x = (a[1]-a[0]).toleft(p-a[0]),y = (a.back()-a[0]).toleft(p-a[0]);
    if(x < 0 || y > 0)return -1;
    if(!x)return onseg(a[0],a[1],p)-1;
    if(!y)return onseg(a[0],a.back(),p)-1;
    int l = 1,r = n-1;
    while(l < r){
        int mid = (l+r+1)>>1;
        if((a[mid]-a[0]).toleft(p-a[0]) >= 0)l = mid;
        else r = mid-1;
    }
    return (a[l+1]-a[l]).toleft(p-a[l]);
}
// O(n)，简单多边形，方向任意；内部 1，边界 0，外部 -1。
// 抄：pt（减法、点积、叉积、toleft）、sgn、onseg、本函数；不用抄其他 relation。
int relation(const vector<pt>& a,pt p){
    int wn = 0,n = a.size();
    for(int i=0;i<n;++i){
        pt u = a[i],v = a[(i+1)%n];
        if(onseg(u,v,p))return 0;
        int z = (v-u).toleft(p-u);
        if(u.y <= p.y && v.y > p.y && z > 0)++wn;
        if(v.y <= p.y && u.y > p.y && z < 0)--wn;
    }
    if(wn)return 1;
    return -1;
}
// 依赖：Cir。整数圆：圆内 1，圆上 0，圆外 -1。
int relation(Cir a,pt p){
    return sgn(sq(a.r)-(p-a.c)*(p-a.c));
}
// 整数点与三点外接圆：圆内 1，圆上 0，圆外 -1；a,b,c 不共线。
// 无额外依赖，不需要构造外接圆。
int relation(pt a,pt b,pt c,pt p){
    int z = (b-a).toleft(c-a);
    assert(z);
    a = a-p,b = b-p,c = c-p;
    return sgn((i128)(a*a)*(b^c)+(i128)(b*b)*(c^a)+(i128)(c*c)*(a^b))*z;
}

// 以下浮点构造用 T=db；与上面的整数判定按题选用。
// inter(两直线)、proj、dis 都只依赖 line，各自独立。
pt inter(line a,line b){ // 不平行
    db cof = (db)((b.p-a.p)^b.v)/(a.v^b.v);
    return a.p+a.v*cof;
}
// O(log n)，T=db；严格逆时针凸包，起点任意，允许空集、点、线段；l.v 非零。
// 返回 0/1/2 个交点，沿 l.v 排列；共边返回该边两端点。
// 依赖：line、extreme、targent_line、inter(两直线)，不用抄其他 inter。
vector<pt> inter(line l,const vector<pt>& a){
    int n = a.size();
    if(!n)return {};
    auto side = [&](int i){return l.v.toleft(a[i]-l.p);};
    auto [x,y] = targent_line(a,l.v);
    if(side(x) < 0 || side(y) > 0)return {};
    vector<pt> ans;
    if(!side(x) || !side(y)){
        int i = x;
        if(side(i))i = y;
        ans.push_back(a[i]);
        for(int j:{(i+1)%n,(i+n-1)%n})if(j != i && !side(j)){
            ans.push_back(a[j]);
            break;
        }
    }else{
        auto find = [&](int u,int v){
            int s = side(u);
            if(v < u)v += n;
            int left = u+1,right = v;
            while(left < right){
                int mid = (left+right)>>1;
                if(side(mid%n) == s)left = mid+1;
                else right = mid;
            }
            int i = left%n,j = (left-1)%n;
            if(!side(i))return a[i];
            return inter(l,line{a[j],a[i]-a[j]});
        };
        ans.push_back(find(x,y));
        if(n > 2)ans.push_back(find(y,x));
    }
    if(ans.size() == 2 && (ans[1]-ans[0])*l.v < 0)swap(ans[0],ans[1]);
    return ans;
}
pt proj(line a,pt p){
    db cof = (db)((p-a.p)*a.v)/(a.v*a.v);
    return a.p+a.v*cof;
}
db dis(line a,pt p){return fabsl((db)(a.v^(p-a.p)))/a.v.len();}
// 依赖：line、dis。
db dis_seg(pt a,pt b,pt p){
    if(a == b || (p-a)*(b-a) <= 0)return (p-a).len();
    if((p-b)*(a-b) <= 0)return (p-b).len();
    return dis({a,b-a},p);
}
// 依赖：line、Cir、proj、eps；不用抄其他 inter。
vector<pt> inter(line a,Cir b){
    pt p = proj(a,b.c);
    db h = sq(b.r)-(p-b.c)*(p-b.c);
    if(h < -eps)return {};
    if(h <= eps)return {p};
    pt v = a.v*(T)(sqrtl(h)/a.v.len());
    return {p-v,p+v};
}
// 重合圆返回空，相切返回一个点；零半径重合圆返回圆心。
// 依赖：Cir、eps。
vector<pt> circle_inter(Cir a,Cir b){
    pt v = b.c-a.c;
    db d = v.len();
    if(d <= eps){
        if(a.r == 0 && b.r == 0)return {a.c};
        return {};
    }
    if(d > a.r+b.r+eps || d < fabsl((db)(a.r-b.r))-eps)return {};
    db x = (sq(d)+sq(a.r)-sq(b.r))/(2*d),h = sq(a.r)-sq(x);
    pt p = a.c+v*(T)(x/d);
    if(h <= eps)return {p};
    v = v.rot()*(T)(sqrtl(h)/d);
    return {p-v,p+v};
}
// 两个 get_circle 各自只依赖 Cir，不互相依赖。
Cir get_circle(pt a,pt b){return {(a+b)*(T)0.5,(T)((a-b).len()/2)};}
Cir get_circle(pt a,pt b,pt c){ // 不共线三点外接圆。
    b = b-a,c = c-a;
    pt p = (b.rot()*(c*c)-c.rot()*(b*b))*(T)(0.5/(b^c));
    return {a+p,(T)p.len()};
}
// 期望 O(n)，T=db；依赖 Cir、两个 get_circle、eps、<random>，空集返回原点半径 0。
Cir min_circle(vector<pt> a){
    static mt19937 rng(712367821);
    shuffle(a.begin(),a.end(),rng);
    Cir c{};
    auto out = [&](pt p)->int{
        return (p-c.c).len() > c.r+eps*max((db)1,(db)c.r);
    };
    for(int i=0;i<(int)a.size();++i)if(i == 0 || out(a[i])){
        c = {a[i],0};
        for(int j=0;j<i;++j)if(out(a[j])){
            c = get_circle(a[i],a[j]);
            for(int k=0;k<j;++k)if(out(a[k]))c = get_circle(a[i],a[j],a[k]);
        }
    }
    return c;
}
// O(n)，T=db，简单多边形与圆交面积，方向任意；圆内弦段 + 两侧扇形。
// 依赖：Cir、sq；pt 只需加减、数乘、点积、叉积、rot，不依赖 sgn/eps/inter。
// 参考 KACTL CirclePolygonIntersection.h（CC0）；用直线距圆心的距离求弦。
db area_inter(const vector<pt>& a,Cir c){
    if(c.r <= 0 || a.size() < 3)return 0;
    db ans = 0,rr = sq(c.r);
    auto arc = [&](pt p,pt q) -> db { // 两倍有向扇形面积；零向量不贡献扇形。
        db x = p^q,y = p*q;
        if(x == 0 && y == 0)return 0;
        return rr*atan2l(x,y);
    };
    for(int i=0,n=a.size();i<n;++i){
        pt p = a[i]-c.c,q = a[(i+1)%n]-c.c,d = q-p;
        db A = d*d;
        if(A == 0)continue;
        db z = d^p,h2 = rr-sq(z)/A;
        if(h2 <= 0){ // 相离或相切：没有正长度的圆内弦段。
            ans += arc(p,q);
            continue;
        }
        db h = sqrtl(h2);
        pt e = d*(1/sqrtl(A)),o = d.rot()*(z/A); // 单位方向、圆心在直线上的投影。
        db l = p*e,r = q*e;
        if(l >= h || r <= -h){ans += arc(p,q);continue;}
        pt u = p,v = q;
        if(l < -h)u = o-e*h;
        if(r > h)v = o+e*h;
        ans += arc(p,u)+(u^v)+arc(v,q);
    }
    return fabsl(ans)/2;
}

int main(){}
