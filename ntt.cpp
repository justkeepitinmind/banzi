#include<bits/stdc++.h>
using namespace std;
using ll = long long;

// 固定 998244353，原根 3；无 modint/泛型多项式封装。
// 缓存单位根写法参考 KACTL NumberTheoreticTransform.h，chilli，CC0。
// https://github.com/kth-competitive-programming/kactl/blob/main/content/numerical/NumberTheoreticTransform.h
namespace NTT{
const int mod = 998244353;
int qpow(int a,int b){
    int ans = 1;
    for(;b;b>>=1,a=(ll)a*a%mod)if(b&1)ans = (ll)ans*a%mod;
    return ans;
}
// a 的长度为 2 的幂且 <= 2^23，元素已归一化到 [0,mod)。
// inv=false 正变换，inv=true 逆变换；O(n log n)，单位根缓存 O(n)。
void ntt(vector<int>& a,bool inv = false){
    int n = a.size();
    assert(n > 0 && !(n&(n-1)) && n <= (1<<23));
    static vector<int> rt(2,1);
    for(int k=rt.size();k<n;k<<=1){
        rt.resize(k*2);
        int z = qpow(3,(mod-1)/(k*2));
        for(int i=k;i<k*2;++i)rt[i] = i&1 ? (ll)rt[i>>1]*z%mod : rt[i>>1];
    }
    for(int i=1,j=0;i<n;++i){
        int bit = n>>1;
        for(;j&bit;bit>>=1)j ^= bit;
        j ^= bit;
        if(i < j)swap(a[i],a[j]);
    }
    for(int k=1;k<n;k<<=1)for(int i=0;i<n;i+=k*2)for(int j=0;j<k;++j){
        int u = a[i+j],v = (ll)a[i+j+k]*rt[k+j]%mod;
        a[i+j] = u+v < mod ? u+v : u+v-mod;
        a[i+j+k] = u-v >= 0 ? u-v : u-v+mod;
    }
    if(inv){
        reverse(a.begin()+1,a.end());
        int z = qpow(n,mod-2);
        for(int& x:a)x = (ll)x*z%mod;
    }
}
// a、b 按常数项在前存；允许负系数，返回长度 a.size()+b.size()-1。
// 值传递使原数组不变；不再使用原数组时可传 move(a)、move(b) 减少复制。
vector<int> multiply(vector<int> a,vector<int> b){
    if(a.empty() || b.empty())return {};
    int need = a.size()+b.size()-1,n = 1;
    assert(need <= (1<<23));
    for(int& x:a){x %= mod;if(x < 0)x += mod;}
    for(int& x:b){x %= mod;if(x < 0)x += mod;}
    if(min(a.size(),b.size()) <= 32){
        vector<int> c(need);
        for(int i=0;i<(int)a.size();++i)for(int j=0;j<(int)b.size();++j)
            c[i+j] = (c[i+j]+(ll)a[i]*b[j])%mod;
        return c;
    }
    while(n < need)n <<= 1;
    a.resize(n),b.resize(n);
    ntt(a),ntt(b);
    for(int i=0;i<n;++i)a[i] = (ll)a[i]*b[i]%mod;
    ntt(a,true);
    a.resize(need);
    return a;
}
}
