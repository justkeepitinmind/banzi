#include<bits/stdc++.h>
#include<ext/pb_ds/assoc_container.hpp>
using namespace std;
using namespace __gnu_pbds;
using ll = long long;
using ull = unsigned long long;

// GNU 扩展，需要 GCC/libstdc++；不是 ISO C++ 的一部分。
// 普通整数键用这个，随机种子防针对固定哈希的构造。
struct chash{
    static ull splitmix64(ull x){
        x += 0x9e3779b97f4a7c15ULL;
        x = (x^(x>>30))*0xbf58476d1ce4e5b9ULL;
        x = (x^(x>>27))*0x94d049bb133111ebULL;
        return x^(x>>31);
    }
    size_t operator()(ull x)const{
        static const ull seed = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x+seed);
    }
};
using hashmap = gp_hash_table<ll,int,chash>; // 值类型按题改；不写泛型封装。

// 用法示例，抄题时不需要本函数。
void hash_example(){
    hashmap h;
    h[7] = 10;
    ++h[8]; // 不存在会先插入 0。
    auto [it,ok] = h.insert({7,20}); // 已存在则不覆盖，ok == false。
    (void)it,(void)ok;
    auto p = h.find(9); // 只查询存在性时用 find，别用 [] 白插入。
    if(p != h.end())cout << p->second << '\n';
    h.erase(7); // 按 key 删；不要照抄 unordered_map 的 erase(iterator)。
    for(auto& [key,val]:h)cout << key << ' ' << val << '\n';
    cout << h.size() << ' ' << h.empty() << '\n';
    h.clear();

    // 可选：直接给初始桶数，必须为 2 的幂；预计插入数的 2~4 倍。
    // 构造参数参考 KACTL HashMap.h（Simon Lindholm/chilli，CC0）。
    hashmap big({},{},{},{},{1<<18});
    big[1] = 2;
}

// 没有 reserve；默认策略不支持直接 resize，照搬 h.resize(...) 会编译失败。
// 自动扩容即可，或使用上面的初始桶数写法。扩容/插入/删除后别保留旧迭代器。
// 二元组键：两个 32 位整数无碰撞打包，配 gp_hash_table<ull,int,chash>。
ull pack(int x,int y){return (ull)(uint32_t)x<<32 | (uint32_t)y;}
