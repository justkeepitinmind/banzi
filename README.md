# 竞赛板子

C++17，按函数/结构体复制使用。新增代码用固定数值类型，不写 `template<class ...>`，不依赖 `modint` 或通用图类。文件没有 `main`，只抄题目需要的部分即可。

| 文件 | 内容 | 主要复杂度 |
| --- | --- | --- |
| [geo.cpp](geo.cpp) | 原有计算几何板子 | 见函数注释 |
| [lichao.cpp](lichao.cpp) | KACTL 动态凸包、离散李超树、区间加线 | 整线插入/查询 `O(log n)`，区间插入 `O(log² n)` |
| [gp_hash_table.cpp](gp_hash_table.cpp) | 随机整数哈希、基本操作、预分配、二元组键 | 通常均摊 `O(1)`，非最坏保证 |
| [scc.cpp](scc.cpp) | Tarjan 强连通分量、缩点 DAG | `O(n+m)` |
| [bcc.cpp](bcc.cpp) | 点双、割点、圆方树；边双、桥、桥树 | `O(n+m)` |
| [flow.cpp](flow.cpp) | Dinic、势能 Dijkstra 最小费用最大流 | 见文件注释 |
| [hlpp.cpp](hlpp.cpp) | 备用最大流：最高标号预流推进 + gap | `O(n² sqrt(m))` |
| [ntt.cpp](ntt.cpp) | `998244353` 卷积、正/逆 NTT | `O(n log n)` |

## 动态凸包与李超树

你给的 [KACTL LineContainer](https://github.com/kth-competitive-programming/kactl/blob/main/content/data-structures/LineContainer.h) 是用 `multiset` 维护的动态凸包。加线、查询顺序都可以任意；插入均摊 `O(log n)`，查询 `O(log n)`，无需限定横坐标范围。它与李超树解决的问题相近，数据结构不同。

```cpp
LineContainer h; // 最大值；查询前至少插一条线。
h.add(2,3);      // y = 2x+3
h.add(-1,8);
i128 ans = h.query(4); // 11
```

交点除法要向下取整；相同斜率只留更优直线。这版将交点差值、乘法和返回值放进 `i128`，包括 `ll` 极值也能比较。返回值确定装得下 `ll` 时再强转输出。

离散李超树适合能预读所有询问横坐标、或者还需要**区间加线**的题。横坐标只做排序去重，不要求等距；比较函数值时仍用真实坐标。

```cpp
LiChao t({-5,0,4,10}); // 最小值；仅能询问已登记坐标。
t.add(2,3);            // 在全部横坐标上加入 y = 2x+3。
t.add_segment(0,4,-1,1); // 在原坐标闭区间 [0,4] 加入 y = -x+1。
i128 ans = t.query(4); // -3
```

李超树允许 `|k|,|b|,|x| <= 1e18`，没有覆盖询问点的直线时返回 `LiChao::inf`。求最大值/最小值切换可以将插入的 `k,b` 和答案一起取负；空值要先判断，且对 `ll` 取负时不能输入 `LLONG_MIN`。两种结构都不支持删线。

## gp_hash_table

```cpp
hashmap h; // gp_hash_table<ll,int,chash>，值类型在别名处改。
++h[123];
if(h.find(123) != h.end())h.erase(123);
```

它是 GNU 扩展，需要 GCC/libstdc++。这台 Mac 可用 `g++-15`；系统自带 `g++` 通常指向 Apple Clang，默认环境不能直接用 PBDS。

`[]` 会插入不存在的键；只判断存在性用 `find`。`insert` 不覆盖旧值；`erase` 按键删除。默认策略没有可用的 `reserve`，也不能直接调用 `resize`。需要预分配时用：

```cpp
hashmap h({},{},{},{},{1<<18}); // 初始桶数为 2 的幂，取预计元素数的 2~4 倍。
```

`find` 返回的迭代器与遍历迭代器类型可能不同，分开声明；容器变化后重新查找。遍历顺序没有保证，不适合依赖有序性的题。`pack(x,y)` 将两个 32 位整数无碰撞编码成 `ull`，对应声明 `gp_hash_table<ull,int,chash>`。随机哈希用于降低被固定构造卡住的风险，实际常数仍取决于题目数据。

## Tarjan、点双、边双和圆方树

点编号统一 `1..n`。无向边只 `add` 一次，函数自动建两个方向；重边按边编号区分。所有图允许不连通；SCC/EBCC 允许自环，VBCC 忽略自环。

```cpp
SCC s(n);
s.add(u,v); // 有向边。
s.run();
// s.cnt：分量数；s.bel[u]：所属分量，编号 1..cnt；s.sz[c]：大小。
auto dag = s.dag(); // 保留缩点重边；需要时自己排序去重。
```

SCC 编号按逆拓扑序产生，跨分量边 `u -> v` 满足 `bel[u] > bel[v]`。需要按拓扑序推 DP 时遍历 `cnt..1`。

```cpp
VBCC v(n);
int id = v.add(u,w); // 无向边；返回输入边编号，从 0 开始。
v.run();
// v.cut[u]：是否割点；v.bcc[i]：第 i 个点双的顶点集合，i 从 0 开始。
auto tr = v.block_tree();
// 原顶点/圆点 1..n；bcc[i] 的方点 n+1+i。
int tot = (int)tr.size()-1;
```

点双以**顶点集合**给出，同一割点可以属于多个点双；桥的两个端点构成一个块。本实现让孤立点也有单点块，因此圆方树中孤立点连接一个单点方点。输入不连通时输出圆方森林，后续 LCA/树 DP 需要遍历每个连通分支。构建过程与 [OI Wiki 的圆方树说明](https://oi-wiki.org/graph/block-forest/) 使用相同的圆点、方点含义。

**可选：输出每条边属于哪个点双。** 默认实现保持原版，相关扩展代码放在 `bcc.cpp` 内，用逐行 `//` 注释并标有「边归属」；有需求时取消各处实现代码行前的 `//`，保留说明注释，输出示例则放到主函数中。

启用后，每条非自环边恰好属于一个点双，`edge_bel[id]` 与 `bcc` 使用同一套从 0 开始的编号，对应方点为 `n+1+edge_bel[id]`。树边和返祖边各入边栈一次，分出点双时弹到当前树边为止并记录归属，仍为 `O(n+m)` 时间、空间。重边分别记录编号；自环沿用忽略约定，归属为 `-1`；孤立点的单点块没有边。再次 `run()` 时会重新计算全部归属。

```cpp
// 仅在取消「边归属」实现注释后使用：
// int c = v.edge_bel[id]; // 这条边所属的点双为 v.bcc[c]；自环返回 -1。
// for(int i=0;i<v.m;++i)cout << v.edge_bel[i] << '\n';
```

```cpp
EBCC e(n);
int id = e.add(u,v); // 无向边编号从 0 开始。
e.run();
// e.bridge[id]：是否桥；e.bel[u]：所属边双，1..cnt；e.sz[c]：大小。
auto tr = e.bridge_tree(); // 以边双为顶点，编号 1..cnt；可能为森林。
```

三个对象都可以再次调用 `run()`，会清空计算状态、保留原图。点双/割点的关键条件是 `low[v] >= dfn[u]`，桥的关键条件是 `low[v] > dfn[u]`；无向图跳过的是父**边**。这些 DFS 使用递归，长链深度可达 `n`，大图需要保证足够的运行栈。

## 网络流

```cpp
Dinic d(n);
int id = d.add(u,v,cap); // 有向边，返回正向残量边下标。
ll f = d.run(s,t);
ll edge_flow = d.e[id^1].cap;
// 完整求到无法增广后，d.dep[u] != -1 表示最小割的源侧。
```

Dinic 使用连续边表、当前弧、失效层次点剪枝，容量为 `ll`。一般图最坏 `O(n²m)`，二分图匹配网络可达 `O(m sqrt(n))`。无向容量边建 `add(u,v,c)` 和 `add(v,u,c)`；互为反向的两条业务边也各自 `add`。支持重边、自环。DFS 同样要注意长链递归深度。

```cpp
MCMF f(n);
f.add(u,v,cap,cost);
auto [flow,cost] = f.run(s,t);      // 最小费用最大流。
// 或 auto [flow,cost] = f.run(s,t,need); // 最多送 need，检查 flow == need。
```

费用流有负费用残量边时先用一次 SPFA 初始化势能，之后每次用 Dijkstra，避免每轮 SPFA。无负边时直接跳过初始化；不支持源点可达的负费用环。复杂度为初始化最坏 `O(nm)` 加 `A` 次 Dijkstra 的 `O(A m log m)`，其中 `A` 是增广次数。

两者要求 `s != t`，`run` 修改残量网络，重复调用返回**本次新增**的流量/费用。`limit` 默认 `LLONG_MAX`。容量和总流量要装进 `ll`；费用流的费用、距离、势能和运算中间值还须满足文件中的数值范围。指定流量提前停止时不能直接把当前可达集当作最小割。

备用最大流单独放在 [hlpp.cpp](hlpp.cpp)，不影响原来的 Dinic、费用流：

```cpp
HLPP f(n);
int id = f.add(u,v,cap); // 自环忽略并返回 -1，其他边返回正向残量边下标。
ll ans = f.run(s,t);
// f.source_side(u)：最小割源侧；id != -1 时 f.e[id^1].cap 为该输入边流量。
```

该实现改写自 [KACTL PushRelabel](https://github.com/kth-competitive-programming/kactl/blob/main/content/graph/PushRelabel.h)，使用高度桶选最高活跃点、当前弧和 gap 优化，无递归。点编号为 `1..n`，允许重边，忽略自环；无向容量边分别加两个方向。**同一对象只能 `run` 一次，运行后不能再加边，也没有 `limit` 参数。** 初始化会推满所有源点出边，因此要求源点出边容量之和不超过 `LLONG_MAX`，即使最终最大流很小也要满足这一点。实际速度取决于图结构，不保证所有题都比 Dinic 快。

## NTT

```cpp
vector<int> a = {1,2,3},b = {4,5}; // 常数项在前。
auto c = NTT::multiply(a,b); // {4,13,22,15}
```

固定模数 `998244353`，原根 `3`，缓存单位根，小数组直接朴素乘法。`multiply` 接受负系数，返回模意义下的卷积，输入为空则返回空；结果长度不能超过 `2^23`。此版本不用于任意模数卷积。

直接调用 `NTT::ntt(a)` / `NTT::ntt(a,true)` 时，数组长度须为不超过 `2^23` 的 2 的幂，元素须已归一化到 `[0,mod)`。卷积会保留原数组；不再用输入时可传 `move(a),move(b)`。

## 验证

新增板子用 GCC 15、C++17 编译；测试固定随机种子以便复现：

```sh
g++-15 -std=c++17 -O2 -Wall -Wextra tests/stress.cpp -o /tmp/banzi-stress
/tmp/banzi-stress
```

测试包含 15100 组随机样例与额外边界检查：直线查询对拍遍历所有直线；SCC 对拍传递闭包；割点/桥对拍删点/删边；点双对拍枚举所有极大无割点连通诱导子图；圆方森林检查分支数和无环性；Dinic 对拍枚举最小割；HLPP 用 4500 组小图对拍枚举最小割并检查流量守恒、另用 500 组较大图对拍 Dinic；费用流对拍枚举所有整数可行流；NTT 对拍朴素卷积。此外有哈希 100000 次操作、大整数边界、重复运行、NTT 正逆变换及大规模卷积检查。原有 `geo.cpp` 不在这次测试范围内。

## 参考

- [KACTL LineContainer](https://github.com/kth-competitive-programming/kactl/blob/main/content/data-structures/LineContainer.h)：动态凸包主体，保留作者与 CC0 来源；本版加宽整数运算。
- [KACTL HashMap](https://github.com/kth-competitive-programming/kactl/blob/main/content/data-structures/HashMap.h)、[GCC PBDS 文档](https://gcc.gnu.org/onlinedocs/libstdc++/manual/policy_data_structures_design.html)：构造、接口与 resize 策略。
- [KACTL SCC](https://github.com/kth-competitive-programming/kactl/blob/main/content/graph/SCC.h)、[KACTL 双连通分量](https://github.com/kth-competitive-programming/kactl/blob/main/content/graph/BiconnectedComponents.h)、[OI Wiki 圆方树](https://oi-wiki.org/graph/block-forest/)：用于核对算法与定义；这里采用普通递归函数，没有回调模板。
- [KACTL Dinic](https://github.com/kth-competitive-programming/kactl/blob/main/content/graph/Dinic.h)、[KACTL 最小费用流](https://github.com/kth-competitive-programming/kactl/blob/main/content/graph/MinCostMaxFlow.h)：网络流参考；这里采用普通 Dinic 和标准堆 Dijkstra。
- [KACTL PushRelabel](https://github.com/kth-competitive-programming/kactl/blob/main/content/graph/PushRelabel.h)：HLPP 备用版本来源，保留作者与 CC0 来源；改用边编号和当前弧下标，源汇不入活跃桶。
- [KACTL NTT](https://github.com/kth-competitive-programming/kactl/blob/main/content/numerical/NumberTheoreticTransform.h)：单位根缓存写法，保留作者与 CC0 来源。
