# Strip Packing Algorithm — 项目说明文档

> Copyright Xiangyi Zhang 2021  
> 学术、非商业用途适用。联系方式：xiangyi.zhang@polymtl.ca

---

## 目录

1. [问题描述](#1-问题描述)
2. [项目结构](#2-项目结构)
3. [核心数据结构](#3-核心数据结构)
4. [模块说明](#4-模块说明)
   - 4.1 [数据读取模块 (datareader)](#41-数据读取模块-datareader)
   - 4.2 [基础工具模块 (spp)](#42-基础工具模块-spp)
   - 4.3 [背包问题模块 (knapsack)](#43-背包问题模块-knapsack)
   - 4.4 [天际线数据结构模块 (skyline)](#44-天际线数据结构模块-skyline)
   - 4.5 [启发式算法模块 (heuristic)](#45-启发式算法模块-heuristic)
   - 4.6 [BLEU 精确算法模块 (BLEU + preprocessing)](#46-bleu-精确算法模块-bleu--preprocessing)
5. [函数调用关系图](#5-函数调用关系图)
6. [整体算法流程](#6-整体算法流程)
   - 6.1 [主流程 (main.cpp)](#61-主流程-maincpp)
   - 6.2 [BLEU 算法初始化流程](#62-bleu-算法初始化流程)
   - 6.3 [预处理流程 (preprocessing)](#63-预处理流程-preprocessing)
   - 6.4 [下界计算流程 (bounds)](#64-下界计算流程-bounds)
   - 6.5 [evaluate 评估主流程](#65-evaluate-评估主流程)
   - 6.6 [分支定界算法 (branchAndBound)](#66-分支定界算法-branchandbound)
   - 6.7 [y-check 算法](#67-y-check-算法)
   - 6.8 [启发式算法流程](#68-启发式算法流程)
7. [测试数据格式](#7-测试数据格式)
8. [关键参数与常量](#8-关键参数与常量)
9. [参考文献](#9-参考文献)

---

## 1. 问题描述

**条带装箱问题（Strip Packing Problem, SPP）**：

给定宽度为 W 的无限高条带，以及 n 个矩形物品，每个物品 i 具有宽度 `w_i` 和高度 `h_i`。要求将所有物品无重叠地放入条带中（物品不允许旋转），目标是最小化所使用的条带总高度。

本项目实现了以下论文中描述的 **BLEU（Branching Lower-bound Exact Update）** 精确算法：

- **主要参考论文**：*"Combinatorial Benders' Cuts for the Strip Packing Problem"*
- **关联论文**：*"An Exact Algorithm for the Two-Dimensional Strip Packing Problem"* (Boschetti & Montaletti, 2010)
- **关联论文**：*"A Branch and Bound Algorithm for the Strip Packing Problem"* (Alvarez-Vales et al., OR Spectrum, 2009)

---

## 2. 项目结构

```
Strip-Packing-Algorithm/
├── src/
│   ├── main.cpp            # 程序入口
│   ├── item.h              # （占位头文件，仅含 #pragma once）
│   ├── spp.h / spp.cpp     # 核心数据结构与基础工具函数
│   ├── datareader.h / datareader.cpp   # 数据读取
│   ├── knapsack.h / knapsack.cpp       # 背包问题 DP
│   ├── skyline.h / skyline.cpp         # 天际线数据结构
│   ├── heuristic.h / heuristic.cpp     # 启发式算法
│   ├── BLEU.h / BLEU.cpp               # BLEU 精确算法主体
│   └── preprocessing.cpp               # y-check 预处理（BLEU 类成员实现）
└── test/
    ├── 2sp/                # 小规模测试用例（HT系列）
    └── all instances/      # 完整测试用例集（BENG、CGCUT、GCUT、HT、NGCUT 系列）
```

---

## 3. 核心数据结构

### 3.1 `item`（矩形物品）

定义于 `spp.h`，实现于 `spp.cpp`。

```cpp
class item {
public:
    int idx;          // 主索引（在 _allItems 中的顺序）
    int width;        // 宽度
    int height;       // 高度
    int idxHelper;    // 辅助索引（在 _processedItems 或局部数组中的位置）
    std::vector<const item*> subItems;  // 合并后的子物品（y-check预处理使用）
    static std::ostringstream ss;       // 静态字符串流（输出工具）
};
```

- 提供两种构造函数：带/不带 `idxHelper`。
- `subItems` 用于 y-check 预处理中的物品合并操作。

### 3.2 `coordinate`（坐标）

```cpp
struct coordinate {
    int x;  // 左下角 x 坐标
    int y;  // 左下角 y 坐标
};
```

- 支持 `==`、`<`、`=` 运算符。
- 用于记录物品的放置位置（左下角坐标）。

### 3.3 `whPair`（宽高对）

```cpp
struct whPair {
    int h;  // 高度
    int w;  // 宽度
};
```

- 支持比较运算符，按高度为主键、宽度为次键排序。

### 3.4 `itemPiece` 和 `itemPieceWidth`

```cpp
class itemPiece {
    std::string id;  // 标识符
    int height;      // 高度（用于并行排程问题）
};

class itemPieceWidth {
    int id;          // 标识符
    int width;       // 宽度（用于1CBP问题）
};
```

- `itemPieceWidth` 用于 `cutItemsAlongHeight()` 将物品切成单行高度的片段。

### 3.5 `Skyline`（天际线节点）

定义于 `skyline.h`。

```cpp
struct Skyline {
    int corX;       // 该段左端 x 坐标
    int corY;       // 该段当前高度（已堆叠物品的顶部 y）
    int length;     // 该段宽度
    Skyline* next;  // 右侧相邻段（双向链表）
    Skyline* prev;  // 左侧相邻段
};
```

- 使用带哑头节点（`leftSide`，corY = 999999）和哑尾节点（`rightSide`，corY = 999999）的双向链表。
- 每个节点代表条带中一段水平可用区域及其当前高度。

### 3.6 `BLEU::BBNode`（分支定界树节点）

定义于 `BLEU.h`（嵌套类）。

```cpp
class BBNode {
public:
    int trialHeight;                       // 当前试验高度
    int leftMostIdx;                       // 当前待填充的最左列索引
    std::vector<int> columnsOccupiedHeight; // 每列已占用高度
    std::vector<const item*> remainingItems;// 尚未放置的物品
    std::vector<const item*> packedItems;  // 已放置的物品（按顺序）
    std::vector<int> maxiItemIdxColumns;   // 每列已放物品的最大索引
    std::vector<coordinate> itemPositions; // 物品的坐标（以 idxHelper 为下标）
};
```

- 提供 4 种构造函数：新根节点、y-check 根节点、拷贝构造、从父节点分支构造。

### 3.7 枚举类型

```cpp
enum solutionStatus { feasible, infeasible, pending, numberStatus };
enum algorithmStatus { approximate, exact, numberAlgStatus };
enum skylineSelectionMode { leftBottom, bestFit };
```

---

## 4. 模块说明

### 4.1 数据读取模块 (datareader)

**文件**：`datareader.h`、`datareader.cpp`

| 函数 | 签名 | 说明 |
|------|------|------|
| `readData` | `int readData(const string& file, vector<const item*>& items)` | 从文件读取实例，填充物品列表，返回条带宽度 W |
| `parseLine` | `vector<int> parseLine(string line)` | 解析空白分隔的一行整数 |
| `convertVec` | `int convertVec(vector<int>& tmp)` | 将数字位向量转换为整数 |

**文件格式**（见第7节）：第1行为 n，第2行为 W，之后每行为 `idx width height`。

---

### 4.2 基础工具模块 (spp)

**文件**：`spp.h`、`spp.cpp`

#### 排序比较函数

| 函数 | 排序规则 |
|------|----------|
| `compareItemByWidth(i, j)` | 宽度降序，同宽按高度降序，再按 idx 降序 |
| `compareItemByWidthLess(i, j)` | 宽度升序，同宽按高度升序，再按 idx 升序 |
| `compareItemByHeight(i, j)` | 高度升序，同高按 idx 降序 |
| `compareItemByIdx(i, j)` | idx 升序 |
| `compareItemByArea(i, j)` | 面积（w×h）升序 |
| `compareItemByWHDifference(H, W)` | `min(W-w, H-h)` 降序（仿函数，需绑定容器尺寸） |
| `compareItemByHeight`（类） | 高度降序，同高按 idx 降序（用于堆） |
| `compareItemByxCords`（类） | 按物品 x 坐标降序（需绑定坐标数组） |

#### 工具函数

| 函数 | 说明 |
|------|------|
| `computeFX(x, idx, items, flag)` | 动态规划：计算从 items 中去除第 idx 个物品后，其余物品的尺寸（flag=true 用宽度，false 用高度）能组成的所有子集和（≤ x），返回可行位置集合 |
| `getMaximalHeight(items)` | 返回所有物品中的最大高度 |
| `getVarName(itemIdx, xPos)` | 生成 CPLEX 变量名字符串 `"itemXassignY"` |
| `subSetSum(v, limit)` | 动态规划：求不超过 limit 的最大子集和 |
| `solve(allItems, mapPosWidth, mapPosHeight, Integer)` | 调用 CPLEX 求解平行机排程 LP/ILP，返回目标值（下界5的核心） |

---

### 4.3 背包问题模块 (knapsack)

**文件**：`knapsack.h`、`knapsack.cpp`

| 函数 | 说明 |
|------|------|
| `dynamicPrg4KnapSack<T>(items, capacity)` | 模板0/1背包 DP，返回最优价值（整数） |
| `dynamicPrg4KnapSack<T>(items, capacity, dummy)` | 同上，额外返回 DP 表最后一行（每个物品数量对应的最优值向量） |
| `dynamicPrg4KnapSack(values, weights, capacity, selected)` | 浮点价值0/1背包 DP，返回最优价值，并回溯填充 selected（选中物品索引） |

- 模板版本要求类型 T 具有 `weight` 和 `value` 整数成员。
- 用于 `LowerBound4`（列生成定价子问题）和 `dynamicCuts`（BLEU 分支定界中的动态割）。

---

### 4.4 天际线数据结构模块 (skyline)

**文件**：`skyline.h`、`skyline.cpp`

| 函数 | 说明 |
|------|------|
| `selectSkyline(head, mode)` | 遍历链表，返回 corY 最小的天际线段（leftBottom 与 bestFit 模式当前实现相同） |
| `addItemOverSkyline(skyline, item)` | 将物品放置在 skyline 段上方，更新坐标和天际线链表；若物品宽度 < 段长，分裂该段 |
| `detectAndMergeSkylines(skyline)` | 向左右合并与当前段等高的相邻段 |
| `removeSkyline(skyline)` | 从链表中删除该段并释放内存 |
| `liftSkyline(skyline)` | 将当前段合并至高度更低的相邻段 |

天际线由哑头（`leftSide`，corY=999999）和哑尾（`rightSide`，corY=999999）夹住，所有实际段从头节点的 `next` 开始。

---

### 4.5 启发式算法模块 (heuristic)

**文件**：`heuristic.h`、`heuristic.cpp`

| 函数 | 说明 |
|------|------|
| `leftBottomHeuristic(items, W)` | 按给定顺序，逐个将物品放在最低、最左的天际线段上；若段宽不够则抬升该段，直到找到合适段 |
| `bestFitHeuristic(items, W)` | 按宽度降序排列物品；每次选最低天际线段，找能放下且宽度最接近（最大宽度优先）的物品放入；无合适物品则抬升该段 |
| `iteratedGreedy(items, W)` | 以 `leftBottomHeuristic` 获得初始解，然后在插入邻域（将物品从位置 i 移至位置 j）中迭代改进，直至无法提升 |
| `generalBestFitHeurisitic(items, bins)` | 变长多箱最优拟合：多个矩形箱，每次选所有箱中最低的天际线段，找最合适物品放入；返回是否所有物品都被放入 |
| `findBestItem(items, skyline)` | 在 items 中找第一个宽度 ≤ skyline->length 的物品（已按宽度降序）并移除 |
| `findBestItem(items, skyline, bin)` | 在 items 中找第一个宽度 ≤ skyline->length 且放入后不超出 bin 高度的物品 |
| `parseSol(skyline)` | 遍历天际线链表，返回最大 corY（条带使用高度），同时释放所有节点内存 |
| `dumpSolution(items)` | 将所有物品的矩形信息和坐标解写入 `Rectangle.output` 和 `Solution.output` |
| `insertItem(v, fromPos, toPos)` | 将 v[fromPos] 移动至 v[toPos]（用于迭代贪婪的邻域探索） |

静态成员 `Heuristic::solutions`：存储所有物品的最终坐标（按 `item::idx` 索引）。

---

### 4.6 BLEU 精确算法模块 (BLEU + preprocessing)

**文件**：`BLEU.h`、`BLEU.cpp`、`preprocessing.cpp`

#### 静态参数

| 静态成员 | 默认值 | 说明 |
|---------|--------|------|
| `tolerance` | 0.0001 | 浮点比较容差 |
| `bigNumber` | 999999 | 无穷大替代值 |
| `BBMaxExplNodesPerPack` | 10,000,000 | 完美装箱时 B&B 最大节点数 |
| `BBMaxExplNodesNonPerPack` | 80,000 | 非完美装箱时 B&B 最大节点数 |
| `ycheckExplNode` | 10,000,000 | y-check 枚举树最大节点数 |
| `nodeLimitFlag` | false | 是否触发节点上限 |
| `algStatus` | exact | 算法状态（exact/approximate） |

#### 构造函数

| 构造函数 | 说明 |
|---------|------|
| `BLEU(items, W, TrialHeight, timeLimit)` | evaluatedMode=true：从给定试验高度开始验证可行性 |
| `BLEU(items, W, timeLimit)` | evaluatedMode=false：自动搜索最优高度 |

两个构造函数均执行：排序 → `reassignItemsIdx` → `preprocessing` → `bounds`。

#### 公共接口

| 函数 | 说明 |
|------|------|
| `preprocessing()` | 执行三步预处理（fixItems + reduceW + modifyWidth） |
| `bounds()` | 计算并记录五个下界的最大值 |
| `evaluate()` | 主评估函数：应用高度预处理后调用 `branchAndBound`，返回 feasible/infeasible/pending |
| `solvePCC()` | 类似 evaluate，但调用 `branchAndBoundYRelax`（不含 y-check，求解平行机排程松弛） |
| `dumpSolution(...)` | 三个重载：将解写入文件 |

---

## 5. 函数调用关系图

```
main()
├── readData()
│   ├── parseLine()
│   └── convertVec()
└── BLEU::BLEU(items, W, TrialHeight, timeLimit)
    ├── std::sort(..., compareItemByWidth)
    ├── reassignItemsIdx()
    ├── preprocessing()
    │   ├── preprocessingFixItems()
    │   │   └── [修改 _processedItems, _processedH]
    │   ├── preprocessingReduceW()
    │   │   └── subSetSum()
    │   └── preprocessingModifyItemWidth()
    │       └── subSetSum()
    └── bounds()
        ├── LowerBound1()                    ← 面积下界
        ├── LowerBound2()                    ← 双可行函数下界
        │   ├── DualFeasibleFunction1()
        │   ├── DualFeasibleFunction2()
        │   └── DualFeasibleFunction3()
        ├── LowerBound3()                    ← Alvarez-Vales 启发式下界
        │   └── compareItemByWHDifference()
        ├── LowerBound4()                    ← 列生成 LP 下界 (CPLEX)
        │   └── dynamicPrg4KnapSack()        ← 定价子问题
        └── LowerBound5()                    ← 平行机排程 LP 下界 (CPLEX)
            ├── computeFX()
            └── solve()

BLEU::evaluate()
├── preprocessItemHeight()
│   └── subSetSum()
└── branchAndBound(Items, binWidth, binHeight)
    ├── BBNode::BBNode(...)                  ← 构建根节点
    ├── bounding(currentNode)
    │   └── dynamicCuts(currentNode)
    │       └── dynamicPrg4KnapSack<oneDimensionItem>()
    ├── makeBranch(currentNode, dfsTree)
    │   └── compareItemByIdx()
    └── yCheckAlgorithm(W, H, positions, items)   ← 叶节点验证
        ├── preprocess4yCheck()
        │   ├── preprocessedFirst4yCheck()
        │   │   ├── getLeftsAndRights()
        │   │   └── mergeItems4yCheck()      ← 循环调用
        │   │       ├── yCheckEnumerationTree()  ← 递归
        │   │       ├── merging()
        │   │       ├── checkSeparable()
        │   │       ├── getFirstColumn()
        │   │       ├── getLastColumn()
        │   │       ├── getMaxWidth()
        │   │       └── transferItemsAndCords4YEnumeration()
        │   ├── preprocessedSecond4yCheck()
        │   │   └── getLeftsAndRights()
        │   └── preprocessedThird4yCheck()
        │       └── getItemsByCol()
        ├── yCheckEnumerationTree(items, cords, H, W)
        │   ├── BBNode::BBNode(items, cords, W, H)
        │   ├── yCheckBounding(currentNode)
        │   └── yCheckMakeBranch(currentNode, yEnTree)
        │       └── getNiche()
        └── releaseTmpItems()

Heuristic::iteratedGreedy()
└── leftBottomHeuristic()
    ├── selectSkyline(..., leftBottom)
    ├── liftSkyline()
    └── addItemOverSkyline()
        └── detectAndMergeSkylines()
            └── removeSkyline()

Heuristic::bestFitHeuristic()
├── std::sort(..., compareItemByWidth)
├── selectSkyline(..., bestFit)
├── findBestItem()
├── addItemOverSkyline()
└── liftSkyline()
```

---

## 6. 整体算法流程

### 6.1 主流程 (main.cpp)

```
1. 遍历 ./2sp/ 目录下的所有文件
2. 对每个文件：
   a. 调用 readData() 读取物品列表和条带宽度 W
   b. 创建 BLEU 实例：BLEU alg(allItems, W, 20, 1000)
      → TrialHeight=20，timeLimit=1000
   c. 调用 alg.evaluate() 执行算法，得到 solutionStatus
   d. 输出状态结果
   e. 释放物品内存
```

> **注意**：当前 `main.cpp` 中 `for (int i = 0; i < 1; ++i)` 只运行一轮，且 `Heuristic` 对象虽然声明但未被使用。

---

### 6.2 BLEU 算法初始化流程

```
BLEU 构造函数
│
├── 1. 按宽度（降序）排序所有物品
├── 2. reassignItemsIdx()：重新分配 idx（从0开始连续）
├── 3. preprocessing()：见 6.3
└── 4. bounds()：见 6.4
```

---

### 6.3 预处理流程 (preprocessing)

根据论文第 2.2 节，共三步：

```
preprocessing()
│
├── Step 1: preprocessingFixItems()
│   目的：找出所有"宽度过大"的物品（不能与任何其他物品并排）
│   条件：item.width + minWidth > W
│   处理：
│     - 这些物品必须独占整个条带宽度，其总高度 _processedH 可直接从试验高度中扣除
│     - 其余物品加入 _processedItems，并分配连续的 idxHelper
│
├── Step 2: preprocessingReduceW()
│   目的：缩减有效条带宽度
│   方法：对 _processedItems 的所有宽度运行 subSetSum(_allWidths, W)
│         → _processedW = 所有物品宽度能凑出的 ≤W 的最大值
│   意义：若物品宽度不能完全填满 W，则实际有效宽度可以更小
│
└── Step 3: preprocessingModifyItemWidth()
    目的：放大部分物品的宽度以减少搜索空间
    方法：对每个物品 i，计算其余物品宽度能凑出的 ≤ (_processedW - w_i) 的最大值 maxFit
          若 w_i + maxFit < _processedW，则 w_i += (_processedW - maxFit - w_i)
    意义：物品的有效宽度可以扩展，而不影响可行性
```

还有一步在 `evaluate()` 中调用的高度预处理：

```
preprocessItemHeight(tmpItems, binHeight, binWidth)
│
├── 类似 Step 3，对高度方向进行同样的宽度修改
│   对每个物品 i，若其余物品的高度之和加 h_i < binHeight，则 h_i 扩展
├── 找最小高度 minHeight
├── 删除"孤立"物品（h_i + minHeight > binHeight，即无法与任何物品同列）
│   同时相应减少 binWidth
└── 重新分配 idxHelper（从0开始）
```

---

### 6.4 下界计算流程 (bounds)

共五个下界，取最大值：

| 下界 | 方法 | 复杂度 |
|------|------|--------|
| **LB1**（面积下界） | `ceil(Σ(w_i × h_i) / W) + _processedH` | O(n) |
| **LB2**（双可行函数） | 应用3种双可行函数变换宽度后计算装箱高度下界；α 从 1 到 W/2 枚举 | O(n × W) |
| **LB3**（Alvarez-Vales） | 贪心模拟：每次选"最紧凑"物品（min(W-w, H-h)最小）放入矩形，缩小矩形直到所有物品放入 | O(n² log n) |
| **LB4**（列生成LP） | 用 CPLEX 求解非连续箱装问题 LP，通过 DP 解定价子问题迭代添加列 | 多项式但较慢 |
| **LB5**（平行机排程LP）| 用 computeFX 计算可行位置集，然后用 CPLEX 求解 LP 松弛 | 依赖 CPLEX |

```
bounds()
├── lb1 = LowerBound1()
├── lb2 = LowerBound2()
├── lb4 = LowerBound4()          ← LB4 依赖 _bestLowerBound，故先计算其他
├── _bestLowerBound = max(lb1, lb2, lb4)
├── lb3 = LowerBound3()          ← LB3 以 _bestLowerBound 为初始搜索点
├── lb5 = LowerBound5()          ← LB5 以 _bestLowerBound 为计算基础
└── _bestLowerBound = max(_bestLowerBound, lb3, lb5)
```

---

### 6.5 evaluate 评估主流程

```
evaluate()
│
├── 快速返回条件：
│   ├── _processedItems 为空 → feasible
│   └── _bestLowerBound > _trialHeight → infeasible
│
├── 设置 binWidth = _processedW, binHeight = _trialHeight - _processedH
├── 复制 _processedItems 为 tmpItems（可修改）
├── preprocessItemHeight(tmpItems, binHeight, binWidth)
│   若 binWidth < 0 → infeasible
│
├── 将 tmpItems 转换为 const vector<const item*> Items
└── branchAndBound(Items, binWidth, binHeight)
    → 返回 feasible / infeasible / pending
```

---

### 6.6 分支定界算法 (branchAndBound)

```
branchAndBound(Items, binWidth, binHeight)
│
├── 初始化：
│   ├── 创建根节点（BBNode）：所有物品未放置，所有列高度为0
│   ├── 计算最大探索节点数（完美装箱 vs 非完美装箱）
│   └── DFS 栈 dfsTree
│
└── DFS 循环（while 栈非空 且 节点数 < 上限）：
    ├── 取出栈顶节点 currentNode
    │
    ├── 若 remainingItems 为空（叶节点）：
    │   └── 调用 yCheckAlgorithm(W, H, positions, Items)
    │       ├── 返回 true → 解可行，return feasible
    │       └── 返回 false → 继续搜索
    │
    ├── 若 bounding(currentNode) 为真 → 剪枝，continue
    │   bounding() 包含：
    │   ├── 对称性剪枝（相同形状物品按 idx 排序）
    │   ├── 面积下界剪枝（剩余面积 > 剩余空间）
    │   └── dynamicCuts()（基于 DP 的精细面积剪枝）
    │
    └── makeBranch(currentNode, dfsTree)：生成子节点
        ├── 选择最左列 selectedColumn = leftMostIdx
        ├── 若该列已有物品：生成"跳过该列"子节点（列高度设为trialHeight，leftMostIdx+1）
        └── 对每个可放入该列的物品 j（按 idx 升序）：
            ├── 检查宽度是否超出条带
            ├── 检查高度是否超出试验高度
            ├── 检查 maxiItemIdxColumns 对称性约束
            └── 创建子节点：记录物品 j 的位置，更新列高度，可能推进 leftMostIdx

若探索节点数 ≥ 上限 → return pending
否则（栈为空）→ return infeasible
```

`branchAndBoundYRelax` 与 `branchAndBound` 的区别：叶节点不调用 `yCheckAlgorithm`，直接返回 feasible（求解 x 坐标松弛版本，用于 `solvePCC`）。

---

### 6.7 y-check 算法

y-check 算法的目的：给定所有物品的 **x 坐标分配**（即每个物品从哪列开始放），判断能否找到满足高度约束的 **y 坐标分配**，使解合法。

```
yCheckAlgorithm(processedW, TrialHeight, itemPositions, processedItems)
│
├── preprocess4yCheck(binWidth, items, coords, Height)
│   ├── preprocessedFirst4yCheck()   ← 合并可互换的相邻物品
│   │   └── 循环：对每个物品，尝试将其左侧/右侧相邻物品合并进来
│   │       → mergeItems4yCheck()（递归检验子问题可行性后合并）
│   ├── preprocessedSecond4yCheck()  ← 将物品宽度扩展至相邻物品之间的间距
│   └── preprocessedThird4yCheck()   ← 合并具有完全相同占用物品集合的相邻列
│
└── yCheckEnumerationTree(items, coords, H, W)
    │
    ├── 创建根节点（BBNode，初始列高度全0，物品坐标已设定）
    └── DFS 枚举树：
        ├── 若 remainingItems 为空 → feasible
        ├── yCheckBounding(currentNode) → 若剪枝条件满足则跳过
        │   包含：
        │   ├── 列高度溢出检查（当前 + 剩余物品高度 > trialHeight）
        │   ├── 对称性剪枝（同形状物品在同列时约束放置顺序）
        │   └── 连续放置对称性约束
        └── yCheckMakeBranch(currentNode, yEnTree)
            ├── getNiche()：找出"凹槽"（高度最低的连续列段 [l, r]）
            ├── 对所有 x 坐标落在凹槽内的剩余物品，生成各自的放置子节点
            │   ├── 每个子节点：物品 y = 当前列高，更新对应列的高度
            │   └── 剪枝条件5：若存在更靠左且更矮的物品可放，跳过当前物品
            └── 若所有放置均会超出凹槽两侧高度，生成"空凹槽"子节点
                （将凹槽列高提升至两侧最小高度）
```

---

### 6.8 启发式算法流程

#### leftBottomHeuristic

```
1. 初始化天际线（单段，宽=W，高=0）
2. 按输入顺序逐个处理物品：
   a. 选最低天际线段
   b. 若段宽 < 物品宽，抬升该段（合并到相邻较低段），回到 a
   c. 将物品放在该段上方，更新天际线
3. 解析天际线，返回最大高度
```

#### bestFitHeuristic

```
1. 按宽度降序排序物品
2. 循环（直到所有物品放完）：
   a. 选最低天际线段
   b. 在剩余物品中找宽度 ≤ 段宽的第一个物品（最宽者）
   c. 若找到：将物品放在该段上方
   d. 若未找到：抬升该段
3. 返回最大高度
```

#### iteratedGreedy

```
1. 调用 leftBottomHeuristic 得到初始解高度 primalBound
2. 重复（直到无改进）：
   for i in [0, n):
     for j in [0, n) where j ≠ i:
       将 items[i] 移到位置 j
       计算新高度
       若高度下降：记录改进，break 内层循环
       否则：撤销移动
3. 返回最优高度
```

---

## 7. 测试数据格式

所有测试文件为纯文本格式：

```
<n>               第1行：物品总数
<W>               第2行：条带宽度
<idx> <w> <h>    第3行起：物品编号、宽度、高度（空白分隔）
...
```

示例（`HT01.TXT`，16个物品，宽度20）：

```
          16
          20
           1           2          12
           2           7          12
           ...
          16          11           2
```

### 测试用例集

| 系列 | 数量 | 说明 |
|------|------|------|
| HT | 09 | Hopper & Turton 经典实例 |
| BENG | 10 | Bengtsson 实例 |
| CGCUT | 03 | Christofides & Whitlock 切割实例 |
| GCUT | 04 | Gilmore & Gomory 切割实例 |
| NGCUT | 12 | Non-guillotine 切割实例 |

---

## 8. 关键参数与常量

| 常量/参数 | 值 | 位置 | 说明 |
|----------|-----|------|------|
| `BigNumber` | 999999 | `spp.h` | 无穷大替代值（高度/宽度上限） |
| `tolerance`（knapsack） | 0.0001 | `knapsack.h` | 浮点精度容差 |
| `BLEU::tolerance` | 0.0001 | `BLEU.cpp` | 浮点精度容差 |
| `BLEU::bigNumber` | 999999 | `BLEU.cpp` | 无穷大替代值 |
| `BLEU::BBMaxExplNodesPerPack` | 10,000,000 | `BLEU.cpp` | 完美装箱时 B&B 节点上限 |
| `BLEU::BBMaxExplNodesNonPerPack` | 80,000 | `BLEU.cpp` | 非完美装箱时 B&B 节点上限 |
| `BLEU::ycheckExplNode` | 10,000,000 | `BLEU.cpp` | y-check 枚举树节点上限 |
| TrialHeight（示例） | 20 | `main.cpp` | 初始试验高度（可调） |
| timeLimit（示例） | 1000 | `main.cpp` | 时间限制（可调） |

---

## 9. 参考文献

1. **Boschetti, M.A., Montaletti, L.** (2010). *An exact algorithm for the two-dimensional strip packing problem*. Operations Research, 58(6), 1774–1791.

2. **Alvarez-Valdes, R., Parreno, F., Tamarit, J.M.** (2009). *A branch and bound algorithm for the strip packing problem*. OR Spectrum, 31(2), 431–459.

3. **Fekete, S.P., Schepers, J.** (2004). *A combinatorial characterization of higher-dimensional orthogonal packing*. Mathematics of Operations Research, 29(2), 353–368.

4. **Clautiaux, F., Carlier, J., Moukrim, A.** (2007). *A new exact method for the two-dimensional orthogonal relaxation problem of the strip packing problem*. European Journal of Operational Research, 183(3), 1196–1211.

5. **Vanderbeck, F.** (1999). *Computational study of a column generation algorithm for bin packing and cutting stock problems*. Mathematical Programming, 86(3), 565–594.

6. **Labbé, M., et al.** *New classes of fast lower bounds for bin packing problems* (双可行函数下界参考).
