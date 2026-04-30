# CMS Onia2MuMu 分析 - MultiLepPAT

## 项目简介

`MultiLepPAT` 是 CMS 实验中用于寻找新型 X 粒子的 EDAnalyzer 模块，通过重建 X(3872) 粒子的衰变过程来探索潜在的新共振态。本分析专注于四μ子加两π子的末态拓扑，利用运动学顶点拟合和质量约束技术精确测量粒子性质。

## 物理目标

本项目的物理目标包括：

1. **寻找 X(3872) 粒子与 J/ψ 粒子的联合产生**
2. **寻找 X 粒子衰变到 J/ψ 和 ψ(2S) 的过程**
3. **寻找 X 粒子衰变到两个 ψ(2S) 的过程**

### 衰变拓扑

```
                    X [候选共振态]
                       │
           ┌───────────┴───────────┐
           │                       │
      ψ(2S) / X(3872)             J/ψ / ψ(2S)
    [mu1+mu2+pi1+pi2]            [mu3+mu4]
           │
     ┌─────┴─────┐
     │           │
   J/ψ₁        π⁺π⁻
 [mu1+mu2]   [pi1+pi2]
```

最终形成的末态为：**4 个 μ 子 + 2 个 π 子**

其中 X 粒子为未知共振态，可能为 X(6900) 或 X(7100)。

## 三种质量假设

为了全面探索潜在的新共振态，程序采用三种质量约束假设：

| 假设编号 | 质量约束配置 | 物理意义 |
|---------|------------|---------|
| **假设 1** | mumupipi 约束到 ψ(2S) (3.6861 GeV)<br>另一个 mumu 约束到 J/ψ (3.0969 GeV) | 寻找 X → ψ(2S) J/ψ 过程 |
| **假设 2** | mumupipi 约束到 X(3872) (3.872 GeV)<br>另一个 mumu 约束到 J/ψ (3.0969 GeV) | 寻找 X(3872) 和 J/ψ 的联合产生过程 |
| **假设 3** | mumupipi 约束到 ψ(2S) (3.6861 GeV)<br>另一个 mumu 约束到 ψ(2S) (3.6861 GeV) | 寻找 X → ψ(2S) ψ(2S) 过程 |

## 核心功能模块

```
MultiLepPAT 主模块
├── 事件预处理层
│   ├── HLT 触发匹配
│   ├── 主顶点获取
│   └── 事件级预筛选
├── 粒子选择层
│   ├── 径迹预选择
│   ├── μ子预选择
│   └── 触发过滤器匹配
├── 粒子组合层
│   ├── μ+μ- 对组合
│   ├── π+π- 对组合
│   ├── J/ψ - ππ 组合
│   └── J/ψ - J/ψ 组合
├── 顶点拟合层
│   ├── 普通顶点拟合
│   └── 质量约束顶点拟合
├── 物理量计算层
│   ├── 四动量计算
│   ├── 顶点概率计算
│   ├── δR 关联计算
│   └── 冲击参数计算
└── 数据存储层
    ├── TTree 分支管理
    ├── 候选数据填充
    └── 变量重置
```

## 粒子筛选条件

### 径迹预选择

| 条件 | 阈值 | 物理意义 |
|-----|------|---------|
| 高纯度标志 | highPurity | 径迹质量标识 |
| 横动量 | > 0.5 GeV | 最低动量要求 |
| 赝快度 | < 2.4 | 探测器接收范围 |
| 动量相对误差 | < 0.1 | 动量分辨率要求 |
| 有效击中数 | ≥ 10 | 径迹探测完整性 |
| 拟合优度 | < 0.18 | 径迹拟合质量 |

### μ子预选择

| 条件 | 要求 | 说明 |
|-----|------|------|
| 软μ子标识 | `isSoftMuon()` | 基于径迹的软μ子识别算法 |
| 分段 p_T-η 阈值 | 见下表 | 不同快度区域的动量要求 |
| 电荷平衡 | 2 个正 μ + 2 个负 μ | 形成两个 J/ψ 的必要条件 |

分段 p_T(η) 函数：

```
p_T^cut(η) =
  3.5 GeV,           when |η| < 1.2
  5.47 - 1.89*|η|,  when 1.2 < |η| < 2.1
  1.5 GeV,           when 2.1 < |η| < 2.4
```

## 程序运行流程

### 生命周期总览

```
程序启动
   │
   ▼
构造函数
├── 读取配置参数
├── 初始化 EDToken
└── 触发器名称排序优化
   │
   ▼
beginJob()
├── 创建 TTree
└── 创建所有分支 (~130 个)
   │
   ▼
事件循环 (每个事件)
├── 数据获取与预筛选
├── HLT 触发匹配
├── 主顶点获取
├── 径迹预选择与电荷分组
├── μ子预选择与电荷分组
├── J/ψ 候选预缓存
├── 嵌套循环组合与拟合
├── 三种质量假设约束拟合
├── 物理量计算
└── TTree Fill()
   │
   ▼
endJob()
   │
   ▼
程序结束
```

### analyze() 函数筛选层级

```
开始
  ↓
[事件级] Muon 数量 >= 4?
  ↓ 否 → 返回
  ↓ 是
[事件级] HLT 触发匹配?
  ↓ 否 → 返回
  ↓ 是
[事件级] 有效主顶点?
  ↓ 否 → 返回
  ↓ 是
[粒子级] 径迹预选择后数量 >= 2?
  ↓ 否 → 返回
  ↓ 是
[粒子级] μ子预选择后数量 >= 4?
  ↓ 否 → 返回
  ↓ 是
[组合级] μ₁+μ₂ 电荷和为 0?
  ↓ 否 → continue
  ↓ 是
[组合级] J/ψ₁ 不变质量在窗口?
  ↓ 否 → continue
  ↓ 是
[拟合级] J/ψ₁ 顶点拟合成功?
  ↓ 否 → continue
  ↓ 是
[拟合级] J/ψ₁ 顶点概率 > 0.01?
  ↓ 否 → continue
  ↓ 是
[组合级] π₁+π₂ 组合筛选
  ↓
[拟合级] ψ(2S) 4 径迹拟合
  ↓
[拟合级] ψ(2S) 顶点概率 > 0.005
  ↓
[拟合级] ψ(2S) p_T > 4 GeV
  ↓
[约束级] 三种质量假设约束拟合
  ↓
[计算级] 所有物理量计算
  ↓
[存储级] TTree 填充
结束
```

## 物理量存储结构

### TTree 分支层级

**TTree 名称**: `X_data`

**总分支数**: ~130 个

```
X_data TTree
├── 事件信息层 (4)
│   ├── runNum, evtNum, lumiNum, nGoodPrimVtx
│
├── X 候选层 (24)
│   ├── 假设 1: X_PJ_mass, X_PJ_VtxProb, X_PJ_massErr, X_PJ_pt, X_PJ_px, ...
│   ├── 假设 2: X_XJ_mass, X_XJ_VtxProb, X_XJ_massErr, X_XJ_pt, X_XJ_px, ...
│   └── 假设 3: X_PP_mass, X_PP_VtxProb, X_PP_massErr, X_PP_pt, X_PP_px, ...
│
├── ψ(2S) 候选层 (9)
│   └── Psi2S_mass, Psi2S_VtxProb, Psi2S_massErr, Psi2S_pt, Psi2S_absEta, ...
│
├── J/ψ 候选层 (16)
│   ├── Jpsi1_mass, Jpsi1_VtxProb, Jpsi1_massErr, Jpsi1_pt, ...
│   └── Jpsi2_mass, Jpsi2_VtxProb, Jpsi2_massErr, Jpsi2_pt, ...
│
├── μ子信息层 (68)
│   ├── μ₁ (16 个变量): pt, eta, charge, trackIso, d0BS, d0PV, dzPV, ...
│   ├── μ₂ (16 个变量)
│   ├── μ₃ (16 个变量)
│   ├── μ₄ (16 个变量)
│   └── μ子 ID 计数: nLooseMuons, nTightMuons, nSoftMuons, nMediumMuons
│
├── π子信息层 (10)
│   ├── π₁: pi1_pt, pi1_px, pi1_py, pi1_pz, pi1_absEta
│   └── π₂: pi2_pt, pi2_px, pi2_py, pi2_pz, pi2_absEta
│
└── δR 关联层 (14)
    ├── dR_mu1_mu2, dR_mu3_mu4, dR_pi1_pi2
    ├── dR_Psi2S_Jpsi1, dR_Psi2S_Jpsi2
    └── dR_Psi2S_pi1, dR_Psi2S_pi2, ...
```

### 无效值处理规范

所有物理量在计算失败或无效时，统一设置为哨兵值：

| 数据类型 | 无效值 | 示例 |
|---------|-------|------|
| 浮点型 (float) | `-999.0` | `mu1_pt = -999.0` |
| 整型 (int/unsigned int) | `0` | `nGoodPrimVtx = 0` |
| 布尔型 (bool) | `false` | `mu1_hasFilterMatch = false` |

## 代码结构

### 主要文件

```
Onia2MuMu/
├── interface/
│   └── MultiLepPAT.h         # 头文件，类定义与 JpsiCandidate 结构体
├── src/
│   └── MultiLepPAT.cc        # 源文件，实现代码
└── test/
    └── runMultiLepPAT_miniAOD.py   # CMSSW 配置文件
```

### 关键数据结构：JpsiCandidate

定义于 `interface/MultiLepPAT.h`，用于缓存 J/ψ 候选的拟合结果，避免在多层循环中重复计算：

```cpp
struct JpsiCandidate {
    // 粒子迭代器
    edm::View<pat::Muon>::const_iterator muPlus;
    edm::View<pat::Muon>::const_iterator muMinus;
    
    // 普通顶点拟合结果
    double mass, vtxProb, massErr;
    ROOT::Math::PxPyPzMVector p4;
    RefCountedKinematicParticle kinematicParticle;
    RefCountedKinematicVertex vertex;
    
    // 质量约束拟合结果
    bool hasConstraintFit = false;
    double constraintMass, constraintVtxProb, constraintMassErr;
    RefCountedKinematicParticle constraintParticle;
    RefCountedKinematicVertex constraintVertex;
    
    // 候选类型标记
    bool isJpsiCandidate = false;
    bool isPsi2SCandidate = false;
    
    // 其他信息
    reco::TransientTrack muonTT1, muonTT2;
    bool filterMatchPlus, filterMatchMinus;
};
```

## 核心优化策略

1. **触发器名称预排序**：按长度升序排列，短模式优先匹配，平均提前退出循环

2. **电荷预分组**：在循环前将径迹和 μ 子按电荷正负分组，减少无效组合

3. **J/ψ 候选预缓存**：在嵌套循环外预先计算所有有效 J/ψ 候选，避免重复拟合

4. **Cheap Cut First**：将计算成本低的筛选（如快速质量计算）放在循环最内层，提前排除无效组合

## 编译与运行

### 环境要求

- CMSSW 版本（根据实际环境调整）
- ROOT 6+
- gcc 7+

### 编译步骤

```bash
# 进入 CMSSW 工作目录
cd $CMSSW_BASE/src

# 克隆代码（如需要）
git clone <repository-url>

# 编译
scram b
```

### 运行

```bash
# 使用配置文件运行
cmsRun test/runMultiLepPAT_miniAOD.py
```

## 文档版本

- 文档版本: v1.0
- 创建日期: 2026.04.30
- 适用代码版本: MultiLepPAT (2024 重构版)
