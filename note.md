# 10.13讨论

## 数据从哪来
http://cmsweb.cern.ch/
DoubleMuonLowMass
一个AOD/MiniAOD->很多.root
http://cmshitinfo.app.cern.ch/
dataset 搜索
在服务器上一个CMSSW环境下，src文件夹中代码
需要有三级文件夹：src->..->TPS..
上次md文件 &#9733; &#9734;
在TPS文件夹中`scram b`

crab 大批量 研究.sh

# 11.5

新分支复制的是王驰原代码

## 问题
456~459行 确有问题，学长有更新
```python
for (unsigned int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++)
		{
			JpsiMatchTrig[JpsiTrig] = 0;
		} // Initiating Jpsi trigger
		for (int itrig = 0; itrig < ntrigs; itrig++)
		{
			string trigName = triggerNames_.triggerName(itrig);
			int hltflag = (*hltresults)[itrig].accept();
			TrigRes->push_back(hltflag);
			TrigNames->push_back(trigName);
			for (unsigned int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++)
			{
				std::regex pattern(".*"+TriggersForJpsi_[JpsiTrig]+".*");
				if (std::regex_search(trigName, pattern))
				{
					JpsiMatchTrig[JpsiTrig] = hltflag;
					bool isDumplicate = false;
					if(hltflag)
					{
						for(unsigned int MatchTrig = 0; MatchTrig < MatchJpsiTrigNames->size(); MatchTrig ++)
						{
							if(trigName == MatchJpsiTrigNames->at(MatchTrig))
							{
								isDumplicate = true;
								break;
							}
						}
					}
					if(!isDumplicate)
					{
						MatchJpsiTrigNames->push_back(trigName);
					}		
					break;
				}
			} // Jpsi Trigger
		}
```

536行似乎多余

601，602也许需要更改

看到607行

## 读代码学到的新概念

### Primary Vtx, Reconstructed Vtx, BeamSpot Vtx三者有什么区别
1. Primary Vertex（主顶点）
定义：主顶点是指粒子碰撞的初始位置，即事件中所有初级粒子（primary particles）产生的地方。

特点：

主顶点是事件的起点，通常是所有重建轨迹的共同交点。
在一个事件中，可能存在多个主顶点，特别是在高亮度碰撞中（pile-up events），多个碰撞事件可能在同一个束交点发生。
主顶点的位置通常通过重建算法从检测器数据中提取。
用途：

用于确定事件的起点，帮助重建粒子的轨迹和动量。
在数据分析中，用于区分不同的碰撞事件。
2. Reconstructed Vertex（重建顶点）
定义：重建顶点是通过重建算法从检测器数据中提取的顶点，表示粒子轨迹的交点。

特点：

重建顶点可以是主顶点，也可以是次级顶点（secondary vertices），如粒子衰变产生的顶点。
重建顶点的位置和误差通过重建算法计算得出，依赖于检测器的分辨率和数据质量。
用途：

用于重建粒子的轨迹和动量。
帮助识别和分析粒子的衰变过程。
3. BeamSpot Vertex（束斑顶点）
定义：束斑顶点是粒子束在碰撞点附近的平均位置，表示粒子束的空间分布特性。

特点：

束斑顶点通常由加速器提供，表示粒子束在碰撞点的中心位置和尺寸。
束斑顶点的位置和尺寸是通过束流监测器测量的，反映了粒子束的横向分布。
用途：

用于校准和优化重建算法，提供粒子束的参考位置。
在重建顶点时，束斑顶点可以作为初始猜测值，提高重建精度。

### 关于ndof
评估拟合质量：

ndof是评估拟合质量的重要指标之一。通常，拟合的质量可以通过卡方（chi-square）或归一化卡方（chi-square per degree of freedom）来评估。归一化卡方定义为：
[ \chi^2/\text{ndof} ]

其中，(\chi^2)是卡方值。理想情况下，(\chi^2/\text{ndof})的值应该接近1，这表明拟合是合理的。如果这个值远大于1，说明拟合可能不佳，数据点与拟合模型之间的偏差较大；如果这个值远小于1，可能表明数据点的误差被高估了。

顶点重建的可靠性：

在顶点重建中，ndof越高，通常意味着拟合使用了更多的测量点，因此重建的顶点位置更可靠。一个高ndof值表明有更多的粒子轨迹参与了顶点的重建，这通常对应于一个更精确和可信的顶点位置。
数据筛选：

在数据分析中，研究人员常常使用ndof来筛选高质量的顶点。例如，在你的代码中，只有那些ndof大于或等于5的顶点才被认为是“好”的顶点。这种筛选可以排除那些基于少量测量点重建的、不可靠的顶点，从而提高分析结果的质量。

### 顶点关于z方向和rho的约束
fabs((*recVtxs)[myi].z()) <= 24：检查顶点在z轴上的位置是否在±24单位以内。z轴通常是沿着粒子束的方向，限制z轴位置可以排除那些远离碰撞点的顶点。
fabs((*recVtxs)[myi].position().rho()) <= 2.0：检查顶点在径向位置（ρ）上的距离是否在2.0单位以内。ρ是顶点在xy平面上的径向距离，限制ρ可以排除那些偏离碰撞点的顶点。

为什么不对拟合出来的顶点做这些限制？

### 不同Muon的区别

#### LooseMuon
LooseMuon工作点是最宽松的选择标准，目的是尽可能多地保留重建的μ子。它通常用于需要高效率的分析，但对背景抑制要求不高的情况。LooseMuon的选择标准可能包括：

基本的轨迹质量要求，例如轨迹的拟合质量（χ²/ndof）。
基本的顶点一致性要求，例如与主要顶点的距离。
基本的探测器响应要求，例如在不同子探测器中的命中数量。

#### TightMuon
TightMuon工作点是最严格的选择标准，目的是最大限度地抑制背景噪声，确保选中的μ子具有高纯度。它通常用于需要高背景抑制的精确测量和分析。TightMuon的选择标准可能包括：

更严格的轨迹质量要求，例如更低的χ²/ndof。
更严格的顶点一致性要求，例如更小的与主要顶点的距离。
更严格的探测器响应要求，例如在所有子探测器中都有足够的命中数量。
额外的隔离要求，例如周围没有过多的其他粒子。

#### SoftMuon
SoftMuon工作点是专门为低动量μ子设计的选择标准，目的是在低动量区域保留更多的μ子。它通常用于研究低动量μ子的物理过程，例如重味物理（heavy flavor physics）中的B介子衰变。SoftMuon的选择标准可能包括：

适中的轨迹质量要求，允许较高的χ²/ndof。
适中的顶点一致性要求，允许较大的与主要顶点的距离。
适中的探测器响应要求，允许在某些子探测器中缺少命中。
低动量区域的特定要求，例如动量小于某个阈值。

#### MediumMuon
MediumMuon工作点是介于LooseMuon和TightMuon之间的选择标准，目的是在效率和纯度之间取得平衡。它通常用于需要适度背景抑制的分析。MediumMuon的选择标准可能包括：

中等严格的轨迹质量要求，例如中等的χ²/ndof。
中等严格的顶点一致性要求，例如中等的与主要顶点的距离。
中等严格的探测器响应要求，例如在大多数子探测器中都有足够的命中数量。

### 关于横向冲击参数（d0）和三维冲击参数(d3)
横向冲击参数（d0）定义为粒子轨迹在最接近光束位置（Beam Spot）或主要顶点（Primary Vertex）处的径向距离。具体来说，d0是粒子轨迹在xy平面上与光束位置的最小距离。这个参数可以用来衡量粒子轨迹的偏移程度。

三维冲击参数（d3）定义为粒子轨迹在三维空间中与参考点的最小距离。具体来说，d3是粒子轨迹在三维空间中与光束位置或主要顶点的最小距离。这个参数可以用来衡量粒子轨迹的偏移程度。

### 关于轨迹隔离(trackIso)
轨迹隔离（trackIso）通常定义为在目标粒子周围一定半径内，所有其他粒子的动量之和。这个半径通常在横向平面（xy平面）上定义。轨迹隔离的计算步骤如下：

选择目标粒子：例如，一个重建的μ子。
定义隔离半径：通常在横向平面上定义一个圆形区域，半径为R。
计算动量之和：在这个圆形区域内，计算所有其他粒子的动量之和。

### 关于$\Delta R$
ΔR定义为粒子在η-φ平面上的距离，计算公式如下：

$ \Delta R = \sqrt{(\Delta \eta)^2 + (\Delta \phi)^2} $

其中：

Δη是两个粒子在伽马角（pseudo-rapidity）方向上的差值。
Δφ是两个粒子在方位角（azimuthal angle）方向上的差值。

# 11.9

## 问题
拟合后对拟合得到粒子加cut和拟合前对mu的和加cut有什么不同

688~715 mass constrain的作用？

最后是在找一个X使得X衰变到X3872和Jpsi

mass constrain可能在计算测量截面时会使用

# 11.11

初步完成改写
未debug

# 11.20
debug完成，可以开始跑，但是拟合得到的phi质量过于离奇（普遍在TeV量级），以及处理一个Event的速度过于之慢。完整程序的有效性需要提交crab来处理。

# 11.21
王驰提出可以先处理muon，先判定一次获得的J/psi有没有两两共顶点，若没有直接跳过整段分析过程。

另外要注意最后三个循环中似乎可以通过把track放外层来避免多次拟合。


