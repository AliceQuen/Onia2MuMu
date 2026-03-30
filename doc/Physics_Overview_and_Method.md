# Physics Overview and Method

This file preserves technical and reference material from the previous `README.md`.

## 整体介绍

在CMS实验上测量质子-质子对撞中一类三介子联合产生过程的反应截面。

这一反应截面部分由三部分子散射（Triple Parton Scattering, TPS）贡献，因而可用于研究质子中部分子动力学分布的关联，并且能为未来超出标准模型物理的搜索提供更为精确的本底描述。

本软件包将用于从`MINIAOD`原始数据中初步进行筛选，将候选事例汇总并存储相关信息，以备后续进一步筛选。

## 反应过程

$pp \rightarrow J/\psi + J/\psi + \Upsilon \rightarrow 6\mu$

$pp \rightarrow J/\psi + J/\psi + \phi \rightarrow 4\mu + 2K$

$pp \rightarrow J/\psi + \Upsilon + \phi \rightarrow 4\mu + 2K$

## 软件基础

使用CERN的`ROOT`框架，以及CMS合作组软件框架`CMSSW_13_0_2`。

从github.com/AliceQuen/Onia2MuMu继承了代码框架。相关代码由清华大学致理书院2021级本科生秦俊凯提供。

## 处理方法总览

### 事件选择

#### 1. 将$\mu^+\mu^-$配对：仅处理几何条件

按照电荷配对并从轨迹拟合顶点，检查是否能够给出一定结果，不限制具体`VtxProb`。

#### 2. 将$K^+K^-$配对：仅处理几何条件

从带电径迹中排除$\mu^{\pm}$， 将余下的部分都视为$K^{\pm}$。按照电荷配对并从轨迹拟合顶点，检查是否能够给出一定结果，不限制具体`VtxProb`。

#### 3. 从$\mu^+\mu^-$对创建候选$J/\psi$和$\Upsilon$：使用部分动力学信息

从通过第1步筛选的$\mu^+\mu^-$对出发，使用不变质量进行粗选：

* 视为$J/\psi$候选：满足 $ 2 \mathrm{GeV/c^2} \leq m_{\mu^+\mu^-} \leq 6 \mathrm{GeV/c^2}$

* 视为$\Upsilon$候选：满足 $ 8 \mathrm{GeV/c^2} \leq m_{\mu^+\mu^-} \leq 12 \mathrm{GeV/c^2}$

#### 4. 从$K^+K^-$对创建候选$\phi$：使用部分动力学信息

从通过第2步筛选的$K^+K^-$对出发，使用不变质量进行粗选：

要求 $ 0.5 \mathrm{GeV/c^2} \leq m_{K^+K^-} \leq 1.5 \mathrm{GeV/c^2}$。

#### 5. 从候选$J/\psi$, $\Upsilon$和$\phi$寻找整体事例候选

仅处理几何条件：从上一步筛选的$J/\psi$, $\Upsilon$和$\phi$出发，在确保末态粒子不相重叠之后，通过顶点拟合，寻找整体事例候选。

### 效率评定

（待完善，或另起repository进行）

### 系统误差

（待完善，或另起repository进行）

## 参考资料

### 1. 如何重建$\phi$

查阅PDG，并参考S. Leontsinis等人的报告（见群聊），选择以$\phi \rightarrow K^+ + K^-$为重建$\phi$的衰变道。

CMS本身不具备直接鉴别出$K^{\pm}$的能力，因而在分析中，我们从带电径迹中选择两条非$\mu$的径迹并直接假定为$K^{\pm}$。

我们将需要挑选出电性相反的这样两条径迹，并进行顶点拟合。凡是能够进行拟合，并且重建不变质量$m_{KK}$满足

$$
0.5 \mathrm{GeV/c^2} \leq m_{KK} \leq 1.5 \mathrm{GeV/c^2}
$$

我们将把这样的径迹对重建结果认为是$\phi$的候选。

我们可以参考$J/\psi + \psi(2S)$过程的分析代码：这一过程需要从一个$J/\psi$和一对$\pi^+\pi^-$重建出$\psi(2S)$，而实际上也只是假定所有的径迹中除了$\mu$以外的径迹都是$\pi^{\pm}$。

因而在我们的分析中，直接处理的代码$K^{\pm}$，相较于原来处理$\pi^{\pm}$的部分，几乎只需要将质量更换为$m_{K^{\pm}}=493.677 \pm 0.015\mathrm{MeV/c^2}$（PDG 2024）。

### 2. 作为参考的代码（来自`Onia2MuMu`代码包）




### 3. 既有代码实现说明

#### 通过`using`声明给出的变量别名

##### (1)  `muon_t`

存储单个$\mu^{\pm}$信息，为`RefCountedKinematicParticle`类型别名。

```cpp
using muon_t  = RefCountedKinematicParticle;
```

##### (2) `muList_t`

存储一系列$\mu^{\pm}$信息，为`std::pair< vector<muon_t>, vector<uint> >`类型别名。

```cpp
using muList_t = std::pair< vector<muon_t>, vector<uint> >;
```

在这一`std::pair`中，第一个元素为`muon_t`列表，第二个元素为`unsigned int`列表，用于存储$\mu^{\pm}$在其列表`thePATMuonHandle`中的下标。**这一下标可以被用于后续排除“重复使用末态粒子”的情况。**

在本分析中，我们将使用`muList_t`来存储$\mu^+\mu^-$对的信息；在这一`std::pair`中，前后两个`vector`的长度都将为2。

#### 实用自定义函数

##### (1) `isOverlapPair`：判断两个$\mu^{\pm}$是否重叠

```cpp
bool MultiLepPAT::isOverlapPair(
    const muList_t& arg_MuonPair1,
    const muList_t& arg_MuonPair2
    )
```

###### 参数

* `arg_MuonPair1`：第一个$\mu^+\mu^-$对，存储在`muList_t`类型的变量中。
* `arg_MuonPair2`：第二个$\mu^+\mu^-$对，存储在`muList_t`类型的变量中。

值得注意的是，这里已经假定了`muList_t`对应的$\mu^{\pm}$的个数为2。

###### 返回值

`bool`类型，为`true`时表示两个$\mu^{\pm}$对重叠，为`false`时表示两个$\mu^{\pm}$对不重叠。

##### (2) `particlesToVtx`：将粒子列表拟合为顶点

这个函数一共引入了三个重载版本，分别用于不同的情况。

代码实现中已经有`try-catch`语句，可以保证不会throw（"no-throw guarantee"）。

```cpp
bool MultiLepPAT::particlesToVtx(
    const vector<RefCountedKinematicParticle>& arg_FromParticles
)
bool MultiLepPAT::particlesToVtx(
    const vector<RefCountedKinematicParticle>& arg_FromParticles,
    const string&                              arg_Message
)
bool MultiLepPAT::particlesToVtx(
    RefCountedKinematicTree&                   arg_VertexFitTree,
    const vector<RefCountedKinematicParticle>& arg_FromParticles,
    const string&                              arg_Message
)
```

###### 参数

* `arg_FromParticles`：待拟合的粒子列表，存储在`vector<RefCountedKinematicParticle>`类型的变量中。
* `arg_VertexFitTree`：拟合结果，存储在相应位置所给定的`RefCountedKinematicTree`（“拟合树”）类型的变量中。
* `arg_Message`：拟合过程中，遇程序错误时输出到标准输出流的描述信息，存储在`string`类型的变量中。

###### 返回值

`bool`类型，为`true`时表示拟合成功，为`false`时表示拟合失败。

##### (3) `extractFitRes`：从`RefCountedKinematicTree`类型的变量中提取信息

这个函数一共引入了三个重载版本，可以按实际需求提取不同信息。

```cpp
bool MultiLepPAT::extractFitRes(
    RefCountedKinematicTree&     arg_VtxTree,
    RefCountedKinematicParticle& res_Part,
    RefCountedKinematicVertex&   res_Vtx,
    KinematicParameters&         res_Param,
    double&                      res_MassErr
)
bool MultiLepPAT::extractFitRes(
    RefCountedKinematicTree&     arg_VtxTree,
    RefCountedKinematicParticle& res_Part,
    RefCountedKinematicVertex&   res_Vtx,
    double&                      res_MassErr
)
bool MultiLepPAT::extractFitRes(
    RefCountedKinematicTree&     arg_VtxTree,
    RefCountedKinematicVertex&   res_Vtx,
    double&                      res_VtxProb
)
```

###### 参数

* `arg_VtxTree`：待提取拟合结果的最初“拟合树”（`RefCountedKinematicTree`类型变量）。
* `res_Part`：提取的初始粒子信息，存储相应位置在所给定的`RefCountedKinematicParticle`类型的变量中。
* `res_Vtx`：提取的顶点信息，存储在相应位置所给定的`RefCountedKinematicVertex`类型的变量中。
* `res_Param`：提取的初始粒子动力学参数信息，存储在相应位置所给定的`KinematicParameters`类型的变量中。
* `res_MassErr`：提取的初始粒子质量不确定度信息，存储在相应位置所给定的`double`类型的变量中。
* `res_VtxProb`：提取的顶点拟合概率信息，存储在相应位置所给定的`double`类型的变量中。

###### 返回值

`bool`类型，为`true`时表示提取成功，为`false`时表示提取失败。

提示：若需要提取`res_MassErr`，则返回的`bool`值实际为这个误差平方值是否非负。除此之外，其他返回值均为`true`。

##### (4) `getDynamics` ：计算$p_T, \eta, \phi$等动力学信息。

这个函数一共引入了两个重载版本，可以按实际需求从不同数据类型提取不同信息。

```cpp
void MultiLepPAT::getDynamics(
    double  arg_mass,
    double  arg_px,  double  arg_py,  double  arg_pz,
    double& res_pt,  double& res_eta, double& res_phi
)
void MultiLepPAT::getDynamics(
    const RefCountedKinematicParticle& arg_Part,
    double& res_pt,  double& res_eta, double& res_phi
)
```

###### 参数

* `arg_mass`：粒子质量，存储在`double`类型的变量中。
* `arg_px`，`arg_py`，`arg_pz`：粒子动量分量，分别存储在`double`类型的变量中。
* `arg_Part`：待提取动力学信息的粒子，存储在`RefCountedKinematicParticle`类型的变量中。
* `res_pt`：提取的**横动量**$p_T$信息，存储在相应位置所给定的`double`类型的变量中。
* `res_eta`：提取的**赝快度**$\eta$信息，存储在相应位置所给定的`double`类型的变量中。
* `res_phi`：提取的**方位角**$\phi$信息，存储在相应位置所给定的`double`类型的变量中。

###### 返回值

无返回值

### 4. 作为参考的代码

#### (1) 初始化`theTrackHandle`，提取出带电粒子径迹

```cpp
edm::Handle<edm::View<pat::PackedCandidate>> theTrackHandle; //  MINIAOD
iEvent.getByToken(trackToken_, theTrackHandle);              //  MINIAOD
std::vector<edm::View<pat::PackedCandidate>::const_iterator> nonMuonPionTrack;
```

#### (2) 初始化`nonMuonPionTrack`：将所有的带电径迹包含在内

```cpp
// Copy tracks iterators
for (edm::View<pat::PackedCandidate>::const_iterator iTrackc = theTrackHandle->begin(); // MINIAOD
iTrackc != theTrackHandle->end(); ++iTrackc)
{
	nonMuonPionTrack.push_back(iTrackc);
}
```

#### (3) 建立`nonMuonPionTrack`：将带电径迹中的$\mu^{\pm}$排除

```cpp
if (thePATMuonHandle->size() >= 4)
{
	// fill muon track block
	for (auto iMuonP  = thePATMuonHandle->begin(); //  MINIAOD
	          iMuonP != thePATMuonHandle->end(); ++iMuonP)
	{
		/*********************************************************
	     * Some muon info processing here. We skip this part.
	    *********************************************************/

		// Find and delete muon Tracks in PionTracks
		for (auto iTrackfID  =  nonMuonPionTrack.begin(); // MINIAOD
		           iTrackfID != nonMuonPionTrack.end(); ++iTrackfID)
		{
			if(iMuonP->track().isNull())
			{
				continue;
			}
			edm::View<pat::PackedCandidate>::const_iterator iTrackf = *(iTrackfID);
			iMuonP->track()->px();

	        /******************************************************************
	         * 通过比较动量确定是否为同一径迹，从而排除muon track。
	         * 注意到CMS通过径迹可以直接给出动量，不需要静止质量或其他参量。
	        ******************************************************************/

			if (   iTrackf->px() == iMuonP->track()->px()
	            && iTrackf->py() == iMuonP->track()->py()
	            && iTrackf->pz() == iMuonP->track()->pz())
			{
				nonMuonPionTrack.erase(iTrackfID);
				iTrackfID = iTrackfID - 1;
			}
		}
	    /******************************************************************
	     * Some trigger matching here. We skip this part.
	    ******************************************************************/
	}
	/******************************************************************
	 * Some trigger matching here. We skip this part.
	******************************************************************/
}
```

#### (4) 枚举`nonMuonPionTrack`：对带电径迹进行配对

注：为了提升可读性，将一些高度复杂的变量类型用`auto`声明。经测试，这是可以直接平移到真实使用代码的。

```cpp
for (auto iTrack1ID  = nonMuonPionTrack.begin(); // MINIAOD
          iTrack1ID != nonMuonPionTrack.end(); ++iTrack1ID){
	edm::View<pat::PackedCandidate>::const_iterator iTrack1 = *(iTrack1ID);
	if (iTrack1->pt() < pionptcut){
		continue;
	}
	for (auto iTrack2ID  = iTrack1ID + 1; // MINIAOD
	          iTrack2ID != nonMuonPionTrack.end(); ++iTrack2ID){
		edm::View<pat::PackedCandidate>::const_iterator iTrack2 = *(iTrack2ID);
		if (iTrack2->pt() < pionptcut)
		{
			continue;
		}
	    /******************************************************************
	     * 校验电荷条件
	    ******************************************************************/
		if ((iTrack1->charge() + iTrack2->charge()) != 0)
			continue;

		// MINIAOD begin
		if (!iTrack1->hasTrackDetails() || iTrack1->charge() == 0)
		{
			continue;
		}
		if (!iTrack2->hasTrackDetails() || iTrack2->charge() == 0)
		{
			continue;
		}
		// MINIAOD end

	    /******************************************************************
	     * 从一个pi+pi-对和一个J/psi（重建自mu+mu-对）重建psi(2S)。
	     * 为此，先进行方向交角筛选。
	    ******************************************************************/

		TLorentzVector P4_Track1, P4_Track2, P4_Jpsipipi;
		P4_Track1.SetPtEtaPhiM(iTrack1->pt(), iTrack1->eta(), iTrack1->phi(), myKMass);
		P4_Track2.SetPtEtaPhiM(iTrack2->pt(), iTrack2->eta(), iTrack2->phi(), myKMass);
		P4_Jpsipipi = P4_mu1 + P4_mu2 + P4_Track1 + P4_Track2;

		if (P4_Track1.DeltaR(P4_Jpsipipi) > pionDRcut)
		{
			continue;
		}
		if (P4_Track2.DeltaR(P4_Jpsipipi) > pionDRcut)
		{
			continue;
		}

	    /******************************************************************
	     * 建立一个pi+pi-对的TransientTrack
	     * 由此构建一个RefCountedKinematicParticle
	     * 我们需要一个KinematicParticleFactoryFromTransientTrack对象
	     * 这个对象将会用于从TransientTrack转换得到RefCountedKinematicParticle
	     * 我们还将需要一个KinematicParticleVertexFitter
	     * 这个对象将会用于拟合得到的粒子
	    ******************************************************************/

		TransientTrack trackTT1(*(iTrack1->bestTrack()), &(bFieldHandle)); // MINIAOD
		TransientTrack trackTT2(*(iTrack2->bestTrack()), &(bFieldHandle)); // MINIAOD
		KinematicParticleFactoryFromTransientTrack JPiPiFactory;
		// The mass of a muon and the insignificant mass sigma
		// to avoid singularities in the covariance matrix.
		ParticleMass pion_mass = myKMass; // pdg mass
		float pion_sigma = myKMasserr;
		// initial chi2 and ndf before kinematic fits.
		float chi = 0.;
		float ndf = 0.;

	    /******************************************************************
	     * 下面开始进行顶点拟合
	    ******************************************************************/

		// mass constrain for psi2s from X6900, now noMC
		// first fit to have a psi2s
		vector<RefCountedKinematicParticle> JPiPiParticles;
		JPiPiParticles.push_back(JPiPiFactory.particle(trackTT1, pion_mass, chi, ndf, pion_sigma));
		JPiPiParticles.push_back(JPiPiFactory.particle(trackTT2, pion_mass, chi, ndf, pion_sigma));
		JPiPiParticles.push_back(pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
		JPiPiParticles.push_back(pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));

	    /******************************************************************
	     * 下面的这个部分相当于我们之前写的代码中的particlesToVtx函数。
	    ******************************************************************/

		KinematicParticleVertexFitter JPiPi_fitter;
		RefCountedKinematicTree JPiPiVertexFitTree;
		Error_t = false;
		try
		{
			JPiPiVertexFitTree = JPiPi_fitter.fit(JPiPiParticles);
		}catch(...)
		{
			Error_t = true;
			std::cout<<"error at JPiPinoMC"<<std::endl;
		}
		if (Error_t || !(JPiPiVertexFitTree->isValid()))
		{
			continue;
		}

	    /******************************************************************
	     * 下面的这个部分相当于我们之前写的代码中的extractFitRes函数。
	    ******************************************************************/

		JPiPiVertexFitTree->movePointerToTheTop();
		RefCountedKinematicParticle JPiPi_vFit_noMC = JPiPiVertexFitTree->currentParticle();
		RefCountedKinematicVertex JPiPi_vFit_vertex_noMC = JPiPiVertexFitTree->currentDecayVertex();

		double JPiPi_vtxprob = ChiSquaredProbability((double)(JPiPi_vFit_vertex_noMC->chiSquared()), (double)(JPiPi_vFit_vertex_noMC->degreesOfFreedom()));

	    /******************************************************************
	     * 之后是一些后续处理，于我们的分析不直接有关，略去。
	    ******************************************************************/
    }
}
```
