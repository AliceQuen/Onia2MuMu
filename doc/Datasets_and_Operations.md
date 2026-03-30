# Datasets and Operations

This file preserves dataset and operational reference material from the previous `README.md`. Some release-specific notes are historical and may not match the current tested CMSSW release.

## 数据集性质

### Run 3：主要数据来源

本次分析将主要使用Run 3的数据，利用较大数据量寻找以上所列各个相对稀有的过程。

结合本课题目标过程，选用Run 3中2022以及2023年`ParkingDoubleMuonLowMass[1-7]`系列数据集。

#### 关于Run2023C-PromptReco-v2的提示

在Run2023C数据采集阶段，曾在RAW与AOD之间出现过一些数据处理问题。

这个问题的影响是，在提交crab job时，可能需要特别指定Run2023C-v2这一部分数据的run range为$[367516, 367620]$。但是重复数据文件应该已经删除，而且从CMS整体分析工作布局来看，这一问题应该已经解决。

具体信息可以参照[Run3 PdmV twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#2023_Era_definition)，以及相关[CMS Talk post](https://cms-talk.web.cern.ch/t/t0-reprocessing-due-to-wrong-raw-version/24528)和[工单](https://its.cern.ch/jira/browse/CMSTZ-1040)。

#### `TODO`

通过Run2023C-v2与Run2023C-v3的数据进行对比，确认是否有出现重复run。可以用处理结果文件的run number做进行对比。

很可能重复数据已经被剔除了，但是仍然需要检查run number确认。

### Run 2：辅助确认

在编写得到$pp \rightarrow J/\psi + J/\psi + \phi \rightarrow 4\mu + 2K$代码之后，将会使用Run 2数据进行分析，用于复现日内瓦大学S. Leontsinis等人结果；通过交叉比对，确认代码可靠性。

## 数据查询实用工具：CMS DAS命令行界面

运行前需要设置好VOMS代理。调用时，需要使用语句：

```bash
/cvmfs/cms.cern.ch/common/dasgoclient -query="YOUR_DAS_QUERY"
```

其中，查询语句`YOUR_DAS_QUERY`与在网页上查看CMS DAS无异。具体使用可以参考 https://cmsweb.cern.ch/das/cli 。

## 所使用数据集

### Run 2

这部分与S. Leontsinis等人的报告一致。

#### 2016

```text
/Charmonium/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD
/Charmonium/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD
/Charmonium/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD
/Charmonium/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD
/Charmonium/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD
/Charmonium/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD
/Charmonium/Run2016G-21Feb2020_UL2016-v1/MINIAOD
```

#### 2017

```text
/Charmonium/Run2017B-UL2017_MiniAODv2-v1/MINIAOD
/Charmonium/Run2017C-UL2017_MiniAODv2-v1/MINIAOD
/Charmonium/Run2017D-UL2017_MiniAODv2-v1/MINIAOD
/Charmonium/Run2017E-UL2017_MiniAODv2-v1/MINIAOD
/Charmonium/Run2017F-UL2017_MiniAODv2-v1/MINIAOD
```

#### 2018

```text
/Charmonium/Run2018A-12Nov2019_UL2018_rsb-v1/MINIAOD
/Charmonium/Run2018B-UL2018_MiniAODv2-v1/MINIAOD
/Charmonium/Run2018C-UL2018_MiniAODv2-v1/MINIAOD
/Charmonium/Run2018D-12Nov2019_UL2018-v1/MINIAOD
```

### Run 3

这部分综合了`test/crabData/Run2023dataList.txt`（cr. 王晰宁）以及PdmV查询结果。

#### 2022 （需使用`CMSSW_12_4_x`）

```text
/ParkingDoubleMuonLowMass0/Run2022C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass0/Run2022D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass0/Run2022D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass0/Run2022E-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass0/Run2022F-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass0/Run2022G-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass1/Run2022C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass1/Run2022D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass1/Run2022D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass1/Run2022E-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass1/Run2022G-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass2/Run2022C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass2/Run2022D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass2/Run2022D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass2/Run2022E-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass2/Run2022F-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass2/Run2022G-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass3/Run2022C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass3/Run2022D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass3/Run2022D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass3/Run2022E-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass3/Run2022F-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass3/Run2022G-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass4/Run2022C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass4/Run2022D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass4/Run2022D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass4/Run2022E-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass4/Run2022F-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass4/Run2022G-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass5/Run2022C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass5/Run2022D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass5/Run2022D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass5/Run2022E-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass5/Run2022F-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass5/Run2022G-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass6/Run2022C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass6/Run2022D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass6/Run2022D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass6/Run2022E-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass6/Run2022F-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass6/Run2022G-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass7/Run2022C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass7/Run2022D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass7/Run2022D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass7/Run2022E-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass7/Run2022F-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass7/Run2022G-PromptReco-v1/MINIAOD
```

#### 2023 （需使用`CMSSW_13_0_x`）

```text
/ParkingDoubleMuonLowMass0/Run2023B-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v3/MINIAOD
/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v4/MINIAOD
/ParkingDoubleMuonLowMass0/Run2023D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass0/Run2023D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass1/Run2023B-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v3/MINIAOD
/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v4/MINIAOD
/ParkingDoubleMuonLowMass1/Run2023D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass1/Run2023D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass2/Run2023B-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v3/MINIAOD
/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v4/MINIAOD
/ParkingDoubleMuonLowMass2/Run2023D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass2/Run2023D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass3/Run2023B-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v3/MINIAOD
/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v4/MINIAOD
/ParkingDoubleMuonLowMass3/Run2023D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass3/Run2023D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass4/Run2023B-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v3/MINIAOD
/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v4/MINIAOD
/ParkingDoubleMuonLowMass4/Run2023D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass4/Run2023D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass5/Run2023B-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v3/MINIAOD
/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v4/MINIAOD
/ParkingDoubleMuonLowMass5/Run2023D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass5/Run2023D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass6/Run2023B-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v3/MINIAOD
/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v4/MINIAOD
/ParkingDoubleMuonLowMass6/Run2023D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass6/Run2023D-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass7/Run2023B-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v2/MINIAOD
/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v3/MINIAOD
/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v4/MINIAOD
/ParkingDoubleMuonLowMass7/Run2023D-PromptReco-v1/MINIAOD
/ParkingDoubleMuonLowMass7/Run2023D-PromptReco-v2/MINIAOD
```
