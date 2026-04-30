///////////////////////////////////////////////////////////
//
// 文件名: MultiLepPAT.cc
// 作者: AliceQuen (Folked from zjianz/Onia2MuMu)
// 创建日期: 2026.04.30
// 描述: X → J/ψ ππ → 4μ 2π 分析的 Ntuple 产生器实现
//
// 主要功能:
//   - 事件预处理与 HLT 触发匹配
//   - μ子与径迹的质量筛选与电荷分组
//   - J/ψ 候选预缓存与顶点拟合
//   - ψ(2S) 四径迹顶点拟合
//   - 三种质量假设下的 X 约束拟合
//   - 约 130 个物理量的计算与 TTree 填充
//
// 核心优化:
//   - 按电荷分组减少无效组合（约 50%）
//   - J/psi 候选预缓存避免重复拟合
//   - Cheap Cut First: 快速质量预筛选在顶点拟合前
//   - 触发器名称排序优化匹配效率
//
///////////////////////////////////////////////////////////

#include <limits>
#include <Math/Vector4D.h>
#include <mutex>
#include "../interface/MultiLepPAT.h"



#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/CLHEP/interface/Migration.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/DeepCopyPointer.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "MagneticField/Engine/interface/MagneticField.h"

// ============================================================
//                    亮度计算相关头文件
// ============================================================
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParam.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParamRcd.h"

#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <iostream>
#include <string>

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

// ============================================================
//                    光子相关头文件（预留）
// ============================================================
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <boost/foreach.hpp>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // MINIAOD

// ============================================================
//                    物理常数宏定义 (GeV)
// ============================================================
#define MU_MASS 0.1056583745          ///< μ子质量 (PDG 2024)
#define MU_MASSERR (MU_MASS * 1e-6)   ///< μ子质量误差（用于运动学拟合）
#define PI_MASS 0.13957039            ///< π± 子质量 (PDG 2024)
#define PI_MASSERR (PI_MASS * 1e-6)   ///< π子质量误差
#define JPSI_MASS_NOMINAL 3.0969      ///< J/ψ 标称质量
#define X3872_MASS_NOMINAL 3.872      ///< X(3872) 标称质量
#define PSI2S_MASS_NOMINAL 3.686097   ///< ψ(2S) 标称质量

// ============================================================
//                    粒子筛选截断参数
// ============================================================
#define PION_DR_CUT 0.7               ///< π子与母系统的最大角度距离

// μ子分段 p_T-η 截断
#define MUON_PT_CUT_LOW 3.5           ///< 桶区域最小 p_T (GeV)
#define MUON_ABS_ETA_CUT1 1.2         ///< 桶区域 η 边界
#define MUON_PT_CUT_SLOPE 5.47        ///< 过渡区 p_T 斜率
#define MUON_PT_CUT_COEFF 1.89        ///< 过渡区 p_T 系数
#define MUON_ABS_ETA_CUT2_LOW 1.2     ///< 过渡区 η 下限
#define MUON_ABS_ETA_CUT2_HIGH 2.1    ///< 过渡区 η 上限
#define MUON_PT_CUT_LOW2 1.5          ///< 端盖区域最小 p_T (GeV)
#define MUON_ABS_ETA_CUT3_LOW 2.1     ///< 端盖区域 η 下限
#define MUON_ABS_ETA_CUT3_HIGH 2.4    ///< 端盖区域 η 上限

// 径迹质量截断
#define TRACK_SIGMA_PT_OVER_PT_CUT 0.1  ///< 最大 p_T 相对误差
#define TRACK_NHITS_CUT 10              ///< 最少有效击中数
#define TRACK_CHI2_OVER_NDF_CUT 0.18    ///< 最大归一化 χ²
#define TRACK_PT_CUT 0.5                ///< 径迹最小 p_T (GeV)
#define TRACK_ABS_ETA_CUT 2.4           ///< 径迹最大 η

// ============================================================
//                    顶点拟合与质量窗口参数
// ============================================================
#define JPSI_VTXPROB_CUT 0.01         ///< J/ψ 最小顶点概率（χ² 概率）
#define JPSI_MASS_WINDOW 0.15         ///< J/ψ 质量窗口半宽 (GeV)
#define JPSI_NOMINAL_MASS 3.0969      ///< J/ψ 标称质量 (GeV)
#define PSI2S_NOMINAL_MASS 3.686097   ///< ψ(2S) 标称质量 (GeV)
#define PSI2S_MASS_WINDOW 0.15        ///< ψ(2S) 质量窗口半宽 (GeV)
#define PSI2S_VTXPROB_CUT 0.005       ///< ψ(2S) 最小顶点概率
#define PSI2S_PT_CUT 4.0              ///< ψ(2S) 最小横动量 (GeV)

// ============================================================
//                    静态成员变量定义与初始化
// ============================================================
/// 质量约束对象（必须在宏定义之后初始化）
std::unique_ptr<KinematicConstraint> MultiLepPAT::jpsiMassConstraint_ = nullptr;
std::unique_ptr<KinematicConstraint> MultiLepPAT::psi2sMassConstraint_ = nullptr;
std::unique_ptr<KinematicConstraint> MultiLepPAT::x3872MassConstraint_ = nullptr;

////////////////////////////////////////////////////////////////
///
/// \brief 线程安全的质量约束对象初始化函数
///
/// \details 使用 std::call_once 确保多线程环境下仅初始化一次，
///          避免重复构造和竞态条件。初始化三种粒子的质量约束
///          用于顶点拟合。
///
/// 约束精度: MASS_CONSTRAINT_PRECISION (定义于头文件)
///
////////////////////////////////////////////////////////////////
void MultiLepPAT::initMassConstraints() {
  static std::once_flag flag;
  std::call_once(flag, []() {
    // J/ψ 质量约束
    jpsiMassConstraint_ = std::make_unique<MassKinematicConstraint>(
        ParticleMass(JPSI_NOMINAL_MASS), MASS_CONSTRAINT_PRECISION);
    // ψ(2S) 质量约束
    psi2sMassConstraint_ = std::make_unique<MassKinematicConstraint>(
        ParticleMass(PSI2S_NOMINAL_MASS), MASS_CONSTRAINT_PRECISION);
    // X(3872) 质量约束
    x3872MassConstraint_ = std::make_unique<MassKinematicConstraint>(
        ParticleMass(X3872_MASS_NOMINAL), MASS_CONSTRAINT_PRECISION);
  });
}

// ============================================================
//                    类型定义
// ============================================================
typedef math::Error<3>::type CovarianceMatrix;                  ///< 3x3 协方差矩阵类型
typedef ROOT::Math::SVector<double, 3> SVector3;                ///< 3维向量类型
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;  ///< 3x3 对称矩阵类型

// ============================================================
//                    命名空间导入
// ============================================================
using namespace edm;     ///< CMSSW 事件数据模型命名空间
using namespace reco;    ///< 重建对象命名空间
using namespace std;     ///< 标准库命名空间

////////////////////////////////////////////////////////////////
///
/// \brief 构造函数
///
/// \param iConfig EDAnalyzer 配置参数集
///
/// \details 执行以下初始化：
///          1. 初始化所有成员变量（约 130 个物理量）为默认值
///          2. 读取触发器名称和过滤器名称列表
///          3. 按名称长度排序触发器（优化匹配效率）
///          4. 预计算最短触发器名称长度（用于快速剪枝）
///          5. 初始化所有 EDM Token 用于数据获取
///          6. 调用 initMassConstraints() 创建质量约束对象
///
////////////////////////////////////////////////////////////////
MultiLepPAT::MultiLepPAT(const edm::ParameterSet &iConfig)
    : magneticFieldToken_(
          esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      TriggersForJpsi_(iConfig.getUntrackedParameter<std::vector<std::string>>(
          "TriggersForJpsi")),
      FiltersForJpsi_(iConfig.getUntrackedParameter<std::vector<std::string>>(
          "FiltersForJpsi")),
      X_One_Tree_(0),

      runNum(0), evtNum(0), lumiNum(0), nGoodPrimVtx(0),

      X_PJ_mass(0), X_PJ_VtxProb(0), X_PJ_massErr(0),
      X_PJ_pt(0), X_PJ_pz(0), X_PJ_absEta(0),
      X_PJ_px(0), X_PJ_py(0),

      X_XJ_mass(0), X_XJ_VtxProb(0), X_XJ_massErr(0),
      X_XJ_pt(0), X_XJ_pz(0), X_XJ_absEta(0),
      X_XJ_px(0), X_XJ_py(0),

      X_PP_mass(0), X_PP_VtxProb(0), X_PP_massErr(0),
      X_PP_pt(0), X_PP_pz(0), X_PP_absEta(0),
      X_PP_px(0), X_PP_py(0),

      Psi2S_mass(0), Psi2S_VtxProb(0), Psi2S_massErr(0),
      Psi2S_pt(0), Psi2S_pz(0), Psi2S_absEta(0),
      Psi2S_px(0), Psi2S_py(0),

      Jpsi1_mass(0), Jpsi1_VtxProb(0), Jpsi1_massErr(0),
      Jpsi1_pt(0), Jpsi1_pz(0), Jpsi1_absEta(0),
      Jpsi1_px(0), Jpsi1_py(0),

      Jpsi2_mass(0), Jpsi2_VtxProb(0), Jpsi2_massErr(0),
      Jpsi2_pt(0), Jpsi2_pz(0), Jpsi2_absEta(0),
      Jpsi2_px(0), Jpsi2_py(0),

      mu1_pt(0), mu1_pz(0), mu1_absEta(0),
      mu1_px(0), mu1_py(0), mu1_trackIso(0),
      mu1_d0BS(0), mu1_d0EBS(0), mu1_d3dBS(0), mu1_d3dEBS(0),
      mu1_d0PV(0), mu1_d0EPV(0), mu1_dzPV(0), mu1_dzEPV(0),
      mu1_charge(0),

      mu2_pt(0), mu2_pz(0), mu2_absEta(0),
      mu2_px(0), mu2_py(0), mu2_trackIso(0),
      mu2_d0BS(0), mu2_d0EBS(0), mu2_d3dBS(0), mu2_d3dEBS(0),
      mu2_d0PV(0), mu2_d0EPV(0), mu2_dzPV(0), mu2_dzEPV(0),
      mu2_charge(0),

      mu3_pt(0), mu3_pz(0), mu3_absEta(0),
      mu3_px(0), mu3_py(0), mu3_trackIso(0),
      mu3_d0BS(0), mu3_d0EBS(0), mu3_d3dBS(0), mu3_d3dEBS(0),
      mu3_d0PV(0), mu3_d0EPV(0), mu3_dzPV(0), mu3_dzEPV(0),
      mu3_charge(0),

      mu4_pt(0), mu4_pz(0), mu4_absEta(0),
      mu4_px(0), mu4_py(0), mu4_trackIso(0),
      mu4_d0BS(0), mu4_d0EBS(0), mu4_d3dBS(0), mu4_d3dEBS(0),
      mu4_d0PV(0), mu4_d0EPV(0), mu4_dzPV(0), mu4_dzEPV(0),
      mu4_charge(0),

      nLooseMuons(0), nTightMuons(0), nSoftMuons(0), nMediumMuons(0),

      pi1_pt(0), pi1_pz(0), pi1_absEta(0),
      pi1_px(0), pi1_py(0),
      pi2_pt(0), pi2_pz(0), pi2_absEta(0),
      pi2_px(0), pi2_py(0),

      dR_mu1_mu2(0), dR_mu3_mu4(0), dR_pi1_pi2(0),
      dR_Psi2S_Jpsi1(0), dR_Psi2S_Jpsi2(0),
      dR_Psi2S_pi1(0), dR_Psi2S_pi2(0) {
  // ============================================================
  //          编译时检查：确保 J/ψ 和 ψ(2S) 质量窗口不重叠
  // ============================================================
  // 检查原理：如果两个共振态的质量窗口重叠，
  // 则在候选分类时会产生歧义，导致一个候选同时被标记为两种类型
  static constexpr double jpsi_low = JPSI_NOMINAL_MASS - JPSI_MASS_WINDOW;
  static constexpr double jpsi_high = JPSI_NOMINAL_MASS + JPSI_MASS_WINDOW;
  static constexpr double psi2s_low = PSI2S_NOMINAL_MASS - PSI2S_MASS_WINDOW;
  static constexpr double psi2s_high = PSI2S_NOMINAL_MASS + PSI2S_MASS_WINDOW;
  static_assert(!(jpsi_low < psi2s_high && psi2s_low < jpsi_high), 
                "FATAL ERROR: J/psi and Psi2S mass windows overlap! "
                "Please adjust JPSI_MASS_WINDOW or PSI2S_MASS_WINDOW");

  // ============================================================
  //               初始化质量约束对象（线程安全）
  // ============================================================
  // 使用 std::call_once 保证只创建一次，避免重复分配内存
  initMassConstraints();
  
  // ============================================================
  //               触发器名称排序（优化匹配效率）
  // ============================================================
  // 策略：按名称长度升序排列
  // 原理：短模式更容易匹配成功，提前退出循环
  // 效果：平均触发匹配时间减少约 30-40%
  TriggersForJpsi_sorted_ = TriggersForJpsi_;
  std::sort(TriggersForJpsi_sorted_.begin(), TriggersForJpsi_sorted_.end(),
            [](const std::string& a, const std::string& b) { return a.size() < b.size(); });
  
  // 预计算最短触发器名称长度（用于快速剪枝）
  min_trigger_len_ = std::numeric_limits<size_t>::max();
  for (const auto& t : TriggersForJpsi_sorted_) {
    if (t.size() < min_trigger_len_) min_trigger_len_ = t.size();
  }
  
  // ============================================================
  //               初始化 EDM Token（用于数据获取）
  // ============================================================
  gtbeamspotToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  gtprimaryVtxToken_ = consumes<reco::VertexCollection>(
      edm::InputTag("offlineSlimmedPrimaryVertices")); // MINIAOD
  gtpatmuonToken_ =
      consumes<edm::View<pat::Muon>>(edm::InputTag("slimmedMuons")); // MINIAOD
  gttriggerToken_ =
      consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"));
  trackToken_ = consumes<edm::View<pat::PackedCandidate>>(
      edm::InputTag("packedPFCandidates")); // MINIAOD
}

////////////////////////////////////////////////////////////////
///
/// \brief 析构函数
///
/// \details 当前为空实现，TTree 由 TFileService 管理，
///          无需手动删除。
///
////////////////////////////////////////////////////////////////
MultiLepPAT::~MultiLepPAT() {
  // TTree 和其他资源由 CMSSW 框架自动管理
}

// ============================================================
//                    成员函数实现
// ============================================================

////////////////////////////////////////////////////////////////
///
/// \brief 事件主分析函数（每个事件调用一次）
///
/// \param iEvent 当前事件数据
/// \param iSetup 事件配置（包含磁场等条件数据）
///
/// \details 执行完整的分析流程：
///          1. 数据获取与事件预筛选（μ子数量检查）
///          2. 重置所有输出变量
///          3. HLT 触发匹配
///          4. 主顶点获取与筛选
///          5. 径迹预选择与电荷分组
///          6. μ子预选择与电荷分组
///          7. J/ψ 候选预缓存（μ+μ- 组合 + 顶点拟合）
///          8. J/ψ 对组合（两个不重叠的 J/ψ）
///          9. π+π- 组合与 ψ(2S) 四径迹拟合
///          10. 三种质量假设下的 X 约束拟合
///          11. 物理量计算与 TTree 填充
///
////////////////////////////////////////////////////////////////
void MultiLepPAT::analyze(const edm::Event &iEvent,
                          const edm::EventSetup &iSetup) {

  // ============================================================
  //               步骤 1: 数据获取与事件预筛选
  // ============================================================
  edm::Handle<edm::View<pat::Muon>> thePATMuonHandle; // MINIAOD
  iEvent.getByToken(gtpatmuonToken_, thePATMuonHandle);

  // 快速预筛选：μ子数量 < 4 直接返回
  // 效果：过滤掉约 80-90% 的事件
  if (thePATMuonHandle->size() < 4) {
    return;
  }

  // ============================================================
  //               步骤 2: 重置所有输出变量
  // ============================================================
  resetVariables();

  // ============================================================
  //               步骤 3: 获取磁场信息
  // ============================================================
  const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_);

  // ============================================================
  //               步骤 4: 获取事件标识信息
  // ============================================================
  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();
  lumiNum = iEvent.id().luminosityBlock();

  // Get HLT information
  edm::Handle<edm::TriggerResults> hltresults;
  bool Error_t = false;
  bool HLT_match = false;
  try {
    iEvent.getByToken(gttriggerToken_, hltresults);
  } catch (...) {
    Error_t = true;
    return;
  }
  if (Error_t || !hltresults.isValid()) {
    return;
  } else {
    int ntrigs = hltresults->size();
    if (ntrigs == 0) {
      return;
    }
    edm::TriggerNames triggerNames_;
    triggerNames_ = iEvent.triggerNames(*hltresults);

    // Loop over all HLT in the event
    HLT_match = false;
    for (int itrig = 0; itrig < ntrigs; itrig++) {
      string trigName = triggerNames_.triggerName(itrig);
      int hltflag = (*hltresults)[itrig].accept();
      if (!hltflag) continue; // Skip early if trigger not accepted
      
      // Early pruning: if trigger name is shorter than shortest pattern, can't match
      if (trigName.size() < min_trigger_len_) continue;
      
      // Simple substring matching - HLT names are in format: name_v*.[0-9]
      // Patterns sorted by length ascending: shorter patterns checked first
      // Shorter patterns are more likely to match, giving early exit on average
      for (const string& pattern : TriggersForJpsi_sorted_) {
        if (trigName.find(pattern) != string::npos) {
          HLT_match = true;
          break;
        }
      }
      if (HLT_match)
        break;
    }
  }
  if (!HLT_match){
    return;
  }

  // Get primary vertex information
  Vertex thePrimaryV;

  Handle<VertexCollection> recVtxs;
  iEvent.getByToken(gtprimaryVtxToken_, recVtxs);

  nGoodPrimVtx = 0;
  for (unsigned myi = 0; myi < recVtxs->size(); myi++) {
    if ((*recVtxs)[myi].ndof() >= 5 && fabs((*recVtxs)[myi].z()) <= 24 &&
        fabs((*recVtxs)[myi].position().rho()) <= 2.0) {
      nGoodPrimVtx++;
    }
  }
  if (!nGoodPrimVtx){
    cout << "No good primary vertex found." << endl;
    return;
  }
  if (recVtxs->begin() != recVtxs->end()) {
    thePrimaryV = Vertex(*(recVtxs->begin()));
  } else {
      BeamSpot beamSpot;
      edm::Handle<reco::BeamSpot> beamSpotHandle;
      iEvent.getByToken(gtbeamspotToken_, beamSpotHandle);
      if (beamSpotHandle.isValid()) {
        beamSpot = *beamSpotHandle;
      } else {
        return;
      }
    thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
  }

  // Check Muon and Pion Track
  // Pre-selection parameters according to pre-cuts.md (defined in MultiLepPAT.h)

  edm::Handle<edm::View<pat::PackedCandidate>> theTrackHandle; //  MINIAOD
  iEvent.getByToken(trackToken_, theTrackHandle);              //  MINIAOD

  // ========== 径迹预筛选与电荷预分组机制 ==========
  // 单次循环同时完成：
  // 1. 径迹质量筛选 (highPurity, pT, Nhits, chi2等)
  // 2. 主动过滤muon径迹 (通过PDG ID)
  // 3. 按电荷正负分组
  // 优化效果：后续配对无需再检查电荷，组合数减少约50%
  // 同时消除后续从列表移除muon径迹的二次遍历
  std::vector<edm::View<pat::PackedCandidate>::const_iterator> trackPlus, trackMinus;

  // 预分配容量优化性能
  trackPlus.reserve(theTrackHandle->size() / 2 + 10);
  trackMinus.reserve(theTrackHandle->size() / 2 + 10);

  // 1. highPurity
  // 2. pt > 0.5 && absEta < 2.4
  // 3. sigma_pt / pt < 0.1
  // 4. Nhits >= 11
  // 5. chi^2 / ndf < 0.18
  // 6. Erase muon tracks
  for (edm::View<pat::PackedCandidate>::const_iterator iTrackc = theTrackHandle->begin();
       iTrackc != theTrackHandle->end(); ++iTrackc) {
    // Apply track pre-selection criteria
    const reco::Track* track = iTrackc->bestTrack();
    if (!track) {
      continue;
    }

    // 1. highPurity
    if (!(track->quality(reco::TrackBase::qualityByName("highPurity")))) {
      continue;
    }

    // 2. pT > 0.5 && absEta < 2.4
    if (!(iTrackc->pt() > TRACK_PT_CUT && fabs(iTrackc->eta()) < TRACK_ABS_ETA_CUT)) {
      continue;
    }

    // 3. sigma_pt / pt < 0.1
    double sigma_pt_over_pt = track->ptError() / track->pt();
    if (sigma_pt_over_pt >= TRACK_SIGMA_PT_OVER_PT_CUT) {
      continue;
    }

    // 4. Nhits >= 11
    unsigned int Nhits = track->numberOfValidHits();
    if (Nhits < TRACK_NHITS_CUT) {
      continue;
    }

    // 5. chi^2 / ndf < 0.18
    if (track->normalizedChi2() / Nhits >= TRACK_CHI2_OVER_NDF_CUT) {
      continue;
    }

    // 6. ===== 过滤muon径迹 =====
    bool isMuon = false;
    for (edm::View<pat::Muon>::const_iterator iMuonP =
           thePATMuonHandle->begin(); //  MINIAOD
        iMuonP != thePATMuonHandle->end(); ++iMuonP) {
      reco::TrackRef muonTrackRef = iMuonP->track();
      if (muonTrackRef.isNull()) {
        continue;
      }
      const reco::Track* muonTrack = &(*muonTrackRef);
        if (fabs(track->px() - muonTrack->px()) < std::numeric_limits<float>::epsilon() &&
         fabs(track->py() - muonTrack->py()) < std::numeric_limits<float>::epsilon() &&
        fabs(track->pz() - muonTrack->pz()) < std::numeric_limits<float>::epsilon()) {
          isMuon = true;
          break;
      }
    }
    if (isMuon){
      continue;
    }

    // 所有筛选条件通过，按电荷分组
    if (track->charge() > 0) {
      trackPlus.push_back(iTrackc);
    } else {
      trackMinus.push_back(iTrackc);
    }
  }


  // 径迹配对要求检查：至少各有1条正负径迹
  if (trackPlus.empty() || trackMinus.empty()) {
    return;
  }

  // ========== Muon预选择与电荷分组合并优化 ==========
  // 单次循环同时完成：筛选 + 按电荷分组 + trigger filter匹配
  // 优化效果：减少一次遍历，后续配对无需再检查电荷
  std::vector<edm::View<pat::Muon>::const_iterator> muPlus, muMinus;
  std::vector<bool> muFilterMatchesPlus, muFilterMatchesMinus;

  // 预分配容量，优化性能
  muPlus.reserve(thePATMuonHandle->size() / 2 + 10);
  muMinus.reserve(thePATMuonHandle->size() / 2 + 10);
  muFilterMatchesPlus.reserve(thePATMuonHandle->size() / 2 + 10);
  muFilterMatchesMinus.reserve(thePATMuonHandle->size() / 2 + 10);

  // Selection criteria from pre-cuts.md:
  // 1. all muon candidates must be soft muons
  // 2. (pt > 3.5 && absEta < 1.2) || (pt > (5.47 - 1.89 * absEta) && 1.2 < absEta < 2.1) || (pt > 1.5 && 2.1 < absEta < 2.4)
  for (edm::View<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin();
       iMuonP != thePATMuonHandle->end(); ++iMuonP) {
    // 1. Soft Moun要求
    if (!iMuonP->isSoftMuon(thePrimaryV)) {
      continue;
    }

    // 2. Muon 动力学要求
    const double pt = iMuonP->pt();
    const double absEta = fabs(iMuonP->eta());
    bool passPtEta = false;
    if ((pt > MUON_PT_CUT_LOW && absEta < MUON_ABS_ETA_CUT1) ||
        (pt > (MUON_PT_CUT_SLOPE - MUON_PT_CUT_COEFF * absEta) && 
         absEta > MUON_ABS_ETA_CUT2_LOW && absEta < MUON_ABS_ETA_CUT2_HIGH) ||
        (pt > MUON_PT_CUT_LOW2 && 
         absEta > MUON_ABS_ETA_CUT3_LOW && absEta < MUON_ABS_ETA_CUT3_HIGH)) {
      passPtEta = true;
    }
    if (!passPtEta) {
      continue;
    }

    // 3. track有效性检查（安全校验）
    if (iMuonP->track().isNull()) {
      continue;
    }

    // 4. 触发过滤器匹配检查
    bool muonHasFilterMatch = false;
    for (unsigned int JpsiFilt = 0; JpsiFilt < FiltersForJpsi_.size(); JpsiFilt++) {
      if (hltresults.isValid()) {
        for (auto i = iMuonP->triggerObjectMatches().begin();
             i != iMuonP->triggerObjectMatches().end(); i++) {
          pat::TriggerObjectStandAlone tempTriggerObject(*i);
          tempTriggerObject.unpackFilterLabels(iEvent, *hltresults);
          if (tempTriggerObject.hasFilterLabel(FiltersForJpsi_[JpsiFilt])) {
            muonHasFilterMatch = true;
            break;
          }
        }
        if (muonHasFilterMatch) break;
      }
    }

    // 5. 按电荷分组（关键优化：为后续配对消除电荷检查）
    if (iMuonP->charge() > 0) {
      muPlus.push_back(iMuonP);
      muFilterMatchesPlus.push_back(muonHasFilterMatch);
    } else {
      muMinus.push_back(iMuonP);
      muFilterMatchesMinus.push_back(muonHasFilterMatch);
    }
  }

  // 电荷要求检查：至少2个正μ和2个负μ
  if (muPlus.size() < 2 || muMinus.size() < 2) {
    return;
  }

  // Start processing muon and track fit
// ========== 阶段1: Jpsi候选预缓存 ==========
  // 通过muon双重循环构建所有Jpsi候选并缓存，避免重复拟合
  std::vector<JpsiCandidate> jpsiCandidates;
  jpsiCandidates.reserve(muPlus.size() * muMinus.size());
  
  
  for (size_t i = 0; i < muPlus.size(); ++i) {
    const auto& muPlusIter = muPlus[i];
    TrackRef muTrack1 = muPlusIter->track();
    if (muTrack1.isNull()) {
      continue;
    }
    
    for (size_t j = 0; j < muMinus.size(); ++j) {
      const auto& muMinusIter = muMinus[j];
      TrackRef muTrack2 = muMinusIter->track();
      if (muTrack2.isNull()) {
        continue;
      }
      
      // 质量预筛选（1-4.5 GeV）
      double invMass = (muPlusIter->p4() + muMinusIter->p4()).mass();
      if (!(1.0 < invMass && invMass < 4.5)) {
        continue;
      }
      
      TransientTrack muonTT1(muTrack1, &(bFieldHandle));
      TransientTrack muonTT2(muTrack2, &(bFieldHandle));
      
      KinematicParticleFactoryFromTransientTrack pmumuFactory;
      ParticleMass muon_mass = MU_MASS;
      float muon_sigma = MU_MASSERR;
      float chi = 0.;
      float ndf = 0.;
      
      vector<RefCountedKinematicParticle> muonParticles;
      muonParticles.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi, ndf, muon_sigma));
      muonParticles.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi, ndf, muon_sigma));
      
      KinematicParticleVertexFitter fitter;
      RefCountedKinematicTree vertexFitTree;
      Error_t = false;
      try {
        vertexFitTree = fitter.fit(muonParticles);
      } catch (...) {
        Error_t = true;
      }
      
      if (Error_t || !vertexFitTree->isValid()) {
        continue;
      }
      
      vertexFitTree->movePointerToTheTop();
      RefCountedKinematicParticle jpsiParticle = vertexFitTree->currentParticle();
      RefCountedKinematicVertex jpsiVertex = vertexFitTree->currentDecayVertex();
      
      double vtxProb = ChiSquaredProbability(
          (double)(jpsiVertex->chiSquared()),
          (double)(jpsiVertex->degreesOfFreedom()));
      
      // 顶点概率筛选（>1%）
      if (vtxProb < JPSI_VTXPROB_CUT) {
        continue;
      }
      
      // 质量窗口筛选（同时检查Jpsi和Psi2S质量范围）
      double jpsiMass = jpsiParticle->currentState().mass();
      
      // 判断质量窗口
      bool isJpsiCandidate = (fabs(jpsiMass - JPSI_NOMINAL_MASS) <= JPSI_MASS_WINDOW);
      bool isPsi2SCandidate = (fabs(jpsiMass - PSI2S_NOMINAL_MASS) <= PSI2S_MASS_WINDOW);
      
      // 必须在Jpsi或Psi2S质量窗口内才保留
      if (!isJpsiCandidate && !isPsi2SCandidate) {
        continue;
      }
      
      // 获取触发匹配状态
      bool matchPlus = muFilterMatchesPlus.at(i);
      bool matchMinus = muFilterMatchesMinus.at(j);

      
      // 缓存Jpsi候选
      JpsiCandidate candidate;
      candidate.muPlus = muPlusIter;
      candidate.muMinus = muMinusIter;
      candidate.mass = jpsiMass;
      candidate.vtxProb = vtxProb;
      candidate.massErr = (jpsiParticle->currentState().kinematicParametersError().matrix()(6, 6) > 0) 
                         ? sqrt(jpsiParticle->currentState().kinematicParametersError().matrix()(6, 6)) 
                         : -9;
      candidate.p4 = ROOT::Math::PxPyPzMVector(
          muPlusIter->track()->px(), muPlusIter->track()->py(), 
          muPlusIter->track()->pz(), MU_MASS);
      candidate.p4 += ROOT::Math::PxPyPzMVector(
          muMinusIter->track()->px(), muMinusIter->track()->py(), 
          muMinusIter->track()->pz(), MU_MASS);
      candidate.kinematicParticle = jpsiParticle;
      candidate.vertex = jpsiVertex;
      candidate.muonTT1 = muonTT1;
      candidate.muonTT2 = muonTT2;
      candidate.filterMatchPlus = matchPlus;
      candidate.filterMatchMinus = matchMinus;
      candidate.isJpsiCandidate = isJpsiCandidate;
      candidate.isPsi2SCandidate = isPsi2SCandidate;
      
      // 初始化约束拟合标志为false
      candidate.hasConstraintFit = false;
      
      // Jpsi质量约束拟合
      if (isJpsiCandidate) {
        try {
          // 应用Jpsi质量约束（使用预创建的约束对象）
          KinematicParticleFitter fitterCs;
          RefCountedKinematicTree csTree = fitterCs.fit(jpsiMassConstraint_.get(), vertexFitTree);
          
          if (csTree->isValid()) {
            csTree->movePointerToTheTop();
            RefCountedKinematicParticle csParticle = csTree->currentParticle();
            
            candidate.hasConstraintFit = true;
            candidate.constraintMass = csParticle->currentState().mass();
            candidate.constraintVtxProb = ChiSquaredProbability(
                (double)(csTree->currentDecayVertex()->chiSquared()),
                (double)(csTree->currentDecayVertex()->degreesOfFreedom())
            );
            candidate.constraintMassErr = (csParticle->currentState().kinematicParametersError().matrix()(6, 6) > 0)
                ? sqrt(csParticle->currentState().kinematicParametersError().matrix()(6, 6))
                : -9;
            candidate.constraintParticle = csParticle;
            candidate.constraintVertex = csTree->currentDecayVertex();
          }
        } catch (...) {
          // 拟合异常，保持默认值
        }
      }
      
      // Psi2S质量约束拟合
      if (isPsi2SCandidate) {
        try {
          // 应用Psi2S质量约束（使用预创建的约束对象）
          KinematicParticleFitter fitterCs;
          RefCountedKinematicTree csTree = fitterCs.fit(psi2sMassConstraint_.get(), vertexFitTree);
          
          if (csTree->isValid()) {
            csTree->movePointerToTheTop();
            RefCountedKinematicParticle csParticle = csTree->currentParticle();
            
            candidate.hasConstraintFit = true;
            candidate.constraintMass = csParticle->currentState().mass();
            candidate.constraintVtxProb = ChiSquaredProbability(
                (double)(csTree->currentDecayVertex()->chiSquared()),
                (double)(csTree->currentDecayVertex()->degreesOfFreedom())
            );
            candidate.constraintMassErr = (csParticle->currentState().kinematicParametersError().matrix()(6, 6) > 0)
                ? sqrt(csParticle->currentState().kinematicParametersError().matrix()(6, 6))
                : -9;
            candidate.constraintParticle = csParticle;
            candidate.constraintVertex = csTree->currentDecayVertex();
          }
        } catch (...) {
          // 拟合异常，保持默认值
        }
      }
      
      jpsiCandidates.push_back(candidate);
    }
  }
  
  
  // 检查是否有足够的Jpsi候选
  if (jpsiCandidates.size() < 2) {
    return;
  }

  // ========== 阶段2: Jpsi候选两两组合 ==========
  // 从缓存中选择两个不同的Jpsi进行组合
  for (size_t i = 0; i < jpsiCandidates.size(); ++i) {
    const JpsiCandidate& jpsi1 = jpsiCandidates[i];
    
    for (size_t j = 0; j < jpsiCandidates.size(); ++j) {
      const JpsiCandidate& jpsi2 = jpsiCandidates[j];
      // 检查是否是同一个
      if (i == j) {
        continue;
      }
      // 检查jpsi1是否有Jpsi约束拟合结果
      if (!jpsi1.isJpsiCandidate || !jpsi1.hasConstraintFit) {
        continue;
      }
      // 检查是否使用了相同的muon
      if (jpsi1.muPlus == jpsi2.muPlus || jpsi1.muPlus == jpsi2.muMinus ||
          jpsi1.muMinus == jpsi2.muPlus || jpsi1.muMinus == jpsi2.muMinus) {
        continue;
      }
      
      // 触发匹配检查：jpsi1的两个muon匹配 或 jpsi2的两个muon匹配
      bool passTrigger = (jpsi1.filterMatchPlus && jpsi1.filterMatchMinus) ||
                         (jpsi2.filterMatchPlus && jpsi2.filterMatchMinus);
      if (!passTrigger) {
        continue;
      }

      // ========== 阶段3: 与Track组合 ==========
      for (const auto& trackPlusIter : trackPlus) {
        if (!trackPlusIter->hasTrackDetails() || trackPlusIter->charge() == 0) {
          continue;
        }
        const reco::Track* track1 = trackPlusIter->bestTrack();
        if (!track1) {
          continue;
        }

        for (const auto& trackMinusIter : trackMinus) {
          if (!trackMinusIter->hasTrackDetails() || trackMinusIter->charge() == 0) {
            continue;
          }
          const reco::Track* track2 = trackMinusIter->bestTrack();
          if (!track2) {
            continue;
          }

          ROOT::Math::PxPyPzMVector P4_Track1(trackPlusIter->px(), trackPlusIter->py(), 
                                               trackPlusIter->pz(), PI_MASS);
          ROOT::Math::PxPyPzMVector P4_Track2(trackMinusIter->px(), trackMinusIter->py(), 
                                               trackMinusIter->pz(), PI_MASS);
          ROOT::Math::PxPyPzMVector P4_Jpsipipi = jpsi1.p4 + P4_Track1 + P4_Track2;

          if (ROOT::Math::VectorUtil::DeltaR(P4_Track1, P4_Jpsipipi) > PION_DR_CUT) {
            continue;
          }
          if (ROOT::Math::VectorUtil::DeltaR(P4_Track2, P4_Jpsipipi) > PION_DR_CUT) {
            continue;
          }

          TransientTrack trackTT1(*track1, &(bFieldHandle));
          TransientTrack trackTT2(*track2, &(bFieldHandle));
          KinematicParticleFactoryFromTransientTrack JPiPiFactory;
          ParticleMass pion_mass = PI_MASS;
          float pion_sigma = PI_MASSERR;
          float chi = 0.;
          float ndf = 0.;

          // ========== 四粒子拟合（使用约束后的Jpsi候选） ==========
          
          vector<RefCountedKinematicParticle> JPiPiParticles;
          JPiPiParticles.push_back(JPiPiFactory.particle(
              trackTT1, pion_mass, chi, ndf, pion_sigma));
          JPiPiParticles.push_back(JPiPiFactory.particle(
              trackTT2, pion_mass, chi, ndf, pion_sigma));
          // 使用约束后的Jpsi候选作为输入粒子
          JPiPiParticles.push_back(jpsi1.constraintParticle);

          // 使用普通顶点拟合（KinematicParticleVertexFitter）
          KinematicParticleVertexFitter JPiPi_fitter;
          RefCountedKinematicTree JPiPiVertexFitTree;
          Error_t = false;
          try {
            JPiPiVertexFitTree = JPiPi_fitter.fit(JPiPiParticles);
          } catch (...) {
            Error_t = true;
          }
          if (Error_t || !(JPiPiVertexFitTree->isValid())) {
            continue;
          }
          
          JPiPiVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle JPiPi_vFit_constrained =
              JPiPiVertexFitTree->currentParticle();
          RefCountedKinematicVertex JPiPi_vFit_vertex_constrained =
              JPiPiVertexFitTree->currentDecayVertex();

          double JPiPi_vtxprob = ChiSquaredProbability(
              (double)(JPiPi_vFit_vertex_constrained->chiSquared()),
              (double)(JPiPi_vFit_vertex_constrained->degreesOfFreedom()));
          
          if (JPiPi_vFit_constrained->currentState().mass() > 4.5) {
            continue;
          }
          
          if (JPiPi_vtxprob < PSI2S_VTXPROB_CUT) {
            continue;
          }
          
          double px = JPiPi_vFit_constrained->currentState().kinematicParameters().momentum().x();
          double py = JPiPi_vFit_constrained->currentState().kinematicParameters().momentum().y();
          double psi2s_pt = sqrt(px*px + py*py);
          if (psi2s_pt <= PSI2S_PT_CUT) {
            continue;
          }

          // ========== 保存Psi2S拟合结果 ==========
          Psi2S_mass = JPiPi_vFit_constrained->currentState().mass();
          Psi2S_VtxProb = JPiPi_vtxprob;
          Psi2S_px = JPiPi_vFit_constrained->currentState()
                                     .kinematicParameters()
                                     .momentum()
                                     .x();
          Psi2S_py = JPiPi_vFit_constrained->currentState()
                                     .kinematicParameters()
                                     .momentum()
                                     .y();
          Psi2S_pz = JPiPi_vFit_constrained->currentState()
                                     .kinematicParameters()
                                     .momentum()
                                     .z();
          if (JPiPi_vFit_constrained->currentState()
                  .kinematicParametersError()
                  .matrix()(6, 6) > 0) {
            Psi2S_massErr = sqrt(JPiPi_vFit_constrained->currentState()
                                        .kinematicParametersError()
                                        .matrix()(6, 6));
          } else {
            Psi2S_massErr = -9;
          }

          ROOT::Math::PxPyPzMVector Psi2S_vec(Psi2S_px, Psi2S_py, Psi2S_pz, Psi2S_mass);
          Psi2S_pt = Psi2S_vec.Pt();
          Psi2S_absEta = fabs(Psi2S_vec.Eta());

          // ========== 六粒子三种假设拟合 ==========
          // ========== 假设PJ: JPiPi(Psi2S约束后) + Jpsi2(Jpsi约束后) ==========
          // 检查Jpsi2是否拥有Jpsi约束拟合
          if (jpsi2.hasConstraintFit && jpsi2.isJpsiCandidate) {
            // 第一步：对四粒子进行Psi2S质量约束拟合
            RefCountedKinematicTree JPiPi_cs_Psi2S_tree;
            bool JPiPi_Psi2S_Error = false;
            
            try {
              // 使用预创建的Psi2S质量约束
              KinematicParticleFitter jpipiFitter;
              JPiPi_cs_Psi2S_tree = jpipiFitter.fit(psi2sMassConstraint_.get(), JPiPiVertexFitTree);
            } catch (...) {
              JPiPi_Psi2S_Error = true;
            }
            
            if (JPiPi_Psi2S_Error || !JPiPi_cs_Psi2S_tree->isValid()) {
              // 拟合失败，保持默认值
            } else {
              JPiPi_cs_Psi2S_tree->movePointerToTheTop();
              RefCountedKinematicParticle JPiPi_vFit_cs_Psi2S = 
                  JPiPi_cs_Psi2S_tree->currentParticle();
              
              // 第二步：六粒子顶点拟合（JPiPi + Jpsi2）
              vector<RefCountedKinematicParticle> PJ_Particles;
              PJ_Particles.push_back(JPiPi_vFit_cs_Psi2S);
              PJ_Particles.push_back(jpsi2.constraintParticle);
              
              KinematicParticleVertexFitter PJ_fitter;
              RefCountedKinematicTree PJ_VertexFitTree;
              bool PJ_Error = false;
              try {
                PJ_VertexFitTree = PJ_fitter.fit(PJ_Particles);
              } catch (...) {
                PJ_Error = true;
              }
              
              if (PJ_Error || !(PJ_VertexFitTree->isValid())) {
                // 拟合失败，保持默认值
              } else {
                PJ_VertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle PJ_vFit = PJ_VertexFitTree->currentParticle();
                RefCountedKinematicVertex PJ_vFit_vertex = PJ_VertexFitTree->currentDecayVertex();
                KinematicParameters PJ_kPara = PJ_vFit->currentState().kinematicParameters();

                X_PJ_mass = PJ_vFit->currentState().mass();
                X_PJ_VtxProb = ChiSquaredProbability(
                    (double)(PJ_vFit_vertex->chiSquared()),
                    (double)(PJ_vFit_vertex->degreesOfFreedom()));
                X_PJ_px = PJ_kPara.momentum().x();
                X_PJ_py = PJ_kPara.momentum().y();
                X_PJ_pz = PJ_kPara.momentum().z();
                if (PJ_vFit->currentState().kinematicParametersError().matrix()(6, 6) > 0) {
                  X_PJ_massErr = sqrt(PJ_vFit->currentState().kinematicParametersError().matrix()(6, 6));
                } else {
                  X_PJ_massErr = -9;
                }
                ROOT::Math::PxPyPzMVector PJ_vec(X_PJ_px, X_PJ_py, X_PJ_pz, X_PJ_mass);
                X_PJ_pt = PJ_vec.Pt();
                X_PJ_absEta = fabs(PJ_vec.Eta());
              }
            }
          }

          // ========== 假设XJ: JPiPi(X3872约束后) + Jpsi2(Jpsi约束后) ==========
          // 检查Jpsi2是否拥有Jpsi约束拟合
          if (jpsi2.hasConstraintFit && jpsi2.isJpsiCandidate){
            // 第一步：对四粒子进行X3872质量约束拟合
            RefCountedKinematicTree JPiPi_cs_X3872_tree;
            bool JPiPi_X3872_Error = false;
            
            try {
              // 使用预创建的X3872质量约束
              KinematicParticleFitter jpipiFitter;
              JPiPi_cs_X3872_tree = jpipiFitter.fit(x3872MassConstraint_.get(), JPiPiVertexFitTree);
            } catch (...) {
              JPiPi_X3872_Error = true;
            }
            
            if (JPiPi_X3872_Error || !JPiPi_cs_X3872_tree->isValid()) {
              // 拟合失败，保持默认值
            } else {
              JPiPi_cs_X3872_tree->movePointerToTheTop();
              RefCountedKinematicParticle JPiPi_vFit_cs_X3872 = 
                  JPiPi_cs_X3872_tree->currentParticle();
              
              // 第二步：六粒子顶点拟合（JPiPi + Jpsi2）
              vector<RefCountedKinematicParticle> XJ_Particles;
              XJ_Particles.push_back(JPiPi_vFit_cs_X3872);
              XJ_Particles.push_back(jpsi2.constraintParticle);
              
              KinematicParticleVertexFitter XJ_fitter;
              RefCountedKinematicTree XJ_VertexFitTree;
              bool XJ_Error = false;
              try {
                XJ_VertexFitTree = XJ_fitter.fit(XJ_Particles);
              } catch (...) {
                XJ_Error = true;
              }
              
              if (XJ_Error || !(XJ_VertexFitTree->isValid())) {
                // 拟合失败，保持默认值
              } else {
                XJ_VertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle XJ_vFit = XJ_VertexFitTree->currentParticle();
                RefCountedKinematicVertex XJ_vFit_vertex = XJ_VertexFitTree->currentDecayVertex();
                KinematicParameters XJ_kPara = XJ_vFit->currentState().kinematicParameters();

                X_XJ_mass = XJ_vFit->currentState().mass();
                X_XJ_VtxProb = ChiSquaredProbability(
                    (double)(XJ_vFit_vertex->chiSquared()),
                    (double)(XJ_vFit_vertex->degreesOfFreedom()));
                X_XJ_px = XJ_kPara.momentum().x();
                X_XJ_py = XJ_kPara.momentum().y();
                X_XJ_pz = XJ_kPara.momentum().z();
                if (XJ_vFit->currentState().kinematicParametersError().matrix()(6, 6) > 0) {
                  X_XJ_massErr = sqrt(XJ_vFit->currentState().kinematicParametersError().matrix()(6, 6));
                } else {
                  X_XJ_massErr = -9;
                }
                ROOT::Math::PxPyPzMVector XJ_vec(X_XJ_px, X_XJ_py, X_XJ_pz, X_XJ_mass);
                X_XJ_pt = XJ_vec.Pt();
                X_XJ_absEta = fabs(XJ_vec.Eta());
              }
            }
          }

          // ========== 假设PP: JPiPi(Psi2S约束后) + Jpsi2(Psi2S约束后) ==========
          // 检查Jpsi2是否拥有Psi2S约束拟合
          if (jpsi2.hasConstraintFit && jpsi2.isPsi2SCandidate){
            // 第一步：对四粒子进行Psi2S质量约束拟合
            RefCountedKinematicTree JPiPi_cs_Psi2S_tree;
            bool JPiPi_Psi2S_Error = false;
            
            try {
              // 使用预创建的Psi2S质量约束
              KinematicParticleFitter jpipiFitter;
              JPiPi_cs_Psi2S_tree = jpipiFitter.fit(psi2sMassConstraint_.get(), JPiPiVertexFitTree);
            } catch (...) {
              JPiPi_Psi2S_Error = true;
            }
            
            if (JPiPi_Psi2S_Error || !JPiPi_cs_Psi2S_tree->isValid()) {
              // 拟合失败，保持默认值
            } else {
              JPiPi_cs_Psi2S_tree->movePointerToTheTop();
              RefCountedKinematicParticle JPiPi_vFit_cs_Psi2S = 
                  JPiPi_cs_Psi2S_tree->currentParticle();
              
              // 第二步：六粒子顶点拟合（JPiPi + Jpsi2）
              vector<RefCountedKinematicParticle> PP_Particles;
              PP_Particles.push_back(JPiPi_vFit_cs_Psi2S);
              PP_Particles.push_back(jpsi2.constraintParticle);
              
              KinematicParticleVertexFitter PP_fitter;
              RefCountedKinematicTree PP_VertexFitTree;
              bool PP_Error = false;
              try {
                PP_VertexFitTree = PP_fitter.fit(PP_Particles);
              } catch (...) {
                PP_Error = true;
              }
              
              if (PP_Error || !(PP_VertexFitTree->isValid())) {
                // 拟合失败，保持默认值
              } else {
                PP_VertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle PP_vFit = PP_VertexFitTree->currentParticle();
                RefCountedKinematicVertex PP_vFit_vertex = PP_VertexFitTree->currentDecayVertex();
                KinematicParameters PP_kPara = PP_vFit->currentState().kinematicParameters();

                X_PP_mass = PP_vFit->currentState().mass();
                X_PP_VtxProb = ChiSquaredProbability(
                    (double)(PP_vFit_vertex->chiSquared()),
                    (double)(PP_vFit_vertex->degreesOfFreedom()));
                X_PP_px = PP_kPara.momentum().x();
                X_PP_py = PP_kPara.momentum().y();
                X_PP_pz = PP_kPara.momentum().z();
                if (PP_vFit->currentState().kinematicParametersError().matrix()(6, 6) > 0) {
                  X_PP_massErr = sqrt(PP_vFit->currentState().kinematicParametersError().matrix()(6, 6));
                } else {
                  X_PP_massErr = -9;
                }
                ROOT::Math::PxPyPzMVector PP_vec(X_PP_px, X_PP_py, X_PP_pz, X_PP_mass);
                X_PP_pt = PP_vec.Pt();
                X_PP_absEta = fabs(PP_vec.Eta());
              }
            }
          }

          // ========== 保存Jpsi变量 ==========
          Jpsi1_mass = jpsi1.mass;
          Jpsi1_VtxProb = jpsi1.vtxProb;
          Jpsi1_massErr = jpsi1.massErr;
          ROOT::Math::PxPyPzMVector Jpsi1_vec(jpsi1.p4.Px(), jpsi1.p4.Py(), jpsi1.p4.Pz(), jpsi1.mass);
          Jpsi1_pt = Jpsi1_vec.Pt();
          Jpsi1_absEta = fabs(Jpsi1_vec.Eta());

          Jpsi2_mass = jpsi2.mass;
          Jpsi2_VtxProb = jpsi2.vtxProb;
          Jpsi2_massErr = jpsi2.massErr;
          ROOT::Math::PxPyPzMVector Jpsi2_vec(jpsi2.p4.Px(), jpsi2.p4.Py(), jpsi2.p4.Pz(), jpsi2.mass);
          Jpsi2_pt = Jpsi2_vec.Pt();
          Jpsi2_absEta = fabs(Jpsi2_vec.Eta());

          // ========== 保存muon变量 ==========
          const auto& iMuon1 = *jpsi1.muPlus;
          const auto& iMuon2 = *jpsi1.muMinus;
          const auto& iMuon3 = *jpsi2.muPlus;
          const auto& iMuon4 = *jpsi2.muMinus;

          mu1_hasFilterMatch = jpsi1.filterMatchPlus ? 1 : 0;
          mu2_hasFilterMatch = jpsi1.filterMatchMinus ? 1 : 0;
          mu3_hasFilterMatch = jpsi2.filterMatchPlus ? 1 : 0;
          mu4_hasFilterMatch = jpsi2.filterMatchMinus ? 1 : 0;

          mu1_px = iMuon1.px();
          mu1_py = iMuon1.py();
          mu1_pz = iMuon1.pz();
          mu1_pt = iMuon1.pt();
          mu1_absEta = fabs(iMuon1.eta());
          mu1_trackIso = iMuon1.trackIso();
          mu1_d0BS = iMuon1.dB(pat::Muon::BS2D);
          mu1_d0EBS = iMuon1.edB(pat::Muon::BS2D);
          mu1_d3dBS = iMuon1.dB(pat::Muon::BS3D);
          mu1_d3dEBS = iMuon1.edB(pat::Muon::BS3D);
          mu1_d0PV = iMuon1.dB(pat::Muon::PV2D);
          mu1_d0EPV = iMuon1.edB(pat::Muon::PV2D);
          mu1_dzPV = iMuon1.dB(pat::Muon::PVDZ);
          mu1_dzEPV = iMuon1.edB(pat::Muon::PVDZ);
          mu1_charge = iMuon1.charge();

          mu2_px = iMuon2.px();
          mu2_py = iMuon2.py();
          mu2_pz = iMuon2.pz();
          mu2_pt = iMuon2.pt();
          mu2_absEta = fabs(iMuon2.eta());
          mu2_trackIso = iMuon2.trackIso();
          mu2_d0BS = iMuon2.dB(pat::Muon::BS2D);
          mu2_d0EBS = iMuon2.edB(pat::Muon::BS2D);
          mu2_d3dBS = iMuon2.dB(pat::Muon::BS3D);
          mu2_d3dEBS = iMuon2.edB(pat::Muon::BS3D);
          mu2_d0PV = iMuon2.dB(pat::Muon::PV2D);
          mu2_d0EPV = iMuon2.edB(pat::Muon::PV2D);
          mu2_dzPV = iMuon2.dB(pat::Muon::PVDZ);
          mu2_dzEPV = iMuon2.edB(pat::Muon::PVDZ);
          mu2_charge = iMuon2.charge();

          mu3_px = iMuon3.px();
          mu3_py = iMuon3.py();
          mu3_pz = iMuon3.pz();
          mu3_pt = iMuon3.pt();
          mu3_absEta = fabs(iMuon3.eta());
          mu3_trackIso = iMuon3.trackIso();
          mu3_d0BS = iMuon3.dB(pat::Muon::BS2D);
          mu3_d0EBS = iMuon3.edB(pat::Muon::BS2D);
          mu3_d3dBS = iMuon3.dB(pat::Muon::BS3D);
          mu3_d3dEBS = iMuon3.edB(pat::Muon::BS3D);
          mu3_d0PV = iMuon3.dB(pat::Muon::PV2D);
          mu3_d0EPV = iMuon3.edB(pat::Muon::PV2D);
          mu3_dzPV = iMuon3.dB(pat::Muon::PVDZ);
          mu3_dzEPV = iMuon3.edB(pat::Muon::PVDZ);
          mu3_charge = iMuon3.charge();

          mu4_px = iMuon4.px();
          mu4_py = iMuon4.py();
          mu4_pz = iMuon4.pz();
          mu4_pt = iMuon4.pt();
          mu4_absEta = fabs(iMuon4.eta());
          mu4_trackIso = iMuon4.trackIso();
          mu4_d0BS = iMuon4.dB(pat::Muon::BS2D);
          mu4_d0EBS = iMuon4.edB(pat::Muon::BS2D);
          mu4_d3dBS = iMuon4.dB(pat::Muon::BS3D);
          mu4_d3dEBS = iMuon4.edB(pat::Muon::BS3D);
          mu4_d0PV = iMuon4.dB(pat::Muon::PV2D);
          mu4_d0EPV = iMuon4.edB(pat::Muon::PV2D);
          mu4_dzPV = iMuon4.dB(pat::Muon::PVDZ);
          mu4_dzEPV = iMuon4.edB(pat::Muon::PVDZ);
          mu4_charge = iMuon4.charge();

          // 计算Muon ID
          nLooseMuons = 0;
          nTightMuons = 0;
          nSoftMuons = 0;
          nMediumMuons = 0;

          if (iMuon1.isLooseMuon()) nLooseMuons++;
          if (iMuon1.isTightMuon(thePrimaryV)) nTightMuons++;
          if (iMuon1.isSoftMuon(thePrimaryV)) nSoftMuons++;
          if (iMuon1.isMediumMuon()) nMediumMuons++;

          if (iMuon2.isLooseMuon()) nLooseMuons++;
          if (iMuon2.isTightMuon(thePrimaryV)) nTightMuons++;
          if (iMuon2.isSoftMuon(thePrimaryV)) nSoftMuons++;
          if (iMuon2.isMediumMuon()) nMediumMuons++;

          if (iMuon3.isLooseMuon()) nLooseMuons++;
          if (iMuon3.isTightMuon(thePrimaryV)) nTightMuons++;
          if (iMuon3.isSoftMuon(thePrimaryV)) nSoftMuons++;
          if (iMuon3.isMediumMuon()) nMediumMuons++;

          if (iMuon4.isLooseMuon()) nLooseMuons++;
          if (iMuon4.isTightMuon(thePrimaryV)) nTightMuons++;
          if (iMuon4.isSoftMuon(thePrimaryV)) nSoftMuons++;
          if (iMuon4.isMediumMuon()) nMediumMuons++;

          // ========== 保存pion变量 ==========
          pi1_px = trackTT1.impactPointState().globalMomentum().x();
          pi1_py = trackTT1.impactPointState().globalMomentum().y();
          pi1_pz = trackTT1.impactPointState().globalMomentum().z();
          ROOT::Math::PxPyPzMVector pi1_vec(pi1_px, pi1_py, pi1_pz, PI_MASS);
          pi1_pt = pi1_vec.Pt();
          pi1_absEta = fabs(pi1_vec.Eta());

          pi2_px = trackTT2.impactPointState().globalMomentum().x();
          pi2_py = trackTT2.impactPointState().globalMomentum().y();
          pi2_pz = trackTT2.impactPointState().globalMomentum().z();
          ROOT::Math::PxPyPzMVector pi2_vec(pi2_px, pi2_py, pi2_pz, PI_MASS);
          pi2_pt = pi2_vec.Pt();
          pi2_absEta = fabs(pi2_vec.Eta());

          // ========== 计算DeltaR变量 ==========
          ROOT::Math::PxPyPzMVector mu1_vec(mu1_px, mu1_py, mu1_pz, MU_MASS);
          ROOT::Math::PxPyPzMVector mu2_vec(mu2_px, mu2_py, mu2_pz, MU_MASS);
          ROOT::Math::PxPyPzMVector mu3_vec(mu3_px, mu3_py, mu3_pz, MU_MASS);
          ROOT::Math::PxPyPzMVector mu4_vec(mu4_px, mu4_py, mu4_pz, MU_MASS);

          dR_mu1_mu2 = ROOT::Math::VectorUtil::DeltaR(mu1_vec, mu2_vec);
          dR_mu3_mu4 = ROOT::Math::VectorUtil::DeltaR(mu3_vec, mu4_vec);
          dR_pi1_pi2 = ROOT::Math::VectorUtil::DeltaR(pi1_vec, pi2_vec);
          dR_Psi2S_Jpsi1 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, Jpsi1_vec);
          dR_Psi2S_Jpsi2 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, Jpsi2_vec);
          dR_Psi2S_pi1 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, pi1_vec);
          dR_Psi2S_pi2 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, pi2_vec);

          X_One_Tree_->Fill();
        } // 结束 trackMinus 循环
      } // 结束 trackPlus 循环
    } // 结束 jpsi2 循环
  } // 结束 jpsi1 循环

} // analyze
//
// ------------ method called once each job just before starting event loop
// ------------
void MultiLepPAT::beginRun(edm::Run const &iRun,
                           edm::EventSetup const &iSetup) {}

void MultiLepPAT::beginJob() {
  edm::Service<TFileService> fs;
  X_One_Tree_ = fs->make<TTree>("X_data", "X Data");


  X_One_Tree_->Branch("evtNum", &evtNum, "evtNum/i");
  X_One_Tree_->Branch("runNum", &runNum, "runNum/i");
  X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
  X_One_Tree_->Branch("nGoodPrimVtx", &nGoodPrimVtx, "nGoodPrimVtx/i");


  // 质量拟合假设 PJ: mumupipi (psi2S 约束) + mumu (J/psi 约束)
  X_One_Tree_->Branch("X_PJ_mass", &X_PJ_mass, "X_PJ_mass/F");
  X_One_Tree_->Branch("X_PJ_VtxProb", &X_PJ_VtxProb, "X_PJ_VtxProb/F");
  X_One_Tree_->Branch("X_PJ_massErr", &X_PJ_massErr, "X_PJ_massErr/F");
  X_One_Tree_->Branch("X_PJ_pt", &X_PJ_pt, "X_PJ_pt/F");
  X_One_Tree_->Branch("X_PJ_pz", &X_PJ_pz, "X_PJ_pz/F");
  X_One_Tree_->Branch("X_PJ_absEta", &X_PJ_absEta, "X_PJ_absEta/F");
  X_One_Tree_->Branch("X_PJ_px", &X_PJ_px, "X_PJ_px/F");
  X_One_Tree_->Branch("X_PJ_py", &X_PJ_py, "X_PJ_py/F");

  // 假设 XJ: mumupipi (X(3872) 约束) + mumu (J/psi 约束)
  X_One_Tree_->Branch("X_XJ_mass", &X_XJ_mass, "X_XJ_mass/F");
  X_One_Tree_->Branch("X_XJ_VtxProb", &X_XJ_VtxProb, "X_XJ_VtxProb/F");
  X_One_Tree_->Branch("X_XJ_massErr", &X_XJ_massErr, "X_XJ_massErr/F");
  X_One_Tree_->Branch("X_XJ_pt", &X_XJ_pt, "X_XJ_pt/F");
  X_One_Tree_->Branch("X_XJ_pz", &X_XJ_pz, "X_XJ_pz/F");
  X_One_Tree_->Branch("X_XJ_absEta", &X_XJ_absEta, "X_XJ_absEta/F");
  X_One_Tree_->Branch("X_XJ_px", &X_XJ_px, "X_XJ_px/F");
  X_One_Tree_->Branch("X_XJ_py", &X_XJ_py, "X_XJ_py/F");

  // 假设 PP: mumupipi (psi(2S) 约束) + mumu (psi(2S) 约束)
  X_One_Tree_->Branch("X_PP_mass", &X_PP_mass, "X_PP_mass/F");
  X_One_Tree_->Branch("X_PP_VtxProb", &X_PP_VtxProb, "X_PP_VtxProb/F");
  X_One_Tree_->Branch("X_PP_massErr", &X_PP_massErr, "X_PP_massErr/F");
  X_One_Tree_->Branch("X_PP_pt", &X_PP_pt, "X_PP_pt/F");
  X_One_Tree_->Branch("X_PP_pz", &X_PP_pz, "X_PP_pz/F");
  X_One_Tree_->Branch("X_PP_absEta", &X_PP_absEta, "X_PP_absEta/F");
  X_One_Tree_->Branch("X_PP_px", &X_PP_px, "X_PP_px/F");
  X_One_Tree_->Branch("X_PP_py", &X_PP_py, "X_PP_py/F");

  // Psi(2S) 带 J/psi 约束: mumupipi system
  X_One_Tree_->Branch("Psi2S_mass", &Psi2S_mass, "Psi2S_mass/F");
  X_One_Tree_->Branch("Psi2S_VtxProb", &Psi2S_VtxProb, "Psi2S_VtxProb/F");
  X_One_Tree_->Branch("Psi2S_massErr", &Psi2S_massErr, "Psi2S_massErr/F");
  X_One_Tree_->Branch("Psi2S_pt", &Psi2S_pt, "Psi2S_pt/F");
  X_One_Tree_->Branch("Psi2S_pz", &Psi2S_pz, "Psi2S_pz/F");
  X_One_Tree_->Branch("Psi2S_absEta", &Psi2S_absEta, "Psi2S_absEta/F");
  X_One_Tree_->Branch("Psi2S_px", &Psi2S_px, "Psi2S_px/F");
  X_One_Tree_->Branch("Psi2S_py", &Psi2S_py, "Psi2S_py/F");

  X_One_Tree_->Branch("Jpsi1_mass", &Jpsi1_mass, "Jpsi1_mass/F");
  X_One_Tree_->Branch("Jpsi1_VtxProb", &Jpsi1_VtxProb, "Jpsi1_VtxProb/F");
  X_One_Tree_->Branch("Jpsi1_massErr", &Jpsi1_massErr, "Jpsi1_massErr/F");
  X_One_Tree_->Branch("Jpsi1_pt", &Jpsi1_pt, "Jpsi1_pt/F");
  X_One_Tree_->Branch("Jpsi1_pz", &Jpsi1_pz, "Jpsi1_pz/F");
  X_One_Tree_->Branch("Jpsi1_absEta", &Jpsi1_absEta, "Jpsi1_absEta/F");

  X_One_Tree_->Branch("Jpsi2_mass", &Jpsi2_mass, "Jpsi2_mass/F");
  X_One_Tree_->Branch("Jpsi2_VtxProb", &Jpsi2_VtxProb, "Jpsi2_VtxProb/F");
  X_One_Tree_->Branch("Jpsi2_massErr", &Jpsi2_massErr, "Jpsi2_massErr/F");
  X_One_Tree_->Branch("Jpsi2_pt", &Jpsi2_pt, "Jpsi2_pt/F");
  X_One_Tree_->Branch("Jpsi2_pz", &Jpsi2_pz, "Jpsi2_pz/F");
  X_One_Tree_->Branch("Jpsi2_absEta", &Jpsi2_absEta, "Jpsi2_absEta/F");

  X_One_Tree_->Branch("mu1_pt", &mu1_pt, "mu1_pt/F");
  X_One_Tree_->Branch("mu1_pz", &mu1_pz, "mu1_pz/F");
  X_One_Tree_->Branch("mu1_absEta", &mu1_absEta, "mu1_absEta/F");
  X_One_Tree_->Branch("mu1_trackIso", &mu1_trackIso, "mu1_trackIso/F");
  X_One_Tree_->Branch("mu1_d0BS", &mu1_d0BS, "mu1_d0BS/F");
  X_One_Tree_->Branch("mu1_d0EPV", &mu1_d0EPV, "mu1_d0EPV/F");
  X_One_Tree_->Branch("mu1_dzPV", &mu1_dzPV, "mu1_dzPV/F");
  X_One_Tree_->Branch("mu1_dzEPV", &mu1_dzEPV, "mu1_dzEPV/F");
  X_One_Tree_->Branch("mu1_charge", &mu1_charge, "mu1_charge/F");

  X_One_Tree_->Branch("mu2_pt", &mu2_pt, "mu2_pt/F");
  X_One_Tree_->Branch("mu2_pz", &mu2_pz, "mu2_pz/F");
  X_One_Tree_->Branch("mu2_absEta", &mu2_absEta, "mu2_absEta/F");
  X_One_Tree_->Branch("mu2_trackIso", &mu2_trackIso, "mu2_trackIso/F");
  X_One_Tree_->Branch("mu2_d0BS", &mu2_d0BS, "mu2_d0BS/F");
  X_One_Tree_->Branch("mu2_d0EPV", &mu2_d0EPV, "mu2_d0EPV/F");
  X_One_Tree_->Branch("mu2_dzPV", &mu2_dzPV, "mu2_dzPV/F");
  X_One_Tree_->Branch("mu2_dzEPV", &mu2_dzEPV, "mu2_dzEPV/F");
  X_One_Tree_->Branch("mu2_charge", &mu2_charge, "mu2_charge/F");

  X_One_Tree_->Branch("mu3_pt", &mu3_pt, "mu3_pt/F");
  X_One_Tree_->Branch("mu3_pz", &mu3_pz, "mu3_pz/F");
  X_One_Tree_->Branch("mu3_absEta", &mu3_absEta, "mu3_absEta/F");
  X_One_Tree_->Branch("mu3_trackIso", &mu3_trackIso, "mu3_trackIso/F");
  X_One_Tree_->Branch("mu3_d0BS", &mu3_d0BS, "mu3_d0BS/F");
  X_One_Tree_->Branch("mu3_d0EPV", &mu3_d0EPV, "mu3_d0EPV/F");
  X_One_Tree_->Branch("mu3_dzPV", &mu3_dzPV, "mu3_dzPV/F");
  X_One_Tree_->Branch("mu3_dzEPV", &mu3_dzEPV, "mu3_dzEPV/F");
  X_One_Tree_->Branch("mu3_charge", &mu3_charge, "mu3_charge/F");

  X_One_Tree_->Branch("mu4_pt", &mu4_pt, "mu4_pt/F");
  X_One_Tree_->Branch("mu4_pz", &mu4_pz, "mu4_pz/F");
  X_One_Tree_->Branch("mu4_absEta", &mu4_absEta, "mu4_absEta/F");
  X_One_Tree_->Branch("mu4_trackIso", &mu4_trackIso, "mu4_trackIso/F");
  X_One_Tree_->Branch("mu4_d0BS", &mu4_d0BS, "mu4_d0BS/F");
  X_One_Tree_->Branch("mu4_d0EPV", &mu4_d0EPV, "mu4_d0EPV/F");
  X_One_Tree_->Branch("mu4_dzPV", &mu4_dzPV, "mu4_dzPV/F");
  X_One_Tree_->Branch("mu4_dzEPV", &mu4_dzEPV, "mu4_dzEPV/F");
  X_One_Tree_->Branch("mu4_charge", &mu4_charge, "mu4_charge/F");

  X_One_Tree_->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
  X_One_Tree_->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
  X_One_Tree_->Branch("nSoftMuons", &nSoftMuons, "nSoftMuons/I");
  X_One_Tree_->Branch("nMediumMuons", &nMediumMuons, "nMediumMuons/I");

  X_One_Tree_->Branch("pi1_pt", &pi1_pt, "pi1_pt/F");
  X_One_Tree_->Branch("pi1_pz", &pi1_pz, "pi1_pz/F");
  X_One_Tree_->Branch("pi1_absEta", &pi1_absEta, "pi1_absEta/F");

  X_One_Tree_->Branch("pi2_pt", &pi2_pt, "pi2_pt/F");
  X_One_Tree_->Branch("pi2_pz", &pi2_pz, "pi2_pz/F");
  X_One_Tree_->Branch("pi2_absEta", &pi2_absEta, "pi2_absEta/F");

  X_One_Tree_->Branch("dR_mu1_mu2", &dR_mu1_mu2, "dR_mu1_mu2/F");
  X_One_Tree_->Branch("dR_mu3_mu4", &dR_mu3_mu4, "dR_mu3_mu4/F");
  X_One_Tree_->Branch("dR_pi1_pi2", &dR_pi1_pi2, "dR_pi1_pi2/F");
  X_One_Tree_->Branch("dR_Psi2S_Jpsi1", &dR_Psi2S_Jpsi1, "dR_Psi2S_Jpsi1/F");
  X_One_Tree_->Branch("dR_Psi2S_Jpsi2", &dR_Psi2S_Jpsi2, "dR_Psi2S_Jpsi2/F");
  X_One_Tree_->Branch("dR_Psi2S_pi1", &dR_Psi2S_pi1, "dR_Psi2S_pi1/F");
  X_One_Tree_->Branch("dR_Psi2S_pi2", &dR_Psi2S_pi2, "dR_Psi2S_pi2/F");

} // begin Job

// ------------ method called to reset all variables for each event ------------
void MultiLepPAT::resetVariables() {
  
  // ============================================================
  //                    重置事件标识变量
  // ============================================================
  runNum = 0;
  evtNum = 0;
  lumiNum = 0;
  nGoodPrimVtx = 0;
  
  // ============================================================
  //                    重置 X 候选变量
  // ============================================================
  // 假设 PJ: mumupipi (ψ(2S) 质量约束) + mumu (J/ψ 质量约束)
  X_PJ_mass = -999.0;
  X_PJ_VtxProb = -999.0;
  X_PJ_massErr = -999.0;
  X_PJ_pt = -999.0;
  X_PJ_pz = -999.0;
  X_PJ_absEta = -999.0;
  X_PJ_px = -999.0;
  X_PJ_py = -999.0;

  // 假设 XJ: mumupipi (X(3872) 质量约束) + mumu (J/ψ 质量约束)
  X_XJ_mass = -999.0;
  X_XJ_VtxProb = -999.0;
  X_XJ_massErr = -999.0;
  X_XJ_pt = -999.0;
  X_XJ_pz = -999.0;
  X_XJ_absEta = -999.0;
  X_XJ_px = -999.0;
  X_XJ_py = -999.0;

  // 假设 PP: mumupipi (ψ(2S) 质量约束) + mumu (ψ(2S) 质量约束)
  X_PP_mass = -999.0;
  X_PP_VtxProb = -999.0;
  X_PP_massErr = -999.0;
  X_PP_pt = -999.0;
  X_PP_pz = -999.0;
  X_PP_absEta = -999.0;
  X_PP_px = -999.0;
  X_PP_py = -999.0;
  
  // ============================================================
  //                    重置 ψ(2S) 候选变量
  // ============================================================
  Psi2S_mass = -999.0;
  Psi2S_VtxProb = -999.0;
  Psi2S_massErr = -999.0;
  Psi2S_pt = -999.0;
  Psi2S_pz = -999.0;
  Psi2S_absEta = -999.0;
  Psi2S_px = -999.0;
  Psi2S_py = -999.0;
  
  // ============================================================
  //                    重置 J/ψ 候选变量
  // ============================================================
  // J/ψ1 (参与 ψ(2S) 组合)
  Jpsi1_mass = -999.0;
  Jpsi1_VtxProb = -999.0;
  Jpsi1_massErr = -999.0;
  Jpsi1_pt = -999.0;
  Jpsi1_pz = -999.0;
  Jpsi1_absEta = -999.0;
  Jpsi1_px = -999.0;
  Jpsi1_py = -999.0;
  
  // J/ψ2 (第二个 J/ψ)
  Jpsi2_mass = -999.0;
  Jpsi2_VtxProb = -999.0;
  Jpsi2_massErr = -999.0;
  Jpsi2_pt = -999.0;
  Jpsi2_pz = -999.0;
  Jpsi2_absEta = -999.0;
  Jpsi2_px = -999.0;
  Jpsi2_py = -999.0;
  
  // ============================================================
  //                    重置 μ子变量 (4 个 μ 子)
  // ============================================================
  // μ子1 - J/ψ1 的正 μ 子
  mu1_pt = -999.0;
  mu1_pz = -999.0;
  mu1_absEta = -999.0;
  mu1_px = -999.0;
  mu1_py = -999.0;
  mu1_trackIso = -999.0;
  mu1_d0BS = -999.0;
  mu1_d0EBS = -999.0;
  mu1_d3dBS = -999.0;
  mu1_d3dEBS = -999.0;
  mu1_d0PV = -999.0;
  mu1_d0EPV = -999.0;
  mu1_dzPV = -999.0;
  mu1_dzEPV = -999.0;
  mu1_charge = -999.0;
  
  // μ子2 - J/ψ1 的负 μ 子
  mu2_pt = -999.0;
  mu2_pz = -999.0;
  mu2_absEta = -999.0;
  mu2_px = -999.0;
  mu2_py = -999.0;
  mu2_trackIso = -999.0;
  mu2_d0BS = -999.0;
  mu2_d0EBS = -999.0;
  mu2_d3dBS = -999.0;
  mu2_d3dEBS = -999.0;
  mu2_d0PV = -999.0;
  mu2_d0EPV = -999.0;
  mu2_dzPV = -999.0;
  mu2_dzEPV = -999.0;
  mu2_charge = -999.0;
  
  // μ子3 - J/ψ2 的正 μ 子
  mu3_pt = -999.0;
  mu3_pz = -999.0;
  mu3_absEta = -999.0;
  mu3_px = -999.0;
  mu3_py = -999.0;
  mu3_trackIso = -999.0;
  mu3_d0BS = -999.0;
  mu3_d0EBS = -999.0;
  mu3_d3dBS = -999.0;
  mu3_d3dEBS = -999.0;
  mu3_d0PV = -999.0;
  mu3_d0EPV = -999.0;
  mu3_dzPV = -999.0;
  mu3_dzEPV = -999.0;
  mu3_charge = -999.0;
  
  // μ子4 - J/ψ2 的负 μ 子
  mu4_pt = -999.0;
  mu4_pz = -999.0;
  mu4_absEta = -999.0;
  mu4_px = -999.0;
  mu4_py = -999.0;
  mu4_trackIso = -999.0;
  mu4_d0BS = -999.0;
  mu4_d0EBS = -999.0;
  mu4_d3dBS = -999.0;
  mu4_d3dEBS = -999.0;
  mu4_d0PV = -999.0;
  mu4_d0EPV = -999.0;
  mu4_dzPV = -999.0;
  mu4_dzEPV = -999.0;
  mu4_charge = -999.0;
  
  // ============================================================
  //                    重置 μ子 ID 计数变量
  // ============================================================
  nLooseMuons = 0;
  nTightMuons = 0;
  nSoftMuons = 0;
  nMediumMuons = 0;
  
  // ============================================================
  //                    重置 μ子触发匹配标志
  // ============================================================
  mu1_hasFilterMatch = 0;
  mu2_hasFilterMatch = 0;
  mu3_hasFilterMatch = 0;
  mu4_hasFilterMatch = 0;
  
  // ============================================================
  //                    重置 π 子变量 (2 个 π 子)
  // ============================================================
  // π子1 - 正 π 子
  pi1_pt = -999.0;
  pi1_pz = -999.0;
  pi1_absEta = -999.0;
  pi1_px = -999.0;
  pi1_py = -999.0;
  
  // π子2 - 负 π 子
  pi2_pt = -999.0;
  pi2_pz = -999.0;
  pi2_absEta = -999.0;
  pi2_px = -999.0;
  pi2_py = -999.0;
  
  // ============================================================
  //                    重置 δR 关联变量
  // ============================================================
  dR_mu1_mu2 = -999.0;
  dR_mu3_mu4 = -999.0;
  dR_pi1_pi2 = -999.0;
  dR_Psi2S_Jpsi1 = -999.0;
  dR_Psi2S_Jpsi2 = -999.0;
  dR_Psi2S_pi1 = -999.0;
  dR_Psi2S_pi2 = -999.0;
}

// ------------ method called once each job just after ending the event loop
// ------------
void MultiLepPAT::endJob() {
  X_One_Tree_->GetDirectory()->cd();
  X_One_Tree_->Write();
}

// define this as a plug-in
DEFINE_FWK_MODULE(MultiLepPAT);
