///////////////////////////////////////////////////////////
//
// 文件名: MultiLepPAT.h
// 作者: AliceQuen (Folked from zjianz/Onia2MuMu)
// 创建日期: 2026.04.30
// 描述: X → J/ψ ππ → 4μ 2π 分析的 Ntuple 产生器实现
//
// 主要功能:
//   - 事件预处理与 HLT 匹配
//   - μ子与径迹的质量筛选
//   - 多粒子组合与顶点拟合
//   - 三种质量假设的约束拟合
//   - 物理量计算与 TTree 存储
//
// 衰变拓扑: X → J/ψ π⁺π⁻ + J/ψ → μ⁺μ⁻μ⁺μ⁻ π⁺π⁻
//
// 依赖: CMSSW 框架, RecoVertex/KinematicFit, DataFormats/PatCandidates
//
///////////////////////////////////////////////////////////

#ifndef _MultiLepPAT_h
#define _MultiLepPAT_h

// ============================================================
//                    物理常数与截断参数
// ============================================================

// 质量约束拟合的精度 (GeV)
// 过小: 拟合不稳定; 过大: 失去约束效果
#define MASS_CONSTRAINT_PRECISION 0.001

// 粒子标称质量 (GeV) - 定义于 MultiLepPAT.cc
// MU_MASS    = 0.1056583745  (μ子)
// PI_MASS    = 0.13957039    (π± 子)
// JPSI_MASS  = 3.0969         (J/ψ)
// PSI2S_MASS = 3.686097       (ψ(2S))

// system include files
#include <memory>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "Math/VectorUtil.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>
#include <map>
#include <string>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// ============================================================
//                    类型定义
// ============================================================

using std::vector;
using namespace edm;
using namespace std;

////////////////////////////////////////////////////////////////
///
/// \struct JpsiCandidate
///
/// \brief J/ψ 候选结构体，存储拟合结果与关联信息
///
/// \details 包含 μ子对的普通顶点拟合结果、质量约束拟合结果、
///          四动量信息、触发匹配状态等。用于缓存预计算的 J/ψ 候选，
///          避免在多层循环中重复拟合。
///
////////////////////////////////////////////////////////////////
struct JpsiCandidate {
    // ============================================================
    //                    粒子引用
    // ============================================================
    edm::View<pat::Muon>::const_iterator muPlus;   ///< 正 μ 子迭代器
    edm::View<pat::Muon>::const_iterator muMinus;  ///< 负 μ 子迭代器

    // ============================================================
    //             普通顶点拟合结果 (无质量约束)
    // ============================================================
    double mass;                                    ///< 拟合不变质量 [GeV]
    double vtxProb;                                 ///< 顶点概率 (χ² 概率)
    double massErr;                                 ///< 质量误差 [GeV]
    ROOT::Math::PxPyPzMVector p4;                   ///< 四动量 (px, py, pz, mass)
    RefCountedKinematicParticle kinematicParticle;  ///< 拟合后的运动学粒子
    RefCountedKinematicVertex vertex;                ///< 拟合后的衰变顶点

    // ============================================================
    //             质量约束拟合结果 (共用一套)
    // ============================================================
    bool hasConstraintFit = false;                   ///< 是否有有效的约束拟合结果
    double constraintMass;                           ///< 约束拟合后的不变质量 [GeV]
    double constraintVtxProb;                        ///< 约束拟合后的顶点概率
    double constraintMassErr;                        ///< 约束拟合后的质量误差 [GeV]
    RefCountedKinematicParticle constraintParticle;  ///< 约束拟合后的运动学粒子
    RefCountedKinematicVertex constraintVertex;      ///< 约束拟合后的衰变顶点

    // ============================================================
    //                    候选类型标记
    // ============================================================
    bool isJpsiCandidate = false;                    ///< 质量在 J/ψ 窗口内
    bool isPsi2SCandidate = false;                   ///< 质量在 ψ(2S) 窗口内

    // ============================================================
    //                    其他关联信息
    // ============================================================
    reco::TransientTrack muonTT1;                   ///< 正 μ 子的瞬态径迹 (含磁场信息)
    reco::TransientTrack muonTT2;                   ///< 负 μ 子的瞬态径迹
    bool filterMatchPlus;                            ///< 正 μ 子触发过滤器匹配标志
    bool filterMatchMinus;                           ///< 负 μ 子触发过滤器匹配标志
};

////////////////////////////////////////////////////////////////
///
/// \class MultiLepPAT
///
/// \brief X 分析的 Ntuple 产生器主类
///
/// \details 继承自 edm::one::EDAnalyzer，实现完整的事例重建流程：
///          1. 事例级预筛选与 HLT 匹配
///          2. 粒子级质量筛选（μ子、径迹）
///          3. 分层组合与顶点拟合（2、4、6 径迹）
///          4. 三种质量假设下的约束拟合
///          5. 约 130 个物理量的计算与 TTree 存储
///
/// 衰变拓扑: X → J/ψ π⁺π⁻ + J/ψ → μ⁺μ⁻μ⁺μ⁻ π⁺π⁻
///
////////////////////////////////////////////////////////////////
class MultiLepPAT : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    ////////////////////////////////////////////////////////////////
    ///
    /// \brief 构造函数
    ///
    /// \param iConfig EDAnalyzer 配置参数集
    ///
    /// \details 执行以下初始化：
    ///          - 读取并存储 HLT 触发器名称列表
    ///          - 读取并存储 HLT 过滤器名称列表
    ///          - 按长度升序排序触发器名称（优化匹配速度）
    ///          - 初始化所有 EDToken
    ///
    ////////////////////////////////////////////////////////////////
    explicit MultiLepPAT(const ParameterSet&);

    ////////////////////////////////////////////////////////////////
    ///
    /// \brief 析构函数
    ///
    ////////////////////////////////////////////////////////////////
    ~MultiLepPAT();

private:
    ////////////////////////////////////////////////////////////////
    ///
    /// \brief 作业开始时调用一次
    ///
    /// \details 创建 TTree 并创建所有约 130 个分支。
    ///          分支按物理类别分组：事件信息、共振态属性、
    ///          μ子属性、π子属性、δR 关联变量。
    ///
    ////////////////////////////////////////////////////////////////
    virtual void beginJob() ;

    ////////////////////////////////////////////////////////////////
    ///
    /// \brief 每个 run 开始时调用
    ///
    /// \param iRun 当前 run 信息
    /// \param iSetup 事件配置
    ///
    /// \details 当前为空实现，预留用于 run 级别的初始化。
    ///
    ////////////////////////////////////////////////////////////////
    virtual void beginRun(Run const & iRun, EventSetup const& iSetup);

    ////////////////////////////////////////////////////////////////
    ///
    /// \brief 每个事件的主分析函数
    ///
    /// \param iEvent 当前事件数据
    /// \param iSetup 事件配置（包含磁场等条件数据）
    ///
    /// \details 执行完整的分析流程：
    ///          1. 数据获取与事件预筛选
    ///          2. HLT 触发匹配
    ///          3. 主顶点获取与筛选
    ///          4. 径迹预选择与电荷分组
    ///          5. μ子预选择与电荷分组
    ///          6. 6 层嵌套循环组合粒子
    ///          7. 顶点拟合与质量约束拟合
    ///          8. 物理量计算与 TTree 填充
    ///
    ////////////////////////////////////////////////////////////////
    virtual void analyze(const Event&, const EventSetup&);

    ////////////////////////////////////////////////////////////////
    ///
    /// \brief 作业结束时调用一次
    ///
    /// \details 当前为空实现，预留用于统计信息输出。
    ///
    ////////////////////////////////////////////////////////////////
    virtual void endJob() ;

    ////////////////////////////////////////////////////////////////
    ///
    /// \brief 重置所有输出变量为默认无效值
    ///
    /// \details 在每个事件开始时调用，将所有 ~130 个物理量
    ///          重置为哨兵值：
    ///          - 浮点型: -999.0
    ///          - 整型: 0
    ///          - 布尔型: false
    ///
    /// 重要: 此函数必须在每个 analyze() 开始时调用，
    ///       确保前一个候选的结果不影响当前候选
    ///
    ////////////////////////////////////////////////////////////////
    void resetVariables();

    ////////////////////////////////////////////////////////////////
    ///
    /// \brief 线程安全的质量约束对象初始化
    ///
    /// \details 使用 std::call_once 保证质量约束对象仅创建一次，
    ///          避免多线程环境下的重复构造。初始化三种粒子的
    ///          质量约束：J/ψ、ψ(2S)、X(3872)。
    ///
    ////////////////////////////////////////////////////////////////
    static void initMassConstraints();

    // ============================================================
    //                    EDM Token 定义
    // ============================================================
    edm::EDGetTokenT<reco::BeamSpot> gtbeamspotToken_;            ///< 束斑数据 Token
    edm::EDGetTokenT<reco::VertexCollection> gtprimaryVtxToken_;  ///< 主顶点集合 Token
    edm::EDGetTokenT<edm::View<pat::Muon>> gtpatmuonToken_;       ///< PAT Muon 视图 Token
    edm::EDGetTokenT<edm::TriggerResults> gttriggerToken_;         ///< 触发结果 Token
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackToken_;  ///< PF 候选径迹 Token (MiniAOD)

    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;  ///< 磁场条件数据 Token

    // ============================================================
    //                    触发配置
    // ============================================================
    vector<string> TriggersForJpsi_;     ///< J/ψ 触发的 HLT 路径名称列表
    vector<string> FiltersForJpsi_;       ///< 对应的 HLT 过滤器名称列表

    vector<string> TriggersForJpsi_sorted_;  ///< 按长度升序排序的触发器名称（优化匹配速度）
    size_t min_trigger_len_;                  ///< 最短触发器名称长度（快速剪枝）

    // ============================================================
    //                    质量约束对象
    // ============================================================
    /// 静态成员，程序启动时初始化一次，避免重复分配
    static std::unique_ptr<KinematicConstraint> jpsiMassConstraint_;    ///< J/ψ 质量约束
    static std::unique_ptr<KinematicConstraint> psi2sMassConstraint_;   ///< ψ(2S) 质量约束
    static std::unique_ptr<KinematicConstraint> x3872MassConstraint_;   ///< X(3872) 质量约束

    // ============================================================
    //                    ROOT 输出
    // ============================================================
    TTree* X_One_Tree_;                    ///< 输出 Ntuple 的 TTree 指针

    // ============================================================
    //                    事件标识变量
    // ============================================================
    unsigned int runNum;                    ///< 运行号
    unsigned int evtNum;                    ///< 事件号
    unsigned int lumiNum;                   ///< 亮度段号
    unsigned int nGoodPrimVtx;              ///< 有效主顶点数量

    // ============================================================
    //              X 候选变量 (假设 1)
    //  约束: mumupipi → ψ(2S), 另一个 mumu → J/ψ
    // ============================================================
    float X_PJ_mass;                        ///< 不变质量 [GeV]
    float X_PJ_VtxProb;                     ///< 顶点概率 (χ² 概率)
    float X_PJ_massErr;                     ///< 质量误差 [GeV]
    float X_PJ_pt;                          ///< 横动量 [GeV]
    float X_PJ_pz;                          ///< Z 方向动量分量 [GeV]
    float X_PJ_absEta;                      ///< 赝快度绝对值
    float X_PJ_px;                          ///< X 方向动量分量 [GeV]
    float X_PJ_py;                          ///< Y 方向动量分量 [GeV]

    // ============================================================
    //              X 候选变量 (假设 2)
    //  约束: mumupipi → X(3872), 另一个 mumu → J/ψ
    // ============================================================
    float X_XJ_mass;                        ///< 不变质量 [GeV]
    float X_XJ_VtxProb;                     ///< 顶点概率
    float X_XJ_massErr;                     ///< 质量误差 [GeV]
    float X_XJ_pt;                          ///< 横动量 [GeV]
    float X_XJ_pz;                          ///< Z 方向动量分量 [GeV]
    float X_XJ_absEta;                      ///< 赝快度绝对值
    float X_XJ_px;                          ///< X 方向动量分量 [GeV]
    float X_XJ_py;                          ///< Y 方向动量分量 [GeV]

    // ============================================================
    //              X 候选变量 (假设 3)
    //  约束: mumupipi → ψ(2S), 另一个 mumu → ψ(2S)
    // ============================================================
    float X_PP_mass;                        ///< 不变质量 [GeV]
    float X_PP_VtxProb;                     ///< 顶点概率
    float X_PP_massErr;                     ///< 质量误差 [GeV]
    float X_PP_pt;                          ///< 横动量 [GeV]
    float X_PP_pz;                          ///< Z 方向动量分量 [GeV]
    float X_PP_absEta;                      ///< 赝快度绝对值
    float X_PP_px;                          ///< X 方向动量分量 [GeV]
    float X_PP_py;                          ///< Y 方向动量分量 [GeV]

    // ============================================================
    //              ψ(2S) 候选变量 (4 径迹拟合)
    // ============================================================
    float Psi2S_mass;                       ///< 不变质量 [GeV]
    float Psi2S_VtxProb;                    ///< 顶点概率
    float Psi2S_massErr;                    ///< 质量误差 [GeV]
    float Psi2S_pt;                         ///< 横动量 [GeV]
    float Psi2S_pz;                         ///< Z 方向动量分量 [GeV]
    float Psi2S_absEta;                     ///< 赝快度绝对值
    float Psi2S_px;                         ///< X 方向动量分量 [GeV]
    float Psi2S_py;                         ///< Y 方向动量分量 [GeV]

    // ============================================================
    //              J/ψ₁ 候选变量 (参与 ψ(2S) 组合)
    // ============================================================
    float Jpsi1_mass;                       ///< 不变质量 [GeV]
    float Jpsi1_VtxProb;                    ///< 顶点概率
    float Jpsi1_massErr;                    ///< 质量误差 [GeV]
    float Jpsi1_pt;                         ///< 横动量 [GeV]
    float Jpsi1_pz;                         ///< Z 方向动量分量 [GeV]
    float Jpsi1_absEta;                     ///< 赝快度绝对值
    float Jpsi1_px;                         ///< X 方向动量分量 [GeV]
    float Jpsi1_py;                         ///< Y 方向动量分量 [GeV]

    // ============================================================
    //              J/ψ₂ 候选变量 (第二个 J/ψ)
    // ============================================================
    float Jpsi2_mass;                       ///< 不变质量 [GeV]
    float Jpsi2_VtxProb;                    ///< 顶点概率
    float Jpsi2_massErr;                    ///< 质量误差 [GeV]
    float Jpsi2_pt;                         ///< 横动量 [GeV]
    float Jpsi2_pz;                         ///< Z 方向动量分量 [GeV]
    float Jpsi2_absEta;                     ///< 赝快度绝对值
    float Jpsi2_px;                         ///< X 方向动量分量 [GeV]
    float Jpsi2_py;                         ///< Y 方向动量分量 [GeV]

    // ============================================================
    //                    μ子变量 (4 个 μ 子)
    // ============================================================

    /// μ₁ - 第一个正 μ子 (参与 J/ψ₁)
    float mu1_pt;                           ///< 横动量 [GeV]
    float mu1_pz;                           ///< Z 方向动量分量 [GeV]
    float mu1_absEta;                       ///< 赝快度绝对值
    float mu1_px;                           ///< X 方向动量分量 [GeV]
    float mu1_py;                           ///< Y 方向动量分量 [GeV]
    float mu1_trackIso;                     ///< 径迹隔离度
    float mu1_d0BS;                         ///< 相对于束斑的横向冲击参数 [cm]
    float mu1_d0EBS;                        ///< d0_BS 的误差 [cm]
    float mu1_d3dBS;                        ///< 相对于束斑的三维冲击参数 [cm]
    float mu1_d3dEBS;                       ///< d3d_BS 的误差 [cm]
    float mu1_d0PV;                         ///< 相对于主顶点的横向冲击参数 [cm]
    float mu1_d0EPV;                        ///< d0_PV 的误差 [cm]
    float mu1_dzPV;                         ///< 相对于主顶点的纵向冲击参数 [cm]
    float mu1_dzEPV;                        ///< dz_PV 的误差 [cm]
    float mu1_charge;                       ///< 电荷

    /// μ₂ - 第一个负 μ子 (参与 J/ψ₁)
    float mu2_pt;                           ///< 横动量 [GeV]
    float mu2_pz;                           ///< Z 方向动量分量 [GeV]
    float mu2_absEta;                       ///< 赝快度绝对值
    float mu2_px;                           ///< X 方向动量分量 [GeV]
    float mu2_py;                           ///< Y 方向动量分量 [GeV]
    float mu2_trackIso;                     ///< 径迹隔离度
    float mu2_d0BS;                         ///< 相对于束斑的横向冲击参数 [cm]
    float mu2_d0EBS;                        ///< d0_BS 的误差 [cm]
    float mu2_d3dBS;                        ///< 相对于束斑的三维冲击参数 [cm]
    float mu2_d3dEBS;                       ///< d3d_BS 的误差 [cm]
    float mu2_d0PV;                         ///< 相对于主顶点的横向冲击参数 [cm]
    float mu2_d0EPV;                        ///< d0_PV 的误差 [cm]
    float mu2_dzPV;                         ///< 相对于主顶点的纵向冲击参数 [cm]
    float mu2_dzEPV;                        ///< dz_PV 的误差 [cm]
    float mu2_charge;                       ///< 电荷

    /// μ₃ - 第二个正 μ子 (参与 J/ψ₂)
    float mu3_pt;                           ///< 横动量 [GeV]
    float mu3_pz;                           ///< Z 方向动量分量 [GeV]
    float mu3_absEta;                       ///< 赝快度绝对值
    float mu3_px;                           ///< X 方向动量分量 [GeV]
    float mu3_py;                           ///< Y 方向动量分量 [GeV]
    float mu3_trackIso;                     ///< 径迹隔离度
    float mu3_d0BS;                         ///< 相对于束斑的横向冲击参数 [cm]
    float mu3_d0EBS;                        ///< d0_BS 的误差 [cm]
    float mu3_d3dBS;                        ///< 相对于束斑的三维冲击参数 [cm]
    float mu3_d3dEBS;                       ///< d3d_BS 的误差 [cm]
    float mu3_d0PV;                         ///< 相对于主顶点的横向冲击参数 [cm]
    float mu3_d0EPV;                        ///< d0_PV 的误差 [cm]
    float mu3_dzPV;                         ///< 相对于主顶点的纵向冲击参数 [cm]
    float mu3_dzEPV;                        ///< dz_PV 的误差 [cm]
    float mu3_charge;                       ///< 电荷

    /// μ₄ - 第二个负 μ子 (参与 J/ψ₂)
    float mu4_pt;                           ///< 横动量 [GeV]
    float mu4_pz;                           ///< Z 方向动量分量 [GeV]
    float mu4_absEta;                       ///< 赝快度绝对值
    float mu4_px;                           ///< X 方向动量分量 [GeV]
    float mu4_py;                           ///< Y 方向动量分量 [GeV]
    float mu4_trackIso;                     ///< 径迹隔离度
    float mu4_d0BS;                         ///< 相对于束斑的横向冲击参数 [cm]
    float mu4_d0EBS;                        ///< d0_BS 的误差 [cm]
    float mu4_d3dBS;                        ///< 相对于束斑的三维冲击参数 [cm]
    float mu4_d3dEBS;                       ///< d3d_BS 的误差 [cm]
    float mu4_d0PV;                         ///< 相对于主顶点的横向冲击参数 [cm]
    float mu4_d0EPV;                        ///< d0_PV 的误差 [cm]
    float mu4_dzPV;                         ///< 相对于主顶点的纵向冲击参数 [cm]
    float mu4_dzEPV;                        ///< dz_PV 的误差 [cm]
    float mu4_charge;                       ///< 电荷

    /// μ子 ID 计数变量
    int nLooseMuons;                         ///< 满足 Loose Muon ID 的数量
    int nTightMuons;                         ///< 满足 Tight Muon ID 的数量
    int nSoftMuons;                          ///< 满足 Soft Muon ID 的数量
    int nMediumMuons;                        ///< 满足 Medium Muon ID 的数量

    /// μ子触发匹配标志
    bool mu1_hasFilterMatch;                 ///< μ₁ 是否通过触发过滤器匹配
    bool mu2_hasFilterMatch;                 ///< μ₂ 是否通过触发过滤器匹配
    bool mu3_hasFilterMatch;                 ///< μ₃ 是否通过触发过滤器匹配
    bool mu4_hasFilterMatch;                 ///< μ₄ 是否通过触发过滤器匹配

    // ============================================================
    //                    π子变量 (2 个 π 子)
    // ============================================================

    /// π₁ - 正 π 子 (参与 ψ(2S))
    float pi1_pt;                           ///< 横动量 [GeV]
    float pi1_pz;                           ///< Z 方向动量分量 [GeV]
    float pi1_absEta;                       ///< 赝快度绝对值
    float pi1_px;                           ///< X 方向动量分量 [GeV]
    float pi1_py;                           ///< Y 方向动量分量 [GeV]

    /// π₂ - 负 π 子 (参与 ψ(2S))
    float pi2_pt;                           ///< 横动量 [GeV]
    float pi2_pz;                           ///< Z 方向动量分量 [GeV]
    float pi2_absEta;                       ///< 赝快度绝对值
    float pi2_px;                           ///< X 方向动量分量 [GeV]
    float pi2_py;                           ///< Y 方向动量分量 [GeV]

    // ============================================================
    //                    δR 关联变量
    // ============================================================
    float dR_mu1_mu2;                        ///< μ₁ 与 μ₂ 的角度距离
    float dR_mu3_mu4;                        ///< μ₃ 与 μ₄ 的角度距离
    float dR_pi1_pi2;                        ///< π₁ 与 π₂ 的角度距离
    float dR_Psi2S_Jpsi1;                    ///< ψ(2S) 与 J/ψ₁ 的角度距离
    float dR_Psi2S_Jpsi2;                    ///< ψ(2S) 与 J/ψ₂ 的角度距离
    float dR_Psi2S_pi1;                       ///< ψ(2S) 与 π₁ 的角度距离
    float dR_Psi2S_pi2;                       ///< ψ(2S) 与 π₂ 的角度距离
};

#endif