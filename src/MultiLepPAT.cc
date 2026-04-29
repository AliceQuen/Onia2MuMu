// system include files

#include "TLorentzVector.h"
#include <limits>
#if DEBUG == 1
#include <chrono>
#endif
#if DEBUG == 2
#include <iostream>
#include <iomanip>
#endif
#include <Math/Vector4D.h>
// user include files
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
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "MagneticField/Engine/interface/MagneticField.h"

//////////This is necessary for lumicalc///////
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

// about photon
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <boost/foreach.hpp>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // MINIAOD

#define MU_MASS 0.1056583745
#define MU_MASSERR (MU_MASS * 1e-6)
#define PI_MASS 0.13957039
#define PI_MASSERR (PI_MASS * 1e-6)

// Particle mass nominal values for constraints
#define JPSI_MASS_NOMINAL 3.0969
#define X3872_MASS_NOMINAL 3.872
#define PSI2S_MASS_NOMINAL 3.686097

typedef math::Error<3>::type CovarianceMatrix;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>>
    SMatrixSym3D;
using namespace edm;
using namespace reco;
using namespace std;

// constructors and destructor
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
  // Sort triggers by length ascending - shorter patterns checked first
  // This allows early exit on average since shorter patterns are more likely to match
  TriggersForJpsi_sorted_ = TriggersForJpsi_;
  std::sort(TriggersForJpsi_sorted_.begin(), TriggersForJpsi_sorted_.end(),
            [](const std::string& a, const std::string& b) { return a.size() < b.size(); });
  
  // Find minimum trigger length for early pruning
  min_trigger_len_ = std::numeric_limits<size_t>::max();
  for (const auto& t : TriggersForJpsi_sorted_) {
    if (t.size() < min_trigger_len_) min_trigger_len_ = t.size();
  }
  
  // get token here for four-muon;
  gtbeamspotToken_ = consumes<BeamSpot>(edm::InputTag("offlineBeamSpot"));
  gtprimaryVtxToken_ = consumes<VertexCollection>(
      edm::InputTag("offlineSlimmedPrimaryVertices")); // MINIAOD
  gtpatmuonToken_ =
      consumes<edm::View<pat::Muon>>(edm::InputTag("slimmedMuons")); // MINIAOD
  gttriggerToken_ =
      consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"));
  trackToken_ = consumes<edm::View<pat::PackedCandidate>>(
      edm::InputTag("packedPFCandidates")); // MINIAOD
}

MultiLepPAT::~MultiLepPAT() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}
// member functions

// ------------ method called to for each event  ------------
void MultiLepPAT::analyze(const edm::Event &iEvent,
                          const edm::EventSetup &iSetup) {

  edm::Handle<edm::View<pat::Muon>> thePATMuonHandle; //  MINIAOD
  iEvent.getByToken(gtpatmuonToken_, thePATMuonHandle);
  // Skip events with less than 4 muons
  if (thePATMuonHandle->size() < 4) {
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    total_return_++;
    #endif
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: Less than 4 muons (line " << __LINE__ << ")" << std::endl;
    #endif
    return;
  }
  // Reset all variables at the beginning of each event
  resetVariables();

  const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_);

  // Get run infomation
  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();
  lumiNum = iEvent.id().luminosityBlock();

  // Get HLT information
#if DEBUG == 1
  auto start_hlt = std::chrono::high_resolution_clock::now();
#endif
  edm::Handle<edm::TriggerResults> hltresults;
  bool Error_t = false;
  bool HLT_match = false;
  try {
    iEvent.getByToken(gttriggerToken_, hltresults);
  } catch (...) {
    Error_t = true;
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    total_return_++;
    #endif
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: Couldn't get handle on HLT Trigger (line " << __LINE__ << ")" << std::endl;
    #endif
    return;
  }
  if (Error_t || !hltresults.isValid()) {
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    total_return_++;
    #endif
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: Invalid HLT results (line " << __LINE__ << ")" << std::endl;
    #endif
    return;
  } else {
    int ntrigs = hltresults->size();
    if (ntrigs == 0) {
      #if DEBUG == 2
      return_counts_[__LINE__]++;
      total_return_++;
      #endif
      #if DEBUG == 1
      std::cout << "[DEBUG] Return encountered in analyze(), location: No triggers in TriggerResults (line " << __LINE__ << ")" << std::endl;
      #endif
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
  #if DEBUG == 1
  auto end_hlt = std::chrono::high_resolution_clock::now();
  auto duration_hlt = std::chrono::duration_cast<std::chrono::milliseconds>(end_hlt - start_hlt).count();
  std::cout << "[DEBUG] GetHLTInformation 运行时间: " << duration_hlt << " ms" << std::endl;
  #endif
  if (!HLT_match){
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    total_return_++;
    #endif
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: No HLT match found (line " << __LINE__ << ")" << std::endl;
    #endif
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
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    total_return_++;
    #endif
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: No good primary vertex found (line " << __LINE__ << ")" << std::endl;
    #endif
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
        #if DEBUG == 2
        return_counts_[__LINE__]++;
        total_return_++;
        #endif
        #if DEBUG == 1
        std::cout << "[DEBUG] Return encountered in analyze(), location: No beam spot available (line " << __LINE__ << ")" << std::endl;
        #endif
        return;
      }
    thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
  }

  // Check Muon and Pion Track
#if DEBUG == 1
  auto start_track = std::chrono::high_resolution_clock::now();
#endif
  // Pre-selection parameters according to pre-cuts.md
  const double pionDRcut = 0.7;
  const double muon_pt_cut_low = 3.5;
  const double muon_absEta_cut1 = 1.2;
  const double muon_pt_cut_slope = 5.47;
  const double muon_pt_cut_coeff = 1.89;
  const double muon_absEta_cut2_low = 1.2;
  const double muon_absEta_cut2_high = 2.1;
  const double muon_pt_cut_low2 = 1.5;
  const double muon_absEta_cut3_low = 2.1;
  const double muon_absEta_cut3_high = 2.4;
  
  const double track_sigma_pt_over_pt_cut = 0.1;
  const int track_nhits_cut = 10;
  const double track_chi2_over_ndf_cut = 0.18;
  const double track_pt_cut = 0.5;
  const double track_absEta_cut = 2.4;
  
  const double jpsi_vtxprob_cut = 0.01;
  const double jpsi_mass_window = 0.15;
  const double jpsi_nominal_mass = 3.0969;
  const double psi2s_nominal_mass = 3.686097;
  
  const double psi2s_vtxprob_cut = 0.005;
  const double psi2s_pt_cut = 4.0;

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

  // Selection criteria from pre-cuts.md:
  // 1. pt > 0.5 && absEta < 2.4
  // 2. sigma_pt / pt < 0.1
  // 3. Nhits >= 11
  // 4. chi^2 / ndf < 0.18
  // 5. 主动过滤muon径迹 (PDG ID = ±13)
  for (edm::View<pat::PackedCandidate>::const_iterator iTrackc = theTrackHandle->begin();
       iTrackc != theTrackHandle->end(); ++iTrackc) {
    // Apply track pre-selection criteria
    const reco::Track* track = iTrackc->bestTrack();
    if (!track) {
      continue;
    }

    // 1. highPurity质量要求
    if (!(track->quality(reco::TrackBase::qualityByName("highPurity")))) {
      continue;
    }

    // 2. pT > 0.5 && absEta < 2.4
    if (!(iTrackc->pt() > track_pt_cut && fabs(iTrackc->eta()) < track_absEta_cut)) {
      continue;
    }

    // 3. sigma_pt / pt < 0.1
    double sigma_pt_over_pt = track->ptError() / track->pt();
    if (!(sigma_pt_over_pt < track_sigma_pt_over_pt_cut)) {
      continue;
    }

    // 4. Nhits >= 11
    unsigned int Nhits = track->numberOfValidHits();
    if (Nhits < track_nhits_cut) {
      continue;
    }

    // 5. chi^2 / ndf < 0.18
    if (track->normalizedChi2() / Nhits >= track_chi2_over_ndf_cut) {
      continue;
    }

    // 6. ===== 主动过滤muon径迹 =====
    // 使用PDG ID识别muon并排除，避免后续再移除
    // PDG ID: ±13 对应muon，±211对应pion，±321对应kaon，±2212对应proton
    int pdgId = abs(iTrackc->pdgId());
    if (pdgId == 13) {
      continue;  // 排除所有muon径迹
    }

    // 所有筛选条件通过，按电荷分组
    if (track->charge() > 0) {
      trackPlus.push_back(iTrackc);
    } else {
      trackMinus.push_back(iTrackc);
    }
  }

  #if DEBUG == 1
  std::cout << "[DEBUG] Track pre-selection with charge grouping: total=" 
            << theTrackHandle->size() 
            << ", plus=" << trackPlus.size() << ", minus=" << trackMinus.size() << std::endl;
  #endif

  // 径迹配对要求检查：至少各有1条正负径迹
  if (trackPlus.empty() || trackMinus.empty()) {
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    #endif
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: Insufficient tracks after pre-selection " 
              << "(plus=" << trackPlus.size() << ", minus=" << trackMinus.size() << ") (line " << __LINE__ << ")" << std::endl;
    #endif
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
    // 1. 软μ子要求
    if (!iMuonP->isSoftMuon(thePrimaryV)) {
      #if DEBUG == 2
      continue_counts_[__LINE__]++;
      #endif
      continue;
    }

    // 2. pT-η分段要求
    const double pt = iMuonP->pt();
    const double absEta = fabs(iMuonP->eta());
    bool passPtEta = false;
    if ((pt > muon_pt_cut_low && absEta < muon_absEta_cut1) ||
        (pt > (muon_pt_cut_slope - muon_pt_cut_coeff * absEta) && 
         absEta > muon_absEta_cut2_low && absEta < muon_absEta_cut2_high) ||
        (pt > muon_pt_cut_low2 && 
         absEta > muon_absEta_cut3_low && absEta < muon_absEta_cut3_high)) {
      passPtEta = true;
    }
    if (!passPtEta) {
      #if DEBUG == 2
      continue_counts_[__LINE__]++;
      #endif
      continue;
    }

    // 3. track有效性检查（安全校验）
    if (iMuonP->track().isNull()) {
      #if DEBUG == 2
      continue_counts_[__LINE__]++;
      #endif
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

  #if DEBUG == 1
  std::cout << "[DEBUG] Muon pre-selection with charge grouping: total=" 
            << thePATMuonHandle->size() 
            << ", plus=" << muPlus.size() << ", minus=" << muMinus.size() << std::endl;
  #endif

  #if DEBUG == 1
  auto end_track = std::chrono::high_resolution_clock::now();
  auto duration_track = std::chrono::duration_cast<std::chrono::milliseconds>(end_track - start_track).count();
  std::cout << "[DEBUG] CheckMuonPionTrack 运行时间: " << duration_track << " ms" << std::endl;
  #endif

  // 电荷要求检查：至少2个正μ和2个负μ
  if (muPlus.size() < 2 || muMinus.size() < 2) {
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    #endif
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: Insufficient muons after pre-selection " 
              << "(plus=" << muPlus.size() << ", minus=" << muMinus.size() << ") (line " << __LINE__ << ")" << std::endl;
    #endif
    return;
  }

  // Start processing muon and track fit
#if DEBUG == 1
  auto start_fit = std::chrono::high_resolution_clock::now();
#endif
  // ========== 第一个J/psi候选构建：使用预分组的muon列表 ==========
  // 遍历所有正负μ子组合（已保证电荷相反，无需额外检查）
  // 优化效果：组合数从 N*(N-1)/2 减少到 N_plus * N_minus，减少约50%
  for (const auto& muPlusIter1 : muPlus) {
    TrackRef muTrack1 = muPlusIter1->track();
    
    for (const auto& muMinusIter1 : muMinus) {
      TrackRef muTrack2 = muMinusIter1->track();
      
      // 质量预筛选（1-4.5 GeV）
      double invMass = (muPlusIter1->p4() + muMinusIter1->p4()).mass();
      if (!(1.0 < invMass && invMass < 4.5)) {
        #if DEBUG == 2
        continue_counts_[__LINE__]++;
        #endif
        #if DEBUG == 1
        std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi1 invariant mass out of range (line " << __LINE__ << ")" << std::endl;
        #endif
        continue;
      }
      TransientTrack muon1TT(muTrack1, &(bFieldHandle)); // MINIAOD
      TransientTrack muon2TT(muTrack2, &(bFieldHandle)); // MINIAOD
      KinematicParticleFactoryFromTransientTrack pmumuFactory;
      // The mass of a muon and the insignificant mass sigma
      // to avoid singularities in the covariance matrix.
      ParticleMass muon_mass = MU_MASS; // pdg mass
      float muon_sigma = MU_MASSERR;
      // initial chi2 and ndf before kinematic fits.
      float chi = 0.;
      float ndf = 0.;
      vector<RefCountedKinematicParticle> muonParticles;
      muonParticles.push_back(
          pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
      muonParticles.push_back(
          pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));
      KinematicParticleVertexFitter mu12_fitter;
      RefCountedKinematicTree Jpsi1VertexFitTree;
      Error_t = false;
      try {
        Jpsi1VertexFitTree = mu12_fitter.fit(muonParticles);
      } catch (...) {
        Error_t = true;
      }
      if (Error_t || !Jpsi1VertexFitTree->isValid()) {
#if DEBUG == 2
        continue_counts_[__LINE__]++;
        total_continue_++;
#endif
        #if DEBUG == 1
        std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi1 vertex fit failed (line " << __LINE__ << ")" << std::endl;
        #endif
        continue;
      }
      Jpsi1VertexFitTree->movePointerToTheTop();
      RefCountedKinematicParticle Jpsi1_vFit_noMC =
          Jpsi1VertexFitTree->currentParticle();
      RefCountedKinematicVertex Jpsi1_vFit_vertex_noMC =
          Jpsi1VertexFitTree->currentDecayVertex();
      KinematicParameters mymumupara =
          Jpsi1_vFit_noMC->currentState().kinematicParameters();
      double Jpsi1_vtxprob = ChiSquaredProbability(
          (double)(Jpsi1_vFit_vertex_noMC->chiSquared()),
          (double)(Jpsi1_vFit_vertex_noMC->degreesOfFreedom()));
      
      // Jpsi pre-selection according to pre-cuts.md
      // 1. opposite muon charges (already checked earlier)
      // 2. VtxProb > 1%
      if (Jpsi1_vtxprob < jpsi_vtxprob_cut) {
#if DEBUG == 2
        continue_counts_[__LINE__]++;
        total_continue_++;
#endif
        #if DEBUG == 1
        std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi1 vertex probability below cut (line " << __LINE__ << ")" << std::endl;
        #endif
        continue;
      }
      // 3. |mass - Jpsi_mass| < 0.15
      double jpsi1_mass = Jpsi1_vFit_noMC->currentState().mass();
      if (fabs(jpsi1_mass - jpsi_nominal_mass) > jpsi_mass_window) {
#if DEBUG == 2
        continue_counts_[__LINE__]++;
        total_continue_++;
#endif
        #if DEBUG == 1
        std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi1 mass outside window (line " << __LINE__ << ")" << std::endl;
        #endif
        continue;
      }

      TLorentzVector P4_mu1;
      P4_mu1.SetPtEtaPhiM(muPlusIter1->track()->pt(), muPlusIter1->track()->eta(),
                          muPlusIter1->track()->phi(), MU_MASS);
      TLorentzVector P4_mu2;
      P4_mu2.SetPtEtaPhiM(muMinusIter1->track()->pt(), muMinusIter1->track()->eta(),
                          muMinusIter1->track()->phi(), MU_MASS);

      // ========== 第二个J/psi候选构建：使用预分组的muon列表 ==========
      // 跳过已用于第一个J/psi的muon，直接遍历正负组合
      for (const auto& muPlusIter2 : muPlus) {
        // 跳过已使用的muon
        if (muPlusIter2 == muPlusIter1 || muPlusIter2 == muMinusIter1) {
#if DEBUG == 2
          continue_counts_[__LINE__]++;
#endif
          continue;
        }
        TrackRef muTrack3 = muPlusIter2->track();

        for (const auto& muMinusIter2 : muMinus) {
          // 跳过已使用的muon
          if (muMinusIter2 == muPlusIter1 || muMinusIter2 == muMinusIter1) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
#endif
            continue;
          }
          TrackRef muTrack4 = muMinusIter2->track();
          
          // 质量预筛选（1-4.5 GeV）
          double invMass2 = (muPlusIter2->p4() + muMinusIter2->p4()).mass();
          if (!(1.0 < invMass2 && invMass2 < 4.5)) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
#endif
            #if DEBUG == 1
            std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi2 invariant mass out of range (line " << __LINE__ << ")" << std::endl;
            #endif
            continue;
          }

          TransientTrack muon3TT(muTrack3, &(bFieldHandle)); // MINIAOD
          TransientTrack muon4TT(muTrack4, &(bFieldHandle)); // MINIAOD
          KinematicParticleFactoryFromTransientTrack pmumuFactory34;
          ParticleMass muon_mass = MU_MASS; // pdg mass
          float muon_sigma = MU_MASSERR;
          float chi = 0.;
          float ndf = 0.;
          vector<RefCountedKinematicParticle> muonParticles34;
          muonParticles34.push_back(pmumuFactory34.particle(
              muon3TT, muon_mass, chi, ndf, muon_sigma));
          muonParticles34.push_back(pmumuFactory34.particle(
              muon4TT, muon_mass, chi, ndf, muon_sigma));
          KinematicParticleVertexFitter mu34_fitter;
          RefCountedKinematicTree Jpsi2VertexFitTree;
          Error_t = false;
          try {
            Jpsi2VertexFitTree = mu34_fitter.fit(muonParticles34);
          } catch (...) {
            Error_t = true;
          }
          if (Error_t || !Jpsi2VertexFitTree->isValid()) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
#endif
            #if DEBUG == 1
            std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi2 vertex fit failed (line " << __LINE__ << ")" << std::endl;
            #endif
            continue;
          }
          Jpsi2VertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle Jpsi2_vFit_noMC =
              Jpsi2VertexFitTree->currentParticle();
          RefCountedKinematicVertex Jpsi2_vFit_vertex_noMC =
              Jpsi2VertexFitTree->currentDecayVertex();
          KinematicParameters mymumupara2 =
              Jpsi2_vFit_noMC->currentState().kinematicParameters();
          double Jpsi2_vtxprob = ChiSquaredProbability(
              (double)(Jpsi2_vFit_vertex_noMC->chiSquared()),
              (double)(Jpsi2_vFit_vertex_noMC->degreesOfFreedom()));
          
          // Jpsi pre-selection according to pre-cuts.md
          // 1. opposite muon charges (already checked earlier)
          // 2. VtxProb > 1%
          if (Jpsi2_vtxprob < jpsi_vtxprob_cut) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
#endif
            #if DEBUG == 1
            std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi2 vertex probability below cut (line " << __LINE__ << ")" << std::endl;
            #endif
            continue;
          }
          // 3. Allow two mass windows:
          //    - |mass - Jpsi_mass (3.1)| < 0.15
          //    - |mass - Psi2S_mass (3.686)| < 0.15 (for Jpsi2 being Psi2S)
          double jpsi2_mass = Jpsi2_vFit_noMC->currentState().mass();
          bool passJpsi2Mass = (fabs(jpsi2_mass - jpsi_nominal_mass) <= jpsi_mass_window) || 
                               (fabs(jpsi2_mass - psi2s_nominal_mass) <= jpsi_mass_window);
          if (!passJpsi2Mass) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
#endif
            #if DEBUG == 1
            std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi2 mass outside window (line " << __LINE__ << ")" << std::endl;
            #endif
            continue;
          }
          
          // Trigger pre-selection according to pre-cuts.md
          // mu1 mu2 filter matched OR mu3 mu4 filter matched
          size_t mu1_idx = std::distance(validMuons.cbegin(), iMuon1ID);
          size_t mu2_idx = std::distance(validMuons.cbegin(), iMuon2ID);
          size_t mu3_idx = std::distance(validMuons.cbegin(), iMuon3ID);
          size_t mu4_idx = std::distance(validMuons.cbegin(), iMuon4ID);
          
          bool mu1_matched = (mu1_idx < muonFilterMatches.size() && muonFilterMatches[mu1_idx]);
          bool mu2_matched = (mu2_idx < muonFilterMatches.size() && muonFilterMatches[mu2_idx]);
          bool mu3_matched = (mu3_idx < muonFilterMatches.size() && muonFilterMatches[mu3_idx]);
          bool mu4_matched = (mu4_idx < muonFilterMatches.size() && muonFilterMatches[mu4_idx]);
          
          bool passTrigger = ((mu1_matched && mu2_matched) || (mu3_matched && mu4_matched));
          if (!passTrigger) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
#endif
            #if DEBUG == 1
            std::cout << "[DEBUG] Continue encountered in analyze(), location: Trigger matching requirement not satisfied (line " << __LINE__ << ")" << std::endl;
            #endif
            continue;
          }

          // ========== Psi2S or X3872 (4mu 2pi) candidate ==========
          // 使用预分组的trackPlus/trackMinus遍历径迹对
          // 优化效果：无需电荷检查，组合数减少约50%
          for (const auto& trackPlusIter : trackPlus) {
            if (!trackPlusIter->hasTrackDetails() || trackPlusIter->charge() == 0) {
#if DEBUG == 2
              continue_counts_[__LINE__]++;
#endif
              continue;
            }
            const reco::Track* track1 = trackPlusIter->bestTrack();
            if (!track1) {
#if DEBUG == 2
              continue_counts_[__LINE__]++;
#endif
              continue;
            }

            for (const auto& trackMinusIter : trackMinus) {
              if (!trackMinusIter->hasTrackDetails() || trackMinusIter->charge() == 0) {
#if DEBUG == 2
                continue_counts_[__LINE__]++;
#endif
                continue;
              }
              const reco::Track* track2 = trackMinusIter->bestTrack();
              if (!track2) {
#if DEBUG == 2
                continue_counts_[__LINE__]++;
#endif
                continue;
              }

              TLorentzVector P4_Track1, P4_Track2, P4_Jpsipipi;
              P4_Track1.SetPtEtaPhiM(trackPlusIter->pt(), trackPlusIter->eta(),
                                     trackPlusIter->phi(), PI_MASS);
              P4_Track2.SetPtEtaPhiM(trackMinusIter->pt(), trackMinusIter->eta(),
                                     trackMinusIter->phi(), PI_MASS);
              P4_Jpsipipi = P4_mu1 + P4_mu2 + P4_Track1 + P4_Track2;

              if (P4_Track1.DeltaR(P4_Jpsipipi) > pionDRcut) {
                #if DEBUG == 2
                continue_counts_[__LINE__]++;
                total_continue_++;
                #endif
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Pion track1 DeltaR to Jpsi-pipi exceeds cut (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }
              if (P4_Track2.DeltaR(P4_Jpsipipi) > pionDRcut) {
#if DEBUG == 2
                continue_counts_[__LINE__]++;
                total_continue_++;
#endif
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Pion track2 DeltaR to Jpsi-pipi exceeds cut (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }

              TransientTrack trackTT1(*track1, &(bFieldHandle));
              TransientTrack trackTT2(*track2, &(bFieldHandle));
              KinematicParticleFactoryFromTransientTrack JPiPiFactory;
              // The mass of a muon and the insignificant mass sigma
              // to avoid singularities in the covariance matrix.
              ParticleMass pion_mass = PI_MASS; // pdg mass
              float pion_sigma = PI_MASSERR;
              // initial chi2 and ndf before kinematic fits.
              float chi = 0.;
              float ndf = 0.;

              // ========== 四粒子带J/psi质量约束拟合 ==========
              // mumupipi系统：将muon1和muon2的不变质量约束到J/psi标称质量
              vector<RefCountedKinematicParticle> JPiPiParticles;
              JPiPiParticles.push_back(JPiPiFactory.particle(
                  trackTT1, pion_mass, chi, ndf, pion_sigma));
              JPiPiParticles.push_back(JPiPiFactory.particle(
                  trackTT2, pion_mass, chi, ndf, pion_sigma));
              JPiPiParticles.push_back(pmumuFactory.particle(
                  muon1TT, muon_mass, chi, ndf, muon_sigma));
              JPiPiParticles.push_back(pmumuFactory.particle(
                  muon2TT, muon_mass, chi, ndf, muon_sigma));

              // 创建J/psi质量约束：约束muon1和muon2（索引2和3）的不变质量到3.0969 GeV
              KinematicConstraint* jpsiMassConstraint = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 2, 3);
              vector<KinematicConstraint*> JPiPiConstraints;
              JPiPiConstraints.push_back(jpsiMassConstraint);

              // 使用带约束的顶点拟合器
              KinematicConstrainedVertexFitter JPiPi_fitter;
              RefCountedKinematicTree JPiPiVertexFitTree;
              Error_t = false;
              try {
                JPiPiVertexFitTree = JPiPi_fitter.fit(JPiPiParticles, JPiPiConstraints);
              } catch (...) {
                Error_t = true;
              }
              if (Error_t || !(JPiPiVertexFitTree->isValid())) {
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Psi2S constrained vertex fit failed (line " << __LINE__ << ")" << std::endl;
                #endif
                delete jpsiMassConstraint;
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
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Psi2S mass exceeds 4.5 (line " << __LINE__ << ")" << std::endl;
                #endif
                delete jpsiMassConstraint;
                continue;
              }
              
              // Psi2S pre-selection according to pre-cuts.md
              // 1. VtxProb > 0.5%
              if (JPiPi_vtxprob < psi2s_vtxprob_cut) {
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Psi2S vertex probability below cut (line " << __LINE__ << ")" << std::endl;
                #endif
                delete jpsiMassConstraint;
                continue;
              }
              
              // 2. pt > 4
              double px = JPiPi_vFit_constrained->currentState().kinematicParameters().momentum().x();
              double py = JPiPi_vFit_constrained->currentState().kinematicParameters().momentum().y();
              double psi2s_pt = sqrt(px*px + py*py);
              if (psi2s_pt <= psi2s_pt_cut) {
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Psi2S pT below cut (line " << __LINE__ << ")" << std::endl;
                #endif
                delete jpsiMassConstraint;
                continue;
              }

              // 清理约束对象内存
              delete jpsiMassConstraint;

              // ========== 保存Psi2S拟合结果 ==========
              // 从J/psi质量约束的四粒子拟合中提取Psi2S相关变量
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

              // ========== 六粒子三种假设质量约束拟合 ==========
              // 构建6粒子列表：pi1(0), pi2(1), mu1(2), mu2(3), mu3(4), mu4(5)
              vector<RefCountedKinematicParticle> X_Particles;
              X_Particles.push_back(JPiPiFactory.particle(
                  trackTT1, pion_mass, chi, ndf, pion_sigma));
              X_Particles.push_back(JPiPiFactory.particle(
                  trackTT2, pion_mass, chi, ndf, pion_sigma));
              X_Particles.push_back(pmumuFactory.particle(
                  muon1TT, muon_mass, chi, ndf, muon_sigma));
              X_Particles.push_back(pmumuFactory.particle(
                  muon2TT, muon_mass, chi, ndf, muon_sigma));
              X_Particles.push_back(pmumuFactory.particle(
                  muon3TT, muon_mass, chi, ndf, muon_sigma));
              X_Particles.push_back(pmumuFactory.particle(
                  muon4TT, muon_mass, chi, ndf, muon_sigma));

              // ========== 假设PJ: mumupipi(Jpsi) + mumup(Jpsi) ==========
              // 约束: (mu1,mu2)→J/psi 和 (mu3,mu4)→J/psi
              {
                vector<KinematicConstraint*> PJConstraints;
                KinematicConstraint* JpsiConstraint1 = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 2, 3);
                KinematicConstraint* JpsiConstraint2 = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 4, 5);
                PJConstraints.push_back(JpsiConstraint1);
                PJConstraints.push_back(JpsiConstraint2);

                KinematicConstrainedVertexFitter PJ_fitter;
                RefCountedKinematicTree PJ_VertexFitTree;
                bool PJ_Error = false;
                try {
                  PJ_VertexFitTree = PJ_fitter.fit(X_Particles, PJConstraints);
                } catch (...) {
                  PJ_Error = true;
                }
                if (PJ_Error || !(PJ_VertexFitTree->isValid())) {
                  #if DEBUG == 1
                  std::cout << "[DEBUG] Hypothesis PJ fit failed (line " << __LINE__ << ")" << std::endl;
                  #endif
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
                delete JpsiConstraint1;
                delete JpsiConstraint2;
              }

              // ========== 假设XJ: mumupipi(X3872) + mumup(Jpsi) ==========
              // 约束: (mu1,mu2,pi1,pi2)→X(3872) 和 (mu3,mu4)→J/psi
              {
                vector<KinematicConstraint*> XJConstraints;
                vector<int> X3872_particleList;
                X3872_particleList.push_back(0); // pi1
                X3872_particleList.push_back(1); // pi2
                X3872_particleList.push_back(2); // mu1
                X3872_particleList.push_back(3); // mu2
                KinematicConstraint* X3872Constraint = new MassKinematicConstraint(X3872_MASS_NOMINAL, X3872_particleList);
                KinematicConstraint* JpsiConstraint2 = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 4, 5);
                XJConstraints.push_back(X3872Constraint);
                XJConstraints.push_back(JpsiConstraint2);

                KinematicConstrainedVertexFitter XJ_fitter;
                RefCountedKinematicTree XJ_VertexFitTree;
                bool XJ_Error = false;
                try {
                  XJ_VertexFitTree = XJ_fitter.fit(X_Particles, XJConstraints);
                } catch (...) {
                  XJ_Error = true;
                }
                if (XJ_Error || !(XJ_VertexFitTree->isValid())) {
                  #if DEBUG == 1
                  std::cout << "[DEBUG] Hypothesis XJ fit failed (line " << __LINE__ << ")" << std::endl;
                  #endif
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
                delete X3872Constraint;
                delete JpsiConstraint2;
              }

              // ========== 假设PP: mumupipi(Jpsi) + mumup(Psi2S) ==========
              // 约束: (mu1,mu2)→J/psi 和 (mu3,mu4)→Psi(2S)
              {
                vector<KinematicConstraint*> PPConstraints;
                KinematicConstraint* JpsiConstraint1 = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 2, 3);
                KinematicConstraint* Psi2SConstraint2 = new MassKinematicConstraint(PSI2S_MASS_NOMINAL, 4, 5);
                PPConstraints.push_back(JpsiConstraint1);
                PPConstraints.push_back(Psi2SConstraint2);

                KinematicConstrainedVertexFitter PP_fitter;
                RefCountedKinematicTree PP_VertexFitTree;
                bool PP_Error = false;
                try {
                  PP_VertexFitTree = PP_fitter.fit(X_Particles, PPConstraints);
                } catch (...) {
                  PP_Error = true;
                }
                if (PP_Error || !(PP_VertexFitTree->isValid())) {
                  #if DEBUG == 1
                  std::cout << "[DEBUG] Hypothesis PP fit failed (line " << __LINE__ << ")" << std::endl;
                  #endif
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
                delete JpsiConstraint1;
                delete Psi2SConstraint2;
              }

              Jpsi1_mass = Jpsi1_vFit_noMC->currentState().mass();
              Jpsi1_VtxProb = Jpsi1_vtxprob;
              Jpsi1_px = mymumupara.momentum().x();
              Jpsi1_py = mymumupara.momentum().y();
              Jpsi1_pz = mymumupara.momentum().z();
              if (Jpsi1_vFit_noMC->currentState()
                      .kinematicParametersError()
                      .matrix()(6, 6) > 0) {
                Jpsi1_massErr = sqrt(Jpsi1_vFit_noMC->currentState()
                                                    .kinematicParametersError()
                                                    .matrix()(6, 6));
              } else {
                Jpsi1_massErr = -9;
              }

              ROOT::Math::PxPyPzMVector Jpsi1_vec(Jpsi1_px, Jpsi1_py, Jpsi1_pz, Jpsi1_mass);
              Jpsi1_pt = Jpsi1_vec.Pt();
              Jpsi1_absEta = fabs(Jpsi1_vec.Eta());

              Jpsi2_mass = Jpsi2_vFit_noMC->currentState().mass();
              Jpsi2_VtxProb = Jpsi2_vtxprob;
              Jpsi2_px = mymumupara2.momentum().x();
              Jpsi2_py = mymumupara2.momentum().y();
              Jpsi2_pz = mymumupara2.momentum().z();
              if (Jpsi2_vFit_noMC->currentState()
                      .kinematicParametersError()
                      .matrix()(6, 6) > 0) {
                Jpsi2_massErr = sqrt(Jpsi2_vFit_noMC->currentState()
                                                    .kinematicParametersError()
                                                    .matrix()(6, 6));
              } else {
                Jpsi2_massErr = -9;
              }

              ROOT::Math::PxPyPzMVector Jpsi2_vec(Jpsi2_px, Jpsi2_py, Jpsi2_pz, Jpsi2_mass);
              Jpsi2_pt = Jpsi2_vec.Pt();
              Jpsi2_absEta = fabs(Jpsi2_vec.Eta());
              
              // Get filter match results for each muon
              mu1_hasFilterMatch = (mu1_matched) ? 1 : 0;
              mu2_hasFilterMatch = (mu2_matched) ? 1 : 0;
              mu3_hasFilterMatch = (mu3_matched) ? 1 : 0;
              mu4_hasFilterMatch = (mu4_matched) ? 1 : 0;
              
              mu1_px = iMuon1->px();
              mu1_py = iMuon1->py();
              mu1_pz = iMuon1->pz();
              mu1_pt = iMuon1->pt();
              mu1_absEta = fabs(iMuon1->eta());
              mu1_trackIso = iMuon1->trackIso();
              mu1_d0BS = iMuon1->dB(pat::Muon::BS2D);
              mu1_d0EBS = iMuon1->edB(pat::Muon::BS2D);
              mu1_d3dBS = iMuon1->dB(pat::Muon::BS3D);
              mu1_d3dEBS = iMuon1->edB(pat::Muon::BS3D);
              mu1_d0PV = iMuon1->dB(pat::Muon::PV2D);
              mu1_d0EPV = iMuon1->edB(pat::Muon::PV2D);
              mu1_dzPV = iMuon1->dB(pat::Muon::PVDZ);
              mu1_dzEPV = iMuon1->edB(pat::Muon::PVDZ);
              mu1_charge = iMuon1->charge();

              mu2_px = iMuon2->px();
              mu2_py = iMuon2->py();
              mu2_pz = iMuon2->pz();
              mu2_pt = iMuon2->pt();
              mu2_absEta = fabs(iMuon2->eta());
              mu2_trackIso = iMuon2->trackIso();
              mu2_d0BS = iMuon2->dB(pat::Muon::BS2D);
              mu2_d0EBS = iMuon2->edB(pat::Muon::BS2D);
              mu2_d3dBS = iMuon2->dB(pat::Muon::BS3D);
              mu2_d3dEBS = iMuon2->edB(pat::Muon::BS3D);
              mu2_d0PV = iMuon2->dB(pat::Muon::PV2D);
              mu2_d0EPV = iMuon2->edB(pat::Muon::PV2D);
              mu2_dzPV = iMuon2->dB(pat::Muon::PVDZ);
              mu2_dzEPV = iMuon2->edB(pat::Muon::PVDZ);
              mu2_charge = iMuon2->charge();

              mu3_px = iMuon3->px();
              mu3_py = iMuon3->py();
              mu3_pz = iMuon3->pz();
              mu3_pt = iMuon3->pt();
              mu3_absEta = fabs(iMuon3->eta());
              mu3_trackIso = iMuon3->trackIso();
              mu3_d0BS = iMuon3->dB(pat::Muon::BS2D);
              mu3_d0EBS = iMuon3->edB(pat::Muon::BS2D);
              mu3_d3dBS = iMuon3->dB(pat::Muon::BS3D);
              mu3_d3dEBS = iMuon3->edB(pat::Muon::BS3D);
              mu3_d0PV = iMuon3->dB(pat::Muon::PV2D);
              mu3_d0EPV = iMuon3->edB(pat::Muon::PV2D);
              mu3_dzPV = iMuon3->dB(pat::Muon::PVDZ);
              mu3_dzEPV = iMuon3->edB(pat::Muon::PVDZ);
              mu3_charge = iMuon3->charge();

              mu4_px = iMuon4->px();
              mu4_py = iMuon4->py();
              mu4_pz = iMuon4->pz();
              mu4_pt = iMuon4->pt();
              mu4_absEta = fabs(iMuon4->eta());
              mu4_trackIso = iMuon4->trackIso();
              mu4_d0BS = iMuon4->dB(pat::Muon::BS2D);
              mu4_d0EBS = iMuon4->edB(pat::Muon::BS2D);
              mu4_d3dBS = iMuon4->dB(pat::Muon::BS3D);
              mu4_d3dEBS = iMuon4->edB(pat::Muon::BS3D);
              mu4_d0PV = iMuon4->dB(pat::Muon::PV2D);
              mu4_d0EPV = iMuon4->edB(pat::Muon::PV2D);
              mu4_dzPV = iMuon4->dB(pat::Muon::PVDZ);
              mu4_dzEPV = iMuon4->edB(pat::Muon::PVDZ);
              mu4_charge = iMuon4->charge();

              // Count muon IDs for all 4 muons
              nLooseMuons = 0;
              nTightMuons = 0;
              nSoftMuons = 0;
              nMediumMuons = 0;

              if (iMuon1->isLooseMuon()) nLooseMuons++;
              if (iMuon1->isTightMuon(thePrimaryV)) nTightMuons++;
              if (iMuon1->isSoftMuon(thePrimaryV)) nSoftMuons++;
              if (iMuon1->isMediumMuon()) nMediumMuons++;

              if (iMuon2->isLooseMuon()) nLooseMuons++;
              if (iMuon2->isTightMuon(thePrimaryV)) nTightMuons++;
              if (iMuon2->isSoftMuon(thePrimaryV)) nSoftMuons++;
              if (iMuon2->isMediumMuon()) nMediumMuons++;

              if (iMuon3->isLooseMuon()) nLooseMuons++;
              if (iMuon3->isTightMuon(thePrimaryV)) nTightMuons++;
              if (iMuon3->isSoftMuon(thePrimaryV)) nSoftMuons++;
              if (iMuon3->isMediumMuon()) nMediumMuons++;

              if (iMuon4->isLooseMuon()) nLooseMuons++;
              if (iMuon4->isTightMuon(thePrimaryV)) nTightMuons++;
              if (iMuon4->isSoftMuon(thePrimaryV)) nSoftMuons++;
              if (iMuon4->isMediumMuon()) nMediumMuons++;

              pi1_px = X_pi1_KP.momentum().x();
              pi1_py = X_pi1_KP.momentum().y();
              pi1_pz = X_pi1_KP.momentum().z();
              ROOT::Math::PxPyPzMVector pi1_vec(pi1_px, pi1_py, pi1_pz, PI_MASS);
              pi1_pt = pi1_vec.Pt();
              pi1_absEta = fabs(pi1_vec.Eta());

              pi2_px = X_pi2_KP.momentum().x();
              pi2_py = X_pi2_KP.momentum().y();
              pi2_pz = X_pi2_KP.momentum().z();
              ROOT::Math::PxPyPzMVector pi2_vec(pi2_px, pi2_py, pi2_pz, PI_MASS);
              pi2_pt = pi2_vec.Pt();
              pi2_absEta = fabs(pi2_vec.Eta());

              Psi2S_mass = (Jpsi1_vec + pi1_vec + pi2_vec).M() - Jpsi1_vec.M() + jpsi_nominal_mass;

              ROOT::Math::PxPyPzMVector mu1_vec(mu1_px, mu1_py, mu1_pz, MU_MASS);
              ROOT::Math::PxPyPzMVector mu2_vec(mu2_px, mu2_py, mu2_pz, MU_MASS);
              ROOT::Math::PxPyPzMVector mu3_vec(mu3_px, mu3_py, mu3_pz, MU_MASS);
              ROOT::Math::PxPyPzMVector mu4_vec(mu4_px, mu4_py, mu4_pz, MU_MASS);

              dR_mu1_mu2 = ROOT::Math::VectorUtil::DeltaR(mu1_vec, mu2_vec);
              dR_mu3_mu4 = ROOT::Math::VectorUtil::DeltaR(mu3_vec, mu4_vec);
              dR_pi1_pi2 = ROOT::Math::VectorUtil::DeltaR(pi1_vec, pi2_vec);
              dR_Jpsi1_X6900 = ROOT::Math::VectorUtil::DeltaR(Jpsi1_vec, X6900_vec);
              dR_Jpsi2_X6900 = ROOT::Math::VectorUtil::DeltaR(Jpsi2_vec, X6900_vec);
              dR_X6900_pi1 = ROOT::Math::VectorUtil::DeltaR(X6900_vec, pi1_vec);
              dR_X6900_pi2 = ROOT::Math::VectorUtil::DeltaR(X6900_vec, pi2_vec);
              dR_X6900_mu1 = ROOT::Math::VectorUtil::DeltaR(X6900_vec, mu1_vec);
              dR_X6900_mu2 = ROOT::Math::VectorUtil::DeltaR(X6900_vec, mu2_vec);
              dR_X6900_mu3 = ROOT::Math::VectorUtil::DeltaR(X6900_vec, mu3_vec);
              dR_X6900_mu4 = ROOT::Math::VectorUtil::DeltaR(X6900_vec, mu4_vec);
              dR_Psi2S_X6900 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, X6900_vec);
              dR_Psi2S_Jpsi1 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, Jpsi1_vec);
              dR_Psi2S_Jpsi2 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, Jpsi2_vec);
              dR_Psi2S_pi1 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, pi1_vec);
              dR_Psi2S_pi2 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, pi2_vec);

              X_One_Tree_->Fill();
              #if DEBUG == 1
              std::cout << "[DEBUG] Filled an Entry";
              #endif

            } // itrack2
          } // itrack1
        } // mu4_loop
      } // mu3_loop
    } // for (std::vector < pat::Muon >::const_iterator iMuon2 = iMuon1 + 1;
  } // for (std::vector < pat::Muon >::const_iterator iMuon1 =
    // thePATMuonHandle->begin();

#if DEBUG == 1
  auto end_fit = std::chrono::high_resolution_clock::now();
  auto duration_fit = std::chrono::duration_cast<std::chrono::milliseconds>(end_fit - start_fit).count();
  std::cout << "[DEBUG] MuonTrackFitting 运行时间: " << duration_fit << " ms" << std::endl;
#endif

} // analyze
//
// ------------ method called once each job just before starting event loop
// ------------
void MultiLepPAT::beginRun(edm::Run const &iRun,
                           edm::EventSetup const &iSetup) {}

void MultiLepPAT::beginJob() {
  edm::Service<TFileService> fs;
  X_One_Tree_ = fs->make<TTree>("X_data", "X(3872) Data");


  X_One_Tree_->Branch("evtNum", &evtNum, "evtNum/i");
  X_One_Tree_->Branch("runNum", &runNum, "runNum/i");
  X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
  X_One_Tree_->Branch("nGoodPrimVtx", &nGoodPrimVtx, "nGoodPrimVtx/i");


  // Hypothesis PJ: mumupipi (J/psi constrained) + mumu (J/psi constrained)
  X_One_Tree_->Branch("X_PJ_mass", &X_PJ_mass, "X_PJ_mass/F");
  X_One_Tree_->Branch("X_PJ_VtxProb", &X_PJ_VtxProb, "X_PJ_VtxProb/F");
  X_One_Tree_->Branch("X_PJ_massErr", &X_PJ_massErr, "X_PJ_massErr/F");
  X_One_Tree_->Branch("X_PJ_pt", &X_PJ_pt, "X_PJ_pt/F");
  X_One_Tree_->Branch("X_PJ_pz", &X_PJ_pz, "X_PJ_pz/F");
  X_One_Tree_->Branch("X_PJ_absEta", &X_PJ_absEta, "X_PJ_absEta/F");
  X_One_Tree_->Branch("X_PJ_px", &X_PJ_px, "X_PJ_px/F");
  X_One_Tree_->Branch("X_PJ_py", &X_PJ_py, "X_PJ_py/F");

  // Hypothesis XJ: mumupipi (X(3872) constrained) + mumu (J/psi constrained)
  X_One_Tree_->Branch("X_XJ_mass", &X_XJ_mass, "X_XJ_mass/F");
  X_One_Tree_->Branch("X_XJ_VtxProb", &X_XJ_VtxProb, "X_XJ_VtxProb/F");
  X_One_Tree_->Branch("X_XJ_massErr", &X_XJ_massErr, "X_XJ_massErr/F");
  X_One_Tree_->Branch("X_XJ_pt", &X_XJ_pt, "X_XJ_pt/F");
  X_One_Tree_->Branch("X_XJ_pz", &X_XJ_pz, "X_XJ_pz/F");
  X_One_Tree_->Branch("X_XJ_absEta", &X_XJ_absEta, "X_XJ_absEta/F");
  X_One_Tree_->Branch("X_XJ_px", &X_XJ_px, "X_XJ_px/F");
  X_One_Tree_->Branch("X_XJ_py", &X_XJ_py, "X_XJ_py/F");

  // Hypothesis PP: mumupipi (J/psi constrained) + mumu (psi(2S) constrained)
  X_One_Tree_->Branch("X_PP_mass", &X_PP_mass, "X_PP_mass/F");
  X_One_Tree_->Branch("X_PP_VtxProb", &X_PP_VtxProb, "X_PP_VtxProb/F");
  X_One_Tree_->Branch("X_PP_massErr", &X_PP_massErr, "X_PP_massErr/F");
  X_One_Tree_->Branch("X_PP_pt", &X_PP_pt, "X_PP_pt/F");
  X_One_Tree_->Branch("X_PP_pz", &X_PP_pz, "X_PP_pz/F");
  X_One_Tree_->Branch("X_PP_absEta", &X_PP_absEta, "X_PP_absEta/F");
  X_One_Tree_->Branch("X_PP_px", &X_PP_px, "X_PP_px/F");
  X_One_Tree_->Branch("X_PP_py", &X_PP_py, "X_PP_py/F");

  // Psi(2S) with J/psi mass constraint: mumupipi system
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
  
  // Clear filter match results vector
  muonFilterMatches.clear();
  
  // Reset event identification variables
  runNum = 0;
  evtNum = 0;
  lumiNum = 0;
  nGoodPrimVtx = 0;
  
  // Reset 3 hypothesis X candidate variables
  // Hypothesis PJ: mumupipi (J/psi constrained) + mumu (J/psi constrained)
  X_PJ_mass = -999.0;
  X_PJ_VtxProb = -999.0;
  X_PJ_massErr = -999.0;
  X_PJ_pt = -999.0;
  X_PJ_pz = -999.0;
  X_PJ_absEta = -999.0;
  X_PJ_px = -999.0;
  X_PJ_py = -999.0;

  // Hypothesis XJ: mumupipi (X(3872) constrained) + mumu (J/psi constrained)
  X_XJ_mass = -999.0;
  X_XJ_VtxProb = -999.0;
  X_XJ_massErr = -999.0;
  X_XJ_pt = -999.0;
  X_XJ_pz = -999.0;
  X_XJ_absEta = -999.0;
  X_XJ_px = -999.0;
  X_XJ_py = -999.0;

  // Hypothesis PP: mumupipi (J/psi constrained) + mumu (psi(2S) constrained)
  X_PP_mass = -999.0;
  X_PP_VtxProb = -999.0;
  X_PP_massErr = -999.0;
  X_PP_pt = -999.0;
  X_PP_pz = -999.0;
  X_PP_absEta = -999.0;
  X_PP_px = -999.0;
  X_PP_py = -999.0;
  
  // Reset Psi(2S) candidate variables with J/psi mass constraint
  Psi2S_mass = -999.0;
  Psi2S_VtxProb = -999.0;
  Psi2S_massErr = -999.0;
  Psi2S_pt = -999.0;
  Psi2S_pz = -999.0;
  Psi2S_absEta = -999.0;
  Psi2S_px = -999.0;
  Psi2S_py = -999.0;
  
  // Reset J/psi1 candidate variables
  Jpsi1_mass = -999.0;
  Jpsi1_VtxProb = -999.0;
  Jpsi1_massErr = -999.0;
  Jpsi1_pt = -999.0;
  Jpsi1_pz = -999.0;
  Jpsi1_absEta = -999.0;
  Jpsi1_px = -999.0;
  Jpsi1_py = -999.0;
  
  // Reset J/psi2 candidate variables
  Jpsi2_mass = -999.0;
  Jpsi2_VtxProb = -999.0;
  Jpsi2_massErr = -999.0;
  Jpsi2_pt = -999.0;
  Jpsi2_pz = -999.0;
  Jpsi2_absEta = -999.0;
  Jpsi2_px = -999.0;
  Jpsi2_py = -999.0;
  
  // Reset muon1 variables
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
  
  // Reset muon2 variables
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
  
  // Reset muon3 variables
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
  
  // Reset muon4 variables
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
  
  // Reset muon ID count variables
  nLooseMuons = 0;
  nTightMuons = 0;
  nSoftMuons = 0;
  nMediumMuons = 0;
  
  // Reset muon filter match variables
  mu1_hasFilterMatch = 0;
  mu2_hasFilterMatch = 0;
  mu3_hasFilterMatch = 0;
  mu4_hasFilterMatch = 0;
  
  // Reset pion variables
  pi1_pt = -999.0;
  pi1_pz = -999.0;
  pi1_absEta = -999.0;
  pi1_px = -999.0;
  pi1_py = -999.0;
  pi2_pt = -999.0;
  pi2_pz = -999.0;
  pi2_absEta = -999.0;
  pi2_px = -999.0;
  pi2_py = -999.0;
  
  // Reset delta R variables
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

#if DEBUG == 2
  std::cout << std::endl;
  std::cout << "================================================================" << std::endl;
  std::cout << "           Enhanced Debug Statistics (DEBUG=2)" << std::endl;
  std::cout << "================================================================" << std::endl;
  std::cout << std::endl;
  std::cout << "Total continue statements triggered: " << total_continue_ << std::endl;
  std::cout << "Total return statements triggered: " << total_return_ << std::endl;
  std::cout << std::endl;

  if (!continue_counts_.empty()) {
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "continue statements statistics (sorted by count descending):" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << std::setw(10) << "Line" << std::setw(15) << "Count" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    std::vector<std::pair<int, unsigned long long>> continue_list;
    for (auto& pair : continue_counts_) {
      continue_list.push_back(pair);
    }

    std::sort(continue_list.begin(), continue_list.end(),
              [](const std::pair<int, unsigned long long>& a,
                 const std::pair<int, unsigned long long>& b) {
                return a.second > b.second;
              });

    for (auto& pair : continue_list) {
      std::cout << std::setw(10) << pair.first << std::setw(15) << pair.second << std::endl;
    }
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  if (!return_counts_.empty()) {
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "return statements statistics (sorted by count descending):" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << std::setw(10) << "Line" << std::setw(15) << "Count" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    std::vector<std::pair<int, unsigned long long>> return_list;
    for (auto& pair : return_counts_) {
      return_list.push_back(pair);
    }

    std::sort(return_list.begin(), return_list.end(),
              [](const std::pair<int, unsigned long long>& a,
                 const std::pair<int, unsigned long long>& b) {
                return a.second > b.second;
              });

    for (auto& pair : return_list) {
      std::cout << std::setw(10) << pair.first << std::setw(15) << pair.second << std::endl;
    }
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  std::cout << "================================================================" << std::endl;
  std::cout << std::endl;
#endif
}

// define this as a plug-in
DEFINE_FWK_MODULE(MultiLepPAT);
