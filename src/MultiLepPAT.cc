/******************************************************************************
 *  [File] 
 *      MultiLepPAT.cc 
 *  [Class]      
 *      MultiLepPAT 
 *  [Directory] 
 *      TPS-Onia2MuMu/src/MultiLepPAT.cc
 *  [Description]
 *      Make rootTuple for quarkonia+quarkonia+meson reconstruction.
 *      Supports J/psi+J/psi+phi, J/psi+J/psi+Upsilon, J/psi+Upsilon+phi
 *      via dynamic config switching (AnalysisMode parameter).
 *  [Implementation]
 *      Refactored from monolithic analyze() into modular helper methods.
 *      All selection cuts externalized to config via StringCutObjectSelector.
 *  [Note]
 *      20240704 [Eric Wang] - Initial version
 *      20241029 [Eric Wang and Shi Zhenpeng] - J/psi+Upsilon+Phi branch
 *      20260306 [Eric Wang - Refactor] - Full modularization
 *          PDG ID: Mu->13, K->321, Pi->211, Jpsi->443, Ups->553, Phi->333
******************************************************************************/

// system include files
#include "TLorentzVector.h"
// user include files
#include "../interface/MultiLepPAT.h"
#include "../interface/VertexReProducer.h"
#include <memory>
#include <regex>
#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <unordered_set>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/DeepCopyPointer.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParam.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParamRcd.h"

#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "../data/TMVAClassification_BDT.class.C"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

typedef math::Error<3>::type CovarianceMatrix;

namespace {

const reco::Track *bestMuonTrackForGenMatch(const pat::Muon &muon) {
  const reco::TrackRef bestTrack = muon.muonBestTrack();
  if (bestTrack.isNonnull() && bestTrack.isAvailable()) {
    return bestTrack.get();
  }

  const reco::TrackRef innerTrack = muon.innerTrack();
  if (innerTrack.isNonnull() && innerTrack.isAvailable()) {
    return innerTrack.get();
  }

  return nullptr;
}

double recoGenMuonChi2(const pat::Muon &muon, const reco::GenParticle &gen) {
  const reco::Track *track = bestMuonTrackForGenMatch(muon);
  if (track == nullptr) {
    return std::numeric_limits<double>::infinity();
  }

  const double sigmaPt = std::max(track->ptError(), 1e-6);
  const double sigmaEta = std::max(track->etaError(), 1e-6);
  const double sigmaPhi = std::max(track->phiError(), 1e-6);

  double dphi = muon.phi() - gen.phi();
  while (dphi > M_PI)
    dphi -= 2. * M_PI;
  while (dphi <= -M_PI)
    dphi += 2. * M_PI;

  return std::pow((muon.pt() - gen.pt()) / sigmaPt, 2) +
         std::pow((muon.eta() - gen.eta()) / sigmaEta, 2) +
         std::pow(dphi / sigmaPhi, 2);
}

}  // namespace

/*****************************************************************************
 * Constructor
 *****************************************************************************/
MultiLepPAT::MultiLepPAT(const edm::ParameterSet &iConfig)
    : hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>(
          "HLTriggerResults", edm::InputTag("TriggerResults::HLT"))),
      inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>(
          "inputGEN", edm::InputTag("prunedGenParticles"))),
      muonLabel_(iConfig.getUntrackedParameter<edm::InputTag>(
          "MuonLabel", edm::InputTag("slimmedMuons"))),
      trackLabel_(iConfig.getUntrackedParameter<edm::InputTag>(
          "TrackLabel", edm::InputTag("packedPFCandidates"))),
      magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      theTTBuilderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(
          edm::ESInputTag("", "TransientTrackBuilder"))),
      doMC(iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", false)),
      doJPsiMassCost(iConfig.getUntrackedParameter<bool>("DoJPsiMassConstraint", false)),
      Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output", false)),
      // Analysis mode: "JpsiJpsiPhi", "JpsiJpsiUps", "JpsiUpsPhi"
      analysisModeName_(iConfig.getUntrackedParameter<std::string>("AnalysisMode", "JpsiJpsiPhi")),
      // StringCutObjectSelector for muon and track
      muonSelectionStr_(iConfig.getUntrackedParameter<std::string>(
          "MuonSelection", "pt > 2.5 && abs(eta) < 2.4")),
      muonSelector_(muonSelectionStr_),
      trackSelectionStr_(iConfig.getUntrackedParameter<std::string>(
          "TrackSelection",
          "pt > 2.0 && abs(eta) < 2.5 && numberOfHits > 4")),
      trackSelector_(trackSelectionStr_),
      // Primary vertex cuts
      pvNdofMin_(iConfig.getUntrackedParameter<int>("PVNdofMin", 5)),
      pvMaxAbsZ_(iConfig.getUntrackedParameter<double>("PVMaxAbsZ", 24.0)),
      pvMaxRho_(iConfig.getUntrackedParameter<double>("PVMaxRho", 2.0)),
      // Muon pair mass windows
      JpsiMassMin_(iConfig.getUntrackedParameter<double>("JpsiMassMin", 1.0)),
      JpsiMassMax_(iConfig.getUntrackedParameter<double>("JpsiMassMax", 4.0)),
      UpsMassMin_(iConfig.getUntrackedParameter<double>("UpsMassMin", 8.0)),
      UpsMassMax_(iConfig.getUntrackedParameter<double>("UpsMassMax", 12.0)),
      // Track pair mass windows
      PhiMassMin_(iConfig.getUntrackedParameter<double>("PhiMassMin", 0.8)),
      PhiMassMax_(iConfig.getUntrackedParameter<double>("PhiMassMax", 1.2)),
      // Track kinematics
      trackPtMin_(iConfig.getUntrackedParameter<double>("TrackPtMin", 2.0)),
      trackDRMax_(iConfig.getUntrackedParameter<double>("TrackDRMax", 0.7)),
      // Vertex prob cuts
      OniaDecayVtxProbCut_(iConfig.getUntrackedParameter<double>("OniaDecayVtxProbCut", 0.001)),
      PriVtxProbCut_(iConfig.getUntrackedParameter<double>("PriVtxProbCut", 0.0)),
      // Per-resonance candidate pT/eta pre-cuts
      jpsiCandPtMin_(iConfig.getUntrackedParameter<double>("JpsiCandPtMin", 0.0)),
      jpsiCandEtaMax_(iConfig.getUntrackedParameter<double>("JpsiCandEtaMax", 999.0)),
      upsCandPtMin_(iConfig.getUntrackedParameter<double>("UpsCandPtMin", 0.0)),
      upsCandEtaMax_(iConfig.getUntrackedParameter<double>("UpsCandEtaMax", 999.0)),
      phiCandPtMin_(iConfig.getUntrackedParameter<double>("PhiCandPtMin", 0.0)),
      phiCandEtaMax_(iConfig.getUntrackedParameter<double>("PhiCandEtaMax", 999.0)),
      // PV selection mode
      pvSelectionMode_(iConfig.getUntrackedParameter<std::string>("PVSelectionMode", "firstVertex")),
      // Minimum fromPV for tracks
      minTrackFromPV_(iConfig.getUntrackedParameter<int>("MinTrackFromPV", 1)),
      // Minimum muon count
      minMuonCount_(iConfig.getUntrackedParameter<unsigned int>("MinMuonCount", 4)),
      // Muon matching
      MuMatchTrkMomentumRelDiffThr_c(iConfig.getUntrackedParameter<double>(
          "MuMatchTrkMomentumRelDiffThr", 0.5)),
      // Muon-track matching method selection
      muTrkMatchMethod_(iConfig.getUntrackedParameter<std::string>("MuTrkMatchMethod", "first")),
      muTrkMatchDebug_(iConfig.getUntrackedParameter<bool>("MuTrkMatchDebug", false)),
      muTrkMatchDzMax_(iConfig.getUntrackedParameter<double>("MuTrkMatchDzMax", 1.0)),
      muTrkMatchSigmaThr_(iConfig.getUntrackedParameter<double>("MuTrkMatchSigmaThr", 5.0)),
      // Store all primary vertices and muon quantities
      storeAllPVs_(iConfig.getUntrackedParameter<bool>("StoreAllPVs", true)),
      storeMuonMomentumErrors_(iConfig.getUntrackedParameter<bool>("StoreMuonMomentumErrors", true)),
      storeMuonPVAssoc_(iConfig.getUntrackedParameter<bool>("StoreMuonPVAssoc", true)),
      recoGenMuonMatchChi2Max_(iConfig.getUntrackedParameter<double>("RecoGenMuonMatchChi2Max", 25.0)),
      // Final fitted mass window check
      checkFinalMass_(iConfig.getUntrackedParameter<bool>("CheckFinalMass", true)),
      // Trigger info
      resolveAmbiguity_(iConfig.getUntrackedParameter<bool>("resolvePileUpAmbiguity", true)),
      addXlessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addXlessPrimaryVertex", true)),
      TriggersForJpsi_(iConfig.getUntrackedParameter<std::vector<std::string>>("TriggersForJpsi")),
      FiltersForJpsi_(iConfig.getUntrackedParameter<std::vector<std::string>>("FiltersForJpsi")),
      TriggersForUpsilon_(iConfig.getUntrackedParameter<std::vector<std::string>>("TriggersForUpsilon")),
      FiltersForUpsilon_(iConfig.getUntrackedParameter<std::vector<std::string>>("FiltersForUpsilon")),
      X_One_Tree_(nullptr),
      runNum(0), evtNum(0), lumiNum(0), nGoodPrimVtx(0),
      trigRes(nullptr), trigNames(nullptr), L1TT(nullptr), MatchJpsiTrigNames(nullptr), MatchUpsTrigNames(nullptr),
      priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0),
      priVtxChiNorm(0), priVtxChi(0), priVtxCL(0),
      PriVtxXCorrX(nullptr), PriVtxXCorrY(nullptr), PriVtxXCorrZ(nullptr),
      PriVtxXCorrEX(nullptr), PriVtxXCorrEY(nullptr), PriVtxXCorrEZ(nullptr),
      PriVtxXCorrC2(nullptr), PriVtxXCorrCL(nullptr),
      // All primary vertices (not just selected)
      nRecVtx(0),
      RecVtx_x(nullptr), RecVtx_y(nullptr), RecVtx_z(nullptr),
      RecVtx_xErr(nullptr), RecVtx_yErr(nullptr), RecVtx_zErr(nullptr),
      RecVtx_chi2(nullptr), RecVtx_ndof(nullptr), RecVtx_vtxProb(nullptr),
      RecVtx_nTracks(nullptr),
      nMu(0),
      muPx(nullptr), muPy(nullptr), muPz(nullptr),
      muD0(nullptr), muD0E(nullptr), muDz(nullptr),
      muChi2(nullptr), muGlChi2(nullptr), mufHits(nullptr),
      muFirstBarrel(nullptr), muFirstEndCap(nullptr), muDzVtx(nullptr), muDxyVtx(nullptr),
      muNDF(nullptr), muGlNDF(nullptr), muPhits(nullptr), muShits(nullptr),
      muGlMuHits(nullptr), muType(nullptr), muQual(nullptr),
      muTrack(nullptr), muCharge(nullptr), muIsoratio(nullptr),
      muIsGoodLooseMuon(nullptr), muIsGoodLooseMuonNew(nullptr),
      muIsGoodSoftMuonNewIlse(nullptr), muIsGoodSoftMuonNewIlseMod(nullptr),
      muIsGlobalMuon(nullptr), muIsGoodTightMuon(nullptr),
      muIsJpsiTrigMatch(nullptr), muIsUpsTrigMatch(nullptr), munMatchedSeg(nullptr),
      muIsJpsiFilterMatch(nullptr), muIsUpsFilterMatch(nullptr),
      muIsPatLooseMuon(nullptr), muIsPatTightMuon(nullptr),
      muIsPatSoftMuon(nullptr), muIsPatMediumMuon(nullptr),
      muFromPV(nullptr), muPVAssocQuality(nullptr),
      // Muon momentum errors
      muPxErr(nullptr), muPyErr(nullptr), muPzErr(nullptr), muPtErr(nullptr),
      // Muon PV association (surplus from sourceCandidatePtr)
      muVertexId(nullptr), muDzAssocPV(nullptr), muDxyAssocPV(nullptr),
      muFromPVAssocPV(nullptr), muPdgId(nullptr),
      muGenMatchIdx(nullptr), muGenMatchSource(nullptr), muGenMatchChi2(nullptr),
      muMVAMuonID(nullptr), musegmentCompatibility(nullptr),
      mupulldXdZ_pos_noArb(nullptr), mupulldYdZ_pos_noArb(nullptr),
      mupulldXdZ_pos_ArbDef(nullptr), mupulldYdZ_pos_ArbDef(nullptr),
      mupulldXdZ_pos_ArbST(nullptr), mupulldYdZ_pos_ArbST(nullptr),
      mupulldXdZ_pos_noArb_any(nullptr), mupulldYdZ_pos_noArb_any(nullptr),
      // Muon-track matching debug (only if muTrkMatchDebug_)
      muMatch_nCandidates(nullptr), muMatch_bestRelDiff(nullptr),
      muMatch_bestDz(nullptr), muMatch_methodUsed(nullptr),
      Jpsi_1_mu_1_Idx(nullptr), Jpsi_1_mu_2_Idx(nullptr),
      Jpsi_2_mu_1_Idx(nullptr), Jpsi_2_mu_2_Idx(nullptr),
      Phi_K_1_Idx(nullptr), Phi_K_2_Idx(nullptr),
      Jpsi_1_mass(nullptr), Jpsi_1_massErr(nullptr), Jpsi_1_massDiff(nullptr),
      Jpsi_2_mass(nullptr), Jpsi_2_massErr(nullptr), Jpsi_2_massDiff(nullptr),
      Phi_mass(nullptr), Phi_massErr(nullptr), Phi_massDiff(nullptr),
      Jpsi_1_ctau(nullptr), Jpsi_1_ctauErr(nullptr), Jpsi_1_Chi2(nullptr),
      Jpsi_1_ndof(nullptr), Jpsi_1_VtxProb(nullptr),
      Jpsi_2_ctau(nullptr), Jpsi_2_ctauErr(nullptr), Jpsi_2_Chi2(nullptr),
      Jpsi_2_ndof(nullptr), Jpsi_2_VtxProb(nullptr),
      Phi_ctau(nullptr), Phi_ctauErr(nullptr), Phi_Chi2(nullptr),
      Phi_ndof(nullptr), Phi_VtxProb(nullptr),
      Jpsi_1_phi(nullptr), Jpsi_1_eta(nullptr), Jpsi_1_pt(nullptr),
      Jpsi_2_phi(nullptr), Jpsi_2_eta(nullptr), Jpsi_2_pt(nullptr),
      Phi_phi(nullptr), Phi_eta(nullptr), Phi_pt(nullptr),
      Jpsi_1_px(nullptr), Jpsi_1_py(nullptr), Jpsi_1_pz(nullptr),
      Jpsi_2_px(nullptr), Jpsi_2_py(nullptr), Jpsi_2_pz(nullptr),
      Phi_px(nullptr), Phi_py(nullptr), Phi_pz(nullptr),
      // Momentum errors
      Jpsi_1_pxErr(nullptr), Jpsi_1_pyErr(nullptr), Jpsi_1_pzErr(nullptr), Jpsi_1_ptErr(nullptr),
      Jpsi_2_pxErr(nullptr), Jpsi_2_pyErr(nullptr), Jpsi_2_pzErr(nullptr), Jpsi_2_ptErr(nullptr),
      Phi_pxErr(nullptr), Phi_pyErr(nullptr), Phi_pzErr(nullptr), Phi_ptErr(nullptr),
      Pri_mass(nullptr), Pri_massErr(nullptr),
      Pri_ctau(nullptr), Pri_ctauErr(nullptr), Pri_Chi2(nullptr),
      Pri_ndof(nullptr), Pri_VtxProb(nullptr),
      Pri_px(nullptr), Pri_py(nullptr), Pri_pz(nullptr),
      Pri_phi(nullptr), Pri_eta(nullptr), Pri_pt(nullptr),
      Pri_pxErr(nullptr), Pri_pyErr(nullptr), Pri_pzErr(nullptr), Pri_ptErr(nullptr),
      Phi_K_1_px(nullptr), Phi_K_1_py(nullptr), Phi_K_1_pz(nullptr),
      Phi_K_2_px(nullptr), Phi_K_2_py(nullptr), Phi_K_2_pz(nullptr),
      Phi_K_1_eta(nullptr), Phi_K_1_phi(nullptr), Phi_K_1_pt(nullptr),
      Phi_K_2_eta(nullptr), Phi_K_2_phi(nullptr), Phi_K_2_pt(nullptr),
      Phi_K_1_fromPV(nullptr), Phi_K_2_fromPV(nullptr),
      Phi_K_1_pvAssocQuality(nullptr), Phi_K_2_pvAssocQuality(nullptr),
      // Upsilon branches
      Ups_mu_1_Idx(nullptr), Ups_mu_2_Idx(nullptr),
      Ups_mass(nullptr), Ups_massErr(nullptr), Ups_massDiff(nullptr),
      Ups_ctau(nullptr), Ups_ctauErr(nullptr), Ups_Chi2(nullptr),
      Ups_ndof(nullptr), Ups_VtxProb(nullptr),
      Ups_px(nullptr), Ups_py(nullptr), Ups_pz(nullptr),
      Ups_phi(nullptr), Ups_eta(nullptr), Ups_pt(nullptr),
      Ups_pxErr(nullptr), Ups_pyErr(nullptr), Ups_pzErr(nullptr), Ups_ptErr(nullptr),
      // MC gen-level (new)
      MC_GenPart_pdgId(nullptr), MC_GenPart_status(nullptr), MC_GenPart_motherPdgId(nullptr),
      MC_GenPart_handleIndex(nullptr),
      MC_GenPart_px(nullptr), MC_GenPart_py(nullptr), MC_GenPart_pz(nullptr), MC_GenPart_mass(nullptr),
      MC_GenPart_pt(nullptr), MC_GenPart_eta(nullptr), MC_GenPart_phi(nullptr),
      // MC gen-level (legacy)
      MC_X_px(nullptr), MC_X_py(nullptr), MC_X_pz(nullptr), MC_X_mass(nullptr),
      MC_Dau_Jpsipx(nullptr), MC_Dau_Jpsipy(nullptr), MC_Dau_Jpsipz(nullptr), MC_Dau_Jpsimass(nullptr),
      MC_Dau_psi2spx(nullptr), MC_Dau_psi2spy(nullptr), MC_Dau_psi2spz(nullptr), MC_Dau_psi2smass(nullptr),
      MC_Granddau_mu1px(nullptr), MC_Granddau_mu1py(nullptr), MC_Granddau_mu1pz(nullptr),
      MC_Granddau_mu2px(nullptr), MC_Granddau_mu2py(nullptr), MC_Granddau_mu2pz(nullptr),
      MC_Granddau_Jpsipx(nullptr), MC_Granddau_Jpsipy(nullptr), MC_Granddau_Jpsipz(nullptr), MC_Granddau_Jpsimass(nullptr),
      MC_Granddau_pi1px(nullptr), MC_Granddau_pi1py(nullptr), MC_Granddau_pi1pz(nullptr),
      MC_Granddau_pi2px(nullptr), MC_Granddau_pi2py(nullptr), MC_Granddau_pi2pz(nullptr),
      MC_Grandgranddau_mu3px(nullptr), MC_Grandgranddau_mu3py(nullptr), MC_Grandgranddau_mu3pz(nullptr),
      MC_Grandgranddau_mu4px(nullptr), MC_Grandgranddau_mu4py(nullptr), MC_Grandgranddau_mu4pz(nullptr),
      MC_X_chg(nullptr),
      MC_Dau_JpsipdgId(nullptr), MC_Dau_psi2spdgId(nullptr),
      MC_Granddau_mu1pdgId(nullptr), MC_Granddau_mu2pdgId(nullptr),
      MC_Granddau_JpsipdgId(nullptr),
      MC_Granddau_pi1pdgId(nullptr), MC_Granddau_pi2pdgId(nullptr),
      MC_Grandgranddau_mu3pdgId(nullptr), MC_Grandgranddau_mu4pdgId(nullptr),
      Match_mu1px(nullptr), Match_mu1py(nullptr), Match_mu1pz(nullptr),
      Match_mu2px(nullptr), Match_mu2py(nullptr), Match_mu2pz(nullptr),
      Match_mu3px(nullptr), Match_mu3py(nullptr), Match_mu3pz(nullptr),
      Match_mu4px(nullptr), Match_mu4py(nullptr), Match_mu4pz(nullptr),
      Match_pi1px(nullptr), Match_pi1py(nullptr), Match_pi1pz(nullptr),
      Match_pi2px(nullptr), Match_pi2py(nullptr), Match_pi2pz(nullptr)
{
    // Parse analysis mode
    if (analysisModeName_ == "JpsiJpsiPhi") {
        analysisChannel_ = AnalysisChannel::JpsiJpsiPhi;
    } else if (analysisModeName_ == "JpsiJpsiUps") {
        analysisChannel_ = AnalysisChannel::JpsiJpsiUps;
    } else if (analysisModeName_ == "JpsiUpsPhi") {
        analysisChannel_ = AnalysisChannel::JpsiUpsPhi;
    } else {
        throw cms::Exception("Configuration")
            << "Unknown AnalysisMode: " << analysisModeName_
            << ". Must be one of: JpsiJpsiPhi, JpsiJpsiUps, JpsiUpsPhi";
    }
    
    // For JpsiJpsiUps, enforce minimum 6 muons
    if (analysisChannel_ == AnalysisChannel::JpsiJpsiUps) {
        minMuonCount_ = std::max(minMuonCount_, 6u);
    }
    
    // EDM tokens
    gtRecordToken_     = consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag("gtDigis"));
    gtbeamspotToken_   = consumes<BeamSpot>(edm::InputTag("offlineBeamSpot"));
    gtprimaryVtxToken_ = consumes<VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    gtpatmuonToken_    = consumes<edm::View<pat::Muon>>(muonLabel_);
    gttriggerToken_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"));
    trackToken_        = consumes<edm::View<pat::PackedCandidate>>(trackLabel_);
    genParticlesToken_ = consumes<reco::GenParticleCollection>(inputGEN_);
}

MultiLepPAT::~MultiLepPAT() {}

/*****************************************************************************
 * Main analyze() - delegates to modular steps
 *****************************************************************************/
void MultiLepPAT::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    using namespace edm;
    using namespace reco;
    using namespace std;

    runNum  = iEvent.id().run();
    evtNum  = iEvent.id().event();
    lumiNum = iEvent.id().luminosityBlock();

    const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_);

    // Step 1: MC gen-level info
    if (doMC) {
        processMCGenInfo(iEvent);
    }
    
    // Step 2: HLT trigger info
    processHLTInfo(iEvent);

    // Step 3: L1 trigger
    edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
    iEvent.getByToken(gtRecordToken_, gtRecord);
    if (gtRecord.isValid()) {
        const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
        for (unsigned int l1i = 0; l1i != ttWord.size(); ++l1i) {
            L1TT->push_back(ttWord.at(l1i));
        }
    }

    // Step 4: Primary vertex reconstruction
    reconstructPrimaryVertex(iEvent);

    // Step 5: Get muon and track handles
    iEvent.getByToken(gtpatmuonToken_, thePATMuonHandle_);
    iEvent.getByToken(trackToken_, theTrackHandle_);

    // Build non-muon track list
    nonMuonTrack_.clear();
    if (theTrackHandle_.isValid()) {
        for (auto iTrack = theTrackHandle_->begin(); iTrack != theTrackHandle_->end(); ++iTrack) {
            nonMuonTrack_.push_back(iTrack);
        }
    }

    // Step 6: Fill the muon block for all available muons.
    if (thePATMuonHandle_.isValid()) {
        fillMuonBlock(iEvent, iSetup, thePrimaryV_);
    }

    // Step 7: MC matching of tracks
    if (doMC && thePATMuonHandle_.isValid() && theTrackHandle_.isValid()) {
        doMCGenMatching(thePATMuonHandle_, theTrackHandle_);
    }

    // Step 8: Pair muons
    if (!thePATMuonHandle_.isValid() || thePATMuonHandle_->size() < minMuonCount_) {
        // Not enough muons: still fill tree if MC
        if (doMC) {
            X_One_Tree_->Fill();
        }
        clearEventData();
        return;
    }

    pairMuons(thePATMuonHandle_, bFieldHandle);

    // Step 9: Pair tracks into meson candidates
    bool needsTrackPairs = (analysisChannel_ == AnalysisChannel::JpsiJpsiPhi ||
                            analysisChannel_ == AnalysisChannel::JpsiUpsPhi);
    if (needsTrackPairs && !theTrackHandle_.isValid()) {
        if (doMC) {
            X_One_Tree_->Fill();
        }
        clearEventData();
        return;
    }
    if (needsTrackPairs) {
        pairTracks(nonMuonTrack_, bFieldHandle);
    }

    // Step 10: Combine candidates and fill branches
    combineCandidates(theBeamSpotV_);

    // Fill tree if candidates found (or if MC)
    if (Pri_VtxProb->size() > 0 || doMC) {
        X_One_Tree_->Fill();
    }

    // Clear everything
    clearEventData();
}

/*****************************************************************************
 * Step 1: Process MC gen-level truth
 *****************************************************************************/
void MultiLepPAT::processMCGenInfo(const edm::Event &iEvent)
{
    handleToNtupleIndex_.clear();

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);
    if (!genParticles.isValid())
        return;

    // Flat gen-particle storage: store all relevant particles
    // (J/psi=443, Upsilon=553, phi=333, mu=13, K=321)
    static const std::unordered_set<int> interestingPdgIds = {443, 553, 333, 13, 321};

    for (size_t i = 0; i < genParticles->size(); ++i) {
        const auto &particle = genParticles->at(i);
        int absPdgId = std::abs(particle.pdgId());
        if (interestingPdgIds.count(absPdgId) || particle.numberOfDaughters() >= 2) {
            // Store if it's a relevant particle or a mother with daughters
            bool isInteresting = interestingPdgIds.count(absPdgId);
            bool hasMuDaughters = false;
            for (unsigned int d = 0; d < particle.numberOfDaughters(); ++d) {
                if (std::abs(particle.daughter(d)->pdgId()) == 13)
                    hasMuDaughters = true;
            }
            if (isInteresting || hasMuDaughters) {
                const int ntupleIndex = MC_GenPart_pdgId->size();
                MC_GenPart_pdgId->push_back(particle.pdgId());
                MC_GenPart_status->push_back(particle.status());
                int motherPdgId = 0;
                if (particle.numberOfMothers() > 0) {
                    motherPdgId = particle.mother(0)->pdgId();
                }
                MC_GenPart_motherPdgId->push_back(motherPdgId);
                MC_GenPart_handleIndex->push_back(static_cast<int>(i));
                MC_GenPart_px->push_back(particle.px());
                MC_GenPart_py->push_back(particle.py());
                MC_GenPart_pz->push_back(particle.pz());
                MC_GenPart_mass->push_back(particle.mass());
                MC_GenPart_pt->push_back(particle.pt());
                MC_GenPart_eta->push_back(particle.eta());
                MC_GenPart_phi->push_back(particle.phi());
                handleToNtupleIndex_[i] = ntupleIndex;
            }
        }
    }
}

/*****************************************************************************
 * Step 2: HLT trigger processing
 *****************************************************************************/
void MultiLepPAT::processHLTInfo(const edm::Event &iEvent)
{
    edm::Handle<edm::TriggerResults> hltresults;
    bool Error_t = false;
    unsigned int nJpsitrigger = TriggersForJpsi_.size();
    unsigned int nUpstrigger  = TriggersForUpsilon_.size();
    try {
        iEvent.getByToken(gttriggerToken_, hltresults);
    } catch (...) {
        Error_t = true;
    }
    if (Error_t || !hltresults.isValid()) {
        return;
    }

    int ntrigs = hltresults->size();
    if (ntrigs == 0) return;

    edm::TriggerNames triggerNames_ = iEvent.triggerNames(*hltresults);

    for (unsigned int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++) {
        JpsiMatchTrig[JpsiTrig] = 0;
    }
    for (unsigned int UpsTrig = 0; UpsTrig < nUpstrigger; UpsTrig++) {
        UpsilonMatchTrig[UpsTrig] = 0;
    }

    for (int itrig = 0; itrig < ntrigs; itrig++) {
        std::string trigName = triggerNames_.triggerName(itrig);
        int hltflag = (*hltresults)[itrig].accept();
        trigRes->push_back(hltflag);
        trigNames->push_back(trigName);

        for (unsigned int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++) {
            std::regex pattern(".*" + TriggersForJpsi_[JpsiTrig] + ".*");
            if (std::regex_search(trigName, pattern)) {
                JpsiMatchTrig[JpsiTrig] = hltflag;
                if (hltflag) {
                    bool isDuplicate = false;
                    for (unsigned int mt = 0; mt < MatchJpsiTrigNames->size(); mt++) {
                        if (trigName == MatchJpsiTrigNames->at(mt)) {
                            isDuplicate = true;
                            break;
                        }
                    }
                    if (!isDuplicate) {
                        MatchJpsiTrigNames->push_back(trigName);
                    }
                }
                break;
            }
        }

        // Upsilon trigger matching
        for (unsigned int UpsTrig = 0; UpsTrig < nUpstrigger; UpsTrig++) {
            std::regex pattern(".*" + TriggersForUpsilon_[UpsTrig] + ".*");
            if (std::regex_search(trigName, pattern)) {
                UpsilonMatchTrig[UpsTrig] = hltflag;
                if (hltflag) {
                    bool isDuplicate = false;
                    for (unsigned int mt = 0; mt < MatchUpsTrigNames->size(); mt++) {
                        if (trigName == MatchUpsTrigNames->at(mt)) {
                            isDuplicate = true;
                            break;
                        }
                    }
                    if (!isDuplicate) {
                        MatchUpsTrigNames->push_back(trigName);
                    }
                }
                break;
            }
        }
    }
}

/*****************************************************************************
 * Step 3: Primary vertex reconstruction
 *****************************************************************************/
void MultiLepPAT::reconstructPrimaryVertex(const edm::Event &iEvent)
{
    BeamSpot beamSpot;
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(gtbeamspotToken_, beamSpotHandle);
    if (beamSpotHandle.isValid()) {
        beamSpot = *beamSpotHandle;
        theBeamSpotV_ = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }

    edm::Handle<VertexCollection> recVtxs;
    iEvent.getByToken(gtprimaryVtxToken_, recVtxs);

    // Count good vertices using configurable cuts
    int myNGoodPrimVtx = 0;
    for (unsigned myi = 0; myi < recVtxs->size(); myi++) {
        if ((*recVtxs)[myi].ndof() >= pvNdofMin_ &&
            fabs((*recVtxs)[myi].z()) <= pvMaxAbsZ_ &&
            fabs((*recVtxs)[myi].position().rho()) <= pvMaxRho_) {
            myNGoodPrimVtx++;
        }
    }
    nGoodPrimVtx = myNGoodPrimVtx;

    // Store ALL PVs information
    if (storeAllPVs_) {
        nRecVtx = recVtxs->size();
        for (unsigned ipv = 0; ipv < recVtxs->size(); ++ipv) {
            const auto& vtx = (*recVtxs)[ipv];
            RecVtx_x->push_back(vtx.x());
            RecVtx_y->push_back(vtx.y());
            RecVtx_z->push_back(vtx.z());
            RecVtx_xErr->push_back(vtx.xError());
            RecVtx_yErr->push_back(vtx.yError());
            RecVtx_zErr->push_back(vtx.zError());
            RecVtx_chi2->push_back(vtx.chi2());
            RecVtx_ndof->push_back(vtx.ndof());
            RecVtx_nTracks->push_back(vtx.nTracks());

            double vtxProb = -9;
            try {
                vtxProb = ChiSquaredProbability(vtx.chi2(), vtx.ndof());
            } catch (...) { vtxProb = -9; }
            RecVtx_vtxProb->push_back(vtxProb);
        }
    } else {
        nRecVtx = 0;
    }

    if (recVtxs->begin() != recVtxs->end()) {
        if (pvSelectionMode_ == "mostTracks") {
            unsigned int nVtxTrks = 0;
            for (auto vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx) {
                if (nVtxTrks < vtx->tracksSize()) {
                    nVtxTrks = vtx->tracksSize();
                    thePrimaryV_ = Vertex(*vtx);
                }
            }
        } else if (pvSelectionMode_ == "highestSumPt2") {
            double maxSumPt2 = -1.0;
            for (auto vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx) {
                double sumPt2 = 0.0;
                for (auto trk = vtx->tracks_begin(); trk != vtx->tracks_end(); ++trk) {
                    sumPt2 += (*trk)->pt() * (*trk)->pt();
                }
                if (sumPt2 > maxSumPt2) {
                    maxSumPt2 = sumPt2;
                    thePrimaryV_ = Vertex(*vtx);
                }
            }
        } else {
            // Default: "firstVertex" (equivalent to legacy addXlessPrimaryVertex=true)
            thePrimaryV_ = Vertex(*(recVtxs->begin()));
        }
    } else {
        thePrimaryV_ = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }

    priVtxX = thePrimaryV_.position().x();
    priVtxY = thePrimaryV_.position().y();
    priVtxZ = thePrimaryV_.position().z();
    priVtxXE = thePrimaryV_.xError();
    priVtxYE = thePrimaryV_.yError();
    priVtxZE = thePrimaryV_.zError();
    priVtxChiNorm = thePrimaryV_.normalizedChi2();
    priVtxChi = thePrimaryV_.chi2();
    try {
        priVtxCL = ChiSquaredProbability((double)(thePrimaryV_.chi2()),
                                          (double)(thePrimaryV_.ndof()));
    } catch (...) {
        priVtxCL = -9;
    }
}

/*****************************************************************************
 * Step 4: Fill muon block
 *****************************************************************************/
void MultiLepPAT::fillMuonBlock(const edm::Event& iEvent,
                                 const edm::EventSetup& iSetup,
                                 const reco::Vertex& thePrimaryV)
{
    edm::Handle<edm::TriggerResults> hltresults;
    try { iEvent.getByToken(gttriggerToken_, hltresults); } catch (...) {}

    edm::Handle<reco::GenParticleCollection> genParticles;
    bool haveGenParticles = false;
    edm::ProductID genProductId;
    if (doMC) {
        iEvent.getByToken(genParticlesToken_, genParticles);
        haveGenParticles = genParticles.isValid();
        if (haveGenParticles) {
            genProductId = genParticles.id();
        }
    }

    auto chargeCompatible = [](const pat::Muon &muon, const reco::GenParticle &gen) {
        return gen.charge() == 0 || muon.charge() == gen.charge();
    };

    auto matchMuonFromPatRefs = [&](const pat::Muon &muon, float &bestChi2) {
        int bestIdx = -1;
        int firstValidIdx = -1;
        double bestChi2Value = std::numeric_limits<double>::infinity();
        bestChi2 = -1.f;

        for (size_t genIdx = 0; genIdx < muon.genParticlesSize(); ++genIdx) {
            const auto genRef = muon.genParticleRef(genIdx);
            if (!genRef.isNonnull() || !genRef.isAvailable())
                continue;
            if (genRef.id() != genProductId || genRef.key() >= genParticles->size())
                continue;

            const auto mapIt = handleToNtupleIndex_.find(genRef.key());
            if (mapIt == handleToNtupleIndex_.end())
                continue;

            const auto &gen = genParticles->at(genRef.key());
            if (std::abs(gen.pdgId()) != 13 || !chargeCompatible(muon, gen))
                continue;

            if (firstValidIdx < 0)
                firstValidIdx = mapIt->second;

            const double chi2 = recoGenMuonChi2(muon, gen);
            if (!std::isfinite(chi2))
                continue;

            if (chi2 < bestChi2Value) {
                bestChi2Value = chi2;
                bestIdx = mapIt->second;
            }
        }

        if (bestIdx >= 0) {
            bestChi2 = static_cast<float>(bestChi2Value);
            return bestIdx;
        }

        return firstValidIdx;
    };

    auto matchMuonByChi2 = [&](const pat::Muon &muon, float &bestChi2) {
        int bestIdx = -1;
        double bestChi2Value = recoGenMuonMatchChi2Max_;
        bestChi2 = -1.f;

        for (size_t genIdx = 0; genIdx < genParticles->size(); ++genIdx) {
            const auto &gen = genParticles->at(genIdx);
            if (std::abs(gen.pdgId()) != 13 || !chargeCompatible(muon, gen))
                continue;

            const auto mapIt = handleToNtupleIndex_.find(genIdx);
            if (mapIt == handleToNtupleIndex_.end())
                continue;

            const double chi2 = recoGenMuonChi2(muon, gen);
            if (!std::isfinite(chi2) || chi2 >= bestChi2Value)
                continue;

            bestChi2Value = chi2;
            bestIdx = mapIt->second;
        }

        if (bestIdx >= 0)
            bestChi2 = static_cast<float>(bestChi2Value);

        return bestIdx;
    };

    for (auto iMuonP = thePATMuonHandle_->begin();
         iMuonP != thePATMuonHandle_->end(); ++iMuonP) {
        ++nMu;
        muIsPatLooseMuon->push_back(iMuonP->isLooseMuon());
        muIsPatTightMuon->push_back(iMuonP->isTightMuon(thePrimaryV));
        muIsPatSoftMuon->push_back(iMuonP->isSoftMuon(thePrimaryV));
        muIsPatMediumMuon->push_back(iMuonP->isMediumMuon());

        muPx->push_back(iMuonP->px());
        muPy->push_back(iMuonP->py());
        muPz->push_back(iMuonP->pz());
        muCharge->push_back(iMuonP->charge());

        // Muon momentum errors
        if (storeMuonMomentumErrors_) {
            float pxErr = -9, pyErr = -9, pzErr = -9, ptErr = -9;
            reco::TrackRef muTrack = iMuonP->track();

            if (muTrack.isNull() || !muTrack.isAvailable()) {
                // Try bestTrack() as fallback
                try {
                    auto bestTrk = iMuonP->muonBestTrack();
                    if (bestTrk.isAvailable()) {
                        auto cov = bestTrk->covariance();
                        pxErr = std::sqrt(std::max(cov(0,0), 0.0));
                        pyErr = std::sqrt(std::max(cov(1,1), 0.0));
                        pzErr = std::sqrt(std::max(cov(2,2), 0.0));

                        float px = bestTrk->momentum().x();
                        float py = bestTrk->momentum().y();
                        float pt = std::sqrt(px*px + py*py);
                        if (pt > 0) {
                            float dpt_dpx = px/pt;
                            float dpt_dpy = py/pt;
                            ptErr = std::sqrt(dpt_dpx*dpt_dpx*cov(0,0) +
                                              dpt_dpy*dpt_dpy*cov(1,1) +
                                              2.0*dpt_dpx*dpt_dpy*cov(0,1));
                        }
                    }
                } catch (...) { /* keep -9 */ }
            } else {
                auto cov = muTrack->covariance();
                pxErr = std::sqrt(std::max(cov(0,0), 0.0));
                pyErr = std::sqrt(std::max(cov(1,1), 0.0));
                pzErr = std::sqrt(std::max(cov(2,2), 0.0));

                float px = muTrack->momentum().x();
                float py = muTrack->momentum().y();
                float pt = std::sqrt(px*px + py*py);
                if (pt > 0) {
                    float dpt_dpx = px/pt;
                    float dpt_dpy = py/pt;
                    ptErr = std::sqrt(dpt_dpx*dpt_dpx*cov(0,0) +
                                      dpt_dpy*dpt_dpy*cov(1,1) +
                                      2.0*dpt_dpx*dpt_dpy*cov(0,1));
                }
            }
            muPxErr->push_back(pxErr);
            muPyErr->push_back(pyErr);
            muPzErr->push_back(pzErr);
            muPtErr->push_back(ptErr);
        } else {
            muPxErr->push_back(-9);
            muPyErr->push_back(-9);
            muPzErr->push_back(-9);
            muPtErr->push_back(-9);
        }

        // Match muon to packed candidate for fromPV info
        unsigned int curFromPV = 0;
        unsigned int curPVAssocQuality = 0;
        double curMuMomentum = std::sqrt(
            iMuonP->px()*iMuonP->px() + iMuonP->py()*iMuonP->py() + iMuonP->pz()*iMuonP->pz());
        
        for (auto iTrackfID = nonMuonTrack_.begin();
             iTrackfID != nonMuonTrack_.end(); ++iTrackfID) {
            try {
                if (iMuonP->track().isNull()) continue;
                auto iTrackf = *iTrackfID;
                if (fabs(iTrackf->px() - iMuonP->px()) < MuMatchTrkMomentumRelDiffThr_c * curMuMomentum &&
                    fabs(iTrackf->py() - iMuonP->py()) < MuMatchTrkMomentumRelDiffThr_c * curMuMomentum &&
                    fabs(iTrackf->pz() - iMuonP->pz()) < MuMatchTrkMomentumRelDiffThr_c * curMuMomentum &&
                    iTrackf->charge() == iMuonP->charge()) {
                    nonMuonTrack_.erase(iTrackfID);
                    curFromPV = iTrackf->fromPV();
                    curPVAssocQuality = iTrackf->pvAssociationQuality();
                    break;
                }
            } catch (...) { continue; }
        }
        
        muFromPV->push_back(curFromPV);
        muPVAssocQuality->push_back(curPVAssocQuality);

        // Surplus: Try to get PV association via sourceCandidatePtr
        if (storeMuonPVAssoc_) {
            int vertexId = -1;
            float dzAssocPV = -9, dxyAssocPV = -9;
            int fromPVAssocPV = 0;
            int pdgId = 0;

            try {
                for (size_t i = 0; i < iMuonP->numberOfSourceCandidatePtrs(); ++i) {
                    reco::CandidatePtr sourcePtr = iMuonP->sourceCandidatePtr(i);
                    if (sourcePtr.isAvailable()) {
                        const pat::PackedCandidate* pfcand =
                            dynamic_cast<const pat::PackedCandidate*>(sourcePtr.get());
                        if (pfcand) {
                            vertexId = pfcand->vertexRef().key();
                            // Get associated PV position for dz/dxy calculation
                            edm::Handle<VertexCollection> recVtxs;
                            iEvent.getByToken(gtprimaryVtxToken_, recVtxs);
                            if (vertexId >= 0 && vertexId < (int)recVtxs->size()) {
                                const auto& assocPV = (*recVtxs)[vertexId];
                                auto bestTrk = iMuonP->muonBestTrack();
                                if (bestTrk.isAvailable()) {
                                    dzAssocPV = bestTrk->dz(assocPV.position());
                                    dxyAssocPV = bestTrk->dxy(assocPV.position());
                                }
                            }
                            fromPVAssocPV = pfcand->fromPV();
                            pdgId = pfcand->pdgId();
                            break;
                        }
                    }
                }
            } catch (...) {
                // Keep default values
            }

            muVertexId->push_back(vertexId);
            muDzAssocPV->push_back(dzAssocPV);
            muDxyAssocPV->push_back(dxyAssocPV);
            muFromPVAssocPV->push_back(fromPVAssocPV);
            muPdgId->push_back(pdgId);
        } else {
            muVertexId->push_back(-1);
            muDzAssocPV->push_back(-9);
            muDxyAssocPV->push_back(-9);
            muFromPVAssocPV->push_back(0);
            muPdgId->push_back(0);
        }

        int genMatchIdx = -1;
        int genMatchSource = 0;
        float genMatchChi2 = -1.f;
        if (haveGenParticles) {
            genMatchIdx = matchMuonFromPatRefs(*iMuonP, genMatchChi2);
            if (genMatchIdx >= 0) {
                genMatchSource = 1;
            } else {
                genMatchIdx = matchMuonByChi2(*iMuonP, genMatchChi2);
                if (genMatchIdx >= 0)
                    genMatchSource = 2;
            }
        }

        muGenMatchIdx->push_back(genMatchIdx);
        muGenMatchSource->push_back(genMatchSource);
        muGenMatchChi2->push_back(genMatchChi2);

        // Trigger matching
        bool isJpsiTrigMatch = false;
        bool isJpsiFilterMatch = false;
        for (unsigned int JpsiTrig = 0; JpsiTrig < TriggersForJpsi_.size(); JpsiTrig++) {
            if (JpsiMatchTrig[JpsiTrig] != 0) isJpsiTrigMatch = true;
        }
        muIsJpsiTrigMatch->push_back(isJpsiTrigMatch);

        for (unsigned int JpsiFilter = 0; JpsiFilter < FiltersForJpsi_.size(); JpsiFilter++) {
            if (isJpsiTrigMatch && hltresults.isValid()) {
                for (auto it = iMuonP->triggerObjectMatches().begin();
                     it != iMuonP->triggerObjectMatches().end(); ++it) {
                    pat::TriggerObjectStandAlone tempObj(*it);
                    tempObj.unpackFilterLabels(iEvent, *hltresults);
                    if (tempObj.hasFilterLabel(FiltersForJpsi_[JpsiFilter])) {
                        isJpsiFilterMatch = true;
                    }
                }
            }
        }
        muIsJpsiFilterMatch->push_back(isJpsiFilterMatch);

        // Upsilon trigger matching
        bool isUpsTrigMatch = false;
        bool isUpsFilterMatch = false;
        for (unsigned int UpsTrig = 0; UpsTrig < TriggersForUpsilon_.size(); UpsTrig++) {
            if (UpsilonMatchTrig[UpsTrig] != 0) isUpsTrigMatch = true;
        }
        muIsUpsTrigMatch->push_back(isUpsTrigMatch);

        for (unsigned int UpsFilter = 0; UpsFilter < FiltersForUpsilon_.size(); UpsFilter++) {
            if (isUpsTrigMatch && hltresults.isValid()) {
                for (auto it = iMuonP->triggerObjectMatches().begin();
                     it != iMuonP->triggerObjectMatches().end(); ++it) {
                    pat::TriggerObjectStandAlone tempObj(*it);
                    tempObj.unpackFilterLabels(iEvent, *hltresults);
                    if (tempObj.hasFilterLabel(FiltersForUpsilon_[UpsFilter])) {
                        isUpsFilterMatch = true;
                    }
                }
            }
        }
        muIsUpsFilterMatch->push_back(isUpsFilterMatch);

        munMatchedSeg->push_back(-1);
    }
}

/*****************************************************************************
 * Step 5: Pair muons into onia candidates
 *****************************************************************************/
void MultiLepPAT::pairMuons(const edm::Handle<edm::View<pat::Muon>>& muonHandle,
                             const MagneticField& bField)
{
    KinematicParticleFactoryFromTransientTrack muPairFactory;
    ParticleMass muMass = myMuMass;
    float muMassSigma = myMuMassErr;
    float chi2 = 0., ndof = 0.;

    std::vector<RefCountedKinematicParticle> transMuonPair;
    std::vector<uint> transMuPairId;

    muPairCand_Onia1_.clear();
    muPairCand_Onia2_.clear();
    muQuad_Onia_.clear();

    // Determine mass windows based on analysis mode
    double onia1MassMin = JpsiMassMin_, onia1MassMax = JpsiMassMax_;
    double onia2MassMin = JpsiMassMin_, onia2MassMax = JpsiMassMax_;
    
    if (analysisChannel_ == AnalysisChannel::JpsiUpsPhi) {
        onia2MassMin = UpsMassMin_;
        onia2MassMax = UpsMassMax_;
    } else if (analysisChannel_ == AnalysisChannel::JpsiJpsiUps) {
        onia2MassMin = UpsMassMin_;
        onia2MassMax = UpsMassMax_;
    }

    for (auto iMuon1 = muonHandle->begin(); iMuon1 != muonHandle->end(); ++iMuon1) {
        TrackRef muTrack1 = iMuon1->track();
        if (muTrack1.isNull()) continue;
        // Use StringCutObjectSelector for muon cuts
        if (!muonSelector_(*iMuon1)) continue;

        TransientTrack transTrk1(muTrack1, &bField);
        transMuonPair.push_back(muPairFactory.particle(transTrk1, muMass, chi2, ndof, muMassSigma));
        transMuPairId.push_back(iMuon1 - muonHandle->begin());
        
        for (auto iMuon2 = iMuon1 + 1; iMuon2 != muonHandle->end(); ++iMuon2) {
            TrackRef muTrack2 = iMuon2->track();
            if (muTrack2.isNull()) continue;
            if (!muonSelector_(*iMuon2)) continue;

            // Charge requirement
            if ((iMuon1->charge() + iMuon2->charge()) != 0) continue;

            double muPairMass = (iMuon1->p4() + iMuon2->p4()).mass();
            
            // Check against both J/psi and Upsilon windows
            bool isOnia1 = (onia1MassMin < muPairMass && muPairMass < onia1MassMax);
            bool isOnia2 = (onia2MassMin < muPairMass && muPairMass < onia2MassMax);
            
            if (!isOnia1 && !isOnia2) continue;

            // Per-resonance pair pT and eta pre-cuts
            auto pairP4 = iMuon1->p4() + iMuon2->p4();
            if (isOnia1) {
                if (pairP4.Pt() < jpsiCandPtMin_ || std::abs(pairP4.Eta()) > jpsiCandEtaMax_)
                    isOnia1 = false;
            }
            if (isOnia2) {
                double ptCut  = (analysisChannel_ == AnalysisChannel::JpsiJpsiPhi) ? jpsiCandPtMin_  : upsCandPtMin_;
                double etaCut = (analysisChannel_ == AnalysisChannel::JpsiJpsiPhi) ? jpsiCandEtaMax_ : upsCandEtaMax_;
                if (pairP4.Pt() < ptCut || std::abs(pairP4.Eta()) > etaCut)
                    isOnia2 = false;
            }
            if (!isOnia1 && !isOnia2) continue;

            TransientTrack transTrk2(muTrack2, &bField);
            transMuonPair.push_back(muPairFactory.particle(transTrk2, muMass, chi2, ndof, muMassSigma));
            transMuPairId.push_back(iMuon2 - muonHandle->begin());

            if (!particlesToVtx(transMuonPair, OniaDecayVtxProbCut_)) {
                transMuonPair.pop_back();
                transMuPairId.pop_back();
                continue;
            }

            RefCountedKinematicTree muVtxFitTree;
            particlesToVtx(muVtxFitTree, transMuonPair, "muon pair", OniaDecayVtxProbCut_);

            // Store in appropriate container based on mass window and analysis mode
            if (analysisChannel_ == AnalysisChannel::JpsiJpsiPhi) {
                // JpsiJpsiPhi: Onia1==Onia2==J/psi, store all in Onia1_ (upper-triangle quartets)
                muPairCand_Onia1_.push_back(std::make_pair(transMuonPair, transMuPairId));
            } else {
                // JpsiJpsiUps or JpsiUpsPhi: separate Jpsi→Onia1, Ups→Onia2
                if (isOnia1) {
                    muPairCand_Onia1_.push_back(std::make_pair(transMuonPair, transMuPairId));
                }
                if (isOnia2) {
                    muPairCand_Onia2_.push_back(std::make_pair(transMuonPair, transMuPairId));
                }
            }

            transMuonPair.pop_back();
            transMuPairId.pop_back();
        }
        transMuonPair.pop_back();
        transMuPairId.pop_back();
    }

    // Build muon quartets from non-overlapping pairs
    RefCountedKinematicTree vtxFitTree_1, vtxFitTree_2;
    RefCountedKinematicParticle fit1, fit2;
    RefCountedKinematicVertex vtx1, vtx2;
    double massErr1, massErr2;
    std::vector<RefCountedKinematicParticle> interOnia;

    if (analysisChannel_ == AnalysisChannel::JpsiUpsPhi) {
        // JpsiUpsPhi: cross-product Onia1(J/psi) × Onia2(Upsilon)
        for (auto& pair1 : muPairCand_Onia1_) {
            for (auto& pair2 : muPairCand_Onia2_) {
                if (isOverlapPair(pair1, pair2)) continue;

                bool valid1 = particlesToVtx(vtxFitTree_1, pair1.first, "Onia_1", OniaDecayVtxProbCut_);
                bool valid2 = particlesToVtx(vtxFitTree_2, pair2.first, "Onia_2", OniaDecayVtxProbCut_);

                if (valid1 && valid2) {
                    extractFitRes(vtxFitTree_1, fit1, vtx1, massErr1);
                    extractFitRes(vtxFitTree_2, fit2, vtx2, massErr2);
                    if (massErr1 >= 0.0 && massErr2 >= 0.0) {
                        interOnia.push_back(fit1);
                        interOnia.push_back(fit2);
                        if (particlesToVtx(interOnia, OniaDecayVtxProbCut_)) {
                            muQuad_Onia_.push_back(std::make_pair(pair1, pair2));
                        }
                        interOnia.clear();
                    }
                }
            }
        }
    } else {
        // JpsiJpsiPhi or JpsiJpsiUps: upper-triangle of Onia1_ for same-species
        for (auto pair1 = muPairCand_Onia1_.begin(); pair1 != muPairCand_Onia1_.end(); ++pair1) {
            for (auto pair2 = pair1 + 1; pair2 != muPairCand_Onia1_.end(); ++pair2) {
                if (isOverlapPair(*pair1, *pair2)) continue;

                bool valid1 = particlesToVtx(vtxFitTree_1, pair1->first, "Onia_1", OniaDecayVtxProbCut_);
                bool valid2 = particlesToVtx(vtxFitTree_2, pair2->first, "Onia_2", OniaDecayVtxProbCut_);

                if (valid1 && valid2) {
                    extractFitRes(vtxFitTree_1, fit1, vtx1, massErr1);
                    extractFitRes(vtxFitTree_2, fit2, vtx2, massErr2);
                    if (massErr1 >= 0.0 && massErr2 >= 0.0) {
                        interOnia.push_back(fit1);
                        interOnia.push_back(fit2);
                        if (particlesToVtx(interOnia, OniaDecayVtxProbCut_)) {
                            muQuad_Onia_.push_back(std::make_pair(*pair1, *pair2));
                        }
                        interOnia.clear();
                    }
                }
            }
        }
    }
}

/*****************************************************************************
 * Step 6: Pair tracks into meson candidates
 *****************************************************************************/
void MultiLepPAT::pairTracks(
    const std::vector<edm::View<pat::PackedCandidate>::const_iterator>& tracks,
    const MagneticField& bField)
{
    KPairCand_Meson_.clear();
    KinematicParticleFactoryFromTransientTrack PhiFactory;
    ParticleMass KMass = myKMass;
    float KMassSigma = myKMassErr;
    float chi2 = 0., ndof = 0.;

    std::vector<RefCountedKinematicParticle> transTrackPair;
    std::vector<uint> transTrackPairId;

    for (auto iTrack1ID = tracks.begin(); iTrack1ID != tracks.end(); ++iTrack1ID) {
        auto iTrack1 = *iTrack1ID;
        if (!iTrack1->hasTrackDetails() || iTrack1->charge() == 0) continue;
        if (static_cast<int>(iTrack1->fromPV()) < minTrackFromPV_) continue;
        // Use StringCutObjectSelector for track cuts
        if (!trackSelector_(*iTrack1)) continue;
        if (iTrack1->bestTrack()->normalizedChi2() > 8 ||
            !iTrack1->bestTrack()->quality(reco::Track::highPurity)) continue;

        TransientTrack trackTT1(*(iTrack1->bestTrack()), &bField);
        transTrackPair.push_back(PhiFactory.particle(trackTT1, KMass, chi2, ndof, KMassSigma));
        transTrackPairId.push_back(iTrack1ID - tracks.begin());

        for (auto iTrack2ID = iTrack1ID + 1; iTrack2ID != tracks.end(); ++iTrack2ID) {
            auto iTrack2 = *iTrack2ID;
            if (!iTrack2->hasTrackDetails() || iTrack2->charge() == 0) continue;
            if (static_cast<int>(iTrack2->fromPV()) < minTrackFromPV_) continue;
            if (!trackSelector_(*iTrack2)) continue;
            if (iTrack2->bestTrack()->normalizedChi2() > 8 ||
                !iTrack2->bestTrack()->quality(reco::Track::highPurity)) continue;
            if ((iTrack1->charge() + iTrack2->charge()) != 0) continue;

            TLorentzVector P4_Track1, P4_Track2, P4_Meson;
            P4_Track1.SetPtEtaPhiM(iTrack1->pt(), iTrack1->eta(), iTrack1->phi(), myKMass);
            P4_Track2.SetPtEtaPhiM(iTrack2->pt(), iTrack2->eta(), iTrack2->phi(), myKMass);
            P4_Meson = P4_Track1 + P4_Track2;

            if (P4_Track1.DeltaR(P4_Meson) > trackDRMax_) continue;
            if (P4_Track2.DeltaR(P4_Meson) > trackDRMax_) continue;

            bool inMassWindow = (PhiMassMin_ < P4_Meson.M() && P4_Meson.M() < PhiMassMax_);
            if (!inMassWindow) continue;

            // Per-resonance Phi pair pT/eta pre-cuts
            if (P4_Meson.Pt() < phiCandPtMin_ || std::abs(P4_Meson.Eta()) > phiCandEtaMax_) continue;

            TransientTrack trackTT2(*(iTrack2->bestTrack()), &bField);
            transTrackPair.push_back(PhiFactory.particle(trackTT2, KMass, chi2, ndof, KMassSigma));
            transTrackPairId.push_back(iTrack2ID - tracks.begin());

            if (particlesToVtx(transTrackPair, OniaDecayVtxProbCut_)) {
                KPairCand_Meson_.push_back(std::make_pair(transTrackPair, transTrackPairId));
            }
            transTrackPair.pop_back();
            transTrackPairId.pop_back();
        }
        transTrackPair.pop_back();
        transTrackPairId.pop_back();
    }
}

/*****************************************************************************
 * Step 7: Combine candidates and fill final branches
 *****************************************************************************/
void MultiLepPAT::combineCandidates(const reco::Vertex& beamSpotV)
{
    RefCountedKinematicTree vtxFitTree_1, vtxFitTree_2, vtxFitTree_Meson, vtxFitTree_Pri;
    RefCountedKinematicParticle fit1, fit2, fitMeson, fitPri;
    RefCountedKinematicVertex vtx1, vtx2, vtxMeson, vtxPri;
    double massErr1, massErr2, massErrMeson, massErrPri;
    std::vector<RefCountedKinematicParticle> interOnia;

    // Determine the PDG mass for the 3rd resonance
    double mesonPdgMass = myPhiMass;
    if (analysisChannel_ == AnalysisChannel::JpsiJpsiUps) {
        mesonPdgMass = myUpsMass;
    }

    reco::Vertex bsV = const_cast<reco::Vertex&>(beamSpotV);

    // For JpsiJpsiPhi and JpsiUpsPhi: loop over track pairs + muon quartets
    if (analysisChannel_ == AnalysisChannel::JpsiJpsiPhi ||
        analysisChannel_ == AnalysisChannel::JpsiUpsPhi) {
        
        for (auto& KPair : KPairCand_Meson_) {
            bool validMeson = particlesToVtx(vtxFitTree_Meson, KPair.first, "Meson", OniaDecayVtxProbCut_);
            if (!validMeson) continue;
            extractFitRes(vtxFitTree_Meson, fitMeson, vtxMeson, massErrMeson);
            if (massErrMeson < 0) continue;

            for (auto& muQuad : muQuad_Onia_) {
                bool valid1 = particlesToVtx(vtxFitTree_1, muQuad.first.first, "Onia_1", OniaDecayVtxProbCut_);
                bool valid2 = particlesToVtx(vtxFitTree_2, muQuad.second.first, "Onia_2", OniaDecayVtxProbCut_);

                if (!valid1 || !valid2) continue;
                extractFitRes(vtxFitTree_1, fit1, vtx1, massErr1);
                extractFitRes(vtxFitTree_2, fit2, vtx2, massErr2);
                if (massErr1 < 0 || massErr2 < 0) continue;

                interOnia.clear();
                interOnia.push_back(fit1);
                interOnia.push_back(fit2);
                interOnia.push_back(fitMeson);

                bool validPri = particlesToVtx(vtxFitTree_Pri, interOnia, "primary vertex");

                // Check FINAL fitted mass against mass windows (not just pre-fit)
                if (checkFinalMass_) {
                    double fittedMass1 = fit1->currentState().mass();
                    double fittedMass2 = fit2->currentState().mass();
                    double fittedMassMeson = fitMeson->currentState().mass();

                    // Determine mass windows based on analysis mode
                    double onia1MassMin = JpsiMassMin_, onia1MassMax = JpsiMassMax_;
                    double onia2MassMin = JpsiMassMin_, onia2MassMax = JpsiMassMax_;
                    if (analysisChannel_ == AnalysisChannel::JpsiUpsPhi) {
                        onia2MassMin = UpsMassMin_;
                        onia2MassMax = UpsMassMax_;
                    }

                    bool fittedMass1InRange = (onia1MassMin < fittedMass1 && fittedMass1 < onia1MassMax);
                    bool fittedMass2InRange = (onia2MassMin < fittedMass2 && fittedMass2 < onia2MassMax);
                    bool fittedMassMesonInRange = (PhiMassMin_ < fittedMassMeson && fittedMassMeson < PhiMassMax_);

                    if (!fittedMass1InRange || !fittedMass2InRange || !fittedMassMesonInRange) {
                        continue;
                    }
                }

                if (!validPri) {
                    // 3-body fit failed: store sentinel for Pri and keep the rest
                    storeSentinelPri();
                } else {
                    extractFitRes(vtxFitTree_Pri, fitPri, vtxPri, massErrPri);

                    // Store primary vertex
                    storeResonanceBranches(fitPri, vtxPri, massErrPri, 0, bsV,
                        Pri_mass, Pri_massErr, nullptr, Pri_ctau, Pri_ctauErr,
                        Pri_Chi2, Pri_ndof, Pri_VtxProb,
                        Pri_px, Pri_py, Pri_pz, Pri_phi, Pri_eta, Pri_pt,
                        Pri_pxErr, Pri_pyErr, Pri_pzErr, Pri_ptErr);
                }

                // Store Onia 1 (J/psi)
                storeResonanceBranches(fit1, vtx1, massErr1, myJpsiMass, bsV,
                    Jpsi_1_mass, Jpsi_1_massErr, Jpsi_1_massDiff,
                    Jpsi_1_ctau, Jpsi_1_ctauErr,
                    Jpsi_1_Chi2, Jpsi_1_ndof, Jpsi_1_VtxProb,
                    Jpsi_1_px, Jpsi_1_py, Jpsi_1_pz,
                    Jpsi_1_phi, Jpsi_1_eta, Jpsi_1_pt,
                    Jpsi_1_pxErr, Jpsi_1_pyErr, Jpsi_1_pzErr, Jpsi_1_ptErr);
                // Store muon indices
                Jpsi_1_mu_1_Idx->push_back(muQuad.first.second[0]);
                Jpsi_1_mu_2_Idx->push_back(muQuad.first.second[1]);
                // Store Onia 2 (J/psi for JpsiJpsiPhi, Upsilon for JpsiUpsPhi)
                double onia2PdgMass = myJpsiMass;
                if (analysisChannel_ == AnalysisChannel::JpsiJpsiPhi) {
                    storeResonanceBranches(fit2, vtx2, massErr2, onia2PdgMass, bsV,
                        Jpsi_2_mass, Jpsi_2_massErr, Jpsi_2_massDiff,
                        Jpsi_2_ctau, Jpsi_2_ctauErr,
                        Jpsi_2_Chi2, Jpsi_2_ndof, Jpsi_2_VtxProb,
                        Jpsi_2_px, Jpsi_2_py, Jpsi_2_pz,
                        Jpsi_2_phi, Jpsi_2_eta, Jpsi_2_pt,
                        Jpsi_2_pxErr, Jpsi_2_pyErr, Jpsi_2_pzErr, Jpsi_2_ptErr);
                    Jpsi_2_mu_1_Idx->push_back(muQuad.second.second[0]);
                    Jpsi_2_mu_2_Idx->push_back(muQuad.second.second[1]);
                    }

                // For JpsiUpsPhi: also store Upsilon info in Ups_* branches
                if (analysisChannel_ == AnalysisChannel::JpsiUpsPhi) {
                    onia2PdgMass = myUpsMass;
                    storeResonanceBranches(fit2, vtx2, massErr2, myUpsMass, bsV,
                        Ups_mass, Ups_massErr, Ups_massDiff,
                        Ups_ctau, Ups_ctauErr,
                        Ups_Chi2, Ups_ndof, Ups_VtxProb,
                        Ups_px, Ups_py, Ups_pz,
                        Ups_phi, Ups_eta, Ups_pt,
                        Ups_pxErr, Ups_pyErr, Ups_pzErr, Ups_ptErr);
                    Ups_mu_1_Idx->push_back(muQuad.second.second[0]);
                    Ups_mu_2_Idx->push_back(muQuad.second.second[1]);
                }

                Phi_K_1_Idx->push_back(KPair.second[0]);
                Phi_K_2_Idx->push_back(KPair.second[1]);

                // Store meson (Phi)
                storeResonanceBranches(fitMeson, vtxMeson, massErrMeson, mesonPdgMass, bsV,
                    Phi_mass, Phi_massErr, Phi_massDiff,
                    Phi_ctau, Phi_ctauErr,
                    Phi_Chi2, Phi_ndof, Phi_VtxProb,
                    Phi_px, Phi_py, Phi_pz,
                    Phi_phi, Phi_eta, Phi_pt,
                    Phi_pxErr, Phi_pyErr, Phi_pzErr, Phi_ptErr);

                // Store kaon kinematics
                Phi_K_1_px->push_back(nonMuonTrack_[KPair.second[0]]->px());
                Phi_K_1_py->push_back(nonMuonTrack_[KPair.second[0]]->py());
                Phi_K_1_pz->push_back(nonMuonTrack_[KPair.second[0]]->pz());
                Phi_K_1_pt->push_back(nonMuonTrack_[KPair.second[0]]->pt());
                Phi_K_1_eta->push_back(nonMuonTrack_[KPair.second[0]]->eta());
                Phi_K_1_phi->push_back(nonMuonTrack_[KPair.second[0]]->phi());
                Phi_K_1_fromPV->push_back(nonMuonTrack_[KPair.second[0]]->fromPV());
                Phi_K_1_pvAssocQuality->push_back(nonMuonTrack_[KPair.second[0]]->pvAssociationQuality());

                Phi_K_2_px->push_back(nonMuonTrack_[KPair.second[1]]->px());
                Phi_K_2_py->push_back(nonMuonTrack_[KPair.second[1]]->py());
                Phi_K_2_pz->push_back(nonMuonTrack_[KPair.second[1]]->pz());
                Phi_K_2_pt->push_back(nonMuonTrack_[KPair.second[1]]->pt());
                Phi_K_2_eta->push_back(nonMuonTrack_[KPair.second[1]]->eta());
                Phi_K_2_phi->push_back(nonMuonTrack_[KPair.second[1]]->phi());
                Phi_K_2_fromPV->push_back(nonMuonTrack_[KPair.second[1]]->fromPV());
                Phi_K_2_pvAssocQuality->push_back(nonMuonTrack_[KPair.second[1]]->pvAssociationQuality());
            }
        }
    }
    // For JpsiJpsiUps: loop over J/psi quartets + Upsilon pair (6 muons total)
    if (analysisChannel_ == AnalysisChannel::JpsiJpsiUps) {
        for (auto& muQuad : muQuad_Onia_) {
            for (auto& upsPair : muPairCand_Onia2_) {
                // Check no overlap between quartet muons and Ups pair muons
                if (isOverlapPair(muQuad.first, upsPair) ||
                    isOverlapPair(muQuad.second, upsPair)) continue;

                bool valid1 = particlesToVtx(vtxFitTree_1, muQuad.first.first, "Jpsi_1", OniaDecayVtxProbCut_);
                bool valid2 = particlesToVtx(vtxFitTree_2, muQuad.second.first, "Jpsi_2", OniaDecayVtxProbCut_);

                RefCountedKinematicTree vtxFitTree_Ups;
                bool validUps = particlesToVtx(vtxFitTree_Ups, upsPair.first, "Ups", OniaDecayVtxProbCut_);

                if (!valid1 || !valid2 || !validUps) continue;

                RefCountedKinematicParticle fitUps;
                RefCountedKinematicVertex vtxUps;
                double massErrUps;

                extractFitRes(vtxFitTree_1, fit1, vtx1, massErr1);
                extractFitRes(vtxFitTree_2, fit2, vtx2, massErr2);
                extractFitRes(vtxFitTree_Ups, fitUps, vtxUps, massErrUps);
                if (massErr1 < 0 || massErr2 < 0 || massErrUps < 0) continue;

                interOnia.clear();
                interOnia.push_back(fit1);
                interOnia.push_back(fit2);
                interOnia.push_back(fitUps);

                bool validPri = particlesToVtx(vtxFitTree_Pri, interOnia, "primary vertex");
                
                // Check FINAL fitted masses
                if (checkFinalMass_) {
                    double fittedMass1 = fit1->currentState().mass();
                    double fittedMass2 = fit2->currentState().mass();
                    double fittedMassUps = fitUps->currentState().mass();

                    bool fittedMass1InRange = (JpsiMassMin_ < fittedMass1 && fittedMass1 < JpsiMassMax_);
                    bool fittedMass2InRange = (JpsiMassMin_ < fittedMass2 && fittedMass2 < JpsiMassMax_);
                    bool fittedMassUpsInRange = (UpsMassMin_ < fittedMassUps && fittedMassUps < UpsMassMax_);

                    if (!fittedMass1InRange || !fittedMass2InRange || !fittedMassUpsInRange) {
                        continue;
                    }
                }
                
                if (!validPri) {
                    // 3-body fit failed: store sentinel for Pri and keep the rest
                    storeSentinelPri();
                } else {
                    extractFitRes(vtxFitTree_Pri, fitPri, vtxPri, massErrPri);

                    // Store primary vertex
                    storeResonanceBranches(fitPri, vtxPri, massErrPri, 0, bsV,
                        Pri_mass, Pri_massErr, nullptr, Pri_ctau, Pri_ctauErr,
                        Pri_Chi2, Pri_ndof, Pri_VtxProb,
                        Pri_px, Pri_py, Pri_pz, Pri_phi, Pri_eta, Pri_pt,
                        Pri_pxErr, Pri_pyErr, Pri_pzErr, Pri_ptErr);
                }

                // Store J/psi 1
                storeResonanceBranches(fit1, vtx1, massErr1, myJpsiMass, bsV,
                    Jpsi_1_mass, Jpsi_1_massErr, Jpsi_1_massDiff,
                    Jpsi_1_ctau, Jpsi_1_ctauErr,
                    Jpsi_1_Chi2, Jpsi_1_ndof, Jpsi_1_VtxProb,
                    Jpsi_1_px, Jpsi_1_py, Jpsi_1_pz,
                    Jpsi_1_phi, Jpsi_1_eta, Jpsi_1_pt,
                    Jpsi_1_pxErr, Jpsi_1_pyErr, Jpsi_1_pzErr, Jpsi_1_ptErr);

                // Store J/psi 2
                storeResonanceBranches(fit2, vtx2, massErr2, myJpsiMass, bsV,
                    Jpsi_2_mass, Jpsi_2_massErr, Jpsi_2_massDiff,
                    Jpsi_2_ctau, Jpsi_2_ctauErr,
                    Jpsi_2_Chi2, Jpsi_2_ndof, Jpsi_2_VtxProb,
                    Jpsi_2_px, Jpsi_2_py, Jpsi_2_pz,
                    Jpsi_2_phi, Jpsi_2_eta, Jpsi_2_pt,
                    Jpsi_2_pxErr, Jpsi_2_pyErr, Jpsi_2_pzErr, Jpsi_2_ptErr);

                // Store Upsilon
                storeResonanceBranches(fitUps, vtxUps, massErrUps, myUpsMass, bsV,
                    Ups_mass, Ups_massErr, Ups_massDiff,
                    Ups_ctau, Ups_ctauErr,
                    Ups_Chi2, Ups_ndof, Ups_VtxProb,
                    Ups_px, Ups_py, Ups_pz,
                    Ups_phi, Ups_eta, Ups_pt,
                    Ups_pxErr, Ups_pyErr, Ups_pzErr, Ups_ptErr);

                // Store muon indices
                Jpsi_1_mu_1_Idx->push_back(muQuad.first.second[0]);
                Jpsi_1_mu_2_Idx->push_back(muQuad.first.second[1]);
                Jpsi_2_mu_1_Idx->push_back(muQuad.second.second[0]);
                Jpsi_2_mu_2_Idx->push_back(muQuad.second.second[1]);
                Ups_mu_1_Idx->push_back(upsPair.second[0]);
                Ups_mu_2_Idx->push_back(upsPair.second[1]);
            }
        }
    }
}

/*****************************************************************************
 * Helper: Store sentinel values for failed 3-body vertex fit
 *****************************************************************************/
void MultiLepPAT::storeSentinelPri()
{
    const float sentinel = -999999.0f;
    Pri_mass->push_back(sentinel);
    Pri_massErr->push_back(sentinel);
    Pri_ctau->push_back(sentinel);
    Pri_ctauErr->push_back(sentinel);
    Pri_Chi2->push_back(sentinel);
    Pri_ndof->push_back(sentinel);
    Pri_VtxProb->push_back(sentinel);
    Pri_px->push_back(sentinel);
    Pri_py->push_back(sentinel);
    Pri_pz->push_back(sentinel);
    Pri_phi->push_back(sentinel);
    Pri_eta->push_back(sentinel);
    Pri_pt->push_back(sentinel);
    Pri_pxErr->push_back(sentinel);
    Pri_pyErr->push_back(sentinel);
    Pri_pzErr->push_back(sentinel);
    Pri_ptErr->push_back(sentinel);
}

/*****************************************************************************
 * Helper: Store resonance fit results into branches
 *****************************************************************************/
void MultiLepPAT::storeResonanceBranches(
    const RefCountedKinematicParticle& fitPart,
    const RefCountedKinematicVertex& vtx,
    double massErr, double pdgMass,
    reco::Vertex& beamSpotV,
    vector<float>* br_mass, vector<float>* br_massErr, vector<float>* br_massDiff,
    vector<float>* br_ctau, vector<float>* br_ctauErr,
    vector<float>* br_Chi2, vector<float>* br_ndof, vector<float>* br_VtxProb,
    vector<float>* br_px, vector<float>* br_py, vector<float>* br_pz,
    vector<float>* br_phi, vector<float>* br_eta, vector<float>* br_pt,
    vector<float>* br_pxErr, vector<float>* br_pyErr, vector<float>* br_pzErr,
    vector<float>* br_ptErr)
{
    double tmp_pt, tmp_eta, tmp_phi;
    getDynamics(fitPart, tmp_pt, tmp_eta, tmp_phi);

    double tmp_pxErr, tmp_pyErr, tmp_pzErr, tmp_ptErr;
    getMomentumErrors(fitPart, tmp_pxErr, tmp_pyErr, tmp_pzErr, tmp_ptErr);

    br_mass->push_back(fitPart->currentState().mass());
    br_massErr->push_back(massErr);
    if (br_massDiff) {
        br_massDiff->push_back(fitPart->currentState().mass() - pdgMass);
    }

    RefCountedKinematicVertex vtxCopy = vtx;
    RefCountedKinematicParticle partCopy = fitPart;
    br_ctau->push_back(GetcTau(vtxCopy, partCopy, beamSpotV));
    br_ctauErr->push_back(GetcTauErr(vtxCopy, partCopy, beamSpotV));
    br_Chi2->push_back(double(vtx->chiSquared()));
    br_ndof->push_back(double(vtx->degreesOfFreedom()));
    br_VtxProb->push_back(ChiSquaredProbability(
        (double)(vtx->chiSquared()), (double)(vtx->degreesOfFreedom())));

    br_px->push_back(fitPart->currentState().kinematicParameters().momentum().x());
    br_py->push_back(fitPart->currentState().kinematicParameters().momentum().y());
    br_pz->push_back(fitPart->currentState().kinematicParameters().momentum().z());
    br_phi->push_back(tmp_phi);
    br_eta->push_back(tmp_eta);
    br_pt->push_back(tmp_pt);

    if (br_pxErr) br_pxErr->push_back(tmp_pxErr);
    if (br_pyErr) br_pyErr->push_back(tmp_pyErr);
    if (br_pzErr) br_pzErr->push_back(tmp_pzErr);
    if (br_ptErr) br_ptErr->push_back(tmp_ptErr);
}

/*****************************************************************************
 * Step 8: MC gen-level matching
 *****************************************************************************/
void MultiLepPAT::doMCGenMatching(
    const edm::Handle<edm::View<pat::Muon>>& muonHandle,
    const edm::Handle<edm::View<pat::PackedCandidate>>& trackHandle)
{
    // Basic deltaR-based matching placeholder
    // This stores matched reco tracks to gen-level kaons
    // (Full matching logic to be expanded as needed)
    if (!Match_pi1px || !MC_GenPart_px) return;
    // Matching is done post-hoc in offline analysis
}

/*****************************************************************************
 * Step 9: Clear all event data
 *****************************************************************************/
void MultiLepPAT::clearEventData()
{
    // MC gen-level (new)
    if (doMC) {
        MC_GenPart_pdgId->clear(); MC_GenPart_status->clear();
        MC_GenPart_motherPdgId->clear();
        MC_GenPart_handleIndex->clear();
        MC_GenPart_px->clear(); MC_GenPart_py->clear();
        MC_GenPart_pz->clear(); MC_GenPart_mass->clear();
        MC_GenPart_pt->clear(); MC_GenPart_eta->clear();
        MC_GenPart_phi->clear();
    }

    trigRes->clear(); trigNames->clear();
    L1TT->clear(); MatchJpsiTrigNames->clear(); MatchUpsTrigNames->clear();
    muIsJpsiTrigMatch->clear(); muIsJpsiFilterMatch->clear();
    muIsUpsTrigMatch->clear(); muIsUpsFilterMatch->clear();

    runNum = 0; evtNum = 0; lumiNum = 0; nGoodPrimVtx = 0;
    priVtxX = 0; priVtxY = 0; priVtxZ = 0;
    priVtxXE = 0; priVtxYE = 0; priVtxZE = 0;
    priVtxChiNorm = 0; priVtxChi = 0; priVtxCL = 0;

    PriVtxXCorrX->clear(); PriVtxXCorrY->clear(); PriVtxXCorrZ->clear();
    PriVtxXCorrEX->clear(); PriVtxXCorrEY->clear(); PriVtxXCorrEZ->clear();
    PriVtxXCorrC2->clear(); PriVtxXCorrCL->clear();

    // All primary vertices
    nRecVtx = 0;
    RecVtx_x->clear(); RecVtx_y->clear(); RecVtx_z->clear();
    RecVtx_xErr->clear(); RecVtx_yErr->clear(); RecVtx_zErr->clear();
    RecVtx_chi2->clear(); RecVtx_ndof->clear(); RecVtx_vtxProb->clear();
    RecVtx_nTracks->clear();

    nMu = 0;
    muPx->clear(); muPy->clear(); muPz->clear();
    muD0->clear(); muD0E->clear(); muDz->clear();
    muChi2->clear(); muGlChi2->clear(); mufHits->clear();
    muFirstBarrel->clear(); muFirstEndCap->clear();
    muDzVtx->clear(); muDxyVtx->clear();
    muNDF->clear(); muGlNDF->clear(); muPhits->clear();
    muShits->clear(); muGlMuHits->clear();
    muType->clear(); muQual->clear();
    muTrack->clear(); muCharge->clear(); muIsoratio->clear();
    muIsGoodLooseMuon->clear(); muIsGoodLooseMuonNew->clear();
    muIsGoodSoftMuonNewIlse->clear(); muIsGoodSoftMuonNewIlseMod->clear();
    muIsGlobalMuon->clear(); muIsGoodTightMuon->clear();
    munMatchedSeg->clear();
    muMVAMuonID->clear(); musegmentCompatibility->clear();
    mupulldXdZ_pos_noArb->clear(); mupulldYdZ_pos_noArb->clear();
    mupulldXdZ_pos_ArbDef->clear(); mupulldYdZ_pos_ArbDef->clear();
    mupulldXdZ_pos_ArbST->clear(); mupulldYdZ_pos_ArbST->clear();
    mupulldXdZ_pos_noArb_any->clear(); mupulldYdZ_pos_noArb_any->clear();
    muIsPatLooseMuon->clear(); muIsPatTightMuon->clear();
    muIsPatSoftMuon->clear(); muIsPatMediumMuon->clear();
    muFromPV->clear(); muPVAssocQuality->clear();

    // Muon momentum errors
    muPxErr->clear(); muPyErr->clear(); muPzErr->clear(); muPtErr->clear();

    // Muon PV association (surplus from sourceCandidatePtr)
    muVertexId->clear(); muDzAssocPV->clear(); muDxyAssocPV->clear();
    muFromPVAssocPV->clear(); muPdgId->clear();
    muGenMatchIdx->clear(); muGenMatchSource->clear(); muGenMatchChi2->clear();

    // Muon-track matching debug
    if (muTrkMatchDebug_) {
        muMatch_nCandidates->clear(); muMatch_bestRelDiff->clear();
        muMatch_bestDz->clear(); muMatch_methodUsed->clear();
    }

    // Resonance branches
    Pri_mass->clear(); Pri_massErr->clear();
    Pri_ctau->clear(); Pri_ctauErr->clear();
    Pri_Chi2->clear(); Pri_ndof->clear(); Pri_VtxProb->clear();
    Pri_px->clear(); Pri_py->clear(); Pri_pz->clear();
    Pri_phi->clear(); Pri_eta->clear(); Pri_pt->clear();
    Pri_pxErr->clear(); Pri_pyErr->clear(); Pri_pzErr->clear(); Pri_ptErr->clear();

    Jpsi_1_mass->clear(); Jpsi_1_massErr->clear(); Jpsi_1_massDiff->clear();
    Jpsi_1_ctau->clear(); Jpsi_1_ctauErr->clear();
    Jpsi_1_Chi2->clear(); Jpsi_1_ndof->clear(); Jpsi_1_VtxProb->clear();
    Jpsi_1_px->clear(); Jpsi_1_py->clear(); Jpsi_1_pz->clear();
    Jpsi_1_phi->clear(); Jpsi_1_eta->clear(); Jpsi_1_pt->clear();
    Jpsi_1_mu_1_Idx->clear(); Jpsi_1_mu_2_Idx->clear();
    Jpsi_1_pxErr->clear(); Jpsi_1_pyErr->clear(); Jpsi_1_pzErr->clear(); Jpsi_1_ptErr->clear();

    Jpsi_2_mass->clear(); Jpsi_2_massErr->clear(); Jpsi_2_massDiff->clear();
    Jpsi_2_ctau->clear(); Jpsi_2_ctauErr->clear();
    Jpsi_2_Chi2->clear(); Jpsi_2_ndof->clear(); Jpsi_2_VtxProb->clear();
    Jpsi_2_px->clear(); Jpsi_2_py->clear(); Jpsi_2_pz->clear();
    Jpsi_2_phi->clear(); Jpsi_2_eta->clear(); Jpsi_2_pt->clear();
    Jpsi_2_mu_1_Idx->clear(); Jpsi_2_mu_2_Idx->clear();
    Jpsi_2_pxErr->clear(); Jpsi_2_pyErr->clear(); Jpsi_2_pzErr->clear(); Jpsi_2_ptErr->clear();

    Phi_mass->clear(); Phi_massErr->clear(); Phi_massDiff->clear();
    Phi_ctau->clear(); Phi_ctauErr->clear();
    Phi_Chi2->clear(); Phi_ndof->clear(); Phi_VtxProb->clear();
    Phi_px->clear(); Phi_py->clear(); Phi_pz->clear();
    Phi_phi->clear(); Phi_eta->clear(); Phi_pt->clear();
    Phi_pxErr->clear(); Phi_pyErr->clear(); Phi_pzErr->clear(); Phi_ptErr->clear();

    Phi_K_1_Idx->clear(); Phi_K_1_px->clear(); Phi_K_1_py->clear(); Phi_K_1_pz->clear();
    Phi_K_1_phi->clear(); Phi_K_1_eta->clear(); Phi_K_1_pt->clear();
    Phi_K_1_fromPV->clear(); Phi_K_1_pvAssocQuality->clear();
    Phi_K_2_Idx->clear(); Phi_K_2_px->clear(); Phi_K_2_py->clear(); Phi_K_2_pz->clear();
    Phi_K_2_phi->clear(); Phi_K_2_eta->clear(); Phi_K_2_pt->clear();
    Phi_K_2_fromPV->clear(); Phi_K_2_pvAssocQuality->clear();

    // Upsilon branches
    Ups_mu_1_Idx->clear(); Ups_mu_2_Idx->clear();
    Ups_mass->clear(); Ups_massErr->clear(); Ups_massDiff->clear();
    Ups_ctau->clear(); Ups_ctauErr->clear();
    Ups_Chi2->clear(); Ups_ndof->clear(); Ups_VtxProb->clear();
    Ups_px->clear(); Ups_py->clear(); Ups_pz->clear();
    Ups_phi->clear(); Ups_eta->clear(); Ups_pt->clear();
    Ups_pxErr->clear(); Ups_pyErr->clear(); Ups_pzErr->clear(); Ups_ptErr->clear();

    // Clear intermediate storage
    handleToNtupleIndex_.clear();
    muPairCand_Onia1_.clear();
    muPairCand_Onia2_.clear();
    muQuad_Onia_.clear();
    KPairCand_Meson_.clear();
    nonMuonTrack_.clear();
}

/*****************************************************************************
 * Utility: Momentum uncertainty extraction
 *****************************************************************************/
void MultiLepPAT::getMomentumErrors(const RefCountedKinematicParticle& arg_Part,
                                     double& res_pxErr, double& res_pyErr,
                                     double& res_pzErr, double& res_ptErr)
{
    try {
        // The 7x7 covariance matrix has indices: (px=3, py=4, pz=5, m=6)
        // but KinematicParametersError uses (x=0,y=1,z=2,px=3,py=4,pz=5,m=6)
        auto errMatrix = arg_Part->currentState().kinematicParametersError().matrix();
        res_pxErr = std::sqrt(std::max(errMatrix(3, 3), 0.0));
        res_pyErr = std::sqrt(std::max(errMatrix(4, 4), 0.0));
        res_pzErr = std::sqrt(std::max(errMatrix(5, 5), 0.0));
        
        // pt error via error propagation: pt = sqrt(px^2 + py^2)
        double px = arg_Part->currentState().kinematicParameters().momentum().x();
        double py = arg_Part->currentState().kinematicParameters().momentum().y();
        double pt = std::sqrt(px*px + py*py);
        if (pt > 0) {
            double dpt_dpx = px / pt;
            double dpt_dpy = py / pt;
            res_ptErr = std::sqrt(
                dpt_dpx * dpt_dpx * errMatrix(3, 3) +
                dpt_dpy * dpt_dpy * errMatrix(4, 4) +
                2.0 * dpt_dpx * dpt_dpy * errMatrix(3, 4));
        } else {
            res_ptErr = -9;
        }
    } catch (...) {
        res_pxErr = -9; res_pyErr = -9; res_pzErr = -9; res_ptErr = -9;
    }
}

/*****************************************************************************
 * Utility methods (getDynamics, tracksToMuonPair, particlesToVtx, etc.)
 * These are carried forward from the original with minimal changes.
 *****************************************************************************/

void MultiLepPAT::getDynamics(double arg_mass, double arg_px, double arg_py, double arg_pz,
                              double& res_pt, double& res_eta, double& res_phi) {
    TLorentzVector myParticle;
    myParticle.SetXYZM(arg_px, arg_py, arg_pz, arg_mass);
    res_pt  = myParticle.Pt();
    res_eta = myParticle.Eta();
    res_phi = myParticle.Phi();
}

void MultiLepPAT::getDynamics(const RefCountedKinematicParticle& arg_Part,
                              double& res_pt, double& res_eta, double& res_phi) {
    getDynamics(arg_Part->currentState().mass(),
                arg_Part->currentState().kinematicParameters().momentum().x(),
                arg_Part->currentState().kinematicParameters().momentum().y(),
                arg_Part->currentState().kinematicParameters().momentum().z(),
                res_pt, res_eta, res_phi);
}

void MultiLepPAT::tracksToMuonPair(
    vector<RefCountedKinematicParticle>& arg_MuonResults,
    KinematicParticleFactoryFromTransientTrack& arg_MuFactory,
    const MagneticField& arg_bField,
    const TrackRef arg_Trk1, const TrackRef arg_Trk2) {
    TransientTrack transTrk1(arg_Trk1, &(arg_bField));
    TransientTrack transTrk2(arg_Trk2, &(arg_bField));
    ParticleMass muMass = myMuMass;
    float muMassSigma = myMuMassErr;
    float chi2 = 0., ndof = 0.;
    arg_MuonResults.push_back(arg_MuFactory.particle(transTrk1, muMass, chi2, ndof, muMassSigma));
    arg_MuonResults.push_back(arg_MuFactory.particle(transTrk2, muMass, chi2, ndof, muMassSigma));
}

// --- particlesToVtx overloads ---

bool MultiLepPAT::particlesToVtx(const vector<RefCountedKinematicParticle>& arg_FromParticles) {
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree vertexFitTree;
    try { vertexFitTree = fitter.fit(arg_FromParticles); } catch (...) { return false; }
    if (!vertexFitTree->isValid()) return false;
    try {
        return (vertexFitTree->currentDecayVertex()->chiSquared() >= 0.0);
    } catch (...) { return false; }
}

bool MultiLepPAT::particlesToVtx(const vector<RefCountedKinematicParticle>& arg_FromParticles,
                                  const string& arg_Message) {
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree vertexFitTree;
    try { vertexFitTree = fitter.fit(arg_FromParticles); }
    catch (...) { return false; }
    if (!vertexFitTree->isValid()) return false;
    try {
        return (vertexFitTree->currentDecayVertex()->chiSquared() >= 0.0);
    } catch (...) { return false; }
}

bool MultiLepPAT::particlesToVtx(RefCountedKinematicTree& arg_VertexFitTree,
                                  const vector<RefCountedKinematicParticle>& arg_FromParticles,
                                  const string& arg_Message) {
    KinematicParticleVertexFitter fitter;
    try { arg_VertexFitTree = fitter.fit(arg_FromParticles); }
    catch (...) { return false; }
    if (!arg_VertexFitTree->isValid()) return false;
    try {
        return (arg_VertexFitTree->currentDecayVertex()->chiSquared() >= 0.0);
    } catch (...) { return false; }
}

bool MultiLepPAT::particlesToVtx(const vector<RefCountedKinematicParticle>& arg_FromParticles,
                                  const double& arg_VtxProbCut) {
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree vertexFitTree;
    try { vertexFitTree = fitter.fit(arg_FromParticles); } catch (...) { return false; }
    if (!vertexFitTree->isValid()) return false;
    auto vtx = vertexFitTree->currentDecayVertex();
    try {
        if (vtx->chiSquared() < 0.0) return false;
    } catch (...) { return false; }
    double vtxprob = ChiSquaredProbability((double)(vtx->chiSquared()),
                                           (double)(vtx->degreesOfFreedom()));
    return (vtxprob >= arg_VtxProbCut);
}

bool MultiLepPAT::particlesToVtx(const vector<RefCountedKinematicParticle>& arg_FromParticles,
                                  const string& arg_Message,
                                  const double& arg_VtxProbCut) {
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree vertexFitTree;
    try { vertexFitTree = fitter.fit(arg_FromParticles); } catch (...) { return false; }
    if (!vertexFitTree->isValid()) return false;
    auto vtx = vertexFitTree->currentDecayVertex();
    try {
        if (vtx->chiSquared() < 0.0) return false;
    } catch (...) { return false; }
    double vtxprob = ChiSquaredProbability((double)(vtx->chiSquared()),
                                           (double)(vtx->degreesOfFreedom()));
    return (vtxprob >= arg_VtxProbCut);
}

bool MultiLepPAT::particlesToVtx(RefCountedKinematicTree& arg_VertexFitTree,
                                  const vector<RefCountedKinematicParticle>& arg_FromParticles,
                                  const string& arg_Message,
                                  const double& arg_VtxProbCut) {
    KinematicParticleVertexFitter fitter;
    try { arg_VertexFitTree = fitter.fit(arg_FromParticles); } catch (...) { return false; }
    if (!arg_VertexFitTree->isValid()) return false;
    auto vtx = arg_VertexFitTree->currentDecayVertex();
    try {
        if (vtx->chiSquared() < 0.0) return false;
    } catch (...) { return false; }
    double vtxprob = ChiSquaredProbability((double)(vtx->chiSquared()),
                                           (double)(vtx->degreesOfFreedom()));
    return (vtxprob >= arg_VtxProbCut);
}

// --- extractFitRes overloads ---

bool MultiLepPAT::extractFitRes(RefCountedKinematicTree& arg_VtxTree,
                                 RefCountedKinematicParticle& res_Part,
                                 RefCountedKinematicVertex& res_Vtx,
                                 KinematicParameters& res_Param,
                                 double& res_MassErr) {
    double tmp_MassErr2 = 0.0;
    arg_VtxTree->movePointerToTheTop();
    res_Part = arg_VtxTree->currentParticle();
    res_Vtx = arg_VtxTree->currentDecayVertex();
    res_Param = res_Part->currentState().kinematicParameters();
    tmp_MassErr2 = res_Part->currentState().kinematicParametersError().matrix()(6, 6);
    res_MassErr = (tmp_MassErr2 < 0.0) ? -9 : std::sqrt(tmp_MassErr2);
    return (res_MassErr >= 0.0);
}

bool MultiLepPAT::extractFitRes(RefCountedKinematicTree& arg_VtxTree,
                                 RefCountedKinematicParticle& res_Part,
                                 RefCountedKinematicVertex& res_Vtx,
                                 double& res_MassErr) {
    double tmp_MassErr2 = 0.0;
    arg_VtxTree->movePointerToTheTop();
    try {
        res_Part = arg_VtxTree->currentParticle();
        res_Vtx = arg_VtxTree->currentDecayVertex();
        tmp_MassErr2 = res_Part->currentState().kinematicParametersError().matrix()(6, 6);
    } catch (...) {
        tmp_MassErr2 = -9;
    }
    res_MassErr = (tmp_MassErr2 < 0.0) ? -9 : std::sqrt(tmp_MassErr2);
    return (res_MassErr >= 0.0);
}

bool MultiLepPAT::extractFitRes(RefCountedKinematicTree& arg_VtxTree,
                                 RefCountedKinematicVertex& res_Vtx,
                                 double& res_VtxProb) {
    try {
        arg_VtxTree->movePointerToTheTop();
        res_Vtx = arg_VtxTree->currentDecayVertex();
        res_VtxProb = ChiSquaredProbability((double)(res_Vtx->chiSquared()),
                                             (double)(res_Vtx->degreesOfFreedom()));
    } catch (...) {
        res_VtxProb = -9.0;
        return false;
    }
    return (res_VtxProb >= 0.0);
}

bool MultiLepPAT::isOverlapPair(const muList_t& arg_MuonPair1,
                                 const muList_t& arg_MuonPair2) {
    return (arg_MuonPair1.second[0] == arg_MuonPair2.second[0] ||
            arg_MuonPair1.second[0] == arg_MuonPair2.second[1] ||
            arg_MuonPair1.second[1] == arg_MuonPair2.second[0] ||
            arg_MuonPair1.second[1] == arg_MuonPair2.second[1]);
}

double MultiLepPAT::fitResEval(double arg_massDiff_1, double arg_massErr_1,
                                double arg_massDiff_2, double arg_massErr_2,
                                double arg_massDiff_3, double arg_massErr_3) {
    return arg_massDiff_1 * arg_massDiff_1 / (arg_massErr_1 * arg_massErr_1) +
           arg_massDiff_2 * arg_massDiff_2 / (arg_massErr_2 * arg_massErr_2) +
           arg_massDiff_3 * arg_massDiff_3 / (arg_massErr_3 * arg_massErr_3);
}

// --- CTau ---
double MultiLepPAT::GetcTau(RefCountedKinematicVertex& decayVrtx,
                              RefCountedKinematicParticle& kinePart,
                              Vertex& bs) {
    TVector3 vtx((*decayVrtx).position().x(), (*decayVrtx).position().y(), 0);
    TVector3 pvtx(bs.position().x(), bs.position().y(), 0);
    TVector3 pperp(kinePart->currentState().globalMomentum().x(),
                   kinePart->currentState().globalMomentum().y(), 0);
    TVector3 vdiff = vtx - pvtx;
    double LxyPV = vdiff.Dot(pperp) / pperp.Mag();
    return LxyPV * kinePart->currentState().mass() / pperp.Perp();
}

double MultiLepPAT::GetcTauErr(RefCountedKinematicVertex& decayVrtx,
                                 RefCountedKinematicParticle& kinePart,
                                 Vertex& bs) {
    TVector3 pperp(kinePart->currentState().globalMomentum().x(),
                   kinePart->currentState().globalMomentum().y(), 0);
    try {
        AlgebraicVector vpperp(3);
        vpperp[0] = pperp.x(); vpperp[1] = pperp.y(); vpperp[2] = 0.;
        GlobalError v1e = (*decayVrtx).error();
        GlobalError v2e = bs.error();
        AlgebraicSymMatrix vXYe = asHepMatrix(v1e.matrix()) + asHepMatrix(v2e.matrix());
        return sqrt(vXYe.similarity(vpperp)) * kinePart->currentState().mass() / (pperp.Perp2());
    } catch (...) {
        return -99999;
    }
}

double MultiLepPAT::deltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > M_PI) dphi -= 2 * M_PI;
    while (dphi <= -M_PI) dphi += 2 * M_PI;
    return sqrt(deta * deta + dphi * dphi);
}

void MultiLepPAT::printKinematics(const RefCountedKinematicParticle& particle,
                                    const std::string& name) {
    const auto& state = particle->currentState();
    std::cout << name << " kinematics:" << std::endl;
    std::cout << "  px=" << state.globalMomentum().x()
              << " py=" << state.globalMomentum().y()
              << " pz=" << state.globalMomentum().z()
              << " pt=" << state.globalMomentum().perp()
              << " eta=" << state.globalMomentum().eta()
              << " phi=" << state.globalMomentum().phi() << std::endl;
}

/*****************************************************************************
 * beginRun / beginJob / endJob
 *****************************************************************************/
void MultiLepPAT::beginJob()
{
    edm::Service<TFileService> fs;
    X_One_Tree_ = fs->make<TTree>("X_data", "Quarkonia Reconstruction Data");

    X_One_Tree_->Branch("TrigRes", &trigRes);
    X_One_Tree_->Branch("TrigNames", &trigNames);
    X_One_Tree_->Branch("MatchJpsiTriggerNames", &MatchJpsiTrigNames);
    X_One_Tree_->Branch("MatchUpsTriggerNames", &MatchUpsTrigNames);
    X_One_Tree_->Branch("L1TrigRes", &L1TT);

    X_One_Tree_->Branch("evtNum", &evtNum, "evtNum/i");
    X_One_Tree_->Branch("runNum", &runNum, "runNum/i");
    X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
    X_One_Tree_->Branch("nGoodPrimVtx", &nGoodPrimVtx, "nGoodPrimVtx/i");

    X_One_Tree_->Branch("priVtxX", &priVtxX, "priVtxX/f");
    X_One_Tree_->Branch("priVtxY", &priVtxY, "priVtxY/f");
    X_One_Tree_->Branch("priVtxZ", &priVtxZ, "priVtxZ/f");
    X_One_Tree_->Branch("priVtxXE", &priVtxXE, "priVtxXE/f");
    X_One_Tree_->Branch("priVtxYE", &priVtxYE, "priVtxYE/f");
    X_One_Tree_->Branch("priVtxZE", &priVtxZE, "priVtxZE/f");
    X_One_Tree_->Branch("priVtxChiNorm", &priVtxChiNorm, "priVtxChiNorm/f");
    X_One_Tree_->Branch("priVtxChi", &priVtxChi, "priVtxChi/f");
    X_One_Tree_->Branch("priVtxCL", &priVtxCL, "priVtxCL/f");
    X_One_Tree_->Branch("PriVtxXCorrX", &PriVtxXCorrX);
    X_One_Tree_->Branch("PriVtxXCorrY", &PriVtxXCorrY);
    X_One_Tree_->Branch("PriVtxXCorrZ", &PriVtxXCorrZ);
    X_One_Tree_->Branch("PriVtxXCorrEX", &PriVtxXCorrEX);
    X_One_Tree_->Branch("PriVtxXCorrEY", &PriVtxXCorrEY);
    X_One_Tree_->Branch("PriVtxXCorrEZ", &PriVtxXCorrEZ);
    X_One_Tree_->Branch("PriVtxXCorrC2", &PriVtxXCorrC2);
    X_One_Tree_->Branch("PriVtxXCorrCL", &PriVtxXCorrCL);

    // All primary vertices (not just selected)
    X_One_Tree_->Branch("nRecVtx", &nRecVtx, "nRecVtx/i");
    X_One_Tree_->Branch("RecVtx_x", &RecVtx_x);
    X_One_Tree_->Branch("RecVtx_y", &RecVtx_y);
    X_One_Tree_->Branch("RecVtx_z", &RecVtx_z);
    X_One_Tree_->Branch("RecVtx_xErr", &RecVtx_xErr);
    X_One_Tree_->Branch("RecVtx_yErr", &RecVtx_yErr);
    X_One_Tree_->Branch("RecVtx_zErr", &RecVtx_zErr);
    X_One_Tree_->Branch("RecVtx_chi2", &RecVtx_chi2);
    X_One_Tree_->Branch("RecVtx_ndof", &RecVtx_ndof);
    X_One_Tree_->Branch("RecVtx_vtxProb", &RecVtx_vtxProb);
    X_One_Tree_->Branch("RecVtx_nTracks", &RecVtx_nTracks);

    // Muon branches
    X_One_Tree_->Branch("nMu", &nMu, "nMu/i");
    X_One_Tree_->Branch("muPx", &muPx); X_One_Tree_->Branch("muPy", &muPy); X_One_Tree_->Branch("muPz", &muPz);
    X_One_Tree_->Branch("muD0", &muD0); X_One_Tree_->Branch("muD0E", &muD0E); X_One_Tree_->Branch("muDz", &muDz);
    X_One_Tree_->Branch("muChi2", &muChi2); X_One_Tree_->Branch("muGlChi2", &muGlChi2);
    X_One_Tree_->Branch("mufHits", &mufHits);
    X_One_Tree_->Branch("muFirstBarrel", &muFirstBarrel); X_One_Tree_->Branch("muFirstEndCap", &muFirstEndCap);
    X_One_Tree_->Branch("muDzVtx", &muDzVtx); X_One_Tree_->Branch("muDxyVtx", &muDxyVtx);
    X_One_Tree_->Branch("muNDF", &muNDF); X_One_Tree_->Branch("muGlNDF", &muGlNDF);
    X_One_Tree_->Branch("muPhits", &muPhits); X_One_Tree_->Branch("muShits", &muShits);
    X_One_Tree_->Branch("muGlMuHits", &muGlMuHits);
    X_One_Tree_->Branch("muType", &muType); X_One_Tree_->Branch("muQual", &muQual);
    X_One_Tree_->Branch("muTrack", &muTrack); X_One_Tree_->Branch("muCharge", &muCharge);
    X_One_Tree_->Branch("muIsoratio", &muIsoratio);
    X_One_Tree_->Branch("munMatchedSeg", &munMatchedSeg);
    X_One_Tree_->Branch("muIsGoodSoftMuonNewIlseMod", &muIsGoodSoftMuonNewIlseMod);
    X_One_Tree_->Branch("muIsGoodSoftMuonNewIlse", &muIsGoodSoftMuonNewIlse);
    X_One_Tree_->Branch("muIsGoodLooseMuonNew", &muIsGoodLooseMuonNew);
    X_One_Tree_->Branch("muIsGoodLooseMuon", &muIsGoodLooseMuon);
    X_One_Tree_->Branch("muIsGoodTightMuon", &muIsGoodTightMuon);
    X_One_Tree_->Branch("muIsGlobalMuon", &muIsGlobalMuon);
    X_One_Tree_->Branch("muIsPatLooseMuon", &muIsPatLooseMuon);
    X_One_Tree_->Branch("muIsPatTightMuon", &muIsPatTightMuon);
    X_One_Tree_->Branch("muIsPatSoftMuon", &muIsPatSoftMuon);
    X_One_Tree_->Branch("muIsPatMediumMuon", &muIsPatMediumMuon);
    X_One_Tree_->Branch("muFromPV", &muFromPV);
    X_One_Tree_->Branch("muPVAssocQuality", &muPVAssocQuality);

    // Muon momentum errors
    X_One_Tree_->Branch("muPxErr", &muPxErr);
    X_One_Tree_->Branch("muPyErr", &muPyErr);
    X_One_Tree_->Branch("muPzErr", &muPzErr);
    X_One_Tree_->Branch("muPtErr", &muPtErr);

    // Muon PV association (surplus from sourceCandidatePtr)
    X_One_Tree_->Branch("muVertexId", &muVertexId);
    X_One_Tree_->Branch("muDzAssocPV", &muDzAssocPV);
    X_One_Tree_->Branch("muDxyAssocPV", &muDxyAssocPV);
    X_One_Tree_->Branch("muFromPVAssocPV", &muFromPVAssocPV);
    X_One_Tree_->Branch("muPdgId", &muPdgId);
    X_One_Tree_->Branch("muGenMatchIdx", &muGenMatchIdx);
    X_One_Tree_->Branch("muGenMatchSource", &muGenMatchSource);
    X_One_Tree_->Branch("muGenMatchChi2", &muGenMatchChi2);

    X_One_Tree_->Branch("muIsJpsiTrigMatch", &muIsJpsiTrigMatch);
    X_One_Tree_->Branch("muIsUpsTrigMatch", &muIsUpsTrigMatch);
    X_One_Tree_->Branch("muIsJpsiFilterMatch", &muIsJpsiFilterMatch);
    X_One_Tree_->Branch("muIsUpsFilterMatch", &muIsUpsFilterMatch);
    X_One_Tree_->Branch("muMVAMuonID", &muMVAMuonID);
    X_One_Tree_->Branch("musegmentCompatibility", &musegmentCompatibility);
    X_One_Tree_->Branch("mupulldXdZ_pos_noArb", &mupulldXdZ_pos_noArb);
    X_One_Tree_->Branch("mupulldYdZ_pos_noArb", &mupulldYdZ_pos_noArb);
    X_One_Tree_->Branch("mupulldXdZ_pos_ArbDef", &mupulldXdZ_pos_ArbDef);
    X_One_Tree_->Branch("mupulldYdZ_pos_ArbDef", &mupulldYdZ_pos_ArbDef);
    X_One_Tree_->Branch("mupulldXdZ_pos_ArbST", &mupulldXdZ_pos_ArbST);
    X_One_Tree_->Branch("mupulldYdZ_pos_ArbST", &mupulldYdZ_pos_ArbST);
    X_One_Tree_->Branch("mupulldXdZ_pos_noArb_any", &mupulldXdZ_pos_noArb_any);
    X_One_Tree_->Branch("mupulldYdZ_pos_noArb_any", &mupulldYdZ_pos_noArb_any);

    // Muon-track matching debug (only if muTrkMatchDebug_)
    if (muTrkMatchDebug_) {
        X_One_Tree_->Branch("muMatch_nCandidates", &muMatch_nCandidates);
        X_One_Tree_->Branch("muMatch_bestRelDiff", &muMatch_bestRelDiff);
        X_One_Tree_->Branch("muMatch_bestDz", &muMatch_bestDz);
        X_One_Tree_->Branch("muMatch_methodUsed", &muMatch_methodUsed);
    }

    // Resonance branches (with momentum errors)
    auto branchReso = [&](const std::string& prefix,
                          vector<float>*& mass, vector<float>*& massErr, vector<float>*& massDiff,
                          vector<float>*& ctau, vector<float>*& ctauErr,
                          vector<float>*& Chi2, vector<float>*& ndof, vector<float>*& VtxProb,
                          vector<float>*& px, vector<float>*& py, vector<float>*& pz,
                          vector<float>*& phi, vector<float>*& eta, vector<float>*& pt,
                          vector<float>*& pxErr, vector<float>*& pyErr, vector<float>*& pzErr,
                          vector<float>*& ptErr) {
        X_One_Tree_->Branch((prefix + "_mass").c_str(), &mass);
        X_One_Tree_->Branch((prefix + "_massErr").c_str(), &massErr);
        if (massDiff) X_One_Tree_->Branch((prefix + "_massDiff").c_str(), &massDiff);
        X_One_Tree_->Branch((prefix + "_ctau").c_str(), &ctau);
        X_One_Tree_->Branch((prefix + "_ctauErr").c_str(), &ctauErr);
        X_One_Tree_->Branch((prefix + "_Chi2").c_str(), &Chi2);
        X_One_Tree_->Branch((prefix + "_ndof").c_str(), &ndof);
        X_One_Tree_->Branch((prefix + "_VtxProb").c_str(), &VtxProb);
        X_One_Tree_->Branch((prefix + "_px").c_str(), &px);
        X_One_Tree_->Branch((prefix + "_py").c_str(), &py);
        X_One_Tree_->Branch((prefix + "_pz").c_str(), &pz);
        X_One_Tree_->Branch((prefix + "_phi").c_str(), &phi);
        X_One_Tree_->Branch((prefix + "_eta").c_str(), &eta);
        X_One_Tree_->Branch((prefix + "_pt").c_str(), &pt);
        X_One_Tree_->Branch((prefix + "_pxErr").c_str(), &pxErr);
        X_One_Tree_->Branch((prefix + "_pyErr").c_str(), &pyErr);
        X_One_Tree_->Branch((prefix + "_pzErr").c_str(), &pzErr);
        X_One_Tree_->Branch((prefix + "_ptErr").c_str(), &ptErr);
    };

    // Pre-allocate massDiff vectors so the if(massDiff) check in branchReso passes.
    // Without this, all *_massDiff pointers remain nullptr (set in the constructor)
    // and clearEventData() would segfault on ->clear().
    Jpsi_1_massDiff = new vector<float>();
    Jpsi_2_massDiff = new vector<float>();
    Phi_massDiff    = new vector<float>();
    Ups_massDiff    = new vector<float>();

    branchReso("Jpsi_1", Jpsi_1_mass, Jpsi_1_massErr, Jpsi_1_massDiff,
               Jpsi_1_ctau, Jpsi_1_ctauErr, Jpsi_1_Chi2, Jpsi_1_ndof, Jpsi_1_VtxProb,
               Jpsi_1_px, Jpsi_1_py, Jpsi_1_pz, Jpsi_1_phi, Jpsi_1_eta, Jpsi_1_pt,
               Jpsi_1_pxErr, Jpsi_1_pyErr, Jpsi_1_pzErr, Jpsi_1_ptErr);
    X_One_Tree_->Branch("Jpsi_1_mu_1_Idx", &Jpsi_1_mu_1_Idx);
    X_One_Tree_->Branch("Jpsi_1_mu_2_Idx", &Jpsi_1_mu_2_Idx);

    branchReso("Jpsi_2", Jpsi_2_mass, Jpsi_2_massErr, Jpsi_2_massDiff,
               Jpsi_2_ctau, Jpsi_2_ctauErr, Jpsi_2_Chi2, Jpsi_2_ndof, Jpsi_2_VtxProb,
               Jpsi_2_px, Jpsi_2_py, Jpsi_2_pz, Jpsi_2_phi, Jpsi_2_eta, Jpsi_2_pt,
               Jpsi_2_pxErr, Jpsi_2_pyErr, Jpsi_2_pzErr, Jpsi_2_ptErr);
    X_One_Tree_->Branch("Jpsi_2_mu_1_Idx", &Jpsi_2_mu_1_Idx);
    X_One_Tree_->Branch("Jpsi_2_mu_2_Idx", &Jpsi_2_mu_2_Idx);

    branchReso("Phi", Phi_mass, Phi_massErr, Phi_massDiff,
               Phi_ctau, Phi_ctauErr, Phi_Chi2, Phi_ndof, Phi_VtxProb,
               Phi_px, Phi_py, Phi_pz, Phi_phi, Phi_eta, Phi_pt,
               Phi_pxErr, Phi_pyErr, Phi_pzErr, Phi_ptErr);
    X_One_Tree_->Branch("Phi_K_1_Idx", &Phi_K_1_Idx);
    X_One_Tree_->Branch("Phi_K_2_Idx", &Phi_K_2_Idx);

    vector<float>* nullDiff = nullptr;
    branchReso("Pri", Pri_mass, Pri_massErr, nullDiff,
               Pri_ctau, Pri_ctauErr, Pri_Chi2, Pri_ndof, Pri_VtxProb,
               Pri_px, Pri_py, Pri_pz, Pri_phi, Pri_eta, Pri_pt,
               Pri_pxErr, Pri_pyErr, Pri_pzErr, Pri_ptErr);

    // Kaon branches
    X_One_Tree_->Branch("Phi_K_1_px", &Phi_K_1_px); X_One_Tree_->Branch("Phi_K_1_py", &Phi_K_1_py);
    X_One_Tree_->Branch("Phi_K_1_pz", &Phi_K_1_pz); X_One_Tree_->Branch("Phi_K_1_phi", &Phi_K_1_phi);
    X_One_Tree_->Branch("Phi_K_1_eta", &Phi_K_1_eta); X_One_Tree_->Branch("Phi_K_1_pt", &Phi_K_1_pt);
    X_One_Tree_->Branch("Phi_K_1_fromPV", &Phi_K_1_fromPV);
    X_One_Tree_->Branch("Phi_K_1_pvAssocQuality", &Phi_K_1_pvAssocQuality);
    X_One_Tree_->Branch("Phi_K_2_px", &Phi_K_2_px); X_One_Tree_->Branch("Phi_K_2_py", &Phi_K_2_py);
    X_One_Tree_->Branch("Phi_K_2_pz", &Phi_K_2_pz); X_One_Tree_->Branch("Phi_K_2_phi", &Phi_K_2_phi);
    X_One_Tree_->Branch("Phi_K_2_eta", &Phi_K_2_eta); X_One_Tree_->Branch("Phi_K_2_pt", &Phi_K_2_pt);
    X_One_Tree_->Branch("Phi_K_2_fromPV", &Phi_K_2_fromPV);
    X_One_Tree_->Branch("Phi_K_2_pvAssocQuality", &Phi_K_2_pvAssocQuality);

    // Upsilon branches
    branchReso("Ups", Ups_mass, Ups_massErr, Ups_massDiff,
               Ups_ctau, Ups_ctauErr, Ups_Chi2, Ups_ndof, Ups_VtxProb,
               Ups_px, Ups_py, Ups_pz, Ups_phi, Ups_eta, Ups_pt,
               Ups_pxErr, Ups_pyErr, Ups_pzErr, Ups_ptErr);
    X_One_Tree_->Branch("Ups_mu_1_Idx", &Ups_mu_1_Idx);
    X_One_Tree_->Branch("Ups_mu_2_Idx", &Ups_mu_2_Idx);

    // MC branches
    if (doMC) {
        // New flat gen-particle branches
        X_One_Tree_->Branch("MC_GenPart_pdgId", &MC_GenPart_pdgId);
        X_One_Tree_->Branch("MC_GenPart_status", &MC_GenPart_status);
        X_One_Tree_->Branch("MC_GenPart_motherPdgId", &MC_GenPart_motherPdgId);
        X_One_Tree_->Branch("MC_GenPart_handleIndex", &MC_GenPart_handleIndex);
        X_One_Tree_->Branch("MC_GenPart_px", &MC_GenPart_px);
        X_One_Tree_->Branch("MC_GenPart_py", &MC_GenPart_py);
        X_One_Tree_->Branch("MC_GenPart_pz", &MC_GenPart_pz);
        X_One_Tree_->Branch("MC_GenPart_mass", &MC_GenPart_mass);
        X_One_Tree_->Branch("MC_GenPart_pt", &MC_GenPart_pt);
        X_One_Tree_->Branch("MC_GenPart_eta", &MC_GenPart_eta);
        X_One_Tree_->Branch("MC_GenPart_phi", &MC_GenPart_phi);
    }
}

/*****************************************************************************
 * Muon-track matching: Method 1 - First pass below threshold (current method)
 *****************************************************************************/
MultiLepPAT::MuonTrackMatchResult
MultiLepPAT::matchMuonToTrack_FirstPass(
    const pat::Muon* muon,
    std::vector<edm::View<pat::PackedCandidate>::const_iterator>& tracks,
    double momentumRelDiffThr)
{
    MuonTrackMatchResult result;
    result.matched = false;
    result.fromPV = 0;
    result.pvAssocQuality = 0;
    result.pdgId = 0;
    result.nCandidates = 0;
    result.relDiff = -9;
    result.dz = -9;

    double muMomentum = std::sqrt(
        muon->px()*muon->px() + muon->py()*muon->py() + muon->pz()*muon->pz());

    for (auto iTrackfID = tracks.begin(); iTrackfID != tracks.end(); ++iTrackfID) {
        try {
            if (muon->track().isNull()) continue;
            auto iTrackf = *iTrackfID;
            result.nCandidates++;

            if (fabs(iTrackf->px() - muon->px()) < momentumRelDiffThr * muMomentum &&
                fabs(iTrackf->py() - muon->py()) < momentumRelDiffThr * muMomentum &&
                fabs(iTrackf->pz() - muon->pz()) < momentumRelDiffThr * muMomentum &&
                iTrackf->charge() == muon->charge()) {

                tracks.erase(iTrackfID);
                result.matched = true;
                result.fromPV = iTrackf->fromPV();
                result.pvAssocQuality = iTrackf->pvAssociationQuality();
                result.pdgId = iTrackf->pdgId();
                result.relDiff = fabs(iTrackf->p() - muMomentum) / muMomentum;
                return result;
            }
        } catch (...) { continue; }
    }
    return result;
}

/*****************************************************************************
 * Muon-track matching: Method 2 - Least momentum difference
 *****************************************************************************/
MultiLepPAT::MuonTrackMatchResult
MultiLepPAT::matchMuonToTrack_LeastDiff(
    const pat::Muon* muon,
    std::vector<edm::View<pat::PackedCandidate>::const_iterator>& tracks)
{
    MuonTrackMatchResult result;
    result.matched = false;
    result.fromPV = 0;
    result.pvAssocQuality = 0;
    result.pdgId = 0;
    result.nCandidates = 0;
    result.relDiff = 999.0;
    result.dz = -9;

    double muMomentum = std::sqrt(
        muon->px()*muon->px() + muon->py()*muon->py() + muon->pz()*muon->pz());
    auto bestIter = tracks.end();

    for (auto iTrackfID = tracks.begin(); iTrackfID != tracks.end(); ++iTrackfID) {
        try {
            if (muon->track().isNull()) continue;
            auto iTrackf = *iTrackfID;
            result.nCandidates++;

            if (iTrackf->charge() != muon->charge()) continue;

            double relDiff = fabs(iTrackf->p() - muMomentum) / muMomentum;

            if (relDiff < result.relDiff) {
                result.relDiff = relDiff;
                bestIter = iTrackfID;
                result.fromPV = iTrackf->fromPV();
                result.pvAssocQuality = iTrackf->pvAssociationQuality();
                result.pdgId = iTrackf->pdgId();
            }
        } catch (...) { continue; }
    }

    if (bestIter != tracks.end()) {
        tracks.erase(bestIter);
        result.matched = true;
    }

    return result;
}

/*****************************************************************************
 * Muon-track matching: Method 3 - Add dZ position consistency
 *****************************************************************************/
MultiLepPAT::MuonTrackMatchResult
MultiLepPAT::matchMuonToTrack_AddDz(
    const pat::Muon* muon,
    std::vector<edm::View<pat::PackedCandidate>::const_iterator>& tracks,
    const reco::Vertex& primaryV,
    double dzMax)
{
    MuonTrackMatchResult result;
    result.matched = false;
    result.fromPV = 0;
    result.pvAssocQuality = 0;
    result.pdgId = 0;
    result.nCandidates = 0;
    result.relDiff = 999.0;
    result.dz = 999.0;

    double muMomentum = std::sqrt(
        muon->px()*muon->px() + muon->py()*muon->py() + muon->pz()*muon->pz());
    auto bestIter = tracks.end();

    double muDz = -999;
    try {
        auto bestTrk = muon->muonBestTrack();
        if (bestTrk.isAvailable()) {
            muDz = bestTrk->dz(primaryV.position());
        }
    } catch (...) { return result; }

    if (muDz < -900) return result;

    for (auto iTrackfID = tracks.begin(); iTrackfID != tracks.end(); ++iTrackfID) {
        try {
            if (muon->track().isNull()) continue;
            auto iTrackf = *iTrackfID;
            result.nCandidates++;

            if (iTrackf->charge() != muon->charge()) continue;

            double trackDz = -999;
            if (iTrackf->hasTrackDetails() && iTrackf->bestTrack() != nullptr) {
                trackDz = iTrackf->bestTrack()->dz(primaryV.position());
            }

            double relDiff = fabs(iTrackf->p() - muMomentum) / muMomentum;

            if (relDiff < 1.0 && fabs(trackDz - muDz) < dzMax) {
                double score = relDiff + fabs(trackDz - muDz) / dzMax;

                if (score < result.relDiff + result.dz) {
                    result.relDiff = relDiff;
                    result.dz = fabs(trackDz - muDz);
                    bestIter = iTrackfID;
                    result.fromPV = iTrackf->fromPV();
                    result.pvAssocQuality = iTrackf->pvAssociationQuality();
                    result.pdgId = iTrackf->pdgId();
                }
            }
        } catch (...) { continue; }
    }

    if (bestIter != tracks.end()) {
        tracks.erase(bestIter);
        result.matched = true;
    }

    return result;
}

/*****************************************************************************
 * Muon-track matching: Method 4 - Sigma-based (using momentum errors)
 *****************************************************************************/
MultiLepPAT::MuonTrackMatchResult
MultiLepPAT::matchMuonToTrack_Sigma(
    const pat::Muon* muon,
    std::vector<edm::View<pat::PackedCandidate>::const_iterator>& tracks,
    double sigmaThr)
{
    MuonTrackMatchResult result;
    result.matched = false;
    result.fromPV = 0;
    result.pvAssocQuality = 0;
    result.pdgId = 0;
    result.nCandidates = 0;
    result.relDiff = 999.0;
    result.dz = -9;

    double muPxErr = -9, muPyErr = -9, muPzErr = -9;
    try {
        reco::TrackRef muTrack = muon->track();
        if (!muTrack.isNull() && muTrack.isAvailable()) {
            auto cov = muTrack->covariance();
            muPxErr = std::sqrt(std::max(cov(0,0), 1e-6));
            muPyErr = std::sqrt(std::max(cov(1,1), 1e-6));
            muPzErr = std::sqrt(std::max(cov(2,2), 1e-6));
        }
    } catch (...) { }

    if (muPxErr < 0) return result;

    double muPx = muon->px(), muPy = muon->py(), muPz = muon->pz();
    auto bestIter = tracks.end();
    double bestChi2 = 999.0;

    for (auto iTrackfID = tracks.begin(); iTrackfID != tracks.end(); ++iTrackfID) {
        try {
            if (muon->track().isNull()) continue;
            auto iTrackf = *iTrackfID;
            result.nCandidates++;

            if (iTrackf->charge() != muon->charge()) continue;

            double dpx = iTrackf->px() - muPx;
            double dpy = iTrackf->py() - muPy;
            double dpz = iTrackf->pz() - muPz;

            double chi2 = (dpx*dpx)/(muPxErr*muPxErr) +
                          (dpy*dpy)/(muPyErr*muPyErr) +
                          (dpz*dpz)/(muPzErr*muPzErr);

            if (chi2 < sigmaThr * sigmaThr && chi2 < bestChi2) {
                bestChi2 = chi2;
                bestIter = iTrackfID;
                result.fromPV = iTrackf->fromPV();
                result.pvAssocQuality = iTrackf->pvAssociationQuality();
                result.pdgId = iTrackf->pdgId();
                result.relDiff = chi2;
            }
        } catch (...) { continue; }
    }

    if (bestIter != tracks.end()) {
        tracks.erase(bestIter);
        result.matched = true;
    }

    return result;
}

void MultiLepPAT::endJob()
{
    X_One_Tree_->GetDirectory()->cd();
    X_One_Tree_->Write();
}

// define this as a plug-in
DEFINE_FWK_MODULE(MultiLepPAT);
