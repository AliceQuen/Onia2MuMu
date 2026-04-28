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

      X6900_mass(0), X6900_VtxProb(0), X6900_massErr(0),
      X6900_pt(0), X6900_pz(0), X6900_absEta(0),
      X6900_px(0), X6900_py(0),

      Psi2S_mass_raw(0), Psi2S_mass(0), Psi2S_VtxProb(0), Psi2S_massErr(0),
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
      dR_Jpsi1_X6900(0), dR_Jpsi2_X6900(0),
      dR_X6900_pi1(0), dR_X6900_pi2(0),
      dR_X6900_mu1(0), dR_X6900_mu2(0), dR_X6900_mu3(0), dR_X6900_mu4(0),
      dR_Psi2S_X6900(0), dR_Psi2S_Jpsi1(0), dR_Psi2S_Jpsi2(0),
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
  std::vector<edm::View<pat::PackedCandidate>::const_iterator> nonMuonPionTrack;

  // Pre-selection for tracks: filter out candidates that don't meet the cuts
  // Selection criteria from pre-cuts.md:
  // 1. pt > 0.5 && absEta < 2.4
  // 2. sigma_pt / pt < 0.1
  // 3. Nhits >= 11
  // 4. chi^2 / ndf < 0.18
  nonMuonPionTrack.clear();
  for (edm::View<pat::PackedCandidate>::const_iterator iTrackc =
           theTrackHandle->begin(); // MINIAOD
       iTrackc != theTrackHandle->end(); ++iTrackc) {
    // Apply track pre-selection criteria
    const reco::Track* track = iTrackc->bestTrack();
    if (!track) {
      continue;
    }
    // 1. Is high purity muon
    if (!(track->quality(reco::TrackBase::qualityByName("highPurity")))){
      continue;
    }
    // 2. pt > 0.5 && absEta < 2.4
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
    if (Nhits < track_nhits_cut){
      continue;
    }
    // 5. chi^2 / ndf < 0.18
   if (track->normalizedChi2() / Nhits >= track_chi2_over_ndf_cut){
      continue;
    }
    // All selection criteria passed
    nonMuonPionTrack.push_back(iTrackc);
  }

  #if DEBUG == 1
  std::cout << "[DEBUG] Track pre-selection: total=" << theTrackHandle->size() << ", selected=" << nonMuonPionTrack.size() << std::endl;
  #endif

  // Check if we have enough tracks after pre-selection
  if (nonMuonPionTrack.size() < 2) {
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: Less than 2 valid non-muon tracks after pre-selection (total found " << nonMuonPionTrack.size() << ")" << " (line " << __LINE__ << ")" << std::endl;
    #endif
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    total_return_++;
    #endif
    return;
  }
  
  // Find and remove muon Tracks from PionTracks, collect filtered muons with pre-selection
  // Clear the muon filter matches vector at the start of each event
  muonFilterMatches.clear();
  std::vector<edm::View<pat::Muon>::const_iterator> validMuons;

  // Pre-selection for muons: filter out candidates that don't meet the cuts
  // Selection criteria from pre-cuts.md:
  // 1. all muon candidates must be soft muons
  // 2. (pt > 3.5 && absEta < 1.2) || (pt > (5.47 - 1.89 * absEta) && 1.2 < absEta < 2.1) || (pt > 1.5 && 2.1 < absEta < 2.4)
  for (edm::View<pat::Muon>::const_iterator iMuonP =
           thePATMuonHandle->begin(); //  MINIAOD
       iMuonP != thePATMuonHandle->end(); ++iMuonP) {
    // Initialize filter match flag for this muon
    bool muonHasFilterMatch = false;

    // Apply muon pre-selection criteria
    // 1. Must be soft muon
    if (!iMuonP->isSoftMuon(thePrimaryV)) {
      #if DEBUG == 2
      continue_counts_[__LINE__]++;
      total_return_++;
      #endif
      continue;
    }
    // 2. pt and eta requirements: (pt > 3.5 && absEta < 1.2) || (pt > (5.47 - 1.89*absEta) && 1.2 < absEta < 2.1) || (pt > 1.5 && 2.1 < absEta < 2.4)
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
      total_return_++;
      #endif
      continue;
    }

    // Remove muon track from pion track list
    if (iMuonP->track().isNull()) {
      #if DEBUG == 2
      continue_counts_[__LINE__]++;
      total_return_++;
      #endif
      #if DEBUG == 1
      std::cout << "[DEBUG] Continue encountered in analyze(), location: Muon track is null during removal (line " << __LINE__ << ")" << std::endl;
      #endif
      continue;
    }
    // All muon selection criteria passed
    validMuons.push_back(iMuonP);
    for (std::vector<edm::View<pat::PackedCandidate>::const_iterator>::
             const_iterator iTrackfID = nonMuonPionTrack.begin(); // MINIAOD
         iTrackfID != nonMuonPionTrack.end(); ++iTrackfID) {
      edm::View<pat::PackedCandidate>::const_iterator iTrackf = *(iTrackfID);
      if (fabs(iTrackf->px() - iMuonP->track()->px()) < std::numeric_limits<float>::epsilon() &&
         fabs(iTrackf->py() - iMuonP->track()->py()) < std::numeric_limits<float>::epsilon() &&
        fabs(iTrackf->pz() - iMuonP->track()->pz()) < std::numeric_limits<float>::epsilon()) {
        nonMuonPionTrack.erase(iTrackfID);
        iTrackfID = iTrackfID - 1;
      }
    }

    // Check filter matching with HLT filters
    for (unsigned int JpsiFilt = 0; JpsiFilt < FiltersForJpsi_.size();
         JpsiFilt++) {
      if (hltresults.isValid()) {
        pat::TriggerObjectStandAlone *tempTriggerObject = nullptr;
        for (auto i = iMuonP->triggerObjectMatches().begin();
             i != iMuonP->triggerObjectMatches().end(); i++) {
          tempTriggerObject = new pat::TriggerObjectStandAlone(*i);
          tempTriggerObject->unpackFilterLabels(iEvent, *hltresults);
          if (tempTriggerObject->hasFilterLabel(FiltersForJpsi_[JpsiFilt])) {
            muonHasFilterMatch = true;
          }
          delete tempTriggerObject;
        }
      }
    }
    
    // Store the filter match result for this muon
    muonFilterMatches.push_back(muonHasFilterMatch);
  }

  #if DEBUG == 1
  std::cout << "[DEBUG] Muon pre-selection: total=" << thePATMuonHandle->size() << ", selected=" << validMuons.size() << std::endl;
  #endif

  #if DEBUG == 1
  auto end_track = std::chrono::high_resolution_clock::now();
  auto duration_track = std::chrono::duration_cast<std::chrono::milliseconds>(end_track - start_track).count();
  std::cout << "[DEBUG] CheckMuonPionTrack 运行时间: " << duration_track << " ms" << std::endl;
  #endif

  // Check if we have enough muons after pre-selection
  // Because muons are always few, we check for the charge as well
  // This maybe replaced by combine candidates first
  unsigned int plus = 0;
  unsigned int minu = 0;
  for (const auto &iMuon :validMuons)
  {
    if (iMuon->charge() > 0) {
      plus++;
    } else {
      minu++;
    }
  }
  if (plus < 2 || minu < 2 || plus + minu < 4) {
    #if DEBUG == 2
    return_counts_[__LINE__]++;
    total_return_++;
    #endif
    #if DEBUG == 1
    std::cout << "[DEBUG] Return encountered in analyze(), location: Less than 4 valid muons after pre-selection (total selected " << validMuons.size() << ")" << " (line " << __LINE__ << ")" << std::endl;
    #endif
    return;
  }

  // Start processing muon and track fit
#if DEBUG == 1
  auto start_fit = std::chrono::high_resolution_clock::now();
#endif
  // Iterate over pre-selected valid muons only (pre-selection already applied above)
  for (std::vector<edm::View<pat::Muon>::const_iterator>::
           const_iterator iMuon1ID = validMuons.begin();
           iMuon1ID != validMuons.end(); ++iMuon1ID) {
    edm::View<pat::Muon>::const_iterator iMuon1 = *(iMuon1ID);
    TrackRef muTrack1 = iMuon1->track();

    // next check for mu2 in pre-selected valid muons
    for (std::vector<edm::View<pat::Muon>::const_iterator>::
             const_iterator iMuon2ID = iMuon1ID + 1;
         iMuon2ID != validMuons.end(); ++iMuon2ID) {
      edm::View<pat::Muon>::const_iterator iMuon2 = *(iMuon2ID);
      TrackRef muTrack2 = iMuon2->track();
      if ((iMuon1->charge() + iMuon2->charge()) != 0) {
#if DEBUG == 2
        continue_counts_[__LINE__]++;
        total_continue_++;
#endif
        #if DEBUG == 1
        std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi1 charge sum not zero (line " << __LINE__ << ")" << std::endl;
        #endif
        continue;
      }
      if (!(1. < (iMuon1->p4() + iMuon2->p4()).mass() &&
            (iMuon1->p4() + iMuon2->p4()).mass() < 4.5)) {
        #if DEBUG == 2
        continue_counts_[__LINE__]++;
        total_continue_++;
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
      P4_mu1.SetPtEtaPhiM(iMuon1->track()->pt(), iMuon1->track()->eta(),
                          iMuon1->track()->phi(), MU_MASS);
      TLorentzVector P4_mu2;
      P4_mu2.SetPtEtaPhiM(iMuon2->track()->pt(), iMuon2->track()->eta(),
                          iMuon2->track()->phi(), MU_MASS);
      // mu3mu4(X6900->Jpsi) - iterate over pre-selected valid muons only
      for (std::vector<edm::View<pat::Muon>::const_iterator>::
               const_iterator iMuon3ID = validMuons.begin();
           iMuon3ID != validMuons.end(); ++iMuon3ID) {
        edm::View<pat::Muon>::const_iterator iMuon3 = *(iMuon3ID);
        if (iMuon3 == iMuon1 || iMuon3 == iMuon2) {
#if DEBUG == 2
          continue_counts_[__LINE__]++;
          total_continue_++;
#endif
          #if DEBUG == 1
          std::cout << "[DEBUG] Continue encountered in analyze(), location: mu3 is same as mu1 or mu2 (line " << __LINE__ << ")" << std::endl;
          #endif
          continue;
        }
        TrackRef muTrack3 = iMuon3->track();

        for (std::vector<edm::View<pat::Muon>::const_iterator>::
                 const_iterator iMuon4ID =
                 iMuon3ID + 1; // MINIAOD
             iMuon4ID != validMuons.end(); ++iMuon4ID) {
          edm::View<pat::Muon>::const_iterator iMuon4 = *(iMuon4ID);
          if (iMuon4 == iMuon1 || iMuon4 == iMuon2) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
#endif
            #if DEBUG == 1
            std::cout << "[DEBUG] Continue encountered in analyze(), location: mu4 is same as mu1 or mu2 (line " << __LINE__ << ")" << std::endl;
            #endif
            continue;
          }
          TrackRef muTrack4 = iMuon4->track();
          if ((iMuon3->charge() + iMuon4->charge()) != 0) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
#endif
            #if DEBUG == 1
            std::cout << "[DEBUG] Continue encountered in analyze(), location: Jpsi2 charge sum not zero (line " << __LINE__ << ")" << std::endl;
            #endif
            continue;
          }
          if (!(1. < (iMuon3->p4() + iMuon4->p4()).mass() &&
                (iMuon3->p4() + iMuon4->p4()).mass() < 4.)) {
#if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
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

          // Psi2S or X3872 (4mu 2pi) candiate
          // Begin processing tracks
          for (std::vector<edm::View<pat::PackedCandidate>::const_iterator>::
                   const_iterator iTrack1ID =
                       nonMuonPionTrack.begin(); // MINIAOD
               iTrack1ID != nonMuonPionTrack.end(); ++iTrack1ID) {
            edm::View<pat::PackedCandidate>::const_iterator iTrack1 =
                *(iTrack1ID);
            for (std::vector<edm::View<pat::PackedCandidate>::const_iterator>::
                     const_iterator iTrack2ID = iTrack1ID + 1; // MINIAOD
                 iTrack2ID != nonMuonPionTrack.end(); ++iTrack2ID) {
              edm::View<pat::PackedCandidate>::const_iterator iTrack2 =
                  *(iTrack2ID);
              if ((iTrack1->charge() + iTrack2->charge()) != 0) {
                #if DEBUG == 2
                continue_counts_[__LINE__]++;
                total_continue_++;
                #endif
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Pion pair same charge (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }
              if (!iTrack1->hasTrackDetails() || iTrack1->charge() == 0) {
#if DEBUG == 2
                continue_counts_[__LINE__]++;
                total_continue_++;
#endif
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Track1 has no valid track details (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }
              if (!iTrack2->hasTrackDetails() || iTrack2->charge() == 0) {
#if DEBUG == 2
                continue_counts_[__LINE__]++;
                total_continue_++;
#endif
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Track2 has no valid track details (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }
              
              // Pre-selection already applied upfront when building nonMuonPionTrack
              const reco::Track* track1 = iTrack1->bestTrack();
              const reco::Track* track2 = iTrack2->bestTrack();
              if (!track1 || !track2) {
#if DEBUG == 2
                continue_counts_[__LINE__]++;
                total_continue_++;
#endif
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Best track pointer null (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }

              TLorentzVector P4_Track1, P4_Track2, P4_Jpsipipi;
              P4_Track1.SetPtEtaPhiM(iTrack1->pt(), iTrack1->eta(),
                                     iTrack1->phi(), PI_MASS);
              P4_Track2.SetPtEtaPhiM(iTrack2->pt(), iTrack2->eta(),
                                     iTrack2->phi(), PI_MASS);
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

              TransientTrack trackTT1(*(iTrack1->bestTrack()),
                                      &(bFieldHandle));
              TransientTrack trackTT2(*(iTrack2->bestTrack()),
                                      &(bFieldHandle));
              KinematicParticleFactoryFromTransientTrack JPiPiFactory;
              // The mass of a muon and the insignificant mass sigma
              // to avoid singularities in the covariance matrix.
              ParticleMass pion_mass = PI_MASS; // pdg mass
              float pion_sigma = PI_MASSERR;
              // initial chi2 and ndf before kinematic fits.
              float chi = 0.;
              float ndf = 0.;

              // first fit to have a psi2s
              vector<RefCountedKinematicParticle> JPiPiParticles;
              JPiPiParticles.push_back(JPiPiFactory.particle(
                  trackTT1, pion_mass, chi, ndf, pion_sigma));
              JPiPiParticles.push_back(JPiPiFactory.particle(
                  trackTT2, pion_mass, chi, ndf, pion_sigma));
              JPiPiParticles.push_back(pmumuFactory.particle(
                  muon1TT, muon_mass, chi, ndf, muon_sigma));
              JPiPiParticles.push_back(pmumuFactory.particle(
                  muon2TT, muon_mass, chi, ndf, muon_sigma));

              KinematicParticleVertexFitter JPiPi_fitter;
              RefCountedKinematicTree JPiPiVertexFitTree;
              Error_t = false;
              try {
                JPiPiVertexFitTree = JPiPi_fitter.fit(JPiPiParticles);
              } catch (...) {
                Error_t = true;
              }
              if (Error_t || !(JPiPiVertexFitTree->isValid())) {
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Psi2S vertex fit failed (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }
              JPiPiVertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle JPiPi_vFit_noMC =
                  JPiPiVertexFitTree->currentParticle();
              RefCountedKinematicVertex JPiPi_vFit_vertex_noMC =
                  JPiPiVertexFitTree->currentDecayVertex();

              double JPiPi_vtxprob = ChiSquaredProbability(
                  (double)(JPiPi_vFit_vertex_noMC->chiSquared()),
                  (double)(JPiPi_vFit_vertex_noMC->degreesOfFreedom()));
              if (JPiPi_vFit_noMC->currentState().mass() > 4.5) {
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Psi2S mass exceeds 4.5 (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }
              
              // Psi2S pre-selection according to pre-cuts.md
              // 1. VtxProb > 0.5%
              if (JPiPi_vtxprob < psi2s_vtxprob_cut) {
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Psi2S vertex probability below cut (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }
              
              // 2. pt > 4
              double px = JPiPi_vFit_noMC->currentState().kinematicParameters().momentum().x();
              double py = JPiPi_vFit_noMC->currentState().kinematicParameters().momentum().y();
              double psi2s_pt = sqrt(px*px + py*py);
              if (psi2s_pt <= psi2s_pt_cut) {
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: Psi2S pT below cut (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }

              // fit the 6 track together to a vertex
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

              KinematicParticleVertexFitter X_fitter;
              RefCountedKinematicTree X_VertexFitTree;
              Error_t = false;
              try {
                X_VertexFitTree = X_fitter.fit(X_Particles);
              } catch (...) {
                Error_t = true;
              }
              if (Error_t || !(X_VertexFitTree->isValid())) {
                #if DEBUG == 1
                std::cout << "[DEBUG] Continue encountered in analyze(), location: X(6900) vertex fit failed (line " << __LINE__ << ")" << std::endl;
                #endif
                continue;
              }
              X_VertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle X_vFit_noMC =
                  X_VertexFitTree->currentParticle();
              RefCountedKinematicVertex X_vFit_vertex_noMC =
                  X_VertexFitTree->currentDecayVertex();
              KinematicParameters X_kPara =
                  X_vFit_noMC->currentState().kinematicParameters();

              double X_vtxprob = ChiSquaredProbability(
                  (double)(X_vFit_vertex_noMC->chiSquared()),
                  (double)(X_vFit_vertex_noMC->degreesOfFreedom()));

              X_VertexFitTree->movePointerToTheFirstChild();
              RefCountedKinematicParticle X_pi1 =
                  X_VertexFitTree->currentParticle();
              X_VertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle X_pi2 =
                  X_VertexFitTree->currentParticle();

              KinematicParameters X_pi1_KP =
                  X_pi1->currentState().kinematicParameters();
              KinematicParameters X_pi2_KP =
                  X_pi2->currentState().kinematicParameters();
              
              // Start calculation variables to be stored
              X6900_mass = X_vFit_noMC->currentState().mass();
              X6900_VtxProb = X_vtxprob;
              X6900_px = X_kPara.momentum().x();
              X6900_py = X_kPara.momentum().y();
              X6900_pz = X_kPara.momentum().z();
              if (X_vFit_noMC->currentState()
                      .kinematicParametersError()
                      .matrix()(6, 6) > 0) {
                X6900_massErr = sqrt(X_vFit_noMC->currentState()
                                              .kinematicParametersError()
                                              .matrix()(6, 6));
              } else {
                X6900_massErr = -9;
              }

              ROOT::Math::PxPyPzMVector X6900_vec(X6900_px, X6900_py, X6900_pz, X6900_mass);
              X6900_pt = X6900_vec.Pt();
              X6900_absEta = fabs(X6900_vec.Eta());

              Psi2S_mass_raw = JPiPi_vFit_noMC->currentState().mass();
              Psi2S_VtxProb = JPiPi_vtxprob;
              Psi2S_px = JPiPi_vFit_noMC->currentState()
                                        .kinematicParameters()
                                        .momentum()
                                        .x();
              Psi2S_py = JPiPi_vFit_noMC->currentState()
                                        .kinematicParameters()
                                        .momentum()
                                        .y();
              Psi2S_pz = JPiPi_vFit_noMC->currentState()
                                        .kinematicParameters()
                                        .momentum()
                                        .z();
              if (JPiPi_vFit_noMC->currentState()
                      .kinematicParametersError()
                      .matrix()(6, 6) > 0) {
                Psi2S_massErr = sqrt(JPiPi_vFit_noMC->currentState()
                                                    .kinematicParametersError()
                                                    .matrix()(6, 6));
              } else {
                Psi2S_massErr = -9;
              }

              ROOT::Math::PxPyPzMVector Psi2S_vec(Psi2S_px, Psi2S_py, Psi2S_pz, Psi2S_mass_raw);
              Psi2S_pt = Psi2S_vec.Pt();
              Psi2S_absEta = fabs(Psi2S_vec.Eta());

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


  X_One_Tree_->Branch("X6900_mass", &X6900_mass, "X6900_mass/F");
  X_One_Tree_->Branch("X6900_VtxProb", &X6900_VtxProb, "X6900_VtxProb/F");
  X_One_Tree_->Branch("X6900_massErr", &X6900_massErr, "X6900_massErr/F");
  X_One_Tree_->Branch("X6900_pt", &X6900_pt, "X6900_pt/F");
  X_One_Tree_->Branch("X6900_pz", &X6900_pz, "X6900_pz/F");
  X_One_Tree_->Branch("X6900_absEta", &X6900_absEta, "X6900_absEta/F");

  X_One_Tree_->Branch("Psi2S_mass_raw", &Psi2S_mass_raw, "Psi2S_mass_raw/F");
  X_One_Tree_->Branch("Psi2S_mass", &Psi2S_mass, "Psi2S_mass/F");
  X_One_Tree_->Branch("Psi2S_VtxProb", &Psi2S_VtxProb, "Psi2S_VtxProb/F");
  X_One_Tree_->Branch("Psi2S_massErr", &Psi2S_massErr, "Psi2S_massErr/F");
  X_One_Tree_->Branch("Psi2S_pt", &Psi2S_pt, "Psi2S_pt/F");
  X_One_Tree_->Branch("Psi2S_pz", &Psi2S_pz, "Psi2S_pz/F");
  X_One_Tree_->Branch("Psi2S_absEta", &Psi2S_absEta, "Psi2S_absEta/F");

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
  X_One_Tree_->Branch("dR_Jpsi1_X6900", &dR_Jpsi1_X6900, "dR_Jpsi1_X6900/F");
  X_One_Tree_->Branch("dR_Jpsi2_X6900", &dR_Jpsi2_X6900, "dR_Jpsi2_X6900/F");
  X_One_Tree_->Branch("dR_X6900_pi1", &dR_X6900_pi1, "dR_X6900_pi1/F");
  X_One_Tree_->Branch("dR_X6900_pi2", &dR_X6900_pi2, "dR_X6900_pi2/F");
  X_One_Tree_->Branch("dR_X6900_mu1", &dR_X6900_mu1, "dR_X6900_mu1/F");
  X_One_Tree_->Branch("dR_X6900_mu2", &dR_X6900_mu2, "dR_X6900_mu2/F");
  X_One_Tree_->Branch("dR_X6900_mu3", &dR_X6900_mu3, "dR_X6900_mu3/F");
  X_One_Tree_->Branch("dR_X6900_mu4", &dR_X6900_mu4, "dR_X6900_mu4/F");
  X_One_Tree_->Branch("dR_Psi2S_X6900", &dR_Psi2S_X6900, "dR_Psi2S_X6900/F");
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
  
  // Reset X(6900) candidate variables
  X6900_mass = -999.0;
  X6900_VtxProb = -999.0;
  X6900_massErr = -999.0;
  X6900_pt = -999.0;
  X6900_pz = -999.0;
  X6900_absEta = -999.0;
  X6900_px = -999.0;
  X6900_py = -999.0;
  
  // Reset Psi(2S) candidate variables
  Psi2S_mass_raw = -999.0;
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
  dR_Jpsi1_X6900 = -999.0;
  dR_Jpsi2_X6900 = -999.0;
  dR_X6900_pi1 = -999.0;
  dR_X6900_pi2 = -999.0;
  dR_X6900_mu1 = -999.0;
  dR_X6900_mu2 = -999.0;
  dR_X6900_mu3 = -999.0;
  dR_X6900_mu4 = -999.0;
  dR_Psi2S_X6900 = -999.0;
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
