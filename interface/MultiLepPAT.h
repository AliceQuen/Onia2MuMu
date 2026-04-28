// -*- C++ -*-
//
// Package:    MultiLepPAT
// Class:      MultiLepPAT
// 
/**\class MultiLepPAT MultiLepPAT.cc myAnalyzers/MultiLepPAT/src/MultiLepPAT.cc

 Description: <one line class summary>
Make rootTuple for JPsiPsi2S reconstruction

 Implementation:
     <Notes on implementation>
*/
//
// Original Author: QIN Junkai 
//
//

#ifndef _MultiLepPAT_h
#define _MultiLepPAT_h

#define DEBUG 0

// system include files
#include <memory>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // xining MINIAODtest

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
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
//#include "RecoVertex/V0Producer/interface/V0Producer.h"
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

//
// class decleration
//

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class MultiLepPAT : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MultiLepPAT(const ParameterSet&);
  ~MultiLepPAT();

  
private:
  virtual void beginJob() ;
  virtual void beginRun(Run const & iRun, EventSetup const& iSetup);
  virtual void analyze(const Event&, const EventSetup&);
  virtual void endJob() ;
  
  // Reset function to initialize/reset all variables for each event
  void resetVariables();

 
//add token here
  edm::EDGetTokenT<BeamSpot> gtbeamspotToken_;
  edm::EDGetTokenT<VertexCollection> gtprimaryVtxToken_;
  edm::EDGetTokenT<edm::View<pat::Muon> > gtpatmuonToken_; // MINIAOD
  edm::EDGetTokenT<edm::TriggerResults> gttriggerToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > trackToken_; // MINIAOD

//  InputTag inputGEN_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>  magneticFieldToken_;
  vector<string>      TriggersForJpsi_;
  vector<string>      FiltersForJpsi_;
  
  // Sorted by length ascending for faster matching - shorter patterns checked first
  vector<string>      TriggersForJpsi_sorted_;
  size_t              min_trigger_len_;

  int JpsiMatchTrig[50];
  vector<bool>        muonFilterMatches;

  TTree* X_One_Tree_;

  unsigned int        runNum, evtNum, lumiNum;
  unsigned int        nGoodPrimVtx;


  float X6900_mass, X6900_VtxProb, X6900_massErr;
  float X6900_pt, X6900_pz, X6900_absEta;
  float X6900_px, X6900_py;

  float Psi2S_mass_raw, Psi2S_mass, Psi2S_VtxProb, Psi2S_massErr;
  float Psi2S_pt, Psi2S_pz, Psi2S_absEta;
  float Psi2S_px, Psi2S_py;

  float Jpsi1_mass, Jpsi1_VtxProb, Jpsi1_massErr;
  float Jpsi1_pt, Jpsi1_pz, Jpsi1_absEta;
  float Jpsi1_px, Jpsi1_py;

  float Jpsi2_mass, Jpsi2_VtxProb, Jpsi2_massErr;
  float Jpsi2_pt, Jpsi2_pz, Jpsi2_absEta;
  float Jpsi2_px, Jpsi2_py;

  float mu1_pt, mu1_pz, mu1_absEta;
  float mu1_px, mu1_py;
  float mu1_trackIso;
  float mu1_d0BS, mu1_d0EBS, mu1_d3dBS, mu1_d3dEBS;
  float mu1_d0PV, mu1_d0EPV, mu1_dzPV, mu1_dzEPV;
  float mu1_charge;

  float mu2_pt, mu2_pz, mu2_absEta;
  float mu2_px, mu2_py;
  float mu2_trackIso;
  float mu2_d0BS, mu2_d0EBS, mu2_d3dBS, mu2_d3dEBS;
  float mu2_d0PV, mu2_d0EPV, mu2_dzPV, mu2_dzEPV;
  float mu2_charge;

  float mu3_pt, mu3_pz, mu3_absEta;
  float mu3_px, mu3_py;
  float mu3_trackIso;
  float mu3_d0BS, mu3_d0EBS, mu3_d3dBS, mu3_d3dEBS;
  float mu3_d0PV, mu3_d0EPV, mu3_dzPV, mu3_dzEPV;
  float mu3_charge;

  float mu4_pt, mu4_pz, mu4_absEta;
  float mu4_px, mu4_py;
  float mu4_trackIso;
  float mu4_d0BS, mu4_d0EBS, mu4_d3dBS, mu4_d3dEBS;
  float mu4_d0PV, mu4_d0EPV, mu4_dzPV, mu4_dzEPV;
  float mu4_charge;

  int nLooseMuons, nTightMuons, nSoftMuons, nMediumMuons;
  
  bool mu1_hasFilterMatch, mu2_hasFilterMatch, mu3_hasFilterMatch, mu4_hasFilterMatch;

  float pi1_pt, pi1_pz, pi1_absEta;
  float pi1_px, pi1_py;
  float pi2_pt, pi2_pz, pi2_absEta;
  float pi2_px, pi2_py;

  float dR_mu1_mu2, dR_mu3_mu4, dR_pi1_pi2;
  float dR_Jpsi1_X6900, dR_Jpsi2_X6900;
  float dR_X6900_pi1, dR_X6900_pi2;
  float dR_X6900_mu1, dR_X6900_mu2, dR_X6900_mu3, dR_X6900_mu4;
  float dR_Psi2S_X6900, dR_Psi2S_Jpsi1, dR_Psi2S_Jpsi2;
  float dR_Psi2S_pi1, dR_Psi2S_pi2;
  #if DEBUG == 2
  private:
  std::map<int, unsigned long long> continue_counts_;
  std::map<int, unsigned long long> return_counts_;
  unsigned long long total_continue_ = 0;
  unsigned long long total_return_ = 0;
  #endif
};

#endif
