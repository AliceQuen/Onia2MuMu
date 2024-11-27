/******************************************************************************
 *  [File] 
 *      MultiLepPAT.cc 
 *  [Class]      
 *      MultiLepPAT 
 *  [Directory] 
 *      TPS-Onia2MuMu/src/MultiLepPAT.cc
 *  [Description]
 *      Make rootTuple for Jpsi+Jpsi+Upsilon reconstruction from 6 muons
 *  [Implementation]
 *     <Notes on implementation>
 *  [Note]
 *      20240704 [Eric Wang]
 *          Upsilon is abbreviated as "Ups"
 *          PDG ID : Mu -> 13,  K -> 321, Pi -> 211, 
 *                   Jpsi -> 443, Ups -> 553, Phi -> 333
 *          Refraction was conducted, especially to the boolean variables.
 *          Room for improvement: (already tagged)
 *              "TO_ENC"            Encapsulation required
 *              "TO_IMPR_CPP11/17"  Use C++11/17 features for better 
 *                                  readability and efficiency
 *                                  (e.g. auto, range-based for loop)
 *          Github copilot is used for code and annotation completion.
 *      
 *      20240811 [Eric Wang]
 *       - The underlying core of edm::View<T> is actually a std::vector<T>.
 *         This fact may help with writing more efficient code.
 *       - It took me some time to realize that I must sort out the muon pairs
 *         before I can proceed with the JPsi and Upsilon reconstruction and 
 *         vertex matching. It would prove too troublesome to cover all 
 *         combinations of muon pairs in the "multi-layer-for-loop" structure.
 *       
 *      20241029 [Eric Wang and Shi Zhenpeng]
 *       - The code is branched from Eric-dev-run-loose.
 *         This will deal with combined production of Jpsi, Upsilon and Phi.
 *       - Relevant changes will be made at spots marked as "[J-U-P]".
******************************************************************************/

// system include files
#include "TLorentzVector.h"
// user include files
#include "../interface/MultiLepPAT.h"
#include "../interface/VertexReProducer.h"

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
#include "DataFormats/Math/interface/Point3D.h"
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

#include "MagneticField/Engine/interface/MagneticField.h"

//////////This is necessary for lumicalc///////
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParam.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParamRcd.h"

#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <algorithm>
#include <memory>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <regex>

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

#include "../data/TMVAClassification_BDT.class.C"

// about photon
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <boost/foreach.hpp>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // MINIAOD

typedef math::Error<3>::type CovarianceMatrix;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;

// constructors and destructor
MultiLepPAT::MultiLepPAT(const edm::ParameterSet &iConfig)
	: hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults", edm::InputTag("TriggerResults::HLT"))),
	  inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN", edm::InputTag("genParticles"))),
	  magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
	  theTTBuilderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
	  vtxSample(iConfig.getUntrackedParameter<std::string>("VtxSample", std::string("offlinePrimaryVertices"))),
	  doMC(iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", false)),
	  MCParticle(iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443)), // 20443 X, 100443 Psi(2S), 9120443  // X from B
	  doJPsiMassCost(iConfig.getUntrackedParameter<bool>("DoJPsiMassConstraint")),
	  MuPixHits_c(iConfig.getUntrackedParameter<int>("MinNumMuPixHits", 0)),
	  MuSiHits_c(iConfig.getUntrackedParameter<int>("MinNumMuSiHits", 0)),
	  MuNormChi_c(iConfig.getUntrackedParameter<double>("MaxMuNormChi2", 1000)),
	  MuD0_c(iConfig.getUntrackedParameter<double>("MaxMuD0", 1000)),
	  JMaxM_c(iConfig.getUntrackedParameter<double>("MaxJPsiMass", 4)),
	  JMinM_c(iConfig.getUntrackedParameter<double>("MinJPsiMass", 2.2)),
	  PiSiHits_c(iConfig.getUntrackedParameter<int>("MinNumTrSiHits", 0)),
	  MuPt_c(iConfig.getUntrackedParameter<double>("MinMuPt", 0)),
	  JPiPiDR_c(iConfig.getUntrackedParameter<double>("JPsiKKKMaxDR", 1)),
	  XPiPiDR_c(iConfig.getUntrackedParameter<double>("XCandPiPiMaxDR", 1.1)),
	  UseXDr_c(iConfig.getUntrackedParameter<bool>("UseXDr", false)),
	  JPiPiMax_c(iConfig.getUntrackedParameter<double>("JPsiKKKMaxMass", 50)),
	  JPiPiMin_c(iConfig.getUntrackedParameter<double>("JPsiKKKMinMass", 0)),
	  resolveAmbiguity_(iConfig.getUntrackedParameter<bool>("resolvePileUpAmbiguity", true)),
	  addXlessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addXlessPrimaryVertex", true)),
	  TriggersForJpsi_(iConfig.getUntrackedParameter<std::vector<std::string>>("TriggersForJpsi")),
	  FiltersForJpsi_(iConfig.getUntrackedParameter<std::vector<std::string>>("FiltersForJpsi")),
	  TriggersForUpsilon_(iConfig.getUntrackedParameter<std::vector<std::string>>("TriggersForUpsilon")),
	  FiltersForUpsilon_(iConfig.getUntrackedParameter<std::vector<std::string>>("FiltersForUpsilon")),
	  Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output", false)),
	  Chi_Track_(iConfig.getUntrackedParameter<double>("Chi2NDF_Track", 10)),
	  X_One_Tree_(0),

	  runNum(0), evtNum(0), lumiNum(0), nGoodPrimVtx(0),
	  trigRes(0), trigNames(0), L1TT(0), MatchTriggerNames(0),

	  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxChiNorm(0), priVtxChi(0), priVtxCL(0),
	  PriVtxXCorrX(0), PriVtxXCorrY(0), PriVtxXCorrZ(0),
	  PriVtxXCorrEX(0), PriVtxXCorrEY(0), PriVtxXCorrEZ(0), PriVtxXCorrC2(0), PriVtxXCorrCL(0),

	  nMu(0),
	  muPx(0), muPy(0), muPz(0), muD0(0), muD0E(0), muDz(0), muChi2(0), muGlChi2(0), mufHits(0),
	  muFirstBarrel(0), muFirstEndCap(0), muDzVtx(0), muDxyVtx(0),
	  muNDF(0), muGlNDF(0), muPhits(0), muShits(0), muGlMuHits(0), muType(0), muQual(0),
	  muTrack(0), muCharge(0), muIsoratio(0), muIsGoodLooseMuon(0), muIsGoodLooseMuonNew(0),
	  muIsGoodSoftMuonNewIlse(0), muIsGoodSoftMuonNewIlseMod(0), muIsGoodTightMuon(0), muIsJpsiTrigMatch(0), muIsUpsTrigMatch(0), munMatchedSeg(0), muJpsiFilterRes(0),

	  muIsPatLooseMuon(0), muIsPatTightMuon(0), muIsPatSoftMuon(0), muIsPatMediumMuon(0),
	  muUpsVrtxMatch(0), muL3TriggerMatch(0),

	  muMVAMuonID(0), musegmentCompatibility(0),
	  mupulldXdZ_pos_noArb(0), mupulldYdZ_pos_noArb(0),
	  mupulldXdZ_pos_ArbDef(0), mupulldYdZ_pos_ArbDef(0),
	  mupulldXdZ_pos_ArbST(0), mupulldYdZ_pos_ArbST(0),
	  mupulldXdZ_pos_noArb_any(0), mupulldYdZ_pos_noArb_any(0),
      
      // [J-U-P] Modify branches to store the original kaon tracks and the reconstructed phi.

      Jpsi_1_mu_1_Idx(0), Jpsi_1_mu_2_Idx(0),
      Jpsi_2_mu_1_Idx(0), Jpsi_2_mu_2_Idx(0),
	     Phi_pi_1_Idx(0),    Phi_pi_2_Idx(0),

      Jpsi_1_mass(0), Jpsi_1_massErr(0), Jpsi_1_massDiff(0),
      Jpsi_2_mass(0), Jpsi_2_massErr(0), Jpsi_2_massDiff(0),
         Phi_mass(0),    Phi_massErr(0),    Phi_massDiff(0),

      Jpsi_1_ctau(0), Jpsi_1_ctauErr(0), Jpsi_1_Chi2(0), Jpsi_1_ndof(0), Jpsi_1_VtxProb(0),
      Jpsi_2_ctau(0), Jpsi_2_ctauErr(0), Jpsi_2_Chi2(0), Jpsi_2_ndof(0), Jpsi_2_VtxProb(0),
                                            Phi_Chi2(0),    Phi_ndof(0),    Phi_VtxProb(0),
      
      Jpsi_1_phi(0), Jpsi_1_eta(0), Jpsi_1_pt(0),
      Jpsi_2_phi(0), Jpsi_2_eta(0), Jpsi_2_pt(0),
         Phi_phi(0),    Phi_eta(0),    Phi_pt(0),

      Jpsi_1_px(0), Jpsi_1_py(0), Jpsi_1_pz(0),
      Jpsi_2_px(0), Jpsi_2_py(0), Jpsi_2_pz(0),
         Phi_px(0),    Phi_py(0),    Phi_pz(0),
      
         Pri_mass(0),  Pri_massErr(0),
         Pri_ctau(0),  Pri_ctauErr(0), Pri_Chi2(0), Pri_ndof(0), Pri_VtxProb(0),
         Pri_px(0),    Pri_py(0),    Pri_pz(0), 
         Pri_phi(0),   Pri_eta(0),   Pri_pt(0),
      
	  // doMC
	  MC_X_py(0),
	  MC_X_pz(0),
	  MC_X_mass(0),
	  MC_Dau_Jpsipx(0),
	  MC_Dau_Jpsipy(0),
	  MC_Dau_Jpsipz(0),
	  MC_Dau_Jpsimass(0),
	  MC_Dau_psi2spx(0),
	  MC_Dau_psi2spy(0),
	  MC_Dau_psi2spz(0),
	  MC_Dau_psi2smass(0),
	  MC_Granddau_mu1px(0),
	  MC_Granddau_mu1py(0),
	  MC_Granddau_mu1pz(0),
	  MC_Granddau_mu2px(0),
	  MC_Granddau_mu2py(0),
	  MC_Granddau_mu2pz(0),
	  MC_Granddau_Jpsipx(0),
	  MC_Granddau_Jpsipy(0),
	  MC_Granddau_Jpsipz(0),
	  MC_Granddau_Jpsimass(0),
	  MC_Granddau_pi1px(0),
	  MC_Granddau_pi1py(0),
	  MC_Granddau_pi1pz(0),
	  MC_Granddau_pi2px(0),
	  MC_Granddau_pi2py(0),
	  MC_Granddau_pi2pz(0),
	  MC_Grandgranddau_mu3px(0),
	  MC_Grandgranddau_mu3py(0),
	  MC_Grandgranddau_mu3pz(0),
	  MC_Grandgranddau_mu4px(0),
	  MC_Grandgranddau_mu4py(0),
	  MC_Grandgranddau_mu4pz(0),

	  MC_X_chg(0),
	  MC_Dau_JpsipdgId(0),
	  MC_Dau_psi2spdgId(0),
	  MC_Granddau_mu1pdgId(0),
	  MC_Granddau_mu2pdgId(0),
	  MC_Granddau_JpsipdgId(0),
	  MC_Granddau_pi1pdgId(0),
	  MC_Granddau_pi2pdgId(0),
	  MC_Grandgranddau_mu3pdgId(0),
	  MC_Grandgranddau_mu4pdgId(0)
	  // mybxlumicorr(0), myrawbxlumi(0)
{
	// get token here for four-muon;
	gtRecordToken_     = consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag("gtDigis"));
	gtbeamspotToken_   = consumes<BeamSpot>(edm::InputTag("offlineBeamSpot"));
	gtprimaryVtxToken_ = consumes<VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices")); //  MINIAOD
	gtpatmuonToken_    = consumes<edm::View<pat::Muon>>(edm::InputTag("slimmedMuons"));				 //  MINIAOD
	gttriggerToken_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"));
	trackToken_        = consumes<edm::View<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates")); //  MINIAOD
	genParticlesToken_ = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
}

MultiLepPAT::~MultiLepPAT()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}
// member functions

//    ofstream myfile("comparison.txt");
// ------------ method called to for each event  ------------
void MultiLepPAT::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	

    // Load the MC results [Annotated by Eric Wang, 20240704]

	TLorentzVector MC_mu1p4, MC_mu2p4, MC_mu3p4, MC_mu4p4, MC_pi1p4, MC_pi2p4;
	if (doMC)
	{
		edm::Handle<reco::GenParticleCollection> genParticles;
		iEvent.getByToken(genParticlesToken_, genParticles);

		for (const auto &particle : *(genParticles.product()))
		{
			if (std::abs(particle.pdgId()) == 35 && particle.numberOfDaughters() >= 2)
			{
				MC_X_px->push_back(particle.px());
				MC_X_py->push_back(particle.py());
				MC_X_pz->push_back(particle.pz());
				MC_X_mass->push_back(particle.mass());
				MC_X_chg->push_back(particle.charge());
				// particle.daughter(0)->pdgId() == 443 Jpsi
				MC_Dau_JpsipdgId->push_back(particle.daughter(0)->pdgId());
				MC_Dau_Jpsipx->push_back(particle.daughter(0)->px());
				MC_Dau_Jpsipy->push_back(particle.daughter(0)->py());
				MC_Dau_Jpsipz->push_back(particle.daughter(0)->pz());
				MC_Dau_Jpsimass->push_back(particle.daughter(0)->mass());
				MC_Dau_psi2spdgId->push_back(particle.daughter(1)->pdgId());
				MC_Dau_psi2spx->push_back(particle.daughter(1)->px());
				MC_Dau_psi2spy->push_back(particle.daughter(1)->py());
				MC_Dau_psi2spz->push_back(particle.daughter(1)->pz());
				MC_Dau_psi2smass->push_back(particle.daughter(1)->mass());
				// particle.daughter(0)->daughter(1)->pdgId() == -13 mu+
				MC_Granddau_mu1pdgId->push_back(particle.daughter(0)->daughter(1)->pdgId());
				MC_Granddau_mu1px->push_back(particle.daughter(0)->daughter(1)->px());
				MC_Granddau_mu1py->push_back(particle.daughter(0)->daughter(1)->py());
				MC_Granddau_mu1pz->push_back(particle.daughter(0)->daughter(1)->pz());
				MC_Granddau_mu2pdgId->push_back(particle.daughter(0)->daughter(0)->pdgId());
				MC_Granddau_mu2px->push_back(particle.daughter(0)->daughter(0)->px());
				MC_Granddau_mu2py->push_back(particle.daughter(0)->daughter(0)->py());
				MC_Granddau_mu2pz->push_back(particle.daughter(0)->daughter(0)->pz());
				// particle.daughter(1)->daughter(0)->pdgId() == 443 Jpsi from psi2s
				MC_Granddau_JpsipdgId->push_back(particle.daughter(1)->daughter(0)->pdgId());
				MC_Granddau_Jpsipx->push_back(particle.daughter(1)->daughter(0)->px());
				MC_Granddau_Jpsipy->push_back(particle.daughter(1)->daughter(0)->py());
				MC_Granddau_Jpsipz->push_back(particle.daughter(1)->daughter(0)->pz());
				MC_Granddau_Jpsimass->push_back(particle.daughter(1)->daughter(0)->mass());
				// particle.daughter(1)->daughter(1)->pdgId() == 211 pi+ from psi2s
				MC_Granddau_pi1pdgId->push_back(particle.daughter(1)->daughter(1)->pdgId());
				MC_Granddau_pi1px->push_back(particle.daughter(1)->daughter(1)->px());
				MC_Granddau_pi1py->push_back(particle.daughter(1)->daughter(1)->py());
				MC_Granddau_pi1pz->push_back(particle.daughter(1)->daughter(1)->pz());
				MC_Granddau_pi2pdgId->push_back(particle.daughter(1)->daughter(2)->pdgId());
				MC_Granddau_pi2px->push_back(particle.daughter(1)->daughter(2)->px());
				MC_Granddau_pi2py->push_back(particle.daughter(1)->daughter(2)->py());
				MC_Granddau_pi2pz->push_back(particle.daughter(1)->daughter(2)->pz());
				// particle.daughter(1)->daughter(0)->daughter(1)->pdgId() == -13 mu+ from Jpsi from psi2s
				MC_Grandgranddau_mu3pdgId->push_back(particle.daughter(1)->daughter(0)->daughter(1)->pdgId());
				MC_Grandgranddau_mu3px->push_back(particle.daughter(1)->daughter(0)->daughter(1)->px());
				MC_Grandgranddau_mu3py->push_back(particle.daughter(1)->daughter(0)->daughter(1)->py());
				MC_Grandgranddau_mu3pz->push_back(particle.daughter(1)->daughter(0)->daughter(1)->pz());
				MC_Grandgranddau_mu4pdgId->push_back(particle.daughter(1)->daughter(0)->daughter(0)->pdgId());
				MC_Grandgranddau_mu4px->push_back(particle.daughter(1)->daughter(0)->daughter(0)->px());
				MC_Grandgranddau_mu4py->push_back(particle.daughter(1)->daughter(0)->daughter(0)->py());
				MC_Grandgranddau_mu4pz->push_back(particle.daughter(1)->daughter(0)->daughter(0)->pz());

				MC_mu1p4.SetXYZM((*MC_Granddau_mu1px)[0], (*MC_Granddau_mu1py)[0], (*MC_Granddau_mu1pz)[0], myMuMass);
				MC_mu2p4.SetXYZM((*MC_Granddau_mu2px)[0], (*MC_Granddau_mu2py)[0], (*MC_Granddau_mu2pz)[0], myMuMass);
				MC_mu3p4.SetXYZM((*MC_Grandgranddau_mu3px)[0], (*MC_Grandgranddau_mu3py)[0], (*MC_Grandgranddau_mu3pz)[0], myMuMass);
				MC_mu4p4.SetXYZM((*MC_Grandgranddau_mu4px)[0], (*MC_Grandgranddau_mu4py)[0], (*MC_Grandgranddau_mu4pz)[0], myMuMass);
				MC_pi1p4.SetXYZM((*MC_Granddau_pi1px)[0], (*MC_Granddau_pi1py)[0], (*MC_Granddau_pi1pz)[0], myPiMass);
				MC_pi2p4.SetXYZM((*MC_Granddau_pi2px)[0], (*MC_Granddau_pi2py)[0], (*MC_Granddau_pi2pz)[0], myPiMass);
			} // for (const auto& particle: *(genParticles.product()))
		}	  // if ( std::abs(particle.pdgId())  == 35 && particle.numberOfDaughters() ==2 )
	}		  // doMC

	// get event content information

	using std::vector;
	using namespace edm;
	using namespace reco;
	using namespace std;

	runNum = iEvent.id().run();
	evtNum = iEvent.id().event();
	lumiNum = iEvent.id().luminosityBlock();

	const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_);

    /**************************************************************************
     * [Section]
     *      HLT Trigger Info
     * [Implementation]
     *      - Call getByToken() to acquire HLT results
     *      - Categorize the 
     * 
    **************************************************************************/

	edm::Handle<edm::TriggerResults> hltresults;
	bool Error_t = false;
	int nJpsitrigger = TriggersForJpsi_.size();
	try
	{
		iEvent.getByToken(gttriggerToken_, hltresults);
	}
	catch (...)
	{	
		Error_t = true;
		cout << "Couldn't get handle on HLT Trigger!" << endl;
	}
	if (Error_t || !hltresults.isValid())
	{
		cout << "No Trigger Results!" << endl;
	}
	else
	{
		int ntrigs = hltresults->size();
		if (ntrigs == 0)
		{
			cout << "No trigger name given in TriggerResults of the input " << endl;
		}

		edm::TriggerNames triggerNames_;
		triggerNames_ = iEvent.triggerNames(*hltresults);   // Get trigger names [Annotated by Eric Wang, 20240704]

		for (int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++)
		{
			JpsiMatchTrig[JpsiTrig] = 0;
		} // Jpsi trigger

		for (int itrig = 0; itrig < ntrigs; itrig++) // Loop over all triggers [Annotated by Eric Wang, 20240704]
		{
			string trigName = triggerNames_.triggerName(itrig);
			int hltflag = (*hltresults)[itrig].accept();  // What is accept()? [Question from Eric Wang, 20240704]
			trigRes->push_back(hltflag);
			trigNames->push_back(trigName);

			for (int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++)
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
	} // end of HLT trigger info

	std::string vrtxFilter("hltVertexmumuFilterUpsilonMuon");
	std::string L3Filter("hltTripleMuL3PreFiltered0");

	edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(theTTBuilderToken_);

	edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
	iEvent.getByToken(gtRecordToken_, gtRecord);
	const DecisionWord dWord = gtRecord->decisionWord();

	const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
	for (unsigned int l1i = 0; l1i != ttWord.size(); ++l1i)
	{
		L1TT->push_back(ttWord.at(l1i));
	}

    // Get the primary vertex information [Annotated by Eric Wang, 20240704]
    //      Initialization
	Vertex thePrimaryV;
	Vertex theRecoVtx;
	Vertex theBeamSpotV;
	BeamSpot beamSpot;
	math::XYZPoint RefVtx;

	// get BeamSplot
	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(gtbeamspotToken_, beamSpotHandle);
	if (beamSpotHandle.isValid())
	{
		beamSpot = *beamSpotHandle;
		theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
	}
	else
	{
		cout << "No beam spot available from EventSetup" << endl;
	}

	Handle<VertexCollection> recVtxs;
	iEvent.getByToken(gtprimaryVtxToken_, recVtxs);
	unsigned int nVtxTrks = 0;

	///////////////////////////////////////////////////////////////////////
	////////////////Check Lines below for Primary Vertex///////////////////
	///////////////////////////////////////////////////////////////////////

    // Determine the good vertices [Annotated by Eric Wang, 20240704]

	int myNGoodPrimVtx = 0;
	for (unsigned myi = 0; myi < recVtxs->size(); myi++)
	{
        //Selection criteria:  ndof    -> "enough particle info"(?)     [Annotated by Eric Wang, 20240704]
        //                      z, rhp  -> "not too far off"            [Annotated by Eric Wang, 20240704]
		if ((*recVtxs)[myi].ndof() >= 5 
                && ( fabs((*recVtxs)[myi].z()) <= 24 && fabs((*recVtxs)[myi].position().rho())) <= 2.0)
		{
			myNGoodPrimVtx++;
		}
	}
	nGoodPrimVtx = myNGoodPrimVtx;

	if (recVtxs->begin() != recVtxs->end())  // At least one vertex [Annotated by Eric Wang, 20240704]
	{
		if (addXlessPrimaryVertex_ || resolveAmbiguity_)  // What criterion? [Question from Eric Wang, 20240704]
		{
			thePrimaryV = Vertex(*(recVtxs->begin()));
		}
		else
		{
            // TO_IMPR_CPP11 (for(auto ...)) [Tagged by Eric Wang, 20240704]
			for (reco::VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx)
			{
                // Choose the vertex with the most tracks [Annotated by Eric Wang, 20240704]
				if (nVtxTrks < vtx->tracksSize())   
				{
					nVtxTrks = vtx->tracksSize();
					thePrimaryV = Vertex(*vtx);
				}
			}
		}
	}
	else
	{
        // If no vertex is found, use the beam spot as the primary vertex [Annotated by Eric Wang, 20240704]
		thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
	}

    // Vertex reconstruction results [Annotated by Eric Wang, 20240704]
	RefVtx = thePrimaryV.position();
	priVtxX = (thePrimaryV.position().x());
	priVtxY = (thePrimaryV.position().y());
	priVtxZ = (thePrimaryV.position().z());
	priVtxXE = (thePrimaryV.xError());
	priVtxYE = (thePrimaryV.yError());
	priVtxZE = (thePrimaryV.zError());
	priVtxChiNorm = (thePrimaryV.normalizedChi2());
	priVtxChi = thePrimaryV.chi2();
	priVtxCL = ChiSquaredProbability((double)(thePrimaryV.chi2()), (double)(thePrimaryV.ndof()));
	// 这里去掉了XYZ两两的协方差

	edm::Handle<edm::View<pat::Muon> > thePATMuonHandle; //  MINIAOD
	iEvent.getByToken(gtpatmuonToken_, thePATMuonHandle);
	edm::Handle<edm::View<pat::PackedCandidate> > theTrackHandle;   //  MINIAOD
	iEvent.getByToken(trackToken_,                theTrackHandle);  //  MINIAOD
	std::vector<edm::View<pat::PackedCandidate>::const_iterator> nonMuonPionTrack;

	// Copy tracks iterators
 	for (edm::View<pat::PackedCandidate>::const_iterator iTrackc = theTrackHandle->begin(); // MINIAOD 
	            iTrackc != theTrackHandle->end(); ++iTrackc                                 )
	{
		nonMuonPionTrack.push_back(iTrackc);
	}

    // Initialize the muon track block [Annotated by Eric Wang, 20240704]

    // Visualize number of muons
    std::cout << "Number of muons: " << thePATMuonHandle->size();
	// 以上非必要

	if (thePATMuonHandle->size() >= 4) // Require at least 6 muons present [Annotated by Eric Wang, 20240704]
	// JJP->4mu+2K, 改回4
	{
		vector<std::string> theInputVariables;
		theInputVariables.push_back("validFrac");
		theInputVariables.push_back("globalChi2");
		theInputVariables.push_back("pt");
		theInputVariables.push_back("eta");
		theInputVariables.push_back("segComp");
		theInputVariables.push_back("chi2LocMom");
		theInputVariables.push_back("chi2LocPos");
		theInputVariables.push_back("glbTrackProb");
		theInputVariables.push_back("NTrkVHits");
		theInputVariables.push_back("NTrkEHitsOut");
		ReadBDT muonID(theInputVariables);
		vector<double> inputValues;
		inputValues.resize(10, 0.);

		// fill muon track block
        // TO_IMPR_CPP11 (for(auto ...)) [Tagged by Eric Wang, 20240704]
		for (edm::View<pat::Muon>::const_iterator iMuonP =  thePATMuonHandle->begin(); //  MINIAOD
			                                      iMuonP != thePATMuonHandle->end(); ++iMuonP)
		{
			// push back all muon information
			++nMu;
			muIsPatLooseMuon->push_back(iMuonP->isLooseMuon());
			muIsPatTightMuon->push_back(iMuonP->isTightMuon(thePrimaryV));
			muIsPatSoftMuon->push_back(iMuonP->isSoftMuon(thePrimaryV));
			muIsPatMediumMuon->push_back(iMuonP->isMediumMuon());

			muPx->push_back(iMuonP->px());
			muPy->push_back(iMuonP->py());
			muPz->push_back(iMuonP->pz());
			muCharge->push_back(iMuonP->charge());
			// 没有引入d0和TrackIso等信息

			int goodSoftMuonNewIlseMod = 0;
			int goodSoftMuonNewIlse = 0;
			int goodLooseMuonNew = 0;
			int goodLooseMuon = 0;
			int goodTightMuon = 0;
			
			
			// Find and delete muon Tracks in PionTracks
			for (std::vector<edm::View<pat::PackedCandidate>::const_iterator>::const_iterator iTrackfID  = nonMuonPionTrack.begin(); // MINIAOD
			                                                                                  iTrackfID != nonMuonPionTrack.end(); 
                                                                                            ++iTrackfID                             )
			{
				if(iMuonP->track().isNull())
				{
					continue;
				}
				edm::View<pat::PackedCandidate>::const_iterator iTrackf = *(iTrackfID);		

                // Why call the function outside? [Question from Eric Wang, 20240704]
				iMuonP->track()->px();

                // Match using the momentum. [Annotated by Eric Wang, 20240704]                  
				if (   iTrackf->px() == iMuonP->track()->px() 
                    && iTrackf->py() == iMuonP->track()->py() 
                    && iTrackf->pz() == iMuonP->track()->pz())
				{
					nonMuonPionTrack.erase(iTrackfID);
					iTrackfID = iTrackfID - 1;
				}
			}
			// float mymuMVABs = -1;

            // Check if match any HLT for Jpsi and Upsilon [Annotated by Eric Wang, 20240704]
            // TO_ENC [Tagged by Eric Wang, 20240704]

			bool isJpsiTrigMatch = false;
			// bool isJpsiFilterMatch = false;

			for (int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++)
			{
				if(JpsiMatchTrig[JpsiTrig] != 0)
				{
					isJpsiTrigMatch = true;
					break;
				}
			}
			// for(unsigned int JpsiTrig = 0; JpsiTrig < FiltersForJpsi_.size();JpsiTrig++)
			// {
			// 	if (isJpsiTrigMatch && hltresults.isValid())
			// 	{
			// 		pat::TriggerObjectStandAlone *tempTriggerObject = nullptr; // 指针初始化
			// 		for (auto i = iMuonP->triggerObjectMatches().begin(); i != iMuonP->triggerObjectMatches().end(); ++i)
			// 		{
			// 			tempTriggerObject = new pat::TriggerObjectStandAlone(*i);
			// 			tempTriggerObject->unpackFilterLabels(iEvent, *hltresults);
			// 			if(tempTriggerObject->hasFilterLabel(FiltersForJpsi_[JpsiTrig]))
			// 			{
			// 				isJpsiFilterMatch = true;
			// 			}
			// 			delete tempTriggerObject;
			// 		}
			// 	}
			// }

			muJpsiFilterRes->push_back(isJpsiTrigMatch); // 这一步内存溢出了？？？
			std::cout << "0" << std::endl;
		}
	} // if two muons

	std::cout << "Finish the part of all muons' trigger match."  << endl;

	if (doMC)
	{
		// pion loop

        // Compare direct MC results and RECO results [Annotated by Eric Wang, 20240704]
        //      For those that matches within the precision limit, store the momentum 

        // TO_IMPR_CPP11 (for(auto ...)) [Tagged by Eric Wang, 20240704]
		for (edm::View<pat::PackedCandidate>::const_iterator iTrack = theTrackHandle->begin(); // MINIAOD
			                                                iTrack != theTrackHandle->end(); ++iTrack)
		{
			TLorentzVector RECO_pip4;
			RECO_pip4.SetXYZM(iTrack->px(), iTrack->py(), iTrack->pz(), myPiMass);
			if (fabs(MC_pi1p4.Pt() - RECO_pip4.Pt()) < 0.08 * MC_pi1p4.Pt() && MC_pi1p4.DeltaR(RECO_pip4) < 0.1)
			{
				Match_pi1px->push_back(RECO_pip4.Px());
				Match_pi1py->push_back(RECO_pip4.Py());
				Match_pi1pz->push_back(RECO_pip4.Pz());
			}
			if ((fabs(MC_pi2p4.Pt() - RECO_pip4.Pt()) < 0.08 * MC_pi2p4.Pt() && MC_pi2p4.DeltaR(RECO_pip4) < 0.1))
			{
				Match_pi2px->push_back(RECO_pip4.Px());
				Match_pi2py->push_back(RECO_pip4.Py());
				Match_pi2pz->push_back(RECO_pip4.Pz());
			}
		}
	} // if(doMC)

    // Set the mass constraints for further reconstruction. [Annotated by Eric Wang, 20240704]

    // [J-U-P] Switch from pions to kaons.

	// It takes a lot less memory out of the loop
	KinematicConstraint *Jpsi_cs = new MassKinematicConstraint(myJpsiMass, myJpsiMassErr);
	KinematicConstraint *Jpsi_cs34 = new MassKinematicConstraint(myJpsiMass, myJpsiMassErr);
	// 要做mass constrain的话这里要加phi的相关变量
	double pionPTcut = 0.25;
	double pionDRcut = 0.7;
	double vtxProbPreCut = 1.0e-7;

    // booleans marking whether the muon pair satisfies the mass constraint
    bool isJpsiMuPair = false;
    bool isGoodVtxFit = false;

    // Muon factory
    KinematicParticleFactoryFromTransientTrack muPairFactory;

    // [J-U-P] Will only need 4 muons.

    // Will be working reco from 3 pairs of muons. 
	if (thePATMuonHandle->size() < 4)
	{
		return;
	}

    // Temporary storage for the muon pair [Annotated by Eric Wang, 20240704]
    std::vector<RefCountedKinematicParticle> transMuonPair;
    std::vector<uint>                        transMuPairId;
    ParticleMass muMass = myMuMass;
    float muMassSigma   = myMuMassErr;
    float chi2 = 0.;
	float ndof = 0.;

    /**************************************************************************
     * [Section]
     *      Muon pairing and preselection with vertex fitting
     * [Implementation]
     *      - Loop over all existing muon pairs.
     *      - Apply kinematic and mass window constraints.
     *      - Fit the vertex to judge.
     *      - Classify muons pairs as Jpsi or Upsilon candidates with mass.
     * [Note]
     *      The intermidiate storage for muon pair stores the muons as 
     *      RefCountedKinematicParticle. This saves repeated reco.
    **************************************************************************/

    // Candidates of muon pairs from Jpsi or Upsilon
    using muon_t   = RefCountedKinematicParticle;
    using muList_t = std::pair< vector<muon_t>, vector<uint> >;
    std::vector< muList_t > muPairCand_Jpsi;

    // Selection for the muon candidates
    for(auto iMuon1 =  thePATMuonHandle->begin(); 
             iMuon1 != thePATMuonHandle->end(); ++iMuon1){
        TrackRef muTrack1 = iMuon1->track();
        if (muTrack1.isNull()){
            continue;
        }
        // Build transient track and store.
        TransientTrack transTrk1(muTrack1, &(bFieldHandle));
        transMuonPair.push_back(muPairFactory.particle(transTrk1, muMass, chi2, ndof, muMassSigma));
        transMuPairId.push_back(iMuon1 - thePATMuonHandle->begin());

        // Next muon candidate.
        for(auto iMuon2  = iMuon1 + 1; 
                 iMuon2 != thePATMuonHandle->end(); ++iMuon2){
            // DEBUG: display current muon pair.
            // Build transient track and store.
            TrackRef muTrack2 = iMuon2->track();
            if (muTrack2.isNull()){
                continue;
            }
            TransientTrack transTrk2(muTrack2, &(bFieldHandle));
            // Charge requirement.
            if ((iMuon1->charge() + iMuon2->charge()) != 0){
				continue;
			}
            // Dynamics selection. A very crude selection.
            // Involves more calculation and is therefore done after kinematics.
            isJpsiMuPair = (1 < (iMuon1->p4() + iMuon2->p4()).mass()
                              && (iMuon1->p4() + iMuon2->p4()).mass() < 4); // 这里mass cut是2~6，与秦俊凯学长定的1~4不一样，Eric说是往宽了放，我先调成1~4吧

            // isJpsiMuPair = true;
            // isUpsMuPair  = true;
            if((!isJpsiMuPair)){
                continue;
            }
            transMuonPair.push_back(muPairFactory.particle(transTrk2,  muMass, 
                                                           chi2, ndof, muMassSigma) );
            transMuPairId.push_back(iMuon2 - thePATMuonHandle->begin());
            // Judging with vertex fitting.
            if(!particlesToVtx(transMuonPair)){
                continue;
            }
            // Passing all the checks, store the muon pair as pairs of RefCountedKinematicParticle.
            if(isJpsiMuPair){
                muPairCand_Jpsi.push_back(
                    std::make_pair(transMuonPair, transMuPairId) );
            }
            // Clear the transient muon pair for the next pair.
            transMuonPair.pop_back();
            transMuPairId.pop_back();
			// pop_back弹掉最后一个元素，不用重复加入iMuon1了
        }
        transMuonPair.pop_back();
        transMuPairId.pop_back();
    }

	std::cout << "Finish the part of muon pair."  << endl;

	// 先判定jpsi是否共顶点
	// Classes for the fitting process.
    RefCountedKinematicTree vtxFitTree_Jpsi_1;
    RefCountedKinematicTree vtxFitTree_Jpsi_2;
    RefCountedKinematicTree vtxFitTree_Phi;
    RefCountedKinematicTree vtxFitTree_Pri;

    // Classes for secondary particles (Jpsi and Upsilon)
    RefCountedKinematicParticle Jpsi_1_Fit_noMC, Jpsi_2_Fit_noMC, Phi_Fit_noMC, Pri_Fit_noMC;
    RefCountedKinematicVertex   Jpsi_1_Vtx_noMC, Jpsi_2_Vtx_noMC, Phi_Vtx_noMC, Pri_Vtx_noMC;
    KinematicParameters         Jpsi_1_Para,     Jpsi_2_Para,     Phi_Para,     Pri_Para;
    std::vector< RefCountedKinematicParticle >  interOnia;

    // Markers for fitting. Only marks if a result is constructed
    bool isValidJpsi_1, isValidJpsi_2, isValidPhi, isValidPri;
    // Fitted mass error is also stricter marker for fitting.
    double tmp_Jpsi_1_massErr, tmp_Jpsi_2_massErr, tmp_Phi_massErr, tmp_Pri_massErr;
    // Temporary storage for particle dynamics.
    double tmp_pt, tmp_eta, tmp_phi;
	bool isValidJpsi_both = false;

	for(auto muPair_Jpsi_1  = muPairCand_Jpsi.begin(); 
             muPair_Jpsi_1 != muPairCand_Jpsi.end();  muPair_Jpsi_1++){
        for(auto muPair_Jpsi_2  = muPair_Jpsi_1 + 1; 
                 muPair_Jpsi_2 != muPairCand_Jpsi.end(); muPair_Jpsi_2++){
            // Check if the muon pairs overlap.
            if(isOverlapPair(*muPair_Jpsi_1, *muPair_Jpsi_2)){
                 continue;
            }
			

			// Initialize the marker for primary vertex
			isValidPri = false;
			// Start constructing the fit tree.
			// Use particlesToVtx() to fit the quarkonia once more.
			// [J-U-P] Add appropriate code for phi reconstruction from kaons.
			isValidJpsi_1 = particlesToVtx(vtxFitTree_Jpsi_1, muPair_Jpsi_1->first, "final Jpsi_1");
			isValidJpsi_2 = particlesToVtx(vtxFitTree_Jpsi_2, muPair_Jpsi_2->first, "final Jpsi_2");

			if(isValidJpsi_1 && isValidJpsi_2){
				// Extract the vertex and the particle parameters from valid results.
				// Here, when an invalid fit is detected, the massErr is set to -9.
				extractFitRes(vtxFitTree_Jpsi_1, Jpsi_1_Fit_noMC, Jpsi_1_Vtx_noMC, tmp_Jpsi_1_massErr);
				extractFitRes(vtxFitTree_Jpsi_2, Jpsi_2_Fit_noMC, Jpsi_2_Vtx_noMC, tmp_Jpsi_2_massErr);
				// Look for "Good Fit". Judge by the massErr.

				if(tmp_Jpsi_1_massErr >= 0.0 && tmp_Jpsi_2_massErr >= 0.0){
					// Initialize the final fitting marker and the secondary particles.
					interOnia.push_back(Jpsi_1_Fit_noMC);
					interOnia.push_back(Jpsi_2_Fit_noMC);
					// Fit the quarkonia to the same vertex
					isValidPri = particlesToVtx(vtxFitTree_Pri, interOnia, "primary vertex");
					interOnia.clear();
				}
			}
			// Work with all fit results above. (Jpsi_1, Jpsi_2, Ups, Pri)
			// Primary vertex fitting comes first.


			if(!isValidPri){
				continue;
			}
			else{
				isValidJpsi_both = true;
			}
        }
    }
	if(!isValidJpsi_both){
		return;
	}

	
    // 看上去需要仿照wc的方式写phi的程序了
	bool isPhiTrackPair = false;
    bool isGoodVtxFit_Phi = false;

    // Muon factory
    KinematicParticleFactoryFromTransientTrack PhiFactory;

    // Temporary storage for the muon pair [Annotated by Eric Wang, 20240704]
    std::vector<RefCountedKinematicParticle> transTrackPair;
    std::vector<uint>                        transTrackPairId;
    ParticleMass piMass = myPiMass;
    float piMassSigma   = myPiMassErr; //见.h文件
    chi2 = 0.;
	ndof = 0.;

    // Candidates of track pairs from Phi
    using pion_t   = RefCountedKinematicParticle;
    using piList_t = std::pair< vector<pion_t>, vector<uint> >;
    std::vector< piList_t > piPairCand_Phi;

	std::cout << "Start the part of track pair."  << endl;
	std::cout << "the number of track" << nonMuonPionTrack.size() << endl;

    // Selection for the muon candidates
    for(std::vector<edm::View<pat::PackedCandidate>::const_iterator>::const_iterator iTrack1ID = nonMuonPionTrack.begin(); // MINIAOD
						 iTrack1ID != nonMuonPionTrack.end(); ++iTrack1ID){
		edm::View<pat::PackedCandidate>::const_iterator iTrack1 = *(iTrack1ID);
		if (iTrack1->pt() < pionPTcut)
		{
			continue;
		}
		if (!iTrack1->hasTrackDetails() || iTrack1->charge() == 0)
		{
			continue;
		}
		if (iTrack1->fromPV() < 1 || abs(iTrack1->eta()) > 2.5)
        {
			continue;
		}
		if (iTrack1->numberOfHits() < 5 || iTrack1->bestTrack()->normalizedChi2() > 8 || !(iTrack1->bestTrack()->quality(reco::Track::highPurity)))
		{
			continue;
		}
        // Build transient track and store.
        TransientTrack trackTT1(*(iTrack1->bestTrack()), &(bFieldHandle));
        transTrackPair.push_back(PhiFactory.particle(trackTT1, piMass, chi2, ndof, piMassSigma));
        transTrackPairId.push_back(iTrack1ID - nonMuonPionTrack.begin());

        // Next muon candidate.
        for(std::vector<edm::View<pat::PackedCandidate>::const_iterator>::const_iterator iTrack2ID = iTrack1ID + 1; // MINIAOD
							 iTrack2ID != nonMuonPionTrack.end(); ++iTrack2ID){
            // DEBUG: display current muon pair.
            // Build transient track and store.
            edm::View<pat::PackedCandidate>::const_iterator iTrack2 = *(iTrack2ID);
			if (iTrack2->pt() < pionPTcut)
			{
				continue;
			}
			if (!iTrack2->hasTrackDetails() || iTrack2->charge() == 0)
			{
				continue;
			}
			if (iTrack2->fromPV() < 1 || abs(iTrack2->eta()) > 2.5)
			{
				continue;
			}
			if (iTrack2->numberOfHits() < 5 || iTrack2->bestTrack()->normalizedChi2() > 8 || !(iTrack2->bestTrack()->quality(reco::Track::highPurity)))
			{
				continue;
			}
            TransientTrack trackTT2(*(iTrack2->bestTrack()), &(bFieldHandle)); 

            // Charge requirement.
            if ((iTrack1->charge() + iTrack2->charge()) != 0)
				continue;

			TLorentzVector P4_Track1, P4_Track2, P4_Phi;
			P4_Track1.SetPtEtaPhiM(iTrack1->pt(), iTrack1->eta(), iTrack1->phi(), myPiMass);
			P4_Track2.SetPtEtaPhiM(iTrack2->pt(), iTrack2->eta(), iTrack2->phi(), myPiMass);
			P4_Phi = P4_Track1 + P4_Track2;

			if (P4_Track1.DeltaR(P4_Phi) > pionDRcut)
			{
				continue;
			}
			if (P4_Track2.DeltaR(P4_Phi) > pionDRcut)
			{
				continue;
			}

            // Dynamics selection. A very crude selection.
            // Involves more calculation and is therefore done after kinematics.
            isPhiTrackPair = (0.8 < (iTrack1->p4() + iTrack2->p4()).mass()
                              && (iTrack1->p4() + iTrack2->p4()).mass() < 1.2); // 要不要加，加多少？

            // isJpsiMuPair = true;
            // isUpsMuPair  = true;
            if((!isPhiTrackPair)){
                continue;
            }
            transTrackPair.push_back(PhiFactory.particle(trackTT2,  piMass, 
                                                           chi2, ndof, piMassSigma) ); //此处之前错为muMass
            transTrackPairId.push_back(iTrack2ID - nonMuonPionTrack.begin());
            // Judging with vertex fitting.
            if(!particlesToVtx(transTrackPair)){
                continue;
            }
            // Passing all the checks, store the muon pair as pairs of RefCountedKinematicParticle.
            if(isPhiTrackPair){
                piPairCand_Phi.push_back(
                    std::make_pair(transTrackPair, transTrackPairId) );
            }
            // Clear the transient muon pair for the next pair.
            transTrackPair.pop_back();
            transTrackPairId.pop_back();
			// pop_back弹掉最后一个元素，不用重复加入iMuon1了
        }
        transTrackPair.pop_back();
		transTrackPairId.pop_back();
    }

	std::cout << "Finish the part of track pair."  << endl;

    /**************************************************************************
     * [Section]
     *      Jpsi and Upsilon reconstruction and fitting.
     * [Implementation]
     *      - Loop over all existing muon pairs.
     *      - Fit the primary vertex for Jpsi and Upsilon candidates.
     *      - Store the fitting results into temporary vectors.
     * [Note]
     *      Mass constriant is not applied to any quarkonia candidates.
     *      TO IMPR: Multiple candidate issue.
     *      - Difficulty: Identify the "multi-candidate" case. (Hashing?)
     *      - Distinction between Jpsi and Upsilon is important!
     *      - Possible selection: massErr ratio; total pT^2;
     *      - Add event number 
    **************************************************************************/


    // // Classes for the fitting process.
    // RefCountedKinematicTree vtxFitTree_Jpsi_1;
    // RefCountedKinematicTree vtxFitTree_Jpsi_2;
    // RefCountedKinematicTree vtxFitTree_Phi;
    // RefCountedKinematicTree vtxFitTree_Pri;

    // // Classes for secondary particles (Jpsi and Upsilon)
    // RefCountedKinematicParticle Jpsi_1_Fit_noMC, Jpsi_2_Fit_noMC, Phi_Fit_noMC, Pri_Fit_noMC;
    // RefCountedKinematicVertex   Jpsi_1_Vtx_noMC, Jpsi_2_Vtx_noMC, Phi_Vtx_noMC, Pri_Vtx_noMC;
    // KinematicParameters         Jpsi_1_Para,     Jpsi_2_Para,     Phi_Para,     Pri_Para;
    // std::vector< RefCountedKinematicParticle >  interOnia;

    // // Markers for fitting. Only marks if a result is constructed
    // bool isValidJpsi_1, isValidJpsi_2, isValidPhi, isValidPri;
    // // Fitted mass error is also stricter marker for fitting.
    // double tmp_Jpsi_1_massErr, tmp_Jpsi_2_massErr, tmp_Phi_massErr, tmp_Pri_massErr;
    // // Temporary storage for particle dynamics.
    // double tmp_pt, tmp_eta, tmp_phi;


	for(auto piPair_Phi  = piPairCand_Phi.begin(); 
                     piPair_Phi != piPairCand_Phi.end(); piPair_Phi++){
		std::cout << "First Kaon Mass: " << piPair_Phi->first[0]->currentState().mass() << std::endl;
		std::cout << "Second Kaon Mass: " << piPair_Phi->first[1]->currentState().mass() << std::endl;
		isValidPhi    = particlesToVtx(vtxFitTree_Phi,    piPair_Phi->first,    "final Phi");
		Phi_pi_1_Idx->push_back(piPair_Phi->second[0]);
        Phi_pi_2_Idx->push_back(piPair_Phi->second[1]);
		if (isValidPhi){
			extractFitRes(vtxFitTree_Phi, Phi_Fit_noMC, Phi_Vtx_noMC, tmp_Phi_massErr);
			std::cout << "Phi mass: " << Phi_Fit_noMC->currentState().mass() << std::endl;
			std::cout << "Phi mass error: " << tmp_Phi_massErr << std::endl;
		}
		else{
			std::cout << "Phi fitting failed." << std::endl;
			continue;
		}
		for(auto muPair_Jpsi_1  = muPairCand_Jpsi.begin(); 
				muPair_Jpsi_1 != muPairCand_Jpsi.end();  muPair_Jpsi_1++){
			for(auto muPair_Jpsi_2  = muPair_Jpsi_1 + 1; 
					muPair_Jpsi_2 != muPairCand_Jpsi.end(); muPair_Jpsi_2++){
				// Check if the muon pairs overlap.
				if(isOverlapPair(*muPair_Jpsi_1, *muPair_Jpsi_2)){
					continue;
				}
                // Initialize the marker for primary vertex
                isValidPri = false;
                // Start constructing the fit tree.
                // Use particlesToVtx() to fit the quarkonia once more.
                // [J-U-P] Add appropriate code for phi reconstruction from kaons.
                isValidJpsi_1 = particlesToVtx(vtxFitTree_Jpsi_1, muPair_Jpsi_1->first, "final Jpsi_1");
                isValidJpsi_2 = particlesToVtx(vtxFitTree_Jpsi_2, muPair_Jpsi_2->first, "final Jpsi_2");
                // isValidPhi    = particlesToVtx(vtxFitTree_Phi,    piPair_Phi->first,    "final Phi");
                // Store the index of the muons.
                Jpsi_1_mu_1_Idx->push_back(muPair_Jpsi_1->second[0]);
                Jpsi_1_mu_2_Idx->push_back(muPair_Jpsi_1->second[1]);
                Jpsi_2_mu_1_Idx->push_back(muPair_Jpsi_2->second[0]);
                Jpsi_2_mu_2_Idx->push_back(muPair_Jpsi_2->second[1]);
                // Phi_pi_1_Idx->push_back(piPair_Phi->second[0]);
                // Phi_pi_2_Idx->push_back(piPair_Phi->second[1]);
                // Check if all fit trees give non-null results.

				if(isValidJpsi_1 && isValidJpsi_2 ){
					extractFitRes(vtxFitTree_Jpsi_1, Jpsi_1_Fit_noMC, Jpsi_1_Vtx_noMC, tmp_Jpsi_1_massErr);
                    extractFitRes(vtxFitTree_Jpsi_2, Jpsi_2_Fit_noMC, Jpsi_2_Vtx_noMC, tmp_Jpsi_2_massErr);
				
                // if(isValidJpsi_1 && isValidJpsi_2 && isValidPhi){
                //     // Extract the vertex and the particle parameters from valid results.
                //     // Here, when an invalid fit is detected, the massErr is set to -9.
                //     extractFitRes(vtxFitTree_Jpsi_1, Jpsi_1_Fit_noMC, Jpsi_1_Vtx_noMC, tmp_Jpsi_1_massErr);
                //     extractFitRes(vtxFitTree_Jpsi_2, Jpsi_2_Fit_noMC, Jpsi_2_Vtx_noMC, tmp_Jpsi_2_massErr);
                //     extractFitRes(vtxFitTree_Phi,       Phi_Fit_noMC,    Phi_Vtx_noMC,    tmp_Phi_massErr);
                //     // Look for "Good Fit". Judge by the massErr.

					if(tmp_Jpsi_1_massErr >= 0.0 && tmp_Jpsi_2_massErr >= 0.0 && tmp_Phi_massErr >= 0.0){
						// Initialize the final fitting marker and the secondary particles.
						interOnia.push_back(Jpsi_1_Fit_noMC);
						interOnia.push_back(Jpsi_2_Fit_noMC);
						interOnia.push_back(Phi_Fit_noMC);
						// Fit the quarkonia to the same vertex
						isValidPri = particlesToVtx(vtxFitTree_Pri, interOnia, "primary vertex");
						interOnia.clear();
					}
				}
                // Work with all fit results above. (Jpsi_1, Jpsi_2, Ups, Pri)
                // Primary vertex fitting comes first.
  

                if(isValidPri){
                    // Extract the vertex and the particle parameters from valid results.
                    extractFitRes(vtxFitTree_Pri, Pri_Fit_noMC, Pri_Vtx_noMC, tmp_Pri_massErr);
                    getDynamics(Pri_Fit_noMC, tmp_pt, tmp_eta, tmp_phi);
                    // Store the fitting results into temporary vectors.
                    Pri_mass->push_back(Pri_Fit_noMC->currentState().mass());
                    Pri_massErr->push_back(tmp_Pri_massErr);
                    Pri_ctau->push_back(   GetcTau(   Pri_Vtx_noMC, Pri_Fit_noMC, theBeamSpotV));
                    Pri_ctauErr->push_back(GetcTauErr(Pri_Vtx_noMC, Pri_Fit_noMC, theBeamSpotV));
                    Pri_VtxProb->push_back(ChiSquaredProbability((double)(Pri_Vtx_noMC->chiSquared()), 
                                                                 (double)(Pri_Vtx_noMC->degreesOfFreedom())));
                    Pri_Chi2->push_back(Pri_Vtx_noMC->chiSquared());
                    Pri_ndof->push_back(Pri_Vtx_noMC->degreesOfFreedom());
                    Pri_px->push_back(Pri_Fit_noMC->currentState().kinematicParameters().momentum().x());
                    Pri_py->push_back(Pri_Fit_noMC->currentState().kinematicParameters().momentum().y()); 
                    Pri_pz->push_back(Pri_Fit_noMC->currentState().kinematicParameters().momentum().z());
                    Pri_phi->push_back(tmp_phi);
                    Pri_eta->push_back(tmp_eta);
                    Pri_pt->push_back(tmp_pt);
					std::cout << "Found candidate" << std::endl;
					std::cout << "Jpsi_1: " << Jpsi_1_Fit_noMC->currentState().mass() << std::endl;
					std::cout << "Jpsi_2: " << Jpsi_2_Fit_noMC->currentState().mass() << std::endl;
					std::cout << "Phi: " << Phi_Fit_noMC->currentState().mass() << std::endl;
					vtxFitTree_Phi->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitKaon1 = vtxFitTree_Phi->currentParticle();
					double Kaon1M_fit_mix = fitKaon1->currentState().mass();
					double Kaon1Px_fit_mix = fitKaon1->currentState().kinematicParameters().momentum().x();
					double Kaon1Py_fit_mix = fitKaon1->currentState().kinematicParameters().momentum().y();
					double Kaon1Pz_fit_mix = fitKaon1->currentState().kinematicParameters().momentum().z();
					std::cout << "Kaon1 mass: " << Kaon1M_fit_mix << std::endl;
					std::cout << "Kaon1 px: " << Kaon1Px_fit_mix << std::endl;
					std::cout << "Kaon1 py: " << Kaon1Py_fit_mix << std::endl;
					std::cout << "Kaon1 pz: " << Kaon1Pz_fit_mix << std::endl;
					vtxFitTree_Phi->movePointerToTheNextChild();
					RefCountedKinematicParticle fitKaon2 = vtxFitTree_Phi->currentParticle();
					double Kaon2M_fit_mix = fitKaon2->currentState().mass();
					double Kaon2Px_fit_mix = fitKaon2->currentState().kinematicParameters().momentum().x();
					double Kaon2Py_fit_mix = fitKaon2->currentState().kinematicParameters().momentum().y();
					double Kaon2Pz_fit_mix = fitKaon2->currentState().kinematicParameters().momentum().z();
					std::cout << "Kaon2 mass: " << Kaon2M_fit_mix << std::endl;
					std::cout << "Kaon2 px: " << Kaon2Px_fit_mix << std::endl;
					std::cout << "Kaon2 py: " << Kaon2Py_fit_mix << std::endl;
					std::cout << "Kaon2 pz: " << Kaon2Pz_fit_mix << std::endl;
                }
                else{
                    // Store "error code" -999 for the primary vertex fitting.
                    Pri_mass->push_back(-999);
                    Pri_massErr->push_back(-999);
                    Pri_ctau->push_back(-999);
                    Pri_ctauErr->push_back(-999);
                    Pri_VtxProb->push_back(-999);
                    Pri_Chi2->push_back(-999);
                    Pri_ndof->push_back(-999);
                    Pri_px->push_back(-999);
                    Pri_py->push_back(-999);
                    Pri_pz->push_back(-999);
                    Pri_phi->push_back(-999);
                    Pri_eta->push_back(-999);
                    Pri_pt->push_back(-999);
                }
                // Then comes the secondary particles (quarkonia).
                if(isValidJpsi_1){
                    getDynamics(Jpsi_1_Fit_noMC, tmp_pt, tmp_eta, tmp_phi);
                    Jpsi_1_mass->push_back(    Jpsi_1_Fit_noMC->currentState().mass());
                    Jpsi_1_massDiff->push_back(Jpsi_1_Fit_noMC->currentState().mass() - myJpsiMass);
                    Jpsi_1_massErr->push_back( tmp_Jpsi_1_massErr);
                    Jpsi_1_ctau->push_back(   GetcTau(   Jpsi_1_Vtx_noMC, Jpsi_1_Fit_noMC, theBeamSpotV));
                    Jpsi_1_ctauErr->push_back(GetcTauErr(Jpsi_1_Vtx_noMC, Jpsi_1_Fit_noMC, theBeamSpotV));
                    Jpsi_1_Chi2->push_back(double(Jpsi_1_Vtx_noMC->chiSquared()));
                    Jpsi_1_ndof->push_back(double(Jpsi_1_Vtx_noMC->degreesOfFreedom()));
                    Jpsi_1_VtxProb->push_back(ChiSquaredProbability((double)(Jpsi_1_Vtx_noMC->chiSquared()), 
                                                                    (double)(Jpsi_1_Vtx_noMC->degreesOfFreedom())));
                    Jpsi_1_px->push_back(Jpsi_1_Fit_noMC->currentState().kinematicParameters().momentum().x());
                    Jpsi_1_py->push_back(Jpsi_1_Fit_noMC->currentState().kinematicParameters().momentum().y());
                    Jpsi_1_pz->push_back(Jpsi_1_Fit_noMC->currentState().kinematicParameters().momentum().z());
                    Jpsi_1_phi->push_back(tmp_pt);
                    Jpsi_1_eta->push_back(tmp_eta);
                    Jpsi_1_pt->push_back(tmp_pt);
                }
                else{
                    
		// Store "error code" -9 for the secondary particles (quarkonia).

                    Jpsi_1_mass->push_back(-9);
                    Jpsi_1_massErr->push_back(-9);
                    Jpsi_1_massDiff->push_back(-9);
                    Jpsi_1_ctau->push_back(-9);
                    Jpsi_1_ctauErr->push_back(-9);
                    Jpsi_1_Chi2->push_back(-9);
                    Jpsi_1_ndof->push_back(-9);
                    Jpsi_1_VtxProb->push_back(-9);
                    Jpsi_1_px->push_back(-9);
                    Jpsi_1_py->push_back(-9);
                    Jpsi_1_pz->push_back(-9);
                    Jpsi_1_phi->push_back(-9);
                    Jpsi_1_eta->push_back(-9);
                    Jpsi_1_pt->push_back(-9);
                }
                if(isValidJpsi_2){
                    getDynamics(Jpsi_2_Fit_noMC, tmp_pt, tmp_eta, tmp_phi);
                    Jpsi_2_mass->push_back(    Jpsi_2_Fit_noMC->currentState().mass());
                    Jpsi_2_massDiff->push_back(Jpsi_2_Fit_noMC->currentState().mass() - myJpsiMass);
                    Jpsi_2_massErr->push_back( tmp_Jpsi_2_massErr);
                    Jpsi_2_ctau->push_back(   GetcTau(   Jpsi_2_Vtx_noMC, Jpsi_2_Fit_noMC, theBeamSpotV));
                    Jpsi_2_ctauErr->push_back(GetcTauErr(Jpsi_2_Vtx_noMC, Jpsi_2_Fit_noMC, theBeamSpotV));
                    Jpsi_2_Chi2->push_back(double(Jpsi_2_Vtx_noMC->chiSquared()));
                    Jpsi_2_ndof->push_back(double(Jpsi_2_Vtx_noMC->degreesOfFreedom()));
                    Jpsi_2_VtxProb->push_back(ChiSquaredProbability((double)(Jpsi_2_Vtx_noMC->chiSquared()), 
                                                                    (double)(Jpsi_2_Vtx_noMC->degreesOfFreedom())));
                    Jpsi_2_px->push_back(Jpsi_2_Fit_noMC->currentState().kinematicParameters().momentum().x());
                    Jpsi_2_py->push_back(Jpsi_2_Fit_noMC->currentState().kinematicParameters().momentum().y());
                    Jpsi_2_pz->push_back(Jpsi_2_Fit_noMC->currentState().kinematicParameters().momentum().z());
                    Jpsi_2_phi->push_back(tmp_pt);
                    Jpsi_2_eta->push_back(tmp_eta);
                    Jpsi_2_pt->push_back(tmp_pt);
                }
                // [TODO] Store the difference between fitted mass with std. mass.
                // [TODO] Store pT eta phi ctau and other kinematic parameters. "As much as possible"
                // [HINT] Only Jpsi ctau required. 
                // [HINT] DR may be useful in BKG suppression. (To deal with pile up. Do it later.)
                else{
                    // Store "error code" -9 for the secondary particles (quarkonia).
                    Jpsi_2_mass->push_back(-9);
                    Jpsi_2_massErr->push_back(-9);
                    Jpsi_2_massDiff->push_back(-9);
                    Jpsi_2_ctau->push_back(-9);
                    Jpsi_2_ctauErr->push_back(-9);
                    Jpsi_2_Chi2->push_back(-9);
                    Jpsi_2_ndof->push_back(-9);
                    Jpsi_2_VtxProb->push_back(-9);
                    Jpsi_2_px->push_back(-9);
                    Jpsi_2_py->push_back(-9);
                    Jpsi_2_pz->push_back(-9);
                    Jpsi_2_phi->push_back(-9);
                    Jpsi_2_eta->push_back(-9);
                    Jpsi_2_pt->push_back(-9);
                }
                if(isValidPhi){
                    getDynamics(Phi_Fit_noMC, tmp_pt, tmp_eta, tmp_phi);
                    Phi_mass->push_back(    Phi_Fit_noMC->currentState().mass());
                    Phi_massDiff->push_back(Phi_Fit_noMC->currentState().mass() - myPhiMass);
                    Phi_massErr->push_back( tmp_Phi_massErr);
                    Phi_Chi2->push_back(double(Phi_Vtx_noMC->chiSquared()));
                    Phi_ndof->push_back(double(Phi_Vtx_noMC->degreesOfFreedom()));
                    Phi_VtxProb->push_back(ChiSquaredProbability((double)(Phi_Vtx_noMC->chiSquared()), 
                                                                 (double)(Phi_Vtx_noMC->degreesOfFreedom())));
                    Phi_px->push_back(Phi_Fit_noMC->currentState().kinematicParameters().momentum().x());
                    Phi_py->push_back(Phi_Fit_noMC->currentState().kinematicParameters().momentum().y());
                    Phi_pz->push_back(Phi_Fit_noMC->currentState().kinematicParameters().momentum().z());
                    Phi_phi->push_back(tmp_pt);
                    Phi_eta->push_back(tmp_eta);
                    Phi_pt->push_back(tmp_pt);
                }
                else{
                    // Store "error code" -9 for the secondary particles (quarkonia).
                    Phi_mass->push_back(-9);
                    Phi_massErr->push_back(-9);
                    Phi_massDiff->push_back(-9);
                    Phi_Chi2->push_back(-9);
                    Phi_ndof->push_back(-9);
                    Phi_VtxProb->push_back(-9);
                    Phi_px->push_back(-9);
                    Phi_py->push_back(-9);
                    Phi_pz->push_back(-9);
                    Phi_phi->push_back(-9);
                    Phi_eta->push_back(-9);
                    Phi_pt->push_back(-9);
                }
            }
        }
    }
    // Currently: Event
	if (Pri_VtxProb->size() > 0 || doMC)
	{
		X_One_Tree_->Fill();
	}

	std::cout << "Finish the part of fitting."  << endl;

	if (Debug_)
	{
	}
    // Reset the vectors [Annotated by Eric Wang, 20240704]
	if (doMC)
	{
		MC_X_px->clear();
		MC_X_py->clear();
		MC_X_pz->clear();
		MC_X_mass->clear();
		MC_X_chg->clear();
		MC_Dau_JpsipdgId->clear();
		MC_Dau_Jpsipx->clear();
		MC_Dau_Jpsipy->clear();
		MC_Dau_Jpsipz->clear();
		MC_Dau_Jpsimass->clear();
		MC_Dau_psi2spdgId->clear();
		MC_Dau_psi2spx->clear();
		MC_Dau_psi2spy->clear();
		MC_Dau_psi2spz->clear();
		MC_Dau_psi2smass->clear();
		MC_Granddau_mu1pdgId->clear();
		MC_Granddau_mu1px->clear();
		MC_Granddau_mu1py->clear();
		MC_Granddau_mu1pz->clear();
		MC_Granddau_mu2pdgId->clear();
		MC_Granddau_mu2px->clear();
		MC_Granddau_mu2py->clear();
		MC_Granddau_mu2pz->clear();
		MC_Granddau_JpsipdgId->clear();
		MC_Granddau_Jpsipx->clear();
		MC_Granddau_Jpsipy->clear();
		MC_Granddau_Jpsipz->clear();
		MC_Granddau_Jpsimass->clear();
		MC_Granddau_pi1pdgId->clear();
		MC_Granddau_pi1px->clear();
		MC_Granddau_pi1py->clear();
		MC_Granddau_pi1pz->clear();
		MC_Granddau_pi2pdgId->clear();
		MC_Granddau_pi2px->clear();
		MC_Granddau_pi2py->clear();
		MC_Granddau_pi2pz->clear();
		MC_Grandgranddau_mu3pdgId->clear();
		MC_Grandgranddau_mu3px->clear();
		MC_Grandgranddau_mu3py->clear();
		MC_Grandgranddau_mu3pz->clear();
		MC_Grandgranddau_mu4pdgId->clear();
		MC_Grandgranddau_mu4px->clear();
		MC_Grandgranddau_mu4py->clear();
		MC_Grandgranddau_mu4pz->clear();

		Match_mu1px->clear();
		Match_mu1py->clear();
		Match_mu1pz->clear();
		Match_mu2px->clear();
		Match_mu2py->clear();
		Match_mu2pz->clear();
		Match_mu3px->clear();
		Match_mu3py->clear();
		Match_mu3pz->clear();
		Match_mu4px->clear();
		Match_mu4py->clear();
		Match_mu4pz->clear();

		Match_pi1px->clear();
		Match_pi1py->clear();
		Match_pi1pz->clear();
		Match_pi2px->clear();
		Match_pi2py->clear();
		Match_pi2pz->clear();
	}

	trigRes->clear();
	trigNames->clear();
	L1TT->clear();
	MatchTriggerNames->clear();
	muIsJpsiTrigMatch->clear();
	muIsUpsTrigMatch->clear();
	runNum = 0;
	evtNum = 0;
	lumiNum = 0;
	nGoodPrimVtx = 0;
	priVtxX = 0;
	priVtxY = 0;
	priVtxZ = 0;
	priVtxXE = 0;
	priVtxYE = 0;
	priVtxZE = 0;
	priVtxChiNorm = 0;
	priVtxChi = 0;
	priVtxCL = 0;
	// mybxlumicorr = 0;
	// myrawbxlumi = 0;
	PriVtxXCorrX->clear();
	PriVtxXCorrY->clear();
	PriVtxXCorrZ->clear();
	PriVtxXCorrEX->clear();
	PriVtxXCorrEY->clear();
	PriVtxXCorrEZ->clear();
	PriVtxXCorrC2->clear();
	PriVtxXCorrCL->clear();

	nMu = 0;
	muPx->clear();
	muPy->clear();
	muPz->clear();
	muD0->clear();
	muD0E->clear();
	muDz->clear();
	muChi2->clear();
	muGlChi2->clear();
	mufHits->clear();
	muFirstBarrel->clear();
	muFirstEndCap->clear();
	muDzVtx->clear();
	muDxyVtx->clear();
	muNDF->clear();
	muGlNDF->clear();
	muPhits->clear();
	muShits->clear();
	muGlMuHits->clear();
	muType->clear();
	muQual->clear();
	muTrack->clear();
	muCharge->clear();
	muIsoratio->clear();
	muIsGoodLooseMuon->clear();
	muIsGoodLooseMuonNew->clear();
	muIsGoodSoftMuonNewIlse->clear();
	muIsGoodSoftMuonNewIlseMod->clear();
	muIsGoodTightMuon->clear();
	munMatchedSeg->clear();
	muMVAMuonID->clear();
	musegmentCompatibility->clear();

	mupulldXdZ_pos_noArb->clear();
	mupulldYdZ_pos_noArb->clear();
	mupulldXdZ_pos_ArbDef->clear();
	mupulldYdZ_pos_ArbDef->clear();
	mupulldXdZ_pos_ArbST->clear();
	mupulldYdZ_pos_ArbST->clear();
	mupulldXdZ_pos_noArb_any->clear();
	mupulldYdZ_pos_noArb_any->clear();

	muIsPatLooseMuon->clear();
	muIsPatTightMuon->clear();
	muIsPatSoftMuon->clear();
	muIsPatMediumMuon->clear();
	muUpsVrtxMatch->clear();
	muL3TriggerMatch->clear();

    Pri_mass->clear();
    Pri_massErr->clear();
    Pri_ctau->clear();
    Pri_ctauErr->clear();
    Pri_Chi2->clear();
    Pri_ndof->clear();
    Pri_VtxProb->clear();
    Pri_px->clear();
    Pri_py->clear();
    Pri_pz->clear();
    Pri_phi->clear();
    Pri_eta->clear();
    Pri_pt->clear();

    Jpsi_1_mass->clear();
    Jpsi_1_massErr->clear();
    Jpsi_1_massDiff->clear();
    Jpsi_1_ctau->clear();
    Jpsi_1_ctauErr->clear();
    Jpsi_1_Chi2->clear();
    Jpsi_1_ndof->clear();
    Jpsi_1_VtxProb->clear();
    Jpsi_1_px->clear();
    Jpsi_1_py->clear();
    Jpsi_1_pz->clear();
    Jpsi_1_phi->clear();
    Jpsi_1_eta->clear();
    Jpsi_1_pt->clear();
    Jpsi_1_mu_1_Idx->clear();
    Jpsi_1_mu_2_Idx->clear();

    Jpsi_2_mass->clear();
    Jpsi_2_massErr->clear();
    Jpsi_2_massDiff->clear();
    Jpsi_2_ctau->clear();
    Jpsi_2_ctauErr->clear();
    Jpsi_2_Chi2->clear();
    Jpsi_2_ndof->clear();
    Jpsi_2_VtxProb->clear();
    Jpsi_2_px->clear();
    Jpsi_2_py->clear();
    Jpsi_2_pz->clear();
    Jpsi_2_phi->clear();
    Jpsi_2_eta->clear();
    Jpsi_2_pt->clear();
    Jpsi_2_mu_1_Idx->clear();
    Jpsi_2_mu_2_Idx->clear();

    Phi_mass->clear();
    Phi_massErr->clear();
    Phi_massDiff->clear();
    Phi_Chi2->clear();
    Phi_ndof->clear();
    Phi_VtxProb->clear();
    Phi_px->clear();
    Phi_py->clear();
    Phi_pz->clear();
    Phi_phi->clear();
    Phi_eta->clear();
    Phi_pt->clear();
    Phi_pi_1_Idx->clear();
    Phi_pi_2_Idx->clear();
} // analyze
// 

/******************************************************************************
 * [Name of function]  
 *      tracksToMuonPair
 * [Description]  
 *      Construct muons from tracks.
 *      Assuming muon mass and mass error as PDG 2023 values.
 *      Adds reconstructed muons to the arg_MuonResults.
 * [Parameters]
 *      vector<RefCountedKinematicParticle>&        arg_MuonResults
 *          - The vector to which reconstructed muons are added.
 *      KinematicParticleFactoryFromTransientTrack& arg_MuFactory
 *          - The class used to reconstruct muons.
 *      const MagneticField&                        arg_bField,
 *          - Magnetic field used in reconstruction.
 *      const TrackRef&                             arg_Trk1, arg_Trk2  
 *          - Tracks identified as muons.      
 * [Return value]
 *      (void)
 * [Note]
 *          
******************************************************************************/

void MultiLepPAT::getDynamics(double  arg_mass, double  arg_px,  double  arg_py, double arg_pz,
                              double& res_pt,   double& res_eta, double& res_phi){
    TLorentzVector myParticle;
    myParticle.SetXYZM(arg_px, arg_py, arg_pz, arg_mass);
    res_pt  = myParticle.Pt();
    res_eta = myParticle.Eta();
    res_phi = myParticle.Phi();
}

/******************************************************************************
 * [Name of function]  
 *      tracksToMuonPair
 * [Description]  
 *      Construct muons from tracks.
 *      Assuming muon mass and mass error as PDG 2023 values.
 *      Adds reconstructed muons to the arg_MuonResults.
 * [Parameters]
 *      vector<RefCountedKinematicParticle>&        arg_MuonResults
 *          - The vector to which reconstructed muons are added.
 *      KinematicParticleFactoryFromTransientTrack& arg_MuFactory
 *          - The class used to reconstruct muons.
 *      const MagneticField&                        arg_bField,
 *          - Magnetic field used in reconstruction.
 *      const TrackRef&                             arg_Trk1, arg_Trk2  
 *          - Tracks identified as muons.      
 * [Return value]
 *      (void)
 * [Note]
 *          
******************************************************************************/
void MultiLepPAT::getDynamics(const RefCountedKinematicParticle& arg_Part,
                              double& res_pt,   double& res_eta, double& res_phi){
    getDynamics(arg_Part->currentState().mass(), 
                arg_Part->currentState().kinematicParameters().momentum().x(),
                arg_Part->currentState().kinematicParameters().momentum().y(),
                arg_Part->currentState().kinematicParameters().momentum().z(),
                res_pt, res_eta, res_phi);
}

/******************************************************************************
 * [Name of function]  
 *      tracksToMuonPair
 * [Description]  
 *      Construct muons from tracks.
 *      Assuming muon mass and mass error as PDG 2023 values.
 *      Adds reconstructed muons to the arg_MuonResults.
 * [Parameters]
 *      vector<RefCountedKinematicParticle>&        arg_MuonResults
 *          - The vector to which reconstructed muons are added.
 *      KinematicParticleFactoryFromTransientTrack& arg_MuFactory
 *          - The class used to reconstruct muons.
 *      const MagneticField&                        arg_bField,
 *          - Magnetic field used in reconstruction.
 *      const TrackRef&                             arg_Trk1, arg_Trk2  
 *          - Tracks identified as muons.      
 * [Return value]
 *      (void)
 * [Note]
 *          
******************************************************************************/
void MultiLepPAT::tracksToMuonPair(vector<RefCountedKinematicParticle>&        arg_MuonResults,
                                   KinematicParticleFactoryFromTransientTrack& arg_MuFactory,
                                   const MagneticField&                        arg_bField,
                                   const TrackRef arg_Trk1,     const TrackRef arg_Trk2        ){
    TransientTrack transTrk1(arg_Trk1, &(arg_bField));
    TransientTrack transTrk2(arg_Trk2, &(arg_bField));
    // Parameters for muon
    ParticleMass muMass = myMuMass;
    float muMassSigma   = myMuMassErr;
    float chi2 = 0.;
	float ndof = 0.;
    arg_MuonResults.push_back(arg_MuFactory.particle(transTrk1, muMass, chi2, ndof, muMassSigma));
    arg_MuonResults.push_back(arg_MuFactory.particle(transTrk2, muMass, chi2, ndof, muMassSigma));
    return ;
}

/******************************************************************************
 * [Name of function]  
 *      particlesToVtx
 * [Description]  
 *      Construct muons from tracks.
 *      Assuming muon mass and mass error as PDG 2023 values.
 *      Adds reconstructed muons to the arg_FromParticles.
 * [Parameters]
 *      vector<RefCountedKinematicParticle>&        arg_FromParticles
 *          - The vector to which reconstructed particles are added.
 * [Return value]
 *      (void)
 * [Note]
 *      A "silent" version of fitting particles to vertex. No error message
 *      will be printed in case of failed fitting.
******************************************************************************/

bool MultiLepPAT::particlesToVtx(const vector<RefCountedKinematicParticle>&  arg_FromParticles){
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree vertexFitTree;
    bool fitError = false;
    try{
        vertexFitTree = fitter.fit(arg_FromParticles);
    }catch(...){
        fitError = true;
    }
	if (fitError || !vertexFitTree->isValid()){
        return false;
    }
	RefCountedKinematicVertex vFit_vertex_noMC = vertexFitTree->currentDecayVertex();
	double vtxprob = ChiSquaredProbability((double)(vFit_vertex_noMC->chiSquared()), (double)(vFit_vertex_noMC->degreesOfFreedom()));
    if (vtxprob < 0.01){
        return false;
    }
    return true;
}

/******************************************************************************
 * [Name of function]  
 *      particlesToVtx
 * [Description]  
 *      Construct muons from tracks.
 *      Assuming muon mass and mass error as PDG 2023 values.
 *      Adds reconstructed muons to the arg_FromParticles.
 * [Parameters]
 *      vector<RefCountedKinematicParticle>&        arg_FromParticles
 *          - The vector to which reconstructed particles are added.
 *      const string&                               arg_Message  
 *          - The message to be displayed in case of error.
 * [Return value]
 *      (void)
 * [Note]
 *      This definition uses an "implicit" VertexFitter and KinematicTree. 
******************************************************************************/

bool MultiLepPAT::particlesToVtx(const vector<RefCountedKinematicParticle>&  arg_FromParticles,
                                 const string&                               arg_Message){
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree vertexFitTree;
    bool fitError = false;
    try{
        vertexFitTree = fitter.fit(arg_FromParticles);
    }catch(...){
        fitError = true;
        std::cout << "[Fit Error] " << arg_Message <<  std::endl;
    }
	if (fitError || !vertexFitTree->isValid()){
        return false;
    }
    RefCountedKinematicVertex vFit_vertex_noMC = vertexFitTree->currentDecayVertex();
	double vtxprob = ChiSquaredProbability((double)(vFit_vertex_noMC->chiSquared()), (double)(vFit_vertex_noMC->degreesOfFreedom()));
    if (vtxprob < 0.01){
        return false;
    }
    return true;
}

/******************************************************************************
 * [Name of function]  
 *      particlesToVtx
 * [Description]  
 *      Construct muons from tracks.
 *      Assuming muon mass and mass error as PDG 2023 values.
 *      Adds reconstructed muons to the arg_FromParticles.
 * [Parameters]
 *      vector<RefCountedKinematicParticle>&        arg_FromParticles
 *          - The vector to which reconstructed particles are added.
 *      const string&                               arg_Message  
 *          - The message to be displayed in case of error.
 *      RefCountedKinematicTree&                    arg_VertexFitTree
 *          - The KinematicTree to which the vertex fit is added.    
 * [Return value]
 *      (void)
 * [Note]
 *      This definition uses an "explicit" KinematicTree.
 *      The KinematicTree is passed as an argument and is modified after call.
******************************************************************************/

bool MultiLepPAT::particlesToVtx(RefCountedKinematicTree&                    arg_VertexFitTree,
                                 const vector<RefCountedKinematicParticle>&  arg_FromParticles,
                                 const string&                               arg_Message){
    KinematicParticleVertexFitter fitter;
    bool fitError = false;
    try{
        arg_VertexFitTree = fitter.fit(arg_FromParticles);
    }catch(...){
        fitError = true;
        std::cout << "[Fit Error] " << arg_Message <<  std::endl;
    }
	if (fitError || !arg_VertexFitTree->isValid()){
        return false;
    }
    RefCountedKinematicVertex vFit_vertex_noMC = arg_VertexFitTree->currentDecayVertex();
	double vtxprob = ChiSquaredProbability((double)(vFit_vertex_noMC->chiSquared()), (double)(vFit_vertex_noMC->degreesOfFreedom()));
    if (vtxprob < 0.01){
        return false;
    }
    return true;
}

/******************************************************************************
 * [Name of function]  
 *      extractFitRes
 * [Description]
 *      Extract kinematic parameters and other results from a KinematicTree.
 *      Calls movePointerToTheTop() .
 * [Parameters]
 *      RefCountedKinematicTree&     arg_VtxTree
 *          - The KinematicTree constructed from fitting.
 *      RefCountedKinematicParticle& res_Part
 *          - The mother particle extracted from arg_VtxTree.
 *      RefCountedKinematicVertex&   res_Vtx
 *          - The primary vertex extracted from arg_VtxTree.
 *      KinematicParameters&         res_Param
 *          - The kinematic parameters of the mother particle.
 *      double&                      res_MassErr
 *          - The mass error of the mother particle.
 * [Return value]
 *      (bool)
 *          - True if the mass error squared is non-negative.
 * [Note]
 *      - Used when the resulting dynamics is important.
 *      - Requires the KinematicTree to be valid.
 *      - Requires "explicit" particle and vertex.
 *      - The mass error is set to -9 if the mass error squared is negative.
******************************************************************************/

bool MultiLepPAT::extractFitRes(RefCountedKinematicTree&     arg_VtxTree,
                                RefCountedKinematicParticle& res_Part,
                                RefCountedKinematicVertex&   res_Vtx,
                                KinematicParameters&         res_Param,
                                double&                      res_MassErr){
    double tmp_MassErr2 = 0.0;
    arg_VtxTree->movePointerToTheTop();
    // Extract particle and vertex.
    res_Part  = arg_VtxTree->currentParticle();
    res_Vtx   = arg_VtxTree->currentDecayVertex();
    // Obtain mass error squared and other parameters for the vertex.
    res_Param    = res_Part->currentState().kinematicParameters();
    tmp_MassErr2 = res_Part->currentState().kinematicParametersError().matrix()(6, 6);
    // Judge if the fit have been a good fit.
    if(tmp_MassErr2 < 0.0){
        res_MassErr = -9;
    }
    else{
        res_MassErr = std::sqrt(tmp_MassErr2);
    }
    return (res_MassErr >= 0.0);
}

/******************************************************************************
 * [Name of function]  
 *      extractFitRes
 * [Description]
 *      Extract kinematic parameters and other results from a KinematicTree.
 *      Calls movePointerToTheTop() .
 * [Parameters]
 *      RefCountedKinematicTree&     arg_VtxTree
 *          - The KinematicTree constructed from fitting.
 *      RefCountedKinematicParticle& res_Part
 *          - The mother particle extracted from arg_VtxTree.
 *      RefCountedKinematicVertex&   res_Vtx
 *          - The primary vertex extracted from arg_VtxTree.
 *      double&                      res_MassErr
 *          - The mass error of the mother particle.
 * [Return value]
 *      (bool)
 *          - True if the mass error squared is non-negative.
 * [Note]
 *      - Used when the resulting dynamics is important.
 *      - Requires the KinematicTree to be valid.
 *      - Requires "explicit" particle and vertex.
 *      - The mass error is set to -9 if the mass error squared is negative.
******************************************************************************/

bool MultiLepPAT::extractFitRes(RefCountedKinematicTree&     arg_VtxTree,
                                RefCountedKinematicParticle& res_Part,
                                RefCountedKinematicVertex&   res_Vtx,
                                double&                      res_MassErr){
    double tmp_MassErr2 = 0.0;
    arg_VtxTree->movePointerToTheTop();
    // Extract particle and vertex.
    res_Part  = arg_VtxTree->currentParticle();
    res_Vtx   = arg_VtxTree->currentDecayVertex();
    // Obtain mass error squared and other parameters for the vertex.
    tmp_MassErr2 = res_Part->currentState().kinematicParametersError().matrix()(6, 6);
    // Judge if the fit have been a good fit.
    if(tmp_MassErr2 < 0.0){
        res_MassErr = -9;
    }
    else{
        res_MassErr = std::sqrt(tmp_MassErr2);
    }
    return (res_MassErr >= 0.0);
}

/******************************************************************************
 * [Name of function]  
 *      extractFitRes
 * [Description]
 *      Extract kinematic parameters and other results from a KinematicTree.
 *      Calls movePointerToTheTop() .
 * [Parameters]
 *      RefCountedKinematicTree&     arg_VtxTree
 *          - The KinematicTree constructed from fitting.
 *      RefCountedKinematicVertex&   res_Vtx
 *          - The primary vertex extracted from arg_VtxTree.
 *      double&                      res_VtxProb
 *          - The vertex probability deduced from res_Vtx parameters.
 * [Return value]
 *      (bool, always true)
 * [Note]
 *      - Used when only the vertex probability is important.
 *      - Requires the KinematicTree to be valid.
 *      - Requires "explicit" vertex.
******************************************************************************/

bool MultiLepPAT::extractFitRes(RefCountedKinematicTree&     arg_VtxTree,
                                RefCountedKinematicVertex&   res_Vtx,
                                double&                      res_VtxProb){
    arg_VtxTree->movePointerToTheTop();
    // Extract particle and vertex.
    res_Vtx   = arg_VtxTree->currentDecayVertex();
    // Obtain mass error squared and other parameters for the vertex.
    res_VtxProb = ChiSquaredProbability((double)(res_Vtx->chiSquared()), 
                                        (double)(res_Vtx->degreesOfFreedom()));
    return true;
}

/******************************************************************************
 * [Name of function]  
 *      isOverlapPair
 * [Description]
 *      Check if two muon pairs overlap from the muon indices.
 * [Parameters]
 *      const muList_t& arg_MuonPair1, arg_MuonPair2
 *          - The muon pairs to be compared.
 * [Return value]
 *      (bool)
 *          - True if the two muon pairs do overlap.
 * [Note]
 *      - Uses muList_t defined in the header.
******************************************************************************/

bool MultiLepPAT::isOverlapPair(const muList_t& arg_MuonPair1, 
                                const muList_t& arg_MuonPair2 ){
    return (arg_MuonPair1.second[0] == arg_MuonPair2.second[0] || 
            arg_MuonPair1.second[0] == arg_MuonPair2.second[1] || 
            arg_MuonPair1.second[1] == arg_MuonPair2.second[0] || 
            arg_MuonPair1.second[1] == arg_MuonPair2.second[1]   );
}

/******************************************************************************
 * [Name of function]  
 *      fitResEval
 * [Description]
 *      Evaluate a fitting result from the mass and mass error.
 * [Parameters]
 *      double arg_mass_Jpsi_1, arg_mass_Jpsi_2, arg_mass_Ups
 *          - The mass of the particles.
 *      double arg_massErr_Jpsi_1, arg_massErr_Jpsi_2, arg_massErr_Ups
 *          - The mass error of the particles.
 * [Return value]
 *      (double)
 *          - Evaluation result.
 * [Note]
 *      To be used in resolving "multi-candidate" events.
******************************************************************************/

double MultiLepPAT::fitResEval(double arg_massDiff_Jpsi_1, double arg_massErr_Jpsi_1,
                               double arg_massDiff_Jpsi_2, double arg_massErr_Jpsi_2,
                               double arg_massDiff_Ups,    double arg_massErr_Ups   ){
    return arg_massDiff_Jpsi_1 * arg_massDiff_Jpsi_1 / (arg_massErr_Jpsi_1 * arg_massErr_Jpsi_1) +
           arg_massDiff_Jpsi_2 * arg_massDiff_Jpsi_2 / (arg_massErr_Jpsi_2 * arg_massErr_Jpsi_2) +
           arg_massDiff_Ups    * arg_massDiff_Ups    / (arg_massErr_Ups    * arg_massErr_Ups   ) ;
}

// ------------ method called once each job just before starting event loop  ------------
void MultiLepPAT::beginRun(edm::Run const &iRun, edm::EventSetup const &iSetup)
{
	// bool changed = true;
	// proccessName_="HLT";
	// hltConfig_.init(iRun,iSetup,proccessName_,changed);
}

void MultiLepPAT::beginJob()
{
	edm::Service<TFileService> fs;

	// estree_ = fs->make<TTree>("eventSummary", "General Event Summary");
	X_One_Tree_ = fs->make<TTree>("X_data", "X(3872) Data");

	X_One_Tree_->Branch("TrigRes", &trigRes);
	X_One_Tree_->Branch("TrigNames", &trigNames);
	X_One_Tree_->Branch("MatchTriggerNames", &MatchTriggerNames);
	X_One_Tree_->Branch("L1TrigRes", &L1TT);

	X_One_Tree_->Branch("evtNum", &evtNum, "evtNum/i");
	X_One_Tree_->Branch("runNum", &runNum, "runNum/i");
	X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
	X_One_Tree_->Branch("nGoodPrimVtx", &nGoodPrimVtx, "nGoodPrimVtx/i");

	// inst. lumi is here
	//X_One_Tree_->Branch("mybxlumicorr", &mybxlumicorr, "mybxlumicorr/f");
	//X_One_Tree_->Branch("myrawbxlumi", &myrawbxlumi, "myrawbxlumi/f");

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

	X_One_Tree_->Branch("nMu", &nMu, "nMu/i");
	X_One_Tree_->Branch("muPx", &muPx);
	X_One_Tree_->Branch("muPy", &muPy);
	X_One_Tree_->Branch("muPz", &muPz);
	X_One_Tree_->Branch("muD0", &muD0);
	X_One_Tree_->Branch("muD0E", &muD0E);
	X_One_Tree_->Branch("muDz", &muDz);
	X_One_Tree_->Branch("muChi2", &muChi2);
	X_One_Tree_->Branch("muGlChi2", &muGlChi2);
	X_One_Tree_->Branch("mufHits", &mufHits);
	X_One_Tree_->Branch("muFirstBarrel", &muFirstBarrel);
	X_One_Tree_->Branch("muFirstEndCap", &muFirstEndCap);
	X_One_Tree_->Branch("muDzVtx", &muDzVtx);
	X_One_Tree_->Branch("muDxyVtx", &muDxyVtx);
	X_One_Tree_->Branch("muNDF", &muNDF);
	X_One_Tree_->Branch("muGlNDF", &muGlNDF);
	X_One_Tree_->Branch("muPhits", &muPhits);
	X_One_Tree_->Branch("muShits", &muShits);
	X_One_Tree_->Branch("muGlMuHits", &muGlMuHits);
	X_One_Tree_->Branch("muType", &muType);
	X_One_Tree_->Branch("muQual", &muQual);
	X_One_Tree_->Branch("muTrack", &muTrack);
	X_One_Tree_->Branch("muCharge", &muCharge);
	X_One_Tree_->Branch("muIsoratio", &muIsoratio);
	X_One_Tree_->Branch("munMatchedSeg", &munMatchedSeg);
	X_One_Tree_->Branch("muIsGoodSoftMuonNewIlseMod", &muIsGoodSoftMuonNewIlseMod);
	X_One_Tree_->Branch("muIsGoodSoftMuonNewIlse", &muIsGoodSoftMuonNewIlse);
	X_One_Tree_->Branch("muIsGoodLooseMuonNew", &muIsGoodLooseMuonNew);
	X_One_Tree_->Branch("muIsGoodLooseMuon", &muIsGoodLooseMuon);
	X_One_Tree_->Branch("muIsGoodTightMuon", &muIsGoodTightMuon);
	X_One_Tree_->Branch("muJpsiFilterRes", &muJpsiFilterRes);

	X_One_Tree_->Branch("muIsPatLooseMuon", &muIsPatLooseMuon);
	X_One_Tree_->Branch("muIsPatTightMuon", &muIsPatTightMuon);
	X_One_Tree_->Branch("muIsPatSoftMuon", &muIsPatSoftMuon);
	X_One_Tree_->Branch("muIsPatMediumMuon", &muIsPatMediumMuon);

	X_One_Tree_->Branch("muIsJpsiTrigMatch", &muIsJpsiTrigMatch);
	X_One_Tree_->Branch("muIsUpsTrigMatch", &muIsUpsTrigMatch);
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

	X_One_Tree_->Branch("muUpsVrtxMatch", &muUpsVrtxMatch);
	X_One_Tree_->Branch("muL3TriggerMatch", &muL3TriggerMatch);

    X_One_Tree_->Branch("Jpsi_1_mass", &Jpsi_1_mass);
    X_One_Tree_->Branch("Jpsi_1_massErr", &Jpsi_1_massErr);
    X_One_Tree_->Branch("Jpsi_1_massDiff", &Jpsi_1_massDiff);
    X_One_Tree_->Branch("Jpsi_1_ctau", &Jpsi_1_ctau);
    X_One_Tree_->Branch("Jpsi_1_ctauErr", &Jpsi_1_ctauErr);
    X_One_Tree_->Branch("Jpsi_1_Chi2", &Jpsi_1_Chi2);
    X_One_Tree_->Branch("Jpsi_1_ndof", &Jpsi_1_ndof);
    X_One_Tree_->Branch("Jpsi_1_VtxProb", &Jpsi_1_VtxProb);
    X_One_Tree_->Branch("Jpsi_1_px", &Jpsi_1_px);
    X_One_Tree_->Branch("Jpsi_1_py", &Jpsi_1_py);
    X_One_Tree_->Branch("Jpsi_1_pz", &Jpsi_1_pz);
    X_One_Tree_->Branch("Jpsi_1_phi", &Jpsi_1_phi);
    X_One_Tree_->Branch("Jpsi_1_eta", &Jpsi_1_eta);
    X_One_Tree_->Branch("Jpsi_1_pt", &Jpsi_1_pt);
    X_One_Tree_->Branch("Jpsi_1_mu_1_Idx", &Jpsi_1_mu_1_Idx);
    X_One_Tree_->Branch("Jpsi_1_mu_2_Idx", &Jpsi_1_mu_2_Idx);

    X_One_Tree_->Branch("Jpsi_2_mass", &Jpsi_2_mass);
    X_One_Tree_->Branch("Jpsi_2_massErr", &Jpsi_2_massErr);
    X_One_Tree_->Branch("Jpsi_2_massDiff", &Jpsi_2_massDiff);
    X_One_Tree_->Branch("Jpsi_2_ctau", &Jpsi_2_ctau);
    X_One_Tree_->Branch("Jpsi_2_ctauErr", &Jpsi_2_ctauErr);
    X_One_Tree_->Branch("Jpsi_2_Chi2", &Jpsi_2_Chi2);
    X_One_Tree_->Branch("Jpsi_2_ndof", &Jpsi_2_ndof);
    X_One_Tree_->Branch("Jpsi_2_VtxProb", &Jpsi_2_VtxProb);
    X_One_Tree_->Branch("Jpsi_2_px", &Jpsi_2_px);
    X_One_Tree_->Branch("Jpsi_2_py", &Jpsi_2_py);
    X_One_Tree_->Branch("Jpsi_2_pz", &Jpsi_2_pz);
    X_One_Tree_->Branch("Jpsi_2_phi", &Jpsi_2_phi);
    X_One_Tree_->Branch("Jpsi_2_eta", &Jpsi_2_eta);
    X_One_Tree_->Branch("Jpsi_2_pt", &Jpsi_2_pt);
    X_One_Tree_->Branch("Jpsi_2_mu_1_Idx", &Jpsi_2_mu_1_Idx);
    X_One_Tree_->Branch("Jpsi_2_mu_2_Idx", &Jpsi_2_mu_2_Idx);

    X_One_Tree_->Branch("Phi_mass", &Phi_mass);
    X_One_Tree_->Branch("Phi_massErr", &Phi_massErr);
    X_One_Tree_->Branch("Phi_massDiff", &Phi_massDiff);
    X_One_Tree_->Branch("Phi_Chi2", &Phi_Chi2);
    X_One_Tree_->Branch("Phi_ndof", &Phi_ndof);
    X_One_Tree_->Branch("Phi_VtxProb", &Phi_VtxProb);
    X_One_Tree_->Branch("Phi_px", &Phi_px);
    X_One_Tree_->Branch("Phi_py", &Phi_py);
    X_One_Tree_->Branch("Phi_pz", &Phi_pz);
    X_One_Tree_->Branch("Phi_phi", &Phi_phi);
    X_One_Tree_->Branch("Phi_eta", &Phi_eta);
    X_One_Tree_->Branch("Phi_pt", &Phi_pt);
    X_One_Tree_->Branch("Phi_pi_1_Idx", &Phi_pi_1_Idx);
    X_One_Tree_->Branch("Phi_pi_2_Idx", &Phi_pi_2_Idx);

    X_One_Tree_->Branch("Pri_mass", &Pri_mass);
    X_One_Tree_->Branch("Pri_massErr", &Pri_massErr);
    X_One_Tree_->Branch("Pri_ctau", &Pri_ctau);
    X_One_Tree_->Branch("Pri_ctauErr", &Pri_ctauErr);
    X_One_Tree_->Branch("Pri_Chi2", &Pri_Chi2);
    X_One_Tree_->Branch("Pri_ndof", &Pri_ndof);
    X_One_Tree_->Branch("Pri_VtxProb", &Pri_VtxProb);
    X_One_Tree_->Branch("Pri_px", &Pri_px);
    X_One_Tree_->Branch("Pri_py", &Pri_py);
    X_One_Tree_->Branch("Pri_pz", &Pri_pz);
    X_One_Tree_->Branch("Pri_phi", &Pri_phi);
    X_One_Tree_->Branch("Pri_eta", &Pri_eta);
    X_One_Tree_->Branch("Pri_pt", &Pri_pt);



	if (doMC)
	{
		X_One_Tree_->Branch("MC_X_px", &MC_X_px);
		X_One_Tree_->Branch("MC_X_py", &MC_X_py);
		X_One_Tree_->Branch("MC_X_pz", &MC_X_pz);
		X_One_Tree_->Branch("MC_X_mass", &MC_X_mass);
		X_One_Tree_->Branch("MC_X_chg", &MC_X_chg);
		X_One_Tree_->Branch("MC_Dau_JpsipdgId", &MC_Dau_JpsipdgId);
		X_One_Tree_->Branch("MC_Dau_Jpsipx", &MC_Dau_Jpsipx);
		X_One_Tree_->Branch("MC_Dau_Jpsipy", &MC_Dau_Jpsipy);
		X_One_Tree_->Branch("MC_Dau_Jpsipz", &MC_Dau_Jpsipz);
		X_One_Tree_->Branch("MC_Dau_Jpsimass", &MC_Dau_Jpsimass);
		X_One_Tree_->Branch("MC_Dau_psi2spdgId", &MC_Dau_psi2spdgId);
		X_One_Tree_->Branch("MC_Dau_psi2spx", &MC_Dau_psi2spx);
		X_One_Tree_->Branch("MC_Dau_psi2spy", &MC_Dau_psi2spy);
		X_One_Tree_->Branch("MC_Dau_psi2spz", &MC_Dau_psi2spz);
		X_One_Tree_->Branch("MC_Dau_psi2smass", &MC_Dau_psi2smass);
		X_One_Tree_->Branch("MC_Granddau_mu1pdgId", &MC_Granddau_mu1pdgId);
		X_One_Tree_->Branch("MC_Granddau_mu1px", &MC_Granddau_mu1px);
		X_One_Tree_->Branch("MC_Granddau_mu1py", &MC_Granddau_mu1py);
		X_One_Tree_->Branch("MC_Granddau_mu1pz", &MC_Granddau_mu1pz);
		X_One_Tree_->Branch("MC_Granddau_mu2pdgId", &MC_Granddau_mu2pdgId);
		X_One_Tree_->Branch("MC_Granddau_mu2px", &MC_Granddau_mu2px);
		X_One_Tree_->Branch("MC_Granddau_mu2py", &MC_Granddau_mu2py);
		X_One_Tree_->Branch("MC_Granddau_mu2pz", &MC_Granddau_mu2pz);
		X_One_Tree_->Branch("MC_Granddau_JpsipdgId", &MC_Granddau_JpsipdgId);
		X_One_Tree_->Branch("MC_Granddau_Jpsipx", &MC_Granddau_Jpsipx);
		X_One_Tree_->Branch("MC_Granddau_Jpsipy", &MC_Granddau_Jpsipy);
		X_One_Tree_->Branch("MC_Granddau_Jpsipz", &MC_Granddau_Jpsipz);
		X_One_Tree_->Branch("MC_Granddau_Jpsimass", &MC_Granddau_Jpsimass);
		X_One_Tree_->Branch("MC_Granddau_pi1pdgId", &MC_Granddau_pi1pdgId);
		X_One_Tree_->Branch("MC_Granddau_pi1px", &MC_Granddau_pi1px);
		X_One_Tree_->Branch("MC_Granddau_pi1py", &MC_Granddau_pi1py);
		X_One_Tree_->Branch("MC_Granddau_pi1pz", &MC_Granddau_pi1pz);
		X_One_Tree_->Branch("MC_Granddau_pi2pdgId", &MC_Granddau_pi2pdgId);
		X_One_Tree_->Branch("MC_Granddau_pi2px", &MC_Granddau_pi2px);
		X_One_Tree_->Branch("MC_Granddau_pi2py", &MC_Granddau_pi2py);
		X_One_Tree_->Branch("MC_Granddau_pi2pz", &MC_Granddau_pi2pz);
		X_One_Tree_->Branch("MC_Grandgranddau_mu3pdgId", &MC_Grandgranddau_mu3pdgId);
		X_One_Tree_->Branch("MC_Grandgranddau_mu3px", &MC_Grandgranddau_mu3px);
		X_One_Tree_->Branch("MC_Grandgranddau_mu3py", &MC_Grandgranddau_mu3py);
		X_One_Tree_->Branch("MC_Grandgranddau_mu3pz", &MC_Grandgranddau_mu3pz);
		X_One_Tree_->Branch("MC_Grandgranddau_mu4pdgId", &MC_Grandgranddau_mu4pdgId);
		X_One_Tree_->Branch("MC_Grandgranddau_mu4px", &MC_Grandgranddau_mu4px);
		X_One_Tree_->Branch("MC_Grandgranddau_mu4py", &MC_Grandgranddau_mu4py);
		X_One_Tree_->Branch("MC_Grandgranddau_mu4pz", &MC_Grandgranddau_mu4pz);
	} // if(doMC)
} // begin Job

// Moved some functions from the header file to the source file [Modified by Eric Wang, 20240705]

// CTau calculation from fitted vertex. [Annotated by Eric Wang, 20240705]
double MultiLepPAT::GetcTau(RefCountedKinematicVertex&   decayVrtx, 
                            RefCountedKinematicParticle& kinePart, 
                            Vertex&                             bs ){	
    TVector3 vtx;
    TVector3 pvtx;
    vtx.SetXYZ((*decayVrtx).position().x(), (*decayVrtx).position().y(), 0);
    pvtx.SetXYZ(bs.position().x(), bs.position().y(), 0);
    TVector3 pperp(kinePart->currentState().globalMomentum().x(),
        kinePart->currentState().globalMomentum().y(), 0);
    TVector3 vdiff = vtx - pvtx;

    //Jinfeng 10.3
    GlobalError DecayErr = (*decayVrtx).error();
    GlobalError PrimaryErr = bs.error();

    AlgebraicSymMatrix33 LxyErrMatrix = DecayErr.matrix() + PrimaryErr.matrix();
    AlgebraicVector3 vdiff_;
    vdiff_[0] = vdiff.x();
    vdiff_[1] = vdiff.y();
    vdiff_[2] = 0;

    double LxyErr2 = ROOT::Math::Similarity(LxyErrMatrix, vdiff_);
    double LxyPVErr = 0;
    // double LxyPVSig = 0;

    double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
    double LxyPV = vdiff.Dot(pperp) / pperp.Mag();
    double Lxy = vdiff.Perp();
    double ctauPV = LxyPV * kinePart->currentState().mass() / pperp.Perp();

    if (Lxy != 0)
        LxyPVErr = TMath::Sqrt(LxyErr2) / Lxy * cosAlpha;

    if (LxyPVErr != 0)
        double LxyPVSig = LxyPV / LxyPVErr;

    //Jinfeng
    // if (MassConstraint) {
    //     temp_Lxy_MC = Lxy;
    //     temp_LxyPV_MC = LxyPV;
    //     temp_CosAlpha_MC = cosAlpha;
    //     temp_LxyPVErr_MC = LxyPVErr;
    //     temp_LxyPVSig_MC = LxyPVSig;
    // }
    // else {
    //     temp_Lxy_noMC = Lxy;
    //     temp_LxyPV_noMC = LxyPV;
    //     temp_CosAlpha_noMC = cosAlpha;
    //     temp_LxyPVErr_noMC = LxyPVErr;
    //     temp_LxyPVSig_noMC = LxyPVSig;
    // }

    return ctauPV;    
}

// CTau error calculation from fitted vertex. [Annotated by Eric Wang, 20240705]
double MultiLepPAT:: GetcTauErr( RefCountedKinematicVertex& decayVrtx, 
                                 RefCountedKinematicParticle& kinePart, 
                                 Vertex& bs                              ){       
    TVector3 pperp(kinePart->currentState().globalMomentum().x(),
	               kinePart->currentState().globalMomentum().y(), 
                   0                                              );
    AlgebraicVector vpperp(3);
    vpperp[0] = pperp.x();
    vpperp[1] = pperp.y();
    vpperp[2] = 0.;

    GlobalError v1e = (*decayVrtx).error();
    GlobalError v2e = bs.error();
    AlgebraicSymMatrix vXYe = asHepMatrix(v1e.matrix()) + asHepMatrix(v2e.matrix());
    double ctauErrPV = sqrt(vXYe.similarity(vpperp)) * kinePart->currentState().mass() / (pperp.Perp2());

    return ctauErrPV;    
}

// deltaR calculation from usual eta and phi [Annotated by Eric Wang, 20240705]
double MultiLepPAT::deltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi >   M_PI) dphi -= 2*M_PI;
    while (dphi <= -M_PI) dphi += 2*M_PI;
    return sqrt(deta*deta + dphi*dphi);
}


/******************************************************************************
 * [Name of function]  
 *      getAllTriggers
 * [Description]   
 *      Get all triggers from the trigger results and categorize them
 * [Parameters]
 *      HLTresult                       trigger results from EDM
 * [Return value]
 *      (void)
 * [Note]
 *      [Eric Wang, 20240705]
 *          
******************************************************************************/
void MultiLepPAT::getAllTriggers(const edm::Handle<edm::TriggerResults>&     HLTresult){
}



/******************************************************************************
 * [Name of function]  
 *      muonMatchTrigType
 * [Description]   
 *      
 * [Parameters]
 *      HLTresult                       trigger results from EDM
 * [Return value]
 *      (void)
 * [Note]
 *      [Eric Wang, 20240705]
 *          
******************************************************************************/
bool MultiLepPAT::muonMatchTrigType(const edm::View<pat::Muon>::const_iterator& muIter,
                                    const vector<string>& trigNames, 
                                          trigType        type                          ){
    return false;
}


// ------------ method called once each job just after ending the event loop  ------------
void MultiLepPAT::endJob()
{
	X_One_Tree_->GetDirectory()->cd();
	X_One_Tree_->Write();
}

// define this as a plug-in
DEFINE_FWK_MODULE(MultiLepPAT);
