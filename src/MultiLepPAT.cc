// system include files
#include <memory>
#include <regex>
#include "TLorentzVector.h"
// user include files
#include "../interface/MultiLepPAT.h"

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

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TrackingTools/IPTools/interface/IPTools.h"


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
	  Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output", false)),
	  Chi_Track_(iConfig.getUntrackedParameter<double>("Chi2NDF_Track", 10)),
	  X_One_Tree_(0),

	  runNum(0), evtNum(0), lumiNum(0), nGoodPrimVtx(0),
	  TrigRes(0), TrigNames(0), L1TrigRes(0), MatchJpsiTrigNames(0),

	  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXErr(0), priVtxYErr(0), priVtxZErr(0), priVtxChiNorm(0), priVtxChi(0), priVtxCL(0),
	  priVtxXCorrY(0), priVtxXCorrZ(0), priVtxYCorrZ(0),

	  nMu(0),
	  muPx(0), muPy(0), muPz(0), muD0BS(0), muD0EBS(0), muD3dBS(0), muD3dEBS(0),
	  muD0PV(0), muD0EPV(0), muDzPV(0), muDzEPV(0),
	  muCharge(0), muTrackIso(0),
	  muIsPatLooseMuon(0), muIsPatTightMuon(0), muIsPatSoftMuon(0), muIsPatMediumMuon(0), muJpsiFilterRes(0),
	  X_mu1Id(0), X_mu2Id(0), X_mu3Id(0), X_mu4Id(0),
	  X_mass(0), X_VtxProb(0), X_Chi2(0), X_ndof(0), X_px(0), X_py(0), X_pz(0), X_massErr(0),
	  X_JPiPi_mass(0), X_JPiPi_VtxProb(0), X_JPiPi_Chi2(0), X_JPiPi_ndof(0), X_JPiPi_px(0), X_JPiPi_py(0), X_JPiPi_pz(0), X_JPiPi_massErr(0),
	  X_Jpsi1_mass(0), X_Jpsi1_VtxProb(0), X_Jpsi1_Chi2(0), X_Jpsi1_ndof(0), X_Jpsi1_px(0), X_Jpsi1_py(0), X_Jpsi1_pz(0), X_Jpsi1_massErr(0),
	  X_Jpsi2_mass(0), X_Jpsi2_VtxProb(0), X_Jpsi2_Chi2(0), X_Jpsi2_ndof(0), X_Jpsi2_px(0), X_Jpsi2_py(0), X_Jpsi2_pz(0), X_Jpsi2_massErr(0),
	  X_JPiPi_Pi1Id(0), X_JPiPi_Pi2Id(0),
	  X_JPiPi_Pi1px(0), X_JPiPi_Pi1py(0), X_JPiPi_Pi1pz(0),
	  X_JPiPi_Pi2px(0), X_JPiPi_Pi2py(0), X_JPiPi_Pi2pz(0),
	  // mass constrain variables on 1208
 	cs_X_Jpsi1_mass(0), cs_X_Jpsi1_VtxProb(0), cs_X_Jpsi1_Chi2(0), cs_X_Jpsi1_ndof(0), cs_X_Jpsi1_px(0), cs_X_Jpsi1_py(0), cs_X_Jpsi1_pz(0), cs_X_Jpsi1_massErr(0),
 	cs_X_Jpsi2_mass(0), cs_X_Jpsi2_VtxProb(0), cs_X_Jpsi2_Chi2(0), cs_X_Jpsi2_ndof(0), cs_X_Jpsi2_px(0),  cs_X_Jpsi2_py(0),  cs_X_Jpsi2_pz(0), cs_X_Jpsi2_massErr(0),
 	cs_X_JPiPi_mass(0), cs_X_JPiPi_VtxProb(0), cs_X_JPiPi_Chi2(0), cs_X_JPiPi_ndof(0), cs_X_JPiPi_px (0), cs_X_JPiPi_py(0), cs_X_JPiPi_pz(0), cs_X_JPiPi_massErr(0),
 	cs_X_mass_Psi2S(0), cs_X_VtxProb_Psi2S(0), cs_X_Chi2_Psi2S(0), cs_X_ndof_Psi2S(0), cs_X_px_Psi2S(0), cs_X_py_Psi2S(0), cs_X_pz_Psi2S(0), cs_X_massErr_Psi2S(0),
 	cs_X_JPiPi_mass_Psi2S(0), cs_X_JPiPi_VtxProb_Psi2S(0), cs_X_JPiPi_Chi2_Psi2S(0), cs_X_JPiPi_ndof_Psi2S(0), cs_X_JPiPi_px_Psi2S(0), cs_X_JPiPi_py_Psi2S(0), cs_X_JPiPi_pz_Psi2S(0), cs_X_JPiPi_massErr_Psi2S(0),
  	cs_X_mass_X3872(0), cs_X_VtxProb_X3872(0), cs_X_Chi2_X3872(0), cs_X_ndof_X3872(0), cs_X_px_X3872(0), cs_X_py_X3872(0), cs_X_pz_X3872(0), cs_X_massErr_X3872(0),
 	cs_X_JPiPi_mass_X3872(0), cs_X_JPiPi_VtxProb_X3872(0), cs_X_JPiPi_Chi2_X3872(0), cs_X_JPiPi_ndof_X3872(0), cs_X_JPiPi_px_X3872(0), cs_X_JPiPi_py_X3872(0), cs_X_JPiPi_pz_X3872(0), cs_X_JPiPi_massErr_X3872(0),
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
	  MC_Grandgranddau_mu4pdgId(0),

	  Match_mu1px(0),
	  Match_mu1py(0),
	  Match_mu1pz(0),
	  Match_mu2px(0),
	  Match_mu2py(0),
	  Match_mu2pz(0),
	  Match_mu3px(0),
	  Match_mu3py(0),
	  Match_mu3pz(0),
	  Match_mu4px(0),
	  Match_mu4py(0),
	  Match_mu4pz(0),

	  Match_pi1px(0),
	  Match_pi1py(0),
	  Match_pi1pz(0),
	  Match_pi2px(0),
	  Match_pi2py(0),
	  Match_pi2pz(0)

	  // mybxlumicorr(0), myrawbxlumi(0)
{
	// get token here for four-muon;
	gtRecordToken_ = consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag("gtDigis"));
	gtbeamspotToken_ = consumes<BeamSpot>(edm::InputTag("offlineBeamSpot"));
	gtprimaryVtxToken_ = consumes<VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices")); //  MINIAOD
	gtpatmuonToken_ = consumes<edm::View<pat::Muon>>(edm::InputTag("slimmedMuons"));				 //  MINIAOD
	gttriggerToken_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"));
	trackToken_ = consumes<edm::View<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates")); //  MINIAOD
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
	// PDG 2020
	const double myJmass = 3.0969, myJmasserr = 0.00004;
	const double mypsi2smass = 3.686097, mypsi2smasserr = 0.00003;
	const double myx3872mass = 3.872690, myx3872masserr = 0.00017;
	const double myMumass = 0.1056583745;
	const double myMumasserr = myMumass * 1e-6;
	const double myPimass = 0.13957039;
	// try
	const double myPimasserr = myPimass * 1e-6;

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

				MC_mu1p4.SetXYZM((*MC_Granddau_mu1px)[0], (*MC_Granddau_mu1py)[0], (*MC_Granddau_mu1pz)[0], myMumass);
				MC_mu2p4.SetXYZM((*MC_Granddau_mu2px)[0], (*MC_Granddau_mu2py)[0], (*MC_Granddau_mu2pz)[0], myMumass);
				MC_mu3p4.SetXYZM((*MC_Grandgranddau_mu3px)[0], (*MC_Grandgranddau_mu3py)[0], (*MC_Grandgranddau_mu3pz)[0], myMumass);
				MC_mu4p4.SetXYZM((*MC_Grandgranddau_mu4px)[0], (*MC_Grandgranddau_mu4py)[0], (*MC_Grandgranddau_mu4pz)[0], myMumass);
				MC_pi1p4.SetXYZM((*MC_Granddau_pi1px)[0], (*MC_Granddau_pi1py)[0], (*MC_Granddau_pi1pz)[0], myPimass);
				MC_pi2p4.SetXYZM((*MC_Granddau_pi2px)[0], (*MC_Granddau_pi2py)[0], (*MC_Granddau_pi2pz)[0], myPimass);
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

	edm::Handle<edm::TriggerResults> hltresults;
	bool Error_t = false;
	unsigned int nJpsitrigger = TriggersForJpsi_.size();
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
		triggerNames_ = iEvent.triggerNames(*hltresults);

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
	} // end of HLT trigger info

	edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
	iEvent.getByToken(gtRecordToken_, gtRecord);

	// L1 Trigger info
	const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
	for (unsigned int l1i = 0; l1i != ttWord.size(); ++l1i)
	{
		L1TrigRes->push_back(ttWord.at(l1i));
	}

	Vertex thePrimaryV;
	Vertex theBeamSpotV;
	BeamSpot beamSpot;

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

	nGoodPrimVtx = 0;
	for (unsigned myi = 0; myi < recVtxs->size(); myi++)
	{
		if ((*recVtxs)[myi].ndof() >= 5 && fabs((*recVtxs)[myi].z()) <= 24 && fabs((*recVtxs)[myi].position().rho()) <= 2.0)
		{
			nGoodPrimVtx++;
		}
	}

	if (recVtxs->begin() != recVtxs->end())
	{
		if (addXlessPrimaryVertex_ || resolveAmbiguity_)
		{
			thePrimaryV = Vertex(*(recVtxs->begin()));
		}
		else
		{
			for (reco::VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx)
			{
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
		thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
	}

	priVtxX = (thePrimaryV.position().x());
	priVtxY = (thePrimaryV.position().y());
	priVtxZ = (thePrimaryV.position().z());
	priVtxXErr = (thePrimaryV.xError());
	priVtxYErr = (thePrimaryV.yError());
	priVtxZErr = (thePrimaryV.zError());
	priVtxChiNorm = (thePrimaryV.normalizedChi2());
	priVtxChi = thePrimaryV.chi2();
	priVtxCL = ChiSquaredProbability((double)(thePrimaryV.chi2()), (double)(thePrimaryV.ndof()));
	priVtxXCorrY->push_back(thePrimaryV.covariance(0,1));
	priVtxXCorrZ->push_back(thePrimaryV.covariance(0,2));
	priVtxYCorrZ->push_back(thePrimaryV.covariance(1,2));

	// Check Moun and Pion Track
	edm::Handle<edm::View<pat::Muon>> thePATMuonHandle; //  MINIAOD
	iEvent.getByToken(gtpatmuonToken_, thePATMuonHandle);
	edm::Handle<edm::View<pat::PackedCandidate>> theTrackHandle; //  MINIAOD
	iEvent.getByToken(trackToken_, theTrackHandle);				 //  MINIAOD
	std::vector<edm::View<pat::PackedCandidate>::const_iterator> nonMuonPionTrack;
	// Copy tracks iterators
 	for (edm::View<pat::PackedCandidate>::const_iterator iTrackc = theTrackHandle->begin(); // MINIAOD 
	iTrackc != theTrackHandle->end(); ++iTrackc)
	{
		nonMuonPionTrack.push_back(iTrackc);
	}

	if (thePATMuonHandle->size() >= 4)
	{
		// fill muon track block
		for (edm::View<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin(); //  MINIAOD
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
			
			muD0BS->push_back(iMuonP->dB(pat::Muon::BS2D));
			muD0EBS->push_back(iMuonP->edB(pat::Muon::BS2D));
			muD3dBS->push_back(iMuonP->dB(pat::Muon::BS3D));
			muD3dEBS->push_back(iMuonP->edB(pat::Muon::BS3D));

			muD0PV->push_back(iMuonP->dB(pat::Muon::PV2D));
			muD0EPV->push_back(iMuonP->edB(pat::Muon::PV2D));
			muDzPV->push_back(iMuonP->dB(pat::Muon::PVDZ));
			muDzEPV->push_back(iMuonP->edB(pat::Muon::PVDZ));

			muCharge->push_back(iMuonP->charge());
			muTrackIso->push_back(iMuonP->trackIso());
			
			// Find and delet muon Tracks in PionTracks
			for (std::vector<edm::View<pat::PackedCandidate>::const_iterator>::const_iterator iTrackfID =  nonMuonPionTrack.begin(); // MINIAOD
			iTrackfID !=  nonMuonPionTrack.end(); ++iTrackfID)
			{
				if(iMuonP->track().isNull())
				{
					continue;
				}
				edm::View<pat::PackedCandidate>::const_iterator iTrackf = *(iTrackfID);		
				iMuonP->track()->px();
				if (iTrackf->px() == iMuonP->track()->px() && iTrackf->py() == iMuonP->track()->py() && iTrackf->pz() == iMuonP->track()->pz())
				{
					nonMuonPionTrack.erase(iTrackfID);
					iTrackfID = iTrackfID - 1;
				}
			}

			bool JpsiFilterMatch = false;
			bool JpsiTriggerMatch = false;

			for (unsigned int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++)
			{
				if(JpsiMatchTrig[JpsiTrig] != 0)
				{
					JpsiTriggerMatch = true;
					break;
				}
			}
			for(unsigned int JpsiTrig = 0; JpsiTrig < FiltersForJpsi_.size();JpsiTrig++)
			{
				if (JpsiTriggerMatch && hltresults.isValid())
				{

					pat::TriggerObjectStandAlone *tempTriggerObject = nullptr;
					for (auto i = iMuonP->triggerObjectMatches().begin(); i != iMuonP->triggerObjectMatches().end(); ++i)
					{
						tempTriggerObject = new pat::TriggerObjectStandAlone(*i);
						tempTriggerObject->unpackFilterLabels(iEvent, *hltresults);
						if(tempTriggerObject->hasFilterLabel(FiltersForJpsi_[JpsiTrig]))
						{
							JpsiFilterMatch = true;
						}
						delete tempTriggerObject;
					}
				}
			}
			muJpsiFilterRes->push_back(JpsiFilterMatch);
		}
	} // if two muons

	if (doMC)
	{
		// pion loop
		for (edm::View<pat::PackedCandidate>::const_iterator iTrack = theTrackHandle->begin(); // MINIAOD
			 iTrack != theTrackHandle->end(); ++iTrack)
		{
			TLorentzVector RECO_pip4;
			RECO_pip4.SetXYZM(iTrack->px(), iTrack->py(), iTrack->pz(), myPimass);
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

	// It takes a lot less memory out of the loop
	KinematicConstraint *Jpsi_cs = new MassKinematicConstraint(myJmass, myJmasserr);
	KinematicConstraint *JPiPi_cs_Psi2S = new MassKinematicConstraint(mypsi2smass, mypsi2smasserr);
	KinematicConstraint *JPiPi_cs_X3872 = new MassKinematicConstraint(myx3872mass, myx3872masserr);
	KinematicConstraint *Jpsi_cs34 = new MassKinematicConstraint(myJmass, myJmasserr);
	double pionptcut = 0.25;
	double pionDRcut = 0.7;
	double vtxprobprecut = 1.0e-7;

	// get the four moun fit, but first also fit dimuon
	if (thePATMuonHandle->size() < 4)
	{
		return;
	}
	//  get X and MyFourMuon cands
	for (edm::View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); // MINIAOD
		 iMuon1 != thePATMuonHandle->end(); ++iMuon1)
	{
		TrackRef muTrack1 = iMuon1->track();
		if (muTrack1.isNull())
		{
			continue;
		}

		reco::Track recoMu1 = *iMuon1->track();

		// next check for mu2
		for (edm::View<pat::Muon>::const_iterator iMuon2 = iMuon1 + 1; // MINIAOD
			 iMuon2 != thePATMuonHandle->end(); ++iMuon2)
		{
			TrackRef muTrack2 = iMuon2->track();
			if (muTrack2.isNull())
			{
				continue;
			}

			reco::Track recoMu2 = *iMuon2->track();

			if (!(1. < (iMuon1->p4() + iMuon2->p4()).mass() && (iMuon1->p4() + iMuon2->p4()).mass() < 4.))
			{
				continue;
			}
			if ((iMuon1->charge() + iMuon2->charge()) != 0)
			{
				continue;
			}

			TransientTrack muon1TT(muTrack1, &(bFieldHandle)); // MINIAOD
			TransientTrack muon2TT(muTrack2, &(bFieldHandle)); // MINIAOD
			KinematicParticleFactoryFromTransientTrack pmumuFactory;
			// The mass of a muon and the insignificant mass sigma
			// to avoid singularities in the covariance matrix.
			ParticleMass muon_mass = myMumass; // pdg mass
			float muon_sigma = myMumasserr;
			// initial chi2 and ndf before kinematic fits.
			float chi = 0.;
			float ndf = 0.;
			vector<RefCountedKinematicParticle> muonParticles;
			muonParticles.push_back(pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
			muonParticles.push_back(pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));
			KinematicParticleVertexFitter mu12_fitter;
			RefCountedKinematicTree Jpsi1VertexFitTree;
			Error_t = false;
			try
			{
				Jpsi1VertexFitTree = mu12_fitter.fit(muonParticles);
			}catch(...)
			{
				Error_t = true;
				std::cout<<"error at Jpsi1noMC"<<std::endl;
			}
			if (Error_t || !Jpsi1VertexFitTree->isValid())
			{
				continue;
			}
			Jpsi1VertexFitTree->movePointerToTheTop();
			RefCountedKinematicParticle Jpsi1_vFit_noMC = Jpsi1VertexFitTree->currentParticle();
			RefCountedKinematicVertex Jpsi1_vFit_vertex_noMC = Jpsi1VertexFitTree->currentDecayVertex();
			KinematicParameters mymumupara = Jpsi1_vFit_noMC->currentState().kinematicParameters();
			double Jpsi1_vtxprob = ChiSquaredProbability((double)(Jpsi1_vFit_vertex_noMC->chiSquared()), (double)(Jpsi1_vFit_vertex_noMC->degreesOfFreedom()));
			if (Jpsi1_vFit_noMC->currentState().mass() > 3.4 || Jpsi1_vFit_noMC->currentState().mass() < 2.8)
			{
				continue;
			}
			if (Jpsi1_vtxprob < vtxprobprecut)
			{
				continue;
			}

			// mass constrain for Jpsi1 from psi2s:
			bool flag_jpsi1 = true;
			bool flag_jpsi2 = true;
			KinematicParticleFitter mu12_fitter_cs;
			RefCountedKinematicParticle Jpsi1_vFit_cs;
			RefCountedKinematicVertex Jpsi1_vFit_vertex_cs;
			KinematicParameters mymumupara_cs;
			double Jpsi1_vtxprob_cs = -9;
			Error_t = false;
			try
			{
				Jpsi1VertexFitTree = mu12_fitter_cs.fit(Jpsi_cs, Jpsi1VertexFitTree);
			}catch(...)
			{
				Error_t = true;
				std::cout<<"error at Jpsi1Mc"<<std::endl;
			}
			if(Error_t || !Jpsi1VertexFitTree->isValid())
			{
				flag_jpsi1 = false;
			}else
			{
				Jpsi1VertexFitTree->movePointerToTheTop();
				Jpsi1_vFit_cs = Jpsi1VertexFitTree->currentParticle();
				Jpsi1_vFit_vertex_cs = Jpsi1VertexFitTree->currentDecayVertex();
				mymumupara_cs = Jpsi1_vFit_cs->currentState().kinematicParameters();
				Jpsi1_vtxprob_cs = ChiSquaredProbability((double)(Jpsi1_vFit_vertex_cs->chiSquared()), (double)(Jpsi1_vFit_vertex_cs->degreesOfFreedom()));
			}
			TLorentzVector P4_mu1;
			P4_mu1.SetPtEtaPhiM(iMuon1->track()->pt(), iMuon1->track()->eta(), iMuon1->track()->phi(), myMumass);
			TLorentzVector P4_mu2;
			P4_mu2.SetPtEtaPhiM(iMuon2->track()->pt(), iMuon2->track()->eta(), iMuon2->track()->phi(), myMumass);
			// mu3mu4(X6900->Jpsi)
			for (edm::View<pat::Muon>::const_iterator iMuon3 = thePATMuonHandle->begin(); // MINIAOD
				 iMuon3 != thePATMuonHandle->end(); ++iMuon3)
			{
				if(iMuon3 == iMuon1 || iMuon3 == iMuon2)
				{
					continue;
				}
				TrackRef muTrack3 = iMuon3->track();
				if (muTrack3.isNull())
				{
					continue;
				}

				reco::Track recoMu3 = *iMuon3->track();
				for (edm::View<pat::Muon>::const_iterator iMuon4 = iMuon3 + 1; // MINIAOD
					 iMuon4 != thePATMuonHandle->end(); ++iMuon4)
				{
					if(iMuon4 == iMuon1 || iMuon4 == iMuon2)
                               		{
                                        	continue;
                                	}
					TrackRef muTrack4 = iMuon4->track();
					if (muTrack4.isNull())
					{
						continue;
					}

					reco::Track recoMu4 = *iMuon4->track(); // MINIAOD

					if (!(1. < (iMuon3->p4() + iMuon4->p4()).mass() && (iMuon3->p4() + iMuon4->p4()).mass() < 4.5))
					{
						continue;
					}
					if ((iMuon3->charge() + iMuon4->charge()) != 0)
					{
						continue;
					}

					TransientTrack muon3TT(muTrack3, &(bFieldHandle)); // MINIAOD
					TransientTrack muon4TT(muTrack4, &(bFieldHandle)); // MINIAOD
					KinematicParticleFactoryFromTransientTrack pmumuFactory34;
					ParticleMass muon_mass = myMumass; // pdg mass
					float muon_sigma = myMumasserr;
					float chi = 0.;
					float ndf = 0.;
					vector<RefCountedKinematicParticle> muonParticles34;
					muonParticles34.push_back(pmumuFactory34.particle(muon3TT, muon_mass, chi, ndf, muon_sigma));
					muonParticles34.push_back(pmumuFactory34.particle(muon4TT, muon_mass, chi, ndf, muon_sigma));
					KinematicParticleVertexFitter mu34_fitter;
					RefCountedKinematicTree Jpsi2VertexFitTree;
					Error_t = false;
					try
					{
						Jpsi2VertexFitTree = mu34_fitter.fit(muonParticles34);
					}catch(...)
					{
						Error_t = true;
						std::cout<<"error at Jpsi2noMC"<<std::endl;
					}
					if (Error_t || !Jpsi2VertexFitTree->isValid())
					{
							continue;
					}
					Jpsi2VertexFitTree->movePointerToTheTop();
					RefCountedKinematicParticle Jpsi2_vFit_noMC = Jpsi2VertexFitTree->currentParticle();
					RefCountedKinematicVertex Jpsi2_vFit_vertex_noMC = Jpsi2VertexFitTree->currentDecayVertex();
					KinematicParameters mymumupara2 = Jpsi2_vFit_noMC->currentState().kinematicParameters();
					double Jpsi2_vtxprob = ChiSquaredProbability((double)(Jpsi2_vFit_vertex_noMC->chiSquared()), (double)(Jpsi2_vFit_vertex_noMC->degreesOfFreedom()));
					if (Jpsi2_vFit_noMC->currentState().mass() > 3.4 || Jpsi2_vFit_noMC->currentState().mass() < 2.8)
					{
						continue;
					} // change 4.0 to 3.4 in 1208
					if (Jpsi2_vtxprob < vtxprobprecut)
					{
						continue;
					}

					// mass constrain for Jpsi from X6900:
					KinematicParticleFitter mu34_fitter_cs;	
					RefCountedKinematicParticle Jpsi2_vFit_cs;
					RefCountedKinematicVertex Jpsi2_vFit_vertex_cs;
					KinematicParameters mymumupara2_cs;
					double Jpsi2_vtxprob_cs = -9;
					Error_t = false;
					try
					{
						Jpsi2VertexFitTree = mu34_fitter_cs.fit(Jpsi_cs34, Jpsi2VertexFitTree);
					}catch(...)
					{
						Error_t = true;
						std::cout<<"error at Jpsi2MC"<<std::endl;
					}
					if(Error_t || !Jpsi2VertexFitTree->isValid())
					{
						flag_jpsi2 = false;
					}else
					{
						Jpsi2VertexFitTree->movePointerToTheTop();
						Jpsi2_vFit_cs = Jpsi2VertexFitTree->currentParticle();
						Jpsi2_vFit_vertex_cs = Jpsi2VertexFitTree->currentDecayVertex();
						mymumupara2_cs = Jpsi2_vFit_cs->currentState().kinematicParameters();
						Jpsi2_vtxprob_cs = ChiSquaredProbability((double)(Jpsi2_vFit_vertex_cs->chiSquared()), (double)(Jpsi2_vFit_vertex_cs->degreesOfFreedom()));
					}
					// psi2s
					for (std::vector<edm::View<pat::PackedCandidate>::const_iterator>::const_iterator iTrack1ID = nonMuonPionTrack.begin(); // MINIAOD
						 iTrack1ID != nonMuonPionTrack.end(); ++iTrack1ID)
					{
						edm::View<pat::PackedCandidate>::const_iterator iTrack1 = *(iTrack1ID);
						if (iTrack1->pt() < pionptcut)
						{
							continue;
						}
						for (std::vector<edm::View<pat::PackedCandidate>::const_iterator>::const_iterator iTrack2ID = iTrack1ID + 1; // MINIAOD
							 iTrack2ID != nonMuonPionTrack.end(); ++iTrack2ID)
						{
							edm::View<pat::PackedCandidate>::const_iterator iTrack2 = *(iTrack2ID);
							if (iTrack2->pt() < pionptcut)
							{
								continue;
							}
							if ((iTrack1->charge() + iTrack2->charge()) != 0)
								continue;

							// MINIAOD begin
							if (!iTrack1->hasTrackDetails() || iTrack1->charge() == 0)
							{
								// cout << "iTrack1->hasTrackDetails() = " << iTrack1->hasTrackDetails() << endl;
								// cout << "iTrack1->charge() = " << iTrack1->charge() << endl;
								continue;
							}
							if (!iTrack2->hasTrackDetails() || iTrack2->charge() == 0)
							{
								// cout << "iTrack2->hasTrackDetails() = " << iTrack2->hasTrackDetails() << endl;
								// cout << "iTrack2->charge() = " << iTrack2->charge() << endl;
								continue;
							}
							// MINIAOD end

							TLorentzVector P4_Track1, P4_Track2, P4_Jpsipipi;
							P4_Track1.SetPtEtaPhiM(iTrack1->pt(), iTrack1->eta(), iTrack1->phi(), myPimass);
							P4_Track2.SetPtEtaPhiM(iTrack2->pt(), iTrack2->eta(), iTrack2->phi(), myPimass);
							P4_Jpsipipi = P4_mu1 + P4_mu2 + P4_Track1 + P4_Track2;

							if (P4_Track1.DeltaR(P4_Jpsipipi) > pionDRcut)
							{
								continue;
							}
							if (P4_Track2.DeltaR(P4_Jpsipipi) > pionDRcut)
							{
								continue;
							}

							// TransientTrack trackTT1(*iTrack1, &(bFieldHandle));
							// TransientTrack trackTT2(*iTrack2, &(bFieldHandle));
							TransientTrack trackTT1(*(iTrack1->bestTrack()), &(bFieldHandle)); // MINIAOD
							TransientTrack trackTT2(*(iTrack2->bestTrack()), &(bFieldHandle)); // MINIAOD
							KinematicParticleFactoryFromTransientTrack JPiPiFactory;
							// The mass of a muon and the insignificant mass sigma
							// to avoid singularities in the covariance matrix.
							ParticleMass pion_mass = myPimass; // pdg mass
							float pion_sigma = myPimasserr;
							// initial chi2 and ndf before kinematic fits.
							float chi = 0.;
							float ndf = 0.;

							// mass constrain for psi2s from X6900, now noMC
							// first fit to have a psi2s
							vector<RefCountedKinematicParticle> JPiPiParticles;
							JPiPiParticles.push_back(JPiPiFactory.particle(trackTT1, pion_mass, chi, ndf, pion_sigma));
							JPiPiParticles.push_back(JPiPiFactory.particle(trackTT2, pion_mass, chi, ndf, pion_sigma));
							JPiPiParticles.push_back(pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
							JPiPiParticles.push_back(pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));

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
							JPiPiVertexFitTree->movePointerToTheTop();
							RefCountedKinematicParticle JPiPi_vFit_noMC = JPiPiVertexFitTree->currentParticle();
							RefCountedKinematicVertex JPiPi_vFit_vertex_noMC = JPiPiVertexFitTree->currentDecayVertex();

							double JPiPi_vtxprob = ChiSquaredProbability((double)(JPiPi_vFit_vertex_noMC->chiSquared()), (double)(JPiPi_vFit_vertex_noMC->degreesOfFreedom()));
							if (JPiPi_vFit_noMC->currentState().mass() > 4.5)
							{
								continue;
							}
							if (JPiPi_vtxprob < vtxprobprecut)
							{
								continue;
							}

							// fit the 6 track together to a vertex
							vector<RefCountedKinematicParticle> X_Particles;
							X_Particles.push_back(JPiPiFactory.particle(trackTT1, pion_mass, chi, ndf, pion_sigma));
							X_Particles.push_back(JPiPiFactory.particle(trackTT2, pion_mass, chi, ndf, pion_sigma));
							X_Particles.push_back(pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
							X_Particles.push_back(pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));
							X_Particles.push_back(pmumuFactory.particle(muon3TT, muon_mass, chi, ndf, muon_sigma));
							X_Particles.push_back(pmumuFactory.particle(muon4TT, muon_mass, chi, ndf, muon_sigma));

							KinematicParticleVertexFitter X_fitter;
							RefCountedKinematicTree X_VertexFitTree;
							Error_t = false;
							try
							{
								X_VertexFitTree = X_fitter.fit(X_Particles);
							}catch(...)
							{
								Error_t = true;	
								std::cout<<"error at XnoMC"<<std::endl;
}
							if (Error_t || !(X_VertexFitTree->isValid()))
							{
								continue;
							}
							X_VertexFitTree->movePointerToTheTop();
							RefCountedKinematicParticle X_vFit_noMC = X_VertexFitTree->currentParticle();
							RefCountedKinematicVertex X_vFit_vertex_noMC = X_VertexFitTree->currentDecayVertex();
							KinematicParameters X_kPara = X_vFit_noMC->currentState().kinematicParameters();

							double X_vtxprob = ChiSquaredProbability((double)(X_vFit_vertex_noMC->chiSquared()), (double)(X_vFit_vertex_noMC->degreesOfFreedom()));

							X_VertexFitTree->movePointerToTheFirstChild();
							RefCountedKinematicParticle X_pi1 = X_VertexFitTree->currentParticle();
							X_VertexFitTree->movePointerToTheNextChild();
							RefCountedKinematicParticle X_pi2 = X_VertexFitTree->currentParticle();

							KinematicParameters X_pi1_KP = X_pi1->currentState().kinematicParameters();
							KinematicParameters X_pi2_KP = X_pi2->currentState().kinematicParameters();

							X_mass->push_back(X_vFit_noMC->currentState().mass());
							X_VtxProb->push_back(X_vtxprob);
							X_Chi2->push_back((double)(X_vFit_vertex_noMC->chiSquared()));
							X_ndof->push_back((double)(X_vFit_vertex_noMC->degreesOfFreedom()));

							X_px->push_back(X_kPara.momentum().x());
							X_py->push_back(X_kPara.momentum().y());
							X_pz->push_back(X_kPara.momentum().z());
 							if (X_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
							{			
								X_massErr->push_back(sqrt(X_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
							}
                                                        else
                                                        {
								X_massErr->push_back(-9); 
							}
	
							X_mu1Id->push_back(std::distance(thePATMuonHandle->begin(), iMuon1));
							X_mu2Id->push_back(std::distance(thePATMuonHandle->begin(), iMuon2));
							X_mu3Id->push_back(std::distance(thePATMuonHandle->begin(), iMuon3));
							X_mu4Id->push_back(std::distance(thePATMuonHandle->begin(), iMuon4));

							X_Jpsi1_mass->push_back(Jpsi1_vFit_noMC->currentState().mass());
							X_Jpsi1_VtxProb->push_back(Jpsi1_vtxprob);
							X_Jpsi1_Chi2->push_back(double(Jpsi1_vFit_vertex_noMC->chiSquared()));
							X_Jpsi1_ndof->push_back(double(Jpsi1_vFit_vertex_noMC->degreesOfFreedom()));
							X_Jpsi1_px->push_back(mymumupara.momentum().x());
							X_Jpsi1_py->push_back(mymumupara.momentum().y());
							X_Jpsi1_pz->push_back(mymumupara.momentum().z());
							if (Jpsi1_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
							{
								X_Jpsi1_massErr->push_back(sqrt(Jpsi1_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
							}
							else
							{
								X_Jpsi1_massErr->push_back(-9);
							}
							X_Jpsi2_mass->push_back(Jpsi2_vFit_noMC->currentState().mass());
							X_Jpsi2_VtxProb->push_back(Jpsi2_vtxprob);
							X_Jpsi2_Chi2->push_back(double(Jpsi2_vFit_vertex_noMC->chiSquared()));
							X_Jpsi2_ndof->push_back(double(Jpsi2_vFit_vertex_noMC->degreesOfFreedom()));
							X_Jpsi2_px->push_back(mymumupara2.momentum().x());
							X_Jpsi2_py->push_back(mymumupara2.momentum().y());
							X_Jpsi2_pz->push_back(mymumupara2.momentum().z());
							if (Jpsi2_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
							{
								X_Jpsi2_massErr->push_back(sqrt(Jpsi2_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
							}
							else
							{
								X_Jpsi2_massErr->push_back(-9);
							}
							X_JPiPi_mass->push_back(JPiPi_vFit_noMC->currentState().mass());
							X_JPiPi_VtxProb->push_back(JPiPi_vtxprob);
							X_JPiPi_Chi2->push_back(double(JPiPi_vFit_vertex_noMC->chiSquared()));
							X_JPiPi_ndof->push_back(double(JPiPi_vFit_vertex_noMC->degreesOfFreedom()));
							X_JPiPi_px->push_back(JPiPi_vFit_noMC->currentState().kinematicParameters().momentum().x());
							X_JPiPi_py->push_back(JPiPi_vFit_noMC->currentState().kinematicParameters().momentum().y());
							X_JPiPi_pz->push_back(JPiPi_vFit_noMC->currentState().kinematicParameters().momentum().z());
							if (JPiPi_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
							{
								X_JPiPi_massErr->push_back(sqrt(JPiPi_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
							}
							else
							{
								X_JPiPi_massErr->push_back(-9);
							}
							X_JPiPi_Pi1Id->push_back(std::distance(theTrackHandle->begin(), iTrack1));
							X_JPiPi_Pi2Id->push_back(std::distance(theTrackHandle->begin(), iTrack2));
							

							X_JPiPi_Pi1px->push_back(X_pi1_KP.momentum().x()); // iTrack1->px()
							X_JPiPi_Pi1py->push_back(X_pi1_KP.momentum().y());
							X_JPiPi_Pi1pz->push_back(X_pi1_KP.momentum().z());

							X_JPiPi_Pi2px->push_back(X_pi2_KP.momentum().x());
							X_JPiPi_Pi2py->push_back(X_pi2_KP.momentum().y());
							X_JPiPi_Pi2pz->push_back(X_pi2_KP.momentum().z());
							// JPiPi Fit with Jpsi1 MassConstraint
							bool flag_vtx = true;
							bool flag_cs_psi2s = true;
							bool flag_cs_x3872 = true;
// flags for fit status, if true fit succeeded, if false fit failed
/* Choice Logic is as below(pseudo code):
 * IF JPiPi fit with Jpsi1 MC == successful
 * THEN
 * 	flag_vtx = true
 * 	DO 
 * 		JPiPi and X fit MC to Psi2S 
 * 		IF successful
 * 		THEN save information & flag_cs_psi2s = true
 * 		ELSE flag_cs_psi2s = false 
 * 	DO 	
 * 		JPiPi and X fit MC to X3872
 * 		IF successful
 * 		THEN save information & flag_cs_x3872 = true
 * 		ELSE flag_cs_x3872 = false
 * ELSE
 * 	flag_vtx = flag_cs_psi2s = flag_cs_x3872 = false
 * IF flag_vtx == false THEN DO save -9 for JPiPi vtx fit and Jpsi1 Jpsi2
 * IF flag_cs_psi2s == true THEN DO save -9 for JPiPi to Psi2S
 * IF flag_cs_x3872 == true THEN DO save -9 for JPiPi to X3872
 */
						
						RefCountedKinematicTree JPiPiVertexFitTree_cs_vtx; 
						if(flag_jpsi1 == false)
						{
							flag_vtx = false;
						}else
						{
							vector<RefCountedKinematicParticle> JPiPiParticles_cs_vtx;
							JPiPiParticles_cs_vtx.push_back(JPiPiFactory.particle(trackTT1, pion_mass, chi, ndf, pion_sigma));
							JPiPiParticles_cs_vtx.push_back(JPiPiFactory.particle(trackTT2, pion_mass, chi, ndf, pion_sigma));
							JPiPiParticles_cs_vtx.push_back(Jpsi1_vFit_cs);

							KinematicParticleVertexFitter JPiPi_fitter_cs_vtx; 
// Variables named with JPiPi_cs_vtx is for JPiPi fitted with 2pi and MassConstraint Jpsi1. And JPiPi_cs_Psi2S/X3872 (no _vtx and MassC target) defined below is for Final MassConstraint fit on JPiPi candidate
							Error_t = false;
							try
							{
								JPiPiVertexFitTree_cs_vtx = JPiPi_fitter_cs_vtx.fit(JPiPiParticles_cs_vtx);
							}catch(...)
							{
								Error_t = true;
								std::cout<<"error at JPiPi_JpsiMC"<<std::endl;
							}	
							if (Error_t || !(JPiPiVertexFitTree_cs_vtx->isValid()))
							{
								flag_vtx   = false;
								flag_cs_psi2s = false;
								flag_cs_x3872 = false; 
							}
							RefCountedKinematicParticle JPiPi_vFit_cs_vtx;
							RefCountedKinematicVertex JPiPi_vFit_vertex_cs_vtx;
							double JPiPi_vtxprob_cs_vtx;
							if (flag_vtx == true)
							{
								JPiPiVertexFitTree_cs_vtx->movePointerToTheTop();
								JPiPi_vFit_cs_vtx = JPiPiVertexFitTree_cs_vtx->currentParticle();
								JPiPi_vFit_vertex_cs_vtx = JPiPiVertexFitTree_cs_vtx->currentDecayVertex();
								JPiPi_vtxprob_cs_vtx = ChiSquaredProbability((double)(JPiPi_vFit_vertex_cs_vtx->chiSquared()), (double)(JPiPi_vFit_vertex_cs_vtx->degreesOfFreedom()));
								if (JPiPi_vFit_cs_vtx->currentState().mass() > 4.5 || JPiPi_vtxprob_cs_vtx < vtxprobprecut)
								{
									flag_vtx   = false;
									flag_cs_psi2s = false;
									flag_cs_x3872 = false;
								}
							}
							if ( flag_vtx == true)
							{
							// Save JPiPi MassConstraint
								cs_X_JPiPi_mass->push_back(JPiPi_vFit_cs_vtx->currentState().mass());
                                                                cs_X_JPiPi_VtxProb->push_back(JPiPi_vtxprob_cs_vtx);
                                                                cs_X_JPiPi_Chi2->push_back((double)(JPiPi_vFit_vertex_cs_vtx->chiSquared()));
                                                                cs_X_JPiPi_ndof->push_back((double)(JPiPi_vFit_vertex_cs_vtx->degreesOfFreedom()));
                                                                cs_X_JPiPi_px->push_back(JPiPi_vFit_cs_vtx->currentState().kinematicParameters().momentum().x());
                                                                cs_X_JPiPi_py->push_back(JPiPi_vFit_cs_vtx->currentState().kinematicParameters().momentum().y());
                                                                cs_X_JPiPi_pz->push_back(JPiPi_vFit_cs_vtx->currentState().kinematicParameters().momentum().z());
								if (JPiPi_vFit_cs_vtx->currentState().kinematicParametersError().matrix()(6, 6) > 0)
                                                                {
                                                                        cs_X_JPiPi_massErr->push_back(sqrt(JPiPi_vFit_cs_vtx->currentState().kinematicParametersError().matrix()(6, 6)));
                                                                }
                                                                else
                                                                {
                                                                        cs_X_JPiPi_massErr->push_back(-9);
                                                                }
								// Save Jpsi and Pion MassConstraint
								RefCountedKinematicParticle JPiPi_pi1_cs, JPiPi_pi2_cs;
                                                                KinematicParameters JPiPi_pi1_KP_cs, JPiPi_pi2_KP_cs;
                                                                JPiPiVertexFitTree_cs_vtx->movePointerToTheFirstChild();
                                                                JPiPi_pi1_cs = JPiPiVertexFitTree_cs_vtx->currentParticle();
                                                                JPiPiVertexFitTree_cs_vtx->movePointerToTheNextChild();
                                                                JPiPi_pi2_cs = JPiPiVertexFitTree_cs_vtx->currentParticle();
                                                                JPiPi_pi1_KP_cs = JPiPi_pi1_cs->currentState().kinematicParameters();
								JPiPi_pi2_KP_cs = JPiPi_pi2_cs->currentState().kinematicParameters();
                                                                cs_X_Jpsi1_mass->push_back(Jpsi1_vFit_cs->currentState().mass());
                                                                cs_X_Jpsi1_VtxProb->push_back(Jpsi1_vtxprob_cs);
								cs_X_Jpsi1_Chi2->push_back(double(Jpsi1_vFit_vertex_cs->chiSquared()));
								cs_X_Jpsi1_ndof->push_back(double(Jpsi1_vFit_vertex_cs->degreesOfFreedom()));
								cs_X_Jpsi1_px->push_back(mymumupara_cs.momentum().x());
								cs_X_Jpsi1_py->push_back(mymumupara_cs.momentum().y());
								cs_X_Jpsi1_pz->push_back(mymumupara_cs.momentum().z());
								if (Jpsi1_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6) > 0)	
								{
									cs_X_Jpsi1_massErr->push_back(sqrt(Jpsi1_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6)));
								}
								else
								{
									cs_X_Jpsi1_massErr->push_back(-9);
								}
								if(flag_jpsi2 == true)
								{
								cs_X_Jpsi2_mass->push_back(Jpsi2_vFit_cs->currentState().mass());
								cs_X_Jpsi2_VtxProb->push_back(Jpsi2_vtxprob_cs);
								cs_X_Jpsi2_Chi2->push_back(double(Jpsi2_vFit_vertex_cs->chiSquared()));
                                                                cs_X_Jpsi2_ndof->push_back(double(Jpsi2_vFit_vertex_cs->degreesOfFreedom())); 
								cs_X_Jpsi2_px->push_back(mymumupara2_cs.momentum().x());
								cs_X_Jpsi2_py->push_back(mymumupara2_cs.momentum().y());
								cs_X_Jpsi2_pz->push_back(mymumupara2_cs.momentum().z());
								if (Jpsi2_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6) > 0)
								{
									cs_X_Jpsi2_massErr->push_back(sqrt(Jpsi2_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6)));
								}
								else
								{
									cs_X_Jpsi2_massErr->push_back(-9);
								}
								}else
								{
								cs_X_Jpsi2_mass->push_back(-9);
								cs_X_Jpsi2_VtxProb->push_back(-9);
								cs_X_Jpsi2_Chi2->push_back(-9);
                                                                cs_X_Jpsi2_ndof->push_back(-9); 
								cs_X_Jpsi2_px->push_back(-99999);
								cs_X_Jpsi2_py->push_back(-99999);
								cs_X_Jpsi2_pz->push_back(-99999);
								cs_X_Jpsi2_massErr->push_back(-9);
								}
							}
// Notice that Jpsi1 and Pion information is actually the same as before, because JPiPi MC fit doesn't change them by default. If you do what to change them, please check JPiPi MC fit and information stored carefully.
						}
							// JPiPi MassConstraint Fit to Psi2S
							if (flag_vtx == true && flag_jpsi2 == true)
							{ 	
								KinematicParticleFitter JPiPi_fitter_cs_Psi2S;
								RefCountedKinematicTree JPiPiVertexFitTree_cs_Psi2S;
								RefCountedKinematicParticle JPiPi_vFit_cs_Psi2S;
								RefCountedKinematicVertex JPiPi_vFit_vertex_cs_Psi2S;
								double JPiPi_vtxprob_cs_Psi2S;
								Error_t = false;
								try
								{
									JPiPiVertexFitTree_cs_Psi2S = JPiPi_fitter_cs_Psi2S.fit(JPiPi_cs_Psi2S, JPiPiVertexFitTree_cs_vtx);
								}catch(...)
								{
									Error_t = true;
									std::cout<<"error at JPiPitoPsi2S"<<std::endl;
								}
								if ( Error_t || !JPiPiVertexFitTree_cs_Psi2S->isValid())	
								{
									flag_cs_psi2s = false;
								}else
								{
									JPiPiVertexFitTree_cs_Psi2S->movePointerToTheTop();
									JPiPi_vFit_cs_Psi2S = JPiPiVertexFitTree_cs_Psi2S->currentParticle();
									JPiPi_vFit_vertex_cs_Psi2S = JPiPiVertexFitTree_cs_Psi2S->currentDecayVertex();
									JPiPi_vtxprob_cs_Psi2S = ChiSquaredProbability((double)(JPiPi_vFit_vertex_cs_Psi2S->chiSquared()), (double)(JPiPi_vFit_vertex_cs_Psi2S->degreesOfFreedom()));
									// X MassConstraint Fit to Psi2S 
									KinematicParticleVertexFitter X_fitter_cs_Psi2S;
									RefCountedKinematicTree X_VertexFitTree_cs_Psi2S;
									vector<RefCountedKinematicParticle> X_Particles_cs_Psi2S;
									X_Particles_cs_Psi2S.push_back(Jpsi2_vFit_cs);
									X_Particles_cs_Psi2S.push_back(JPiPi_vFit_cs_Psi2S);
									Error_t = false;
									try
									{
										X_VertexFitTree_cs_Psi2S = X_fitter_cs_Psi2S.fit(X_Particles_cs_Psi2S);
									}catch(...)
									{
										Error_t = true;
										std::cout<<"error at XtoPsi2S"<<std::endl;
									}
									if (Error_t || !(X_VertexFitTree_cs_Psi2S->isValid()))
									{
										flag_cs_psi2s = false;
									}else
									{
										X_VertexFitTree_cs_Psi2S->movePointerToTheTop();
										RefCountedKinematicParticle X_vFit_cs_Psi2S = X_VertexFitTree_cs_Psi2S->currentParticle();
										RefCountedKinematicVertex X_vFit_vertex_cs_Psi2S = X_VertexFitTree_cs_Psi2S->currentDecayVertex();
										double X_vtxprob_cs_Psi2S = ChiSquaredProbability((double)(X_vFit_vertex_cs_Psi2S->chiSquared()), (double)(X_vFit_vertex_cs_Psi2S->degreesOfFreedom()));
										cs_X_JPiPi_mass_Psi2S->push_back(JPiPi_vFit_cs_Psi2S->currentState().mass());
										cs_X_JPiPi_VtxProb_Psi2S->push_back(JPiPi_vtxprob_cs_Psi2S); 
										cs_X_JPiPi_Chi2_Psi2S->push_back((double)(JPiPi_vFit_vertex_cs_Psi2S->chiSquared()));
										cs_X_JPiPi_ndof_Psi2S->push_back((double)(JPiPi_vFit_vertex_cs_Psi2S->degreesOfFreedom()));
										cs_X_JPiPi_px_Psi2S->push_back(JPiPi_vFit_cs_Psi2S->currentState().kinematicParameters().momentum().x());
										cs_X_JPiPi_py_Psi2S->push_back(JPiPi_vFit_cs_Psi2S->currentState().kinematicParameters().momentum().y());
										cs_X_JPiPi_pz_Psi2S->push_back(JPiPi_vFit_cs_Psi2S->currentState().kinematicParameters().momentum().z());
										if (JPiPi_vFit_cs_Psi2S->currentState().kinematicParametersError().matrix()(6, 6) > 0)
										{
										cs_X_JPiPi_massErr_Psi2S->push_back(sqrt(JPiPi_vFit_cs_Psi2S->currentState().kinematicParametersError().matrix()(6, 6)));
										}		
										else
										{
										cs_X_JPiPi_massErr_Psi2S->push_back(-9);
										}

										cs_X_mass_Psi2S->push_back(X_vFit_cs_Psi2S->currentState().mass());
										cs_X_VtxProb_Psi2S->push_back(X_vtxprob_cs_Psi2S);
										cs_X_Chi2_Psi2S->push_back((double)(X_vFit_vertex_cs_Psi2S->chiSquared()));
										cs_X_ndof_Psi2S->push_back((double)(X_vFit_vertex_cs_Psi2S->degreesOfFreedom()));
										cs_X_px_Psi2S->push_back(X_vFit_cs_Psi2S->currentState().kinematicParameters().momentum().x());
										cs_X_py_Psi2S->push_back(X_vFit_cs_Psi2S->currentState().kinematicParameters().momentum().y());
										cs_X_pz_Psi2S->push_back(X_vFit_cs_Psi2S->currentState().kinematicParameters().momentum().z());

										if (X_vFit_cs_Psi2S->currentState().kinematicParametersError().matrix()(6, 6) > 0)
										{
										cs_X_massErr_Psi2S->push_back(sqrt(X_vFit_cs_Psi2S->currentState().kinematicParametersError().matrix()(6, 6)));
										}		
										else
										{
										cs_X_massErr_Psi2S->push_back(-9);
										}
									}
								}
							}


							// JPiPi MassConstraint Fit to X3872
							if (flag_vtx == true && flag_jpsi2 == true)
							{ 	
								KinematicParticleFitter JPiPi_fitter_cs_X3872;
								RefCountedKinematicTree JPiPiVertexFitTree_cs_X3872;
								RefCountedKinematicParticle JPiPi_vFit_cs_X3872;
								RefCountedKinematicVertex JPiPi_vFit_vertex_cs_X3872;
								double JPiPi_vtxprob_cs_X3872;
								Error_t = false;
								try
								{
									JPiPiVertexFitTree_cs_X3872 = JPiPi_fitter_cs_X3872.fit(JPiPi_cs_X3872, JPiPiVertexFitTree_cs_vtx);
								}catch(...)
								{
									Error_t = true;
									std::cout<<"error at JPiPitoX3872"<<std::endl;
								}
								if ( Error_t || !JPiPiVertexFitTree_cs_X3872->isValid())	
								{
									flag_cs_x3872 = false;
								}else
								{
									JPiPiVertexFitTree_cs_X3872->movePointerToTheTop();
									JPiPi_vFit_cs_X3872 = JPiPiVertexFitTree_cs_X3872->currentParticle();
									JPiPi_vFit_vertex_cs_X3872 = JPiPiVertexFitTree_cs_X3872->currentDecayVertex();
									JPiPi_vtxprob_cs_X3872 = ChiSquaredProbability((double)(JPiPi_vFit_vertex_cs_X3872->chiSquared()), (double)(JPiPi_vFit_vertex_cs_X3872->degreesOfFreedom()));
									// X MassConstraint Fit to X3872 
									KinematicParticleVertexFitter X_fitter_cs_X3872;
									RefCountedKinematicTree X_VertexFitTree_cs_X3872;
									vector<RefCountedKinematicParticle> X_Particles_cs_X3872;
									X_Particles_cs_X3872.push_back(Jpsi2_vFit_cs);
									X_Particles_cs_X3872.push_back(JPiPi_vFit_cs_X3872);
									Error_t = false;
									try
									{
										X_VertexFitTree_cs_X3872 = X_fitter_cs_X3872.fit(X_Particles_cs_X3872);
									}catch(...)
									{
										Error_t = true;
										std::cout<<"error at XtoX3872"<<std::endl;
									}
									if (Error_t || !(X_VertexFitTree_cs_X3872->isValid()))
									{
										flag_cs_x3872 = false;
									}else
									{
										X_VertexFitTree_cs_X3872->movePointerToTheTop();
										RefCountedKinematicParticle X_vFit_cs_X3872 = X_VertexFitTree_cs_X3872->currentParticle();
										RefCountedKinematicVertex X_vFit_vertex_cs_X3872 = X_VertexFitTree_cs_X3872->currentDecayVertex();
										double X_vtxprob_cs_X3872 = ChiSquaredProbability((double)(X_vFit_vertex_cs_X3872->chiSquared()), (double)(X_vFit_vertex_cs_X3872->degreesOfFreedom()));
										cs_X_JPiPi_mass_X3872->push_back(JPiPi_vFit_cs_X3872->currentState().mass());
										cs_X_JPiPi_VtxProb_X3872->push_back(JPiPi_vtxprob_cs_X3872);
										cs_X_JPiPi_Chi2_X3872->push_back((double)(JPiPi_vFit_vertex_cs_X3872->chiSquared()));
										cs_X_JPiPi_ndof_X3872->push_back((double)(JPiPi_vFit_vertex_cs_X3872->degreesOfFreedom()));
										cs_X_JPiPi_px_X3872->push_back(JPiPi_vFit_cs_X3872->currentState().kinematicParameters().momentum().x());
										cs_X_JPiPi_py_X3872->push_back(JPiPi_vFit_cs_X3872->currentState().kinematicParameters().momentum().y());
										cs_X_JPiPi_pz_X3872->push_back(JPiPi_vFit_cs_X3872->currentState().kinematicParameters().momentum().z());
										if (JPiPi_vFit_cs_X3872->currentState().kinematicParametersError().matrix()(6, 6) > 0)
										{
										cs_X_JPiPi_massErr_X3872->push_back(sqrt(JPiPi_vFit_cs_X3872->currentState().kinematicParametersError().matrix()(6, 6)));
										}		
										else
										{
										cs_X_JPiPi_massErr_X3872->push_back(-9);
										}

										cs_X_mass_X3872->push_back(X_vFit_cs_X3872->currentState().mass());
										cs_X_VtxProb_X3872->push_back(X_vtxprob_cs_X3872);
										cs_X_Chi2_X3872->push_back((double)(X_vFit_vertex_cs_X3872->chiSquared()));
										cs_X_ndof_X3872->push_back((double)(X_vFit_vertex_cs_X3872->degreesOfFreedom()));
										cs_X_px_X3872->push_back(X_vFit_cs_X3872->currentState().kinematicParameters().momentum().x());
										cs_X_py_X3872->push_back(X_vFit_cs_X3872->currentState().kinematicParameters().momentum().y());
										cs_X_pz_X3872->push_back(X_vFit_cs_X3872->currentState().kinematicParameters().momentum().z());

										if (X_vFit_cs_X3872->currentState().kinematicParametersError().matrix()(6, 6) > 0)
										{
										cs_X_massErr_X3872->push_back(sqrt(X_vFit_cs_X3872->currentState().kinematicParametersError().matrix()(6, 6)));
										}		
										else
										{
										cs_X_massErr_X3872->push_back(-9);
										}
									}
								}
							}
							if (flag_vtx == false)
							{
								cs_X_Jpsi1_mass->push_back(-9); 
								cs_X_Jpsi1_VtxProb->push_back(-9); 
								cs_X_Jpsi1_Chi2->push_back(-9);
								cs_X_Jpsi1_ndof->push_back(-9);
								cs_X_Jpsi1_px->push_back(-99999); 
								cs_X_Jpsi1_py->push_back(-99999); 
								cs_X_Jpsi1_pz->push_back(-99999); 
								cs_X_Jpsi1_massErr->push_back(-9);	
							 	cs_X_Jpsi2_mass->push_back(-9); 
								cs_X_Jpsi2_VtxProb->push_back(-9);
								cs_X_Jpsi2_Chi2->push_back(-9);
								cs_X_Jpsi2_ndof->push_back(-9); 
								cs_X_Jpsi2_px->push_back(-99999); 
								cs_X_Jpsi2_py->push_back(-99999); 
	 							cs_X_Jpsi2_pz->push_back(-99999); 
								cs_X_Jpsi2_massErr->push_back(-9);
								cs_X_JPiPi_mass->push_back(-9);
								cs_X_JPiPi_VtxProb->push_back(-9);
								cs_X_JPiPi_Chi2->push_back(-9);
								cs_X_JPiPi_ndof->push_back(-9);
								cs_X_JPiPi_px ->push_back(-99999);
								cs_X_JPiPi_py->push_back(-99999);
								cs_X_JPiPi_pz->push_back(-99999);
								cs_X_JPiPi_massErr->push_back(-9);
							}
							if ( flag_cs_psi2s == false)
							{
								cs_X_mass_Psi2S->push_back(-9);
								cs_X_VtxProb_Psi2S->push_back(-9);
								cs_X_Chi2_Psi2S->push_back(-9);
								cs_X_ndof_Psi2S->push_back(-9);
								cs_X_px_Psi2S->push_back(-99999);
								cs_X_py_Psi2S->push_back(-99999);
								cs_X_pz_Psi2S->push_back(-99999);
								cs_X_massErr_Psi2S->push_back(-9);
								cs_X_JPiPi_mass_Psi2S->push_back(-9);
								cs_X_JPiPi_VtxProb_Psi2S->push_back(-9);
								cs_X_JPiPi_Chi2_Psi2S->push_back(-9);
								cs_X_JPiPi_ndof_Psi2S->push_back(-9);
								cs_X_JPiPi_px_Psi2S->push_back(-99999);
								cs_X_JPiPi_py_Psi2S->push_back(-99999);
								cs_X_JPiPi_pz_Psi2S->push_back(-99999);
								cs_X_JPiPi_massErr_Psi2S->push_back(-9);
							}
							if ( flag_cs_x3872 == false)
							{
								cs_X_mass_X3872->push_back(-9);
								cs_X_VtxProb_X3872->push_back(-9);
								cs_X_Chi2_X3872->push_back(-9);
								cs_X_ndof_X3872->push_back(-9);
								cs_X_px_X3872->push_back(-99999);
								cs_X_py_X3872->push_back(-99999);
								cs_X_pz_X3872->push_back(-99999);
								cs_X_massErr_X3872->push_back(-9);
								cs_X_JPiPi_mass_X3872->push_back(-9);
								cs_X_JPiPi_VtxProb_X3872->push_back(-9);
								cs_X_JPiPi_Chi2_X3872->push_back(-9);
								cs_X_JPiPi_ndof_X3872->push_back(-9);
								cs_X_JPiPi_px_X3872->push_back(-99999);
								cs_X_JPiPi_py_X3872->push_back(-99999);
								cs_X_JPiPi_pz_X3872->push_back(-99999);
								cs_X_JPiPi_massErr_X3872->push_back(-9);
							}
						}	  // itrack2
					}		  // itrack1
				}			  // mu4_loop
			}				  // mu3_loop
		}					  // for (std::vector < pat::Muon >::const_iterator iMuon2 = iMuon1 + 1;
	}						  // for (std::vector < pat::Muon >::const_iterator iMuon1 = thePATMuonHandle->begin();
	if (X_VtxProb->size() > 0 || doMC)
	{
		X_One_Tree_->Fill();
	}
	if (Debug_)
	{
	}

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

	TrigRes->clear();
	TrigNames->clear();
	L1TrigRes->clear();
	MatchJpsiTrigNames->clear();
	runNum = 0;
	evtNum = 0;
	lumiNum = 0;
	nGoodPrimVtx = 0;
	priVtxX = 0;
	priVtxY = 0;
	priVtxZ = 0;
	priVtxXErr = 0;
	priVtxYErr = 0;
	priVtxZErr = 0;
	priVtxChiNorm = 0;
	priVtxChi = 0;
	priVtxCL = 0;
	priVtxXCorrY->clear();
	priVtxXCorrZ->clear();
	priVtxYCorrZ->clear();

	nMu = 0;
	muPx->clear();
	muPy->clear();
	muPz->clear();
	muD0BS->clear();
	muD0EBS->clear();
	muD3dBS->clear();
	muD3dEBS->clear();
	muD0PV->clear();
	muD0EPV->clear();
	muDzPV->clear();
	muDzEPV->clear();
	muCharge->clear();
	muTrackIso->clear();

	muIsPatLooseMuon->clear();
	muIsPatTightMuon->clear();
	muIsPatSoftMuon->clear();
	muIsPatMediumMuon->clear();
	muJpsiFilterRes->clear();

	X_mu1Id->clear();
	X_mu2Id->clear();
	X_mu3Id->clear();
	X_mu4Id->clear();
	X_mass->clear();
	X_VtxProb->clear();
	X_Chi2->clear();
	X_ndof->clear();
	X_px->clear();
	X_py->clear();
	X_pz->clear();
	X_massErr->clear();
	X_JPiPi_mass->clear();
	X_JPiPi_VtxProb->clear();
	X_JPiPi_Chi2->clear();
	X_JPiPi_ndof->clear();
	X_JPiPi_px->clear();
	X_JPiPi_py->clear();
	X_JPiPi_pz->clear();
	X_JPiPi_massErr->clear();
	X_Jpsi1_mass->clear();
	X_Jpsi1_VtxProb->clear();
	X_Jpsi1_Chi2->clear();
	X_Jpsi1_ndof->clear();
	X_Jpsi1_px->clear();
	X_Jpsi1_py->clear();
	X_Jpsi1_pz->clear();
	X_Jpsi1_massErr->clear();
	X_Jpsi2_mass->clear();
	X_Jpsi2_VtxProb->clear();
	X_Jpsi2_Chi2->clear();
	X_Jpsi2_ndof->clear();
	X_Jpsi2_px->clear();
	X_Jpsi2_py->clear();
	X_Jpsi2_pz->clear();
	X_Jpsi2_massErr->clear();
	X_JPiPi_Pi1Id->clear();
	X_JPiPi_Pi2Id->clear();
	X_JPiPi_Pi1px->clear();
	X_JPiPi_Pi1py->clear();
	X_JPiPi_Pi1pz->clear();
	X_JPiPi_Pi2px->clear();
	X_JPiPi_Pi2py->clear();
	X_JPiPi_Pi2pz->clear();

	// mass constrain variables
 	cs_X_Jpsi1_mass->clear();
	cs_X_Jpsi1_VtxProb->clear();
	cs_X_Jpsi1_Chi2->clear();
	cs_X_Jpsi1_ndof->clear();
	cs_X_Jpsi1_px->clear();
	cs_X_Jpsi1_py->clear();
	cs_X_Jpsi1_pz->clear();
	cs_X_Jpsi1_massErr->clear();
 	cs_X_Jpsi2_mass->clear();
	cs_X_Jpsi2_VtxProb->clear();
	cs_X_Jpsi2_Chi2->clear();
	cs_X_Jpsi2_ndof->clear();
	cs_X_Jpsi2_px->clear();
	cs_X_Jpsi2_py->clear();
	cs_X_Jpsi2_pz->clear();
	cs_X_Jpsi2_massErr->clear();
 	cs_X_JPiPi_mass->clear();
	cs_X_JPiPi_VtxProb->clear();
	cs_X_JPiPi_Chi2->clear();
	cs_X_JPiPi_ndof->clear();
	cs_X_JPiPi_px ->clear();
	cs_X_JPiPi_py->clear();
	cs_X_JPiPi_pz->clear();
	cs_X_JPiPi_massErr->clear();
 	cs_X_mass_Psi2S->clear();
	cs_X_VtxProb_Psi2S->clear();
	cs_X_Chi2_Psi2S->clear();
	cs_X_ndof_Psi2S->clear();
	cs_X_px_Psi2S->clear();
	cs_X_py_Psi2S->clear();
	cs_X_pz_Psi2S->clear();
	cs_X_massErr_Psi2S->clear();
 	cs_X_JPiPi_mass_Psi2S->clear();
	cs_X_JPiPi_VtxProb_Psi2S->clear();
	cs_X_JPiPi_Chi2_Psi2S->clear();
	cs_X_JPiPi_ndof_Psi2S->clear();
	cs_X_JPiPi_px_Psi2S->clear();
	cs_X_JPiPi_py_Psi2S->clear();
	cs_X_JPiPi_pz_Psi2S->clear();
	cs_X_JPiPi_massErr_Psi2S->clear();
 	cs_X_mass_X3872->clear();
	cs_X_VtxProb_X3872->clear();
	cs_X_Chi2_X3872->clear();
	cs_X_ndof_X3872->clear();
	cs_X_px_X3872->clear();
	cs_X_py_X3872->clear();
	cs_X_pz_X3872->clear();
	cs_X_massErr_X3872->clear();
 	cs_X_JPiPi_mass_X3872->clear();
	cs_X_JPiPi_VtxProb_X3872->clear();
	cs_X_JPiPi_Chi2_X3872->clear();
	cs_X_JPiPi_ndof_X3872->clear();
	cs_X_JPiPi_px_X3872->clear();
	cs_X_JPiPi_py_X3872->clear();
	cs_X_JPiPi_pz_X3872->clear();
	cs_X_JPiPi_massErr_X3872->clear();
} // analyze
// 
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
	X_One_Tree_ = fs->make<TTree>("X_data", "X(3872) Data");

	X_One_Tree_->Branch("TrigRes", &TrigRes);
	X_One_Tree_->Branch("TrigNames", &TrigNames);
	X_One_Tree_->Branch("L1TrigRes", &L1TrigRes);
	X_One_Tree_->Branch("MatchJpsiTrigNames", &MatchJpsiTrigNames);

	X_One_Tree_->Branch("evtNum", &evtNum, "evtNum/i");
	X_One_Tree_->Branch("runNum", &runNum, "runNum/i");
	X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
	X_One_Tree_->Branch("nGoodPrimVtx", &nGoodPrimVtx, "nGoodPrimVtx/i");

	X_One_Tree_->Branch("priVtxX", &priVtxX, "priVtxX/f");
	X_One_Tree_->Branch("priVtxY", &priVtxY, "priVtxY/f");
	X_One_Tree_->Branch("priVtxZ", &priVtxZ, "priVtxZ/f");
	X_One_Tree_->Branch("priVtxXErr", &priVtxXErr, "priVtxXErr/f");
	X_One_Tree_->Branch("priVtxYErr", &priVtxYErr, "priVtxYErr/f");
	X_One_Tree_->Branch("priVtxZErr", &priVtxZErr, "priVtxZErr/f");
	X_One_Tree_->Branch("priVtxChiNorm", &priVtxChiNorm, "priVtxChiNorm/f");
	X_One_Tree_->Branch("priVtxChi", &priVtxChi, "priVtxChi/f");
	X_One_Tree_->Branch("priVtxCL", &priVtxCL, "priVtxCL/f");

	X_One_Tree_->Branch("priVtxXCorrY", &priVtxXCorrY);
	X_One_Tree_->Branch("priVtxXCorrZ", &priVtxXCorrZ);
	X_One_Tree_->Branch("priVtxYCorrZ", &priVtxYCorrZ);

	X_One_Tree_->Branch("nMu", &nMu, "nMu/i");
	X_One_Tree_->Branch("muPx", &muPx);
	X_One_Tree_->Branch("muPy", &muPy);
	X_One_Tree_->Branch("muPz", &muPz);
	X_One_Tree_->Branch("muD0BS", &muD0BS);
	X_One_Tree_->Branch("muD0EBS", &muD0EBS);
	X_One_Tree_->Branch("muD3dBS", &muD3dBS);
	X_One_Tree_->Branch("muD3dEBS", &muD3dEBS);
	X_One_Tree_->Branch("muD0PV", &muD0PV);
	X_One_Tree_->Branch("muD0EPV", &muD0EPV);
	X_One_Tree_->Branch("muDzPV", &muDzPV);
	X_One_Tree_->Branch("muDzEPV", &muDzEPV);
	X_One_Tree_->Branch("muCharge", &muCharge);
	X_One_Tree_->Branch("muTrackIso", &muTrackIso);

	X_One_Tree_->Branch("muIsPatLooseMuon", &muIsPatLooseMuon);
	X_One_Tree_->Branch("muIsPatTightMuon", &muIsPatTightMuon);
	X_One_Tree_->Branch("muIsPatSoftMuon", &muIsPatSoftMuon);
	X_One_Tree_->Branch("muIsPatMediumMuon", &muIsPatMediumMuon);

	X_One_Tree_->Branch("muJpsiFilterRes", &muJpsiFilterRes);

	X_One_Tree_->Branch("X_mu1Id", &X_mu1Id);
	X_One_Tree_->Branch("X_mu2Id", &X_mu2Id);
	X_One_Tree_->Branch("X_mu3Id", &X_mu3Id);
	X_One_Tree_->Branch("X_mu4Id", &X_mu4Id);
	X_One_Tree_->Branch("X_mass", &X_mass);
	X_One_Tree_->Branch("X_VtxProb", &X_VtxProb);
	X_One_Tree_->Branch("X_Chi2", &X_Chi2);
	X_One_Tree_->Branch("X_ndof", &X_ndof);
	X_One_Tree_->Branch("X_px", &X_px);
	X_One_Tree_->Branch("X_py", &X_py);
	X_One_Tree_->Branch("X_pz", &X_pz);
	X_One_Tree_->Branch("X_massErr", &X_massErr);
	X_One_Tree_->Branch("X_JPiPi_mass", &X_JPiPi_mass);
	X_One_Tree_->Branch("X_JPiPi_VtxProb", &X_JPiPi_VtxProb);
	X_One_Tree_->Branch("X_JPiPi_Chi2", &X_JPiPi_Chi2);
	X_One_Tree_->Branch("X_JPiPi_ndof", &X_JPiPi_ndof);
	X_One_Tree_->Branch("X_JPiPi_px", &X_JPiPi_px);
	X_One_Tree_->Branch("X_JPiPi_py", &X_JPiPi_py);
	X_One_Tree_->Branch("X_JPiPi_pz", &X_JPiPi_pz);
	X_One_Tree_->Branch("X_JPiPi_massErr", &X_JPiPi_massErr);
	X_One_Tree_->Branch("X_Jpsi1_mass", &X_Jpsi1_mass);
	X_One_Tree_->Branch("X_Jpsi1_VtxProb", &X_Jpsi1_VtxProb);
	X_One_Tree_->Branch("X_Jpsi1_Chi2", &X_Jpsi1_Chi2);
	X_One_Tree_->Branch("X_Jpsi1_ndof", &X_Jpsi1_ndof);
	X_One_Tree_->Branch("X_Jpsi1_px", &X_Jpsi1_px);
	X_One_Tree_->Branch("X_Jpsi1_py", &X_Jpsi1_py);
	X_One_Tree_->Branch("X_Jpsi1_pz", &X_Jpsi1_pz);
	X_One_Tree_->Branch("X_Jpsi1_massErr", &X_Jpsi1_massErr);
	X_One_Tree_->Branch("X_Jpsi2_mass", &X_Jpsi2_mass);
	X_One_Tree_->Branch("X_Jpsi2_VtxProb", &X_Jpsi2_VtxProb);
	X_One_Tree_->Branch("X_Jpsi2_Chi2", &X_Jpsi2_Chi2);
	X_One_Tree_->Branch("X_Jpsi2_ndof", &X_Jpsi2_ndof);
	X_One_Tree_->Branch("X_Jpsi2_px", &X_Jpsi2_px);
	X_One_Tree_->Branch("X_Jpsi2_py", &X_Jpsi2_py);
	X_One_Tree_->Branch("X_Jpsi2_pz", &X_Jpsi2_pz);
	X_One_Tree_->Branch("X_Jpsi2_massErr", &X_Jpsi2_massErr);
	X_One_Tree_->Branch("X_JPiPi_Pi1Id", &X_JPiPi_Pi1Id);
	X_One_Tree_->Branch("X_JPiPi_Pi2Id", &X_JPiPi_Pi2Id);
	X_One_Tree_->Branch("X_JPiPi_Pi1px", &X_JPiPi_Pi1px);
	X_One_Tree_->Branch("X_JPiPi_Pi1py", &X_JPiPi_Pi1py);
	X_One_Tree_->Branch("X_JPiPi_Pi1pz", &X_JPiPi_Pi1pz);
	X_One_Tree_->Branch("X_JPiPi_Pi2px", &X_JPiPi_Pi2px);
	X_One_Tree_->Branch("X_JPiPi_Pi2py", &X_JPiPi_Pi2py);
	X_One_Tree_->Branch("X_JPiPi_Pi2pz", &X_JPiPi_Pi2pz);

	// mass constrain variables
 	X_One_Tree_->Branch("cs_X_Jpsi1_mass", &cs_X_Jpsi1_mass);
	X_One_Tree_->Branch("cs_X_Jpsi1_VtxProb", &cs_X_Jpsi1_VtxProb);
	X_One_Tree_->Branch("cs_X_Jpsi1_Chi2", &cs_X_Jpsi1_Chi2);
	X_One_Tree_->Branch("cs_X_Jpsi1_ndof", &cs_X_Jpsi1_ndof);
	X_One_Tree_->Branch("cs_X_Jpsi1_px", &cs_X_Jpsi1_px);
	X_One_Tree_->Branch("cs_X_Jpsi1_py", &cs_X_Jpsi1_py);
	X_One_Tree_->Branch("cs_X_Jpsi1_pz", &cs_X_Jpsi1_pz);
	X_One_Tree_->Branch("cs_X_Jpsi1_massErr", &cs_X_Jpsi1_massErr);
 	X_One_Tree_->Branch("cs_X_Jpsi2_mass", &cs_X_Jpsi2_mass);
	X_One_Tree_->Branch("cs_X_Jpsi2_VtxProb", &cs_X_Jpsi2_VtxProb);
	X_One_Tree_->Branch("cs_X_Jpsi2_Chi2", &cs_X_Jpsi2_Chi2);
	X_One_Tree_->Branch("cs_X_Jpsi2_ndof", &cs_X_Jpsi2_ndof);
	X_One_Tree_->Branch("cs_X_Jpsi2_px", &cs_X_Jpsi2_px);
	X_One_Tree_->Branch("cs_X_Jpsi2_py", &cs_X_Jpsi2_py);
	X_One_Tree_->Branch("cs_X_Jpsi2_pz", &cs_X_Jpsi2_pz);
	X_One_Tree_->Branch("cs_X_Jpsi2_massErr", &cs_X_Jpsi2_massErr);
 	X_One_Tree_->Branch("cs_X_JPiPi_mass", &cs_X_JPiPi_mass);
	X_One_Tree_->Branch("cs_X_JPiPi_VtxProb", &cs_X_JPiPi_VtxProb);
	X_One_Tree_->Branch("cs_X_JPiPi_Chi2", &cs_X_JPiPi_Chi2);
	X_One_Tree_->Branch("cs_X_JPiPi_ndof", &cs_X_JPiPi_ndof);
	X_One_Tree_->Branch("cs_X_JPiPi_px ", &cs_X_JPiPi_px );
	X_One_Tree_->Branch("cs_X_JPiPi_py", &cs_X_JPiPi_py);
	X_One_Tree_->Branch("cs_X_JPiPi_pz", &cs_X_JPiPi_pz);
	X_One_Tree_->Branch("cs_X_JPiPi_massErr", &cs_X_JPiPi_massErr);
 	X_One_Tree_->Branch("cs_X_mass_Psi2S", &cs_X_mass_Psi2S);
	X_One_Tree_->Branch("cs_X_VtxProb_Psi2S", &cs_X_VtxProb_Psi2S);
	X_One_Tree_->Branch("cs_X_Chi2_Psi2S", &cs_X_Chi2_Psi2S);
	X_One_Tree_->Branch("cs_X_ndof_Psi2S", &cs_X_ndof_Psi2S);
	X_One_Tree_->Branch("cs_X_px_Psi2S", &cs_X_px_Psi2S);
	X_One_Tree_->Branch("cs_X_py_Psi2S", &cs_X_py_Psi2S);
	X_One_Tree_->Branch("cs_X_pz_Psi2S", &cs_X_pz_Psi2S);
	X_One_Tree_->Branch("cs_X_massErr_Psi2S", &cs_X_massErr_Psi2S);
 	X_One_Tree_->Branch("cs_X_JPiPi_mass_Psi2S", &cs_X_JPiPi_mass_Psi2S);
	X_One_Tree_->Branch("cs_X_JPiPi_VtxProb_Psi2S", &cs_X_JPiPi_VtxProb_Psi2S);
	X_One_Tree_->Branch("cs_X_JPiPi_Chi2_Psi2S", &cs_X_JPiPi_Chi2_Psi2S);
	X_One_Tree_->Branch("cs_X_JPiPi_ndof_Psi2S", &cs_X_JPiPi_ndof_Psi2S);
	X_One_Tree_->Branch("cs_X_JPiPi_px_Psi2S", &cs_X_JPiPi_px_Psi2S);
	X_One_Tree_->Branch("cs_X_JPiPi_py_Psi2S", &cs_X_JPiPi_py_Psi2S);
	X_One_Tree_->Branch("cs_X_JPiPi_pz_Psi2S", &cs_X_JPiPi_pz_Psi2S);
	X_One_Tree_->Branch("cs_X_JPiPi_massErr_Psi2S", &cs_X_JPiPi_massErr_Psi2S);
 	X_One_Tree_->Branch("cs_X_mass_X3872", &cs_X_mass_X3872);
	X_One_Tree_->Branch("cs_X_VtxProb_X3872", &cs_X_VtxProb_X3872);
	X_One_Tree_->Branch("cs_X_Chi2_X3872", &cs_X_Chi2_X3872);
	X_One_Tree_->Branch("cs_X_ndof_X3872", &cs_X_ndof_X3872);
	X_One_Tree_->Branch("cs_X_px_X3872", &cs_X_px_X3872);
	X_One_Tree_->Branch("cs_X_py_X3872", &cs_X_py_X3872);
	X_One_Tree_->Branch("cs_X_pz_X3872", &cs_X_pz_X3872);
	X_One_Tree_->Branch("cs_X_massErr_X3872", &cs_X_massErr_X3872);
 	X_One_Tree_->Branch("cs_X_JPiPi_mass_X3872", &cs_X_JPiPi_mass_X3872);
	X_One_Tree_->Branch("cs_X_JPiPi_VtxProb_X3872", &cs_X_JPiPi_VtxProb_X3872);
	X_One_Tree_->Branch("cs_X_JPiPi_Chi2_X3872", &cs_X_JPiPi_Chi2_X3872);
	X_One_Tree_->Branch("cs_X_JPiPi_ndof_X3872", &cs_X_JPiPi_ndof_X3872);
	X_One_Tree_->Branch("cs_X_JPiPi_px_X3872", &cs_X_JPiPi_px_X3872);
	X_One_Tree_->Branch("cs_X_JPiPi_py_X3872", &cs_X_JPiPi_py_X3872);
	X_One_Tree_->Branch("cs_X_JPiPi_pz_X3872", &cs_X_JPiPi_pz_X3872);
	X_One_Tree_->Branch("cs_X_JPiPi_massErr_X3872", &cs_X_JPiPi_massErr_X3872);
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
		X_One_Tree_->Branch("Match_mu1px", &Match_mu1px);
		X_One_Tree_->Branch("Match_mu1py", &Match_mu1py);
		X_One_Tree_->Branch("Match_mu1pz", &Match_mu1pz);
		X_One_Tree_->Branch("Match_mu2px", &Match_mu2px);
		X_One_Tree_->Branch("Match_mu2py", &Match_mu2py);
		X_One_Tree_->Branch("Match_mu2pz", &Match_mu2pz);
		X_One_Tree_->Branch("Match_mu3px", &Match_mu3px);
		X_One_Tree_->Branch("Match_mu3py", &Match_mu3py);
		X_One_Tree_->Branch("Match_mu3pz", &Match_mu3pz);
		X_One_Tree_->Branch("Match_mu4px", &Match_mu4px);
		X_One_Tree_->Branch("Match_mu4py", &Match_mu4py);
		X_One_Tree_->Branch("Match_mu4pz", &Match_mu4pz);

		X_One_Tree_->Branch("Match_pi1px", &Match_pi1px);
		X_One_Tree_->Branch("Match_pi1py", &Match_pi1py);
		X_One_Tree_->Branch("Match_pi1pz", &Match_pi1pz);
		X_One_Tree_->Branch("Match_pi2px", &Match_pi2px);
		X_One_Tree_->Branch("Match_pi2py", &Match_pi2py);
		X_One_Tree_->Branch("Match_pi2pz", &Match_pi2pz);
	} // if(doMC)
} // begin Job

// ------------ method called once each job just after ending the event loop  ------------
void MultiLepPAT::endJob()
{
	X_One_Tree_->GetDirectory()->cd();
	X_One_Tree_->Write();
}

// define this as a plug-in
DEFINE_FWK_MODULE(MultiLepPAT);
