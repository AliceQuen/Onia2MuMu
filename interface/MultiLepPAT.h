/******************************************************************************
 *  [File] 
 *      MultiLepPAT.h 
 *  [Class]      
 *      MultiLepPAT 
 *  [Directory] 
 *      TPS-Onia2MuMu/interface/MultiLepPAT.h
 *  [Description]
 *      Make rootTuple for JPsi+Jpsi+Upsilon reconstruction
 *  [Implementation]
 *     <Notes on implementation>
 *  [Original Author]
 *  
 *  [Note]
 *      20240626 [Eric Wang]
 *          Upsilon is abbreviated as "Ups"
 *      20240705 [Eric Wang]
 *          Adding new member functions for the following purposes:
 *              Judging if the event 
******************************************************************************/



#ifndef _MultiLepPAT_h
#define _MultiLepPAT_h

// system include files
#include <memory>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // xining MINIAODtest

// user include files
#include "../interface/VertexReProducer.h"

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

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//
// class decleration
//

using std::vector;
using std::string;
using namespace edm;
using namespace reco;

class MultiLepPAT : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit MultiLepPAT(const ParameterSet&);
    ~MultiLepPAT();

    // Note: always ensures that in muList_t, the first muon is the one with the smaller index.
    using muon_t   = RefCountedKinematicParticle;
    using muList_t = std::pair< vector<muon_t>, vector<uint> >;

    // Using enumeration to define trigger types [Annotation by Eric Wang, 20240626] 
    static constexpr unsigned int trigCount = 50;
    static constexpr unsigned int trigTypeCount = 3;

    enum class trigType{JPSI, UPS, PHI};
  
private:
    // Framework methods
    virtual void beginJob() ;
    virtual void beginRun(Run const & iRun, EventSetup const& iSetup);
    virtual void analyze(const Event&, const EventSetup&);
    virtual void endJob() ;

    // Essential methods (ctau)
    static double GetcTau(RefCountedKinematicVertex&   decayVrtx, 
                          RefCountedKinematicParticle& kinePart, 
                          Vertex&                      bs );
    static double GetcTauErr(RefCountedKinematicVertex& decayVrtx, 
                             RefCountedKinematicParticle& kinePart, 
                             Vertex&                      bs );
    
    // Essential methods (deltaR)
    static double deltaR(double eta1, double phi1, double eta2, double phi2);
    
    // Essential methods (pT, eta, phi all in one)
    static void getDynamics(const RefCountedKinematicParticle& arg_Part,
                            double& res_pt,   double& res_eta, double& res_phi);
    static void getDynamics(double  arg_mass, double  arg_px,  double  arg_py, double arg_pz,
                            double& res_pt,   double& res_eta, double& res_phi);


 
    //add token here
    edm::EDGetTokenT< L1GlobalTriggerReadoutRecord    > gtRecordToken_;
    edm::EDGetTokenT< BeamSpot                        > gtbeamspotToken_;
    edm::EDGetTokenT< VertexCollection                > gtprimaryVtxToken_;
    edm::EDGetTokenT< edm::View<pat::Muon>            > gtpatmuonToken_; // MINIAOD
    edm::EDGetTokenT< edm::TriggerResults             > gttriggerToken_;
    edm::EDGetTokenT< edm::View<pat::PackedCandidate> > trackToken_; // MINIAOD
    edm::EDGetTokenT< reco::GenParticleCollection     > genParticlesToken_;

    // ----------member data ---------------------------

    virtual int GetHitsBefore(const EventSetup& setup, const vector<reco::Track>& trks, RefCountedKinematicTree& fitTr){
        return -1;
    }

    virtual int GetMissesAfter(const EventSetup& setup, const vector<reco::Track>& trks, RefCountedKinematicTree& VrtxTree){
        return -1;
    }

    // Some useful methods for maintenance, readability and elegance    
    static void tracksToMuonPair(vector<RefCountedKinematicParticle>&        arg_MuonResults,
                                 KinematicParticleFactoryFromTransientTrack& arg_MuFactory,
                                 const MagneticField&                        arg_bField,
                                 const TrackRef arg_Trk1,     const TrackRef arg_Trk2       );

    static bool particlesToVtx(const vector<RefCountedKinematicParticle>&  arg_MuonResults);
    static bool particlesToVtx(const vector<RefCountedKinematicParticle>&  arg_MuonResults,
                               const string&                               arg_Message);
    static bool particlesToVtx(RefCountedKinematicTree&                    arg_VertexFitTree,
                               const vector<RefCountedKinematicParticle>&  arg_Muons,
                               const string&                               arg_Message);

    static bool extractFitRes(RefCountedKinematicTree&     arg_VtxTree,
                              RefCountedKinematicParticle& res_Part,
                              RefCountedKinematicVertex&   res_Vtx,
                              KinematicParameters&         res_Param,
                              double&                      res_MassErr);
    static bool extractFitRes(RefCountedKinematicTree&     arg_VtxTree,
                              RefCountedKinematicParticle& res_Part,
                              RefCountedKinematicVertex&   res_Vtx,
                              double&                      res_MassErr);
    static bool extractFitRes(RefCountedKinematicTree&     arg_VtxTree,
                              RefCountedKinematicVertex&   res_Vtx,
                              double&                      res_VtxProb);

    // To avoid overlapping muon pairs
    static bool isOverlapPair(const muList_t& arg_MuonPair1, 
                              const muList_t& arg_MuonPair2 );

    // Deal with "multi-candidate" issue
    static double fitResEval(double arg_massDiff_Jpsi_1, double arg_massErr_Jpsi_1,
                             double arg_massDiff_Jpsi_2, double arg_massErr_Jpsi_2,
                             double arg_massDiff_Phi,    double arg_massErr_Phi   );

    
    // Member data

    // Essentials [Annotation by Eric Wang, 20240626]
    string proccessName_;
    HLTConfigProvider hltConfig_;
    
    InputTag hlTriggerResults_;
    InputTag inputGEN_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>  magneticFieldToken_;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord>  theTTBuilderToken_;
    string  vtxSample;
    bool    doMC;
    int     MCParticle;
    bool    doJPsiMassCost;
    int     MuPixHits_c;
    int     MuSiHits_c;
    double  MuNormChi_c;
    double  MuD0_c;

    // Limits for secondary particles [Annotation by Eric Wang, 20240626]
    
    double JMaxM_c;
    double JMinM_c;
    int PiSiHits_c;
    double MuPt_c;
    double JPiPiDR_c;
    double XPiPiDR_c;
    bool UseXDr_c;	
    double JPiPiMax_c;
    double JPiPiMin_c;

    // Trigger info [Annotation by Eric Wang, 20240626]
    
    bool resolveAmbiguity_; 
    bool addXlessPrimaryVertex_;

    // Identifying triggers and filters with their name [Annotation by Eric Wang, 20240626]
    vector<string>      TriggersForJpsi_;
    vector<string>      FiltersForJpsi_;
    vector<string>      TriggersForUpsilon_;
    vector<string>      FiltersForUpsilon_;
    vector<std::string>* MatchJpsiTrigNames;
    int JpsiMatchTrig[50], UpsilonMatchTrig[50];


    virtual void getAllTriggers(   const edm::Handle<edm::TriggerResults>&     HLTresult);
    virtual bool muonMatchTrigType(const edm::View<pat::Muon>::const_iterator& muonIter,
                                   const vector<string>& trigNames, 
                                         trigType        type                          );
    
    vector<string>      TriggersForMatching_;
    vector<string>      FiltersForMatching_;
    int  MatchingTriggerResult[50];
    bool Debug_;
    double Chi_Track_;

    // PDG 2023
	static constexpr double myJpsiMass = 3.0969,   myJpsiMassErr = 0.00004;
	static constexpr double myUpsMass  = 9.4603,   myUpsMassErr  = 0.0003;
	static constexpr double myPhiMass  = 1.019455, myPhiMassErr    = 0.000020;
	static constexpr double myMuMass = 0.1056583745;
	static constexpr double myMuMassErr = 0.0000000023; // From PDG 2024
	static constexpr double myKMass = 0.493677; // Kaon mass and error
	// try
	static constexpr double myKMassErr = 0.000015; // From PDG 2024

    // general restrictions for vtx prob
    static constexpr double VtxProbCut = 0.01;

    // Constructing TTree object [Annotation by Eric Wang, 20240626]
    
    TTree* X_One_Tree_;
    
    unsigned int        runNum, evtNum, lumiNum;
    unsigned int        nGoodPrimVtx;
    
    vector<unsigned int>*   trigRes;
    vector<std::string>*    trigNames;
    vector<unsigned int>*   L1TT;
    vector<std::string>*    MatchTriggerNames;

    // primary vertices [Annotation by Eric Wang, 20240626]
    float               priVtxX,    priVtxY,    priVtxZ, 
                        priVtxXE,   priVtxYE,   priVtxZE, 
                        priVtxChiNorm, priVtxChi, priVtxCL;
    vector<float>       *PriVtxXCorrX, *PriVtxXCorrY, *PriVtxXCorrZ;
    vector<double>      *PriVtxXCorrEX, *PriVtxXCorrEY, *PriVtxXCorrEZ;
    vector<float>	    *PriVtxXCorrC2, *PriVtxXCorrCL;
    
    // all muons: kinematic and other information [Annotation by Eric Wang, 20240626]
    unsigned int        nMu;
    vector<float>       *muPx, *muPy, *muPz, 
                        *muD0, *muD0E, *muDz, 
                        *muChi2, *muGlChi2,  *mufHits;   
    vector<bool>        *muFirstBarrel, *muFirstEndCap;
    vector<float>       *muDzVtx, *muDxyVtx;   
    vector<int>         *muNDF, *muGlNDF, *muPhits, *muShits, *muGlMuHits, *muType, *muQual;
    vector<int>         *muTrack;
    vector<float>       *muCharge;
    vector<float>       *muIsoratio;

    // all muons: selection result [Annotation by Eric Wang, 20240626]
    vector<int>         *muIsGoodLooseMuon,         *muIsGoodLooseMuonNew, 
                        *muIsGoodSoftMuonNewIlse,   *muIsGoodSoftMuonNewIlseMod, 
                        *muIsGoodTightMuon,         *muIsJpsiTrigMatch,         
                        *muIsUpsTrigMatch,          *munMatchedSeg,
                        *muJpsiFilterRes;
    vector<int>         *muIsPatLooseMuon, *muIsPatTightMuon, *muIsPatSoftMuon, *muIsPatMediumMuon;
    
    //for Maksat trigger match [Annotation by Eric Wang, 20240626]
    vector<int> *muUpsVrtxMatch, *muL3TriggerMatch;
    
    //added by zhenhu for MuonID
    vector<float>  *muMVAMuonID, *musegmentCompatibility; 
    //for Stoyan slope pull
    vector<float>  *mupulldXdZ_pos_noArb, *mupulldYdZ_pos_noArb;
    //addition 3,4,5
    vector<float>  *mupulldXdZ_pos_ArbDef, *mupulldYdZ_pos_ArbDef;
    vector<float>  *mupulldXdZ_pos_ArbST, *mupulldYdZ_pos_ArbST;
    vector<float>  *mupulldXdZ_pos_noArb_any, *mupulldYdZ_pos_noArb_any;

    // Muons from Jpsi and Upsilon.
    vector<unsigned int> *Jpsi_1_mu_1_Idx, *Jpsi_1_mu_2_Idx, 
                         *Jpsi_2_mu_1_Idx, *Jpsi_2_mu_2_Idx,
                             *Phi_K_1_Idx,     *Phi_K_2_Idx;
    
    // Branches for the supposed kaon tracks from Phi decay.
    vector<float>        *Phi_K_1_px, *Phi_K_1_py, *Phi_K_1_pz,
                         *Phi_K_2_px, *Phi_K_2_py, *Phi_K_2_pz,
                         *Phi_K_1_eta, *Phi_K_1_phi, *Phi_K_1_pt,
                         *Phi_K_2_eta, *Phi_K_2_phi, *Phi_K_2_pt;

    // Reconstructed Jpsi and Upsilon.
    // Note: Used "vector<T>* a, b" instead of "vector<T> *a, *b"
    vector<float> *Jpsi_1_mass, *Jpsi_1_massErr, *Jpsi_1_massDiff,
                  *Jpsi_2_mass, *Jpsi_2_massErr, *Jpsi_2_massDiff,
                     *Phi_mass,    *Phi_massErr,    *Phi_massDiff ;
               
    vector<float> *Jpsi_1_ctau, *Jpsi_1_ctauErr, *Jpsi_1_Chi2, *Jpsi_1_ndof, *Jpsi_1_VtxProb,
                  *Jpsi_2_ctau, *Jpsi_2_ctauErr, *Jpsi_2_Chi2, *Jpsi_2_ndof, *Jpsi_2_VtxProb,
                                                    *Phi_Chi2,    *Phi_ndof,    *Phi_VtxProb;
                  
    vector<float> *Jpsi_1_phi, *Jpsi_1_eta, *Jpsi_1_pt,
                  *Jpsi_2_phi, *Jpsi_2_eta, *Jpsi_2_pt,
                     *Phi_phi,    *Phi_eta,    *Phi_pt;
               
    vector<float> *Jpsi_1_px, *Jpsi_1_py, *Jpsi_1_pz,
                  *Jpsi_2_px, *Jpsi_2_py, *Jpsi_2_pz,
                     *Phi_px,    *Phi_py,    *Phi_pz;
    // Primary vertex reconstructied from Jpsi and Upsilon.              
    vector<float>    *Pri_mass,  *Pri_massErr,
                     *Pri_ctau,  *Pri_ctauErr, *Pri_Chi2, *Pri_ndof, *Pri_VtxProb,
                     *Pri_px,    *Pri_py,    *Pri_pz, 
                     *Pri_phi,   *Pri_eta,   *Pri_pt;  

    //doMC
    vector<float> 
    *MC_X_px,
    *MC_X_py,
    *MC_X_pz,
    *MC_X_mass,
    *MC_Dau_Jpsipx,
    *MC_Dau_Jpsipy,
    *MC_Dau_Jpsipz,
    *MC_Dau_Jpsimass,
    *MC_Dau_psi2spx,
    *MC_Dau_psi2spy,
    *MC_Dau_psi2spz,
    *MC_Dau_psi2smass,
    *MC_Granddau_mu1px,
    *MC_Granddau_mu1py,
    *MC_Granddau_mu1pz,
    *MC_Granddau_mu2px,
    *MC_Granddau_mu2py,
    *MC_Granddau_mu2pz,
    *MC_Granddau_Jpsipx,
    *MC_Granddau_Jpsipy,
    *MC_Granddau_Jpsipz,
    *MC_Granddau_Jpsimass,
    *MC_Granddau_pi1px,
    *MC_Granddau_pi1py,
    *MC_Granddau_pi1pz,
    *MC_Granddau_pi2px,
    *MC_Granddau_pi2py,
    *MC_Granddau_pi2pz,
    *MC_Grandgranddau_mu3px,
    *MC_Grandgranddau_mu3py,
    *MC_Grandgranddau_mu3pz,
    *MC_Grandgranddau_mu4px,
    *MC_Grandgranddau_mu4py,
    *MC_Grandgranddau_mu4pz;
    
    vector<int>
    *MC_X_chg,
    *MC_Dau_JpsipdgId,
    *MC_Dau_psi2spdgId,
    *MC_Granddau_mu1pdgId,
    *MC_Granddau_mu2pdgId,
    *MC_Granddau_JpsipdgId,
    *MC_Granddau_pi1pdgId,
    *MC_Granddau_pi2pdgId,
    *MC_Grandgranddau_mu3pdgId,
    *MC_Grandgranddau_mu4pdgId;
    
    vector<float> 
    *Match_mu1px,
    *Match_mu1py,
    *Match_mu1pz,
    *Match_mu2px,
    *Match_mu2py,
    *Match_mu2pz,
    *Match_mu3px,
    *Match_mu3py,
    *Match_mu3pz,
    *Match_mu4px,
    *Match_mu4py,
    *Match_mu4pz,
    
    *Match_pi1px,
    *Match_pi1py,
    *Match_pi1pz,
    *Match_pi2px,
    *Match_pi2py,
    *Match_pi2pz; 

    // vector<float> 
    //     *Match_Jpsi_1_mu1px, *Match_Jpsi_1_mu1py, *Match_Jpsi1_mu1pz,
    //     *Match_Jpsi_1_mu2px, *Match_Jpsi_1_mu2py, *Match_Jpsi1_mu2pz,
    //     *Match_Jpsi_2_mu1px, *Match_Jpsi_2_mu1py, *Match_Jpsi2_mu1pz,
    //     *Match_Jpsi_2_mu2px, *Match_Jpsi_2_mu2py, *Match_Jpsi2_mu2pz,
    //     *Match_Ups_mu1px,    *Match_Ups_mu1py,    *Match_Ups_mu1pz,
    //     *Match_Ups_mu2px,    *Match_Ups_mu2py,    *Match_Ups_mu2pz;


  
};

#endif
