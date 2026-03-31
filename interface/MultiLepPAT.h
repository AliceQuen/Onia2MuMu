/******************************************************************************
 *  [File] 
 *      MultiLepPAT.h 
 *  [Class]      
 *      MultiLepPAT 
 *  [Directory] 
 *      TPS-Onia2MuMu/interface/MultiLepPAT.h
 *  [Description]
 *      Make rootTuple for quarkonia+quarkonia+meson reconstruction.
 *      Supports J/psi+J/psi+phi, J/psi+J/psi+Upsilon, J/psi+Upsilon+phi
 *      via dynamic config switching.
 *  [Implementation]
 *     Refactored: modularized analyze(), externalized all cuts to config,
 *     added momentum uncertainties, MC gen-matching branches,
 *     and StringCutObjectSelector for flexible runtime cuts.
 *  [Note]
 *      20240626 [Eric Wang]
 *          Upsilon is abbreviated as "Ups"
 *      20240705 [Eric Wang]
 *          Adding new member functions.
 *      20260306 [Eric Wang - Refactor]
 *          - Split monolithic analyze() into modular helper methods.
 *          - All hardcoded thresholds extracted to config parameters.
 *          - Added momentum uncertainty branches (ptErr, pxErr, pyErr, pzErr).
 *          - Added MC gen-level branches for proper reco-gen matching.
 *          - Parameterized quarkonia states via AnalysisMode config string.
 *          - StringCutObjectSelector for muon and track selection.
 *      [Legacy variables - kept for backward compatibility]:
 *          vtxSample: Originally used to select PV collection name.
 *              Now hardcoded to "offlineSlimmedPrimaryVertices" for MINIAOD.
 *          MCParticle: Originally the PDG ID of the MC mother particle.
 *              Replaced by AnalysisMode-based automatic lookup.
 *          JPiPiDR_c, XPiPiDR_c, UseXDr_c, JPiPiMax_c, JPiPiMin_c:
 *              Legacy from X(3872)->J/psi pi pi analysis. Now replaced by
 *              the generic track pair DR and mass window cuts.
 *          PiSiHits_c: Track silicon hits cut; now in TrackSelection string.
 *          MuPixHits_c, MuSiHits_c, MuNormChi_c, MuD0_c:
 *              Muon quality cuts; now in MuonSelection string.
 *          doJPsiMassCost: Whether to apply J/psi mass constraint.
 *              Still functional but rarely changed at runtime.
******************************************************************************/

#ifndef _MultiLepPAT_h
#define _MultiLepPAT_h

// system include files
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <utility>

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

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

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

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

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

using std::vector;
using std::string;
using namespace edm;
using namespace reco;

class MultiLepPAT : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit MultiLepPAT(const ParameterSet&);
    ~MultiLepPAT();

    // Type aliases for intermediate storage
    using muon_t   = RefCountedKinematicParticle;
    using muList_t = std::pair< vector<muon_t>, vector<uint> >;

    // Trigger type enumeration
    enum class trigType{JPSI, UPS, PHI};

    // Analysis mode enumeration
    enum class AnalysisChannel { JpsiJpsiPhi, JpsiJpsiUps, JpsiUpsPhi };

    // Muon-track matching method enumeration
    enum class MuTrkMatchMethod {
        SourceCandidatePtr = 1,
        Vector = 2,
        Chi2 = 3,
        DzAssoc = 4,
        DzPv = 5
    };

    // Result structure for muon-track matching
    struct MuonTrackMatchResult {
        bool matched;
        int packedHandleIdx;
        int trackPoolIdx;
        int vertexId;
        int fromPV;
        int pvAssocQuality;
        int pdgId;
        int nCandidates;
        int methodUsed;
        float vectorRelP;
        float chi2;
        float dzSelectedPV;
        float dzAssocPV;
        float dxyAssocPV;
    };

    struct PriCandidateDiagnostics {
        bool fitValid;
        bool fitPass;
        bool assocPVPass;
        bool trackPVPass;
        bool passAny;
        int assocPVIdx;
        float maxAbsDzPV;
        float maxAbsDxyPV;
    };

    struct PhiKaonDiagnostics {
        int vertexId;
        int fromPV;
        int pvAssocQuality;
        int genMatchIdx;
        int genMatchSource;
        bool hasAssocPV;
        bool passDzPV;
        bool passDxyPV;
        bool passTrackPV;
        float dzPV;
        float dxyPV;
        float dzAssocPV;
        float dxyAssocPV;
        float genMatchChi2;
    };

    struct PhiVertexDiagnostics {
        bool fitPass;
        bool commonAssocPVPass;
        bool trackPVPass;
        bool vertexCriteriaPass;
        int commonAssocPVIdx;
        float maxAbsDzPV;
        float maxAbsDxyPV;
    };

private:
    // ======================== Framework methods ========================
    virtual void beginJob() override;
    virtual void analyze(const Event&, const EventSetup&) override;
    virtual void endJob() override;

    // ======================== Modularized analysis steps ========================
    // Step 1: Process MC gen-level truth
    void processMCGenInfo(const edm::Event& iEvent);
    // Step 2: Process HLT trigger info
    void processHLTInfo(const edm::Event& iEvent);
    // Step 3: Reconstruct primary vertex
    void reconstructPrimaryVertex(const edm::Event& iEvent);
    // Step 4: Fill muon block
    void fillMuonBlock(const edm::Event& iEvent,
                       const edm::EventSetup& iSetup,
                       const reco::Vertex& thePrimaryV);
    // Step 5: Pair muons into J/psi (or Upsilon) candidates
    void pairMuons(const edm::Handle<edm::View<pat::Muon>>& muonHandle,
                   const MagneticField& bField);
    // Step 6: Pair tracks into phi (or other meson) candidates
    void pairTracks(const std::vector<edm::View<pat::PackedCandidate>::const_iterator>& tracks,
                    const MagneticField& bField);
    // Step 7: Combine resonances and fill final candidate branches
    void combineCandidates(const reco::Vertex& beamSpotV,
                           const edm::Handle<reco::GenParticleCollection>& genParticles);
    // Step 8: MC gen-level matching for efficiency studies
    void doMCGenMatching(const edm::Handle<edm::View<pat::Muon>>& muonHandle,
                         const edm::Handle<edm::View<pat::PackedCandidate>>& trackHandle);
    // Step 9: Clear all vectors at end of event
    void clearEventData();

    // ======================== Utility methods ========================
    static double GetcTau(RefCountedKinematicVertex& decayVrtx, 
                          RefCountedKinematicParticle& kinePart, 
                          Vertex& bs);
    static double GetcTauErr(RefCountedKinematicVertex& decayVrtx, 
                             RefCountedKinematicParticle& kinePart, 
                             Vertex& bs);
    static double deltaR(double eta1, double phi1, double eta2, double phi2);
    
    // Kinematics extraction (with and without uncertainties)
    static void getDynamics(const RefCountedKinematicParticle& arg_Part,
                            double& res_pt, double& res_eta, double& res_phi);
    static void getDynamics(double arg_mass, double arg_px, double arg_py, double arg_pz,
                            double& res_pt, double& res_eta, double& res_phi);
    // Momentum uncertainty extraction from kinematic fit covariance
    static void getMomentumErrors(const RefCountedKinematicParticle& arg_Part,
                                  double& res_pxErr, double& res_pyErr, double& res_pzErr,
                                  double& res_ptErr);

    // Particle reconstruction helpers
    static void tracksToMuonPair(vector<RefCountedKinematicParticle>& arg_MuonResults,
                                 KinematicParticleFactoryFromTransientTrack& arg_MuFactory,
                                 const MagneticField& arg_bField,
                                 const TrackRef arg_Trk1, const TrackRef arg_Trk2);

    // Vertex fitting overloads
    static bool particlesToVtx(const vector<RefCountedKinematicParticle>& arg_FromParticles);
    static bool particlesToVtx(const vector<RefCountedKinematicParticle>& arg_FromParticles,
                               const string& arg_Message);
    static bool particlesToVtx(RefCountedKinematicTree& arg_VertexFitTree,
                               const vector<RefCountedKinematicParticle>& arg_Muons,
                               const string& arg_Message);
    static bool particlesToVtx(const vector<RefCountedKinematicParticle>& arg_FromParticles,
                               const double& arg_VtxProbCut);
    static bool particlesToVtx(const vector<RefCountedKinematicParticle>& arg_FromParticles,
                               const string& arg_Message,
                               const double& arg_VtxProbCut);
    static bool particlesToVtx(RefCountedKinematicTree& arg_VertexFitTree,
                               const vector<RefCountedKinematicParticle>& arg_Muons,
                               const string& arg_Message,
                               const double& arg_VtxProbCut);

    // Fit result extraction
    static bool extractFitRes(RefCountedKinematicTree& arg_VtxTree,
                              RefCountedKinematicParticle& res_Part,
                              RefCountedKinematicVertex& res_Vtx,
                              KinematicParameters& res_Param,
                              double& res_MassErr);
    static bool extractFitRes(RefCountedKinematicTree& arg_VtxTree,
                              RefCountedKinematicParticle& res_Part,
                              RefCountedKinematicVertex& res_Vtx,
                              double& res_MassErr);
    static bool extractFitRes(RefCountedKinematicTree& arg_VtxTree,
                              RefCountedKinematicVertex& res_Vtx,
                              double& res_VtxProb);

    static bool isOverlapPair(const muList_t& arg_MuonPair1, 
                              const muList_t& arg_MuonPair2);

    static double fitResEval(double arg_massDiff_1, double arg_massErr_1,
                             double arg_massDiff_2, double arg_massErr_2,
                             double arg_massDiff_3, double arg_massErr_3);
    
    void printKinematics(const RefCountedKinematicParticle& particle, const std::string& name);

    // ======================== Muon-track matching methods ========================
    MuonTrackMatchResult matchMuonToTrack(
        const pat::Muon& muon,
        const edm::Handle<edm::View<pat::PackedCandidate>>& packedHandle,
        const edm::Handle<VertexCollection>& recVtxs,
        const reco::Vertex& primaryV) const;
    PriCandidateDiagnostics evaluatePriCandidateDiagnostics(
        bool priFitValid,
        float priVtxProb,
        const std::vector<const pat::PackedCandidate*>& daughterPackedCands,
        const std::vector<const reco::Track*>& daughterTracks) const;
    PhiKaonDiagnostics buildPhiKaonDiagnostics(
        const pat::PackedCandidate& cand,
        const reco::Vertex& primaryV,
        const edm::Handle<reco::GenParticleCollection>& genParticles) const;
    PhiVertexDiagnostics buildPhiVertexDiagnostics(
        bool fitPass,
        const PhiKaonDiagnostics& kaon1,
        const PhiKaonDiagnostics& kaon2) const;
    void storePriDiagnostics(const PriCandidateDiagnostics& diagnostics);

    // Store resonance fit results into branches (reduces code duplication)
    void storeResonanceBranches(
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
        vector<float>* br_ptErr);

    // Store sentinel values for failed 3-body vertex fit
    void storeSentinelPri();

    // ======================== EDM tokens ========================
    edm::EDGetTokenT<L1GlobalTriggerReadoutRecord>    gtRecordToken_;
    edm::EDGetTokenT<BeamSpot>                        gtbeamspotToken_;
    edm::EDGetTokenT<VertexCollection>                gtprimaryVtxToken_;
    edm::EDGetTokenT<edm::View<pat::Muon>>            gtpatmuonToken_;
    edm::EDGetTokenT<edm::TriggerResults>             gttriggerToken_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackToken_;
    edm::EDGetTokenT<reco::GenParticleCollection>     genParticlesToken_;

    // ======================== Config parameters ========================
    
    // -- Essentials --
    InputTag hlTriggerResults_;
    InputTag inputGEN_;
    InputTag muonLabel_;
    InputTag trackLabel_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBuilderToken_;
    
    bool    doMC;
    bool    requireAcceptedCandidatesForMonteCarloTree_;
    bool    doJPsiMassCost;
    bool    Debug_;
    
    // -- Analysis mode: "JpsiJpsiPhi", "JpsiJpsiUps", "JpsiUpsPhi" --
    std::string analysisModeName_;
    AnalysisChannel analysisChannel_;
    
    // -- Muon selection (StringCutObjectSelector based) --
    std::string muonSelectionStr_;
    StringCutObjectSelector<pat::Muon> muonSelector_;
    
    // -- Track (kaon/pion) selection (StringCutObjectSelector based) --
    std::string trackSelectionStr_;
    StringCutObjectSelector<pat::PackedCandidate> trackSelector_;
    
    // -- Primary vertex quality cuts (configurable) --
    int    pvNdofMin_;
    double pvMaxAbsZ_;
    double pvMaxRho_;
    
    // -- Muon pair (onia) mass windows --
    double JpsiMassMin_, JpsiMassMax_;
    double UpsMassMin_,  UpsMassMax_;
    
    // -- Track pair (meson) mass windows --
    double PhiMassMin_,  PhiMassMax_;
    
    // -- Track kinematics --
    double trackPtMin_;
    double trackDRMax_;
    
    // -- Vertex probability cuts --
    double OniaDecayVtxProbCut_;
    double PriVtxProbCut_;
    
    // -- Per-resonance candidate pT and eta pre-cuts --
    double jpsiCandPtMin_,  jpsiCandEtaMax_;
    double upsCandPtMin_,   upsCandEtaMax_;
    double phiCandPtMin_,   phiCandEtaMax_;
    
    // -- Primary vertex selection mode ("firstVertex", "mostTracks", "highestSumPt2") --
    std::string pvSelectionMode_;
    
    // -- Minimum fromPV value for tracks (0=none, 1=PVLoose, 2=PVTight, 3=PVUsedInFit) --
    int minTrackFromPV_;
    
    // -- Minimum number of muons required --
    unsigned int minMuonCount_;
    
    // -- Muon-track matching method selection --
    std::string muTrkMatchMethod_;
    MuTrkMatchMethod muTrkMatchMode_;
    bool muTrkMatchDebug_;
    double muonPackedMatchVectorRelPMax_;
    double muonPackedMatchChi2Max_;
    double muonPackedMatchDzPvChi2Max_;
    double muonPackedMatchDzAssocChi2Max_;

    // -- Store all primary vertices --
    bool storeAllPVs_;
    bool storeMuonMomentumErrors_;
    bool storeMuonPVAssoc_;
    double recoGenMuonMatchChi2Max_;
    double recoGenKaonMatchChi2Max_;
    bool priRequireCommonAssocPV_;
    bool priRequireTrackPVCompatibility_;
    double priTrackDzPVMax_;
    double priTrackDxyPVMax_;

    // -- Final fitted mass window check --
    bool checkFinalMass_;

    // -- Trigger info --
    bool resolveAmbiguity_; 
    bool addXlessPrimaryVertex_;
    vector<string> TriggersForJpsi_;
    vector<string> FiltersForJpsi_;
    vector<string> TriggersForUpsilon_;
    vector<string> FiltersForUpsilon_;
    int JpsiMatchTrig[50], UpsilonMatchTrig[50];

    // -- HLT config --
    string proccessName_;
    HLTConfigProvider hltConfig_;

    // PDG 2024 constants
    static constexpr double myJpsiMass = 3.0969,   myJpsiMassErr = 0.00004;
    static constexpr double myUpsMass  = 9.4603,   myUpsMassErr  = 0.0003;
    static constexpr double myPhiMass  = 1.019455, myPhiMassErr  = 0.000020;
    static constexpr double myMuMass   = 0.1056583745;
    static constexpr double myMuMassErr = 0.0000000023;
    static constexpr double myKMass    = 0.493677;
    static constexpr double myKMassErr = 0.000015;

    // ======================== Intermediate storage (event-level) ========================
    // These are populated during analyze() and cleared at end
    reco::Vertex thePrimaryV_;
    reco::Vertex theBeamSpotV_;
    edm::Handle<edm::View<pat::Muon>> thePATMuonHandle_;
    edm::Handle<edm::View<pat::PackedCandidate>> theTrackHandle_;
    std::vector<edm::View<pat::PackedCandidate>::const_iterator> nonMuonTrack_;
    std::unordered_map<unsigned int, int> handleToNtupleIndex_;
    
    // Muon pair candidates
    std::vector<muList_t> muPairCand_Onia1_;  // J/psi (or 1st quarkonium)
    std::vector<muList_t> muPairCand_Onia2_;  // Upsilon (or 2nd quarkonium, if different)
    std::vector<std::pair<muList_t, muList_t>> muQuad_Onia_;
    
    // Track pair candidates
    using Kaon_t  = RefCountedKinematicParticle;
    using KList_t = std::pair<vector<Kaon_t>, vector<uint>>;
    std::vector<KList_t> KPairCand_Meson_;

    // ======================== TTree ========================
    TTree* X_One_Tree_;
    
    // -- Event info --
    unsigned int runNum, evtNum, lumiNum;
    unsigned int nGoodPrimVtx;
    
    vector<unsigned int>* trigRes;
    vector<std::string>*  trigNames;
    vector<unsigned int>* L1TT;
    vector<std::string>*  MatchJpsiTrigNames;
    vector<std::string>*  MatchUpsTrigNames;

    // -- Primary vertex --
    float priVtxX, priVtxY, priVtxZ;
    float priVtxXE, priVtxYE, priVtxZE;
    float priVtxChiNorm, priVtxChi, priVtxCL;
    vector<float>  *PriVtxXCorrX, *PriVtxXCorrY, *PriVtxXCorrZ;
    vector<double> *PriVtxXCorrEX, *PriVtxXCorrEY, *PriVtxXCorrEZ;
    vector<float>  *PriVtxXCorrC2, *PriVtxXCorrCL;

    // -- All primary vertices (not just selected) --
    unsigned int nRecVtx;
    vector<float> *RecVtx_x;
    vector<float> *RecVtx_y;
    vector<float> *RecVtx_z;
    vector<float> *RecVtx_xErr;
    vector<float> *RecVtx_yErr;
    vector<float> *RecVtx_zErr;
    vector<float> *RecVtx_chi2;
    vector<float> *RecVtx_ndof;
    vector<float> *RecVtx_vtxProb;
    vector<int>   *RecVtx_nTracks;

    // -- All muons: kinematics --
    unsigned int nMu;
    vector<float> *muPx, *muPy, *muPz;
    vector<float> *muD0, *muD0E, *muDz;
    vector<float> *muChi2, *muGlChi2, *mufHits;
    vector<bool>  *muFirstBarrel, *muFirstEndCap;
    vector<float> *muDzVtx, *muDxyVtx;
    vector<int>   *muNDF, *muGlNDF, *muPhits, *muShits, *muGlMuHits, *muType, *muQual;
    vector<int>   *muTrack;
    vector<float> *muCharge;
    vector<float> *muIsoratio;

    // -- All muons: selection results --
    vector<int> *muIsGoodLooseMuon, *muIsGoodLooseMuonNew;
    vector<int> *muIsGoodSoftMuonNewIlse, *muIsGoodSoftMuonNewIlseMod;
    vector<int> *muIsGlobalMuon;
    vector<int> *muIsGoodTightMuon, *muIsJpsiTrigMatch;
    vector<int> *muIsUpsTrigMatch, *munMatchedSeg;
    vector<int> *muIsJpsiFilterMatch, *muIsUpsFilterMatch;
    vector<int> *muIsPatLooseMuon, *muIsPatTightMuon, *muIsPatSoftMuon, *muIsPatMediumMuon;
    vector<int> *muFromPV, *muPVAssocQuality;

    // -- Muon momentum errors --
    vector<float> *muPxErr;
    vector<float> *muPyErr;
    vector<float> *muPzErr;
    vector<float> *muPtErr;

    // -- Muon PV association (surplus from sourceCandidatePtr) --
    vector<int>   *muVertexId;          // Which PV this muon is associated to (from sourceCandidatePtr)
    vector<float> *muDzAssocPV;         // dZ w.r.t. associated PV
    vector<float> *muDxyAssocPV;        // dxy w.r.t. associated PV
    vector<int>   *muFromPVAssocPV;     // fromPV from sourceCandidatePtr
    vector<int>   *muPdgId;             // PDG ID from sourceCandidatePtr
    vector<int>   *muPackedMatchIdx;    // Handle index in packedPFCandidates
    vector<int>   *muPackedMatchMethod; // MuTrkMatchMethod enum value
    vector<float> *muPackedMatchVectorRelP;
    vector<float> *muPackedMatchChi2;
    vector<float> *muPackedMatchDzPV;
    vector<float> *muPackedMatchDzAssocPV;
    vector<int>   *muGenMatchIdx;       // Index into MC_GenPart_* branches
    vector<int>   *muGenMatchSource;    // 0=none, 1=PAT ref, 2=chi2 fallback
    vector<float> *muGenMatchChi2;      // Best chi2 used for the assignment

    // -- Muon ID variables --
    vector<float> *muMVAMuonID, *musegmentCompatibility;
    vector<float> *mupulldXdZ_pos_noArb, *mupulldYdZ_pos_noArb;
    vector<float> *mupulldXdZ_pos_ArbDef, *mupulldYdZ_pos_ArbDef;
    vector<float> *mupulldXdZ_pos_ArbST, *mupulldYdZ_pos_ArbST;
    vector<float> *mupulldXdZ_pos_noArb_any, *mupulldYdZ_pos_noArb_any;

    // -- Muon-track matching debug (only if muTrkMatchDebug_) --
    vector<int>   *muMatch_nCandidates;
    vector<float> *muMatch_bestRelDiff;
    vector<float> *muMatch_bestDz;
    vector<int>   *muMatch_methodUsed;

    // -- Resonance candidate indices --
    vector<float> *Jpsi_1_mu_1_Idx, *Jpsi_1_mu_2_Idx;
    vector<float> *Jpsi_2_mu_1_Idx, *Jpsi_2_mu_2_Idx;
    vector<float> *Phi_K_1_Idx, *Phi_K_2_Idx;
    // [Comment] Storing raw dxy and dz values for muons and the fitted candidates is good redundancy.
    // -- Reconstructed resonances: mass, vertex, kinematics --
    vector<float> *Jpsi_1_mass, *Jpsi_1_massErr, *Jpsi_1_massDiff;
    vector<float> *Jpsi_2_mass, *Jpsi_2_massErr, *Jpsi_2_massDiff;
    vector<float> *Phi_mass, *Phi_massErr, *Phi_massDiff;

    vector<float> *Jpsi_1_ctau, *Jpsi_1_ctauErr, *Jpsi_1_Chi2, *Jpsi_1_ndof, *Jpsi_1_VtxProb;
    vector<float> *Jpsi_2_ctau, *Jpsi_2_ctauErr, *Jpsi_2_Chi2, *Jpsi_2_ndof, *Jpsi_2_VtxProb;
    vector<float> *Phi_ctau, *Phi_ctauErr, *Phi_Chi2, *Phi_ndof, *Phi_VtxProb;

    vector<float> *Jpsi_1_phi, *Jpsi_1_eta, *Jpsi_1_pt;
    vector<float> *Jpsi_2_phi, *Jpsi_2_eta, *Jpsi_2_pt;
    vector<float> *Phi_phi, *Phi_eta, *Phi_pt;

    vector<float> *Jpsi_1_px, *Jpsi_1_py, *Jpsi_1_pz;
    vector<float> *Jpsi_2_px, *Jpsi_2_py, *Jpsi_2_pz;
    vector<float> *Phi_px, *Phi_py, *Phi_pz;

    // -- Momentum uncertainties for resonances --
    vector<float> *Jpsi_1_pxErr, *Jpsi_1_pyErr, *Jpsi_1_pzErr, *Jpsi_1_ptErr;
    vector<float> *Jpsi_2_pxErr, *Jpsi_2_pyErr, *Jpsi_2_pzErr, *Jpsi_2_ptErr;
    vector<float> *Phi_pxErr, *Phi_pyErr, *Phi_pzErr, *Phi_ptErr;
    vector<int>   *Phi_fitPass, *Phi_commonAssocPVPass, *Phi_commonAssocPVIdx, *Phi_trackPVPass, *Phi_vertexCriteriaPass;
    vector<float> *Phi_maxAbsDzPV, *Phi_maxAbsDxyPV;
    
    // -- Primary (combined) vertex --
    vector<float> *Pri_mass, *Pri_massErr;
    vector<float> *Pri_ctau, *Pri_ctauErr, *Pri_Chi2, *Pri_ndof, *Pri_VtxProb;
    vector<float> *Pri_px, *Pri_py, *Pri_pz;
    vector<float> *Pri_phi, *Pri_eta, *Pri_pt;
    vector<float> *Pri_pxErr, *Pri_pyErr, *Pri_pzErr, *Pri_ptErr;
    vector<int>   *Pri_fitValid, *Pri_fitPass, *Pri_assocPVPass, *Pri_assocPVIdx, *Pri_trackPVPass, *Pri_passAny;
    vector<float> *Pri_maxAbsDzPV, *Pri_maxAbsDxyPV;

    // -- Kaon tracks from Phi (or other meson) decay --
    vector<float> *Phi_K_1_px, *Phi_K_1_py, *Phi_K_1_pz;
    vector<float> *Phi_K_2_px, *Phi_K_2_py, *Phi_K_2_pz;
    vector<float> *Phi_K_1_eta, *Phi_K_1_phi, *Phi_K_1_pt;
    vector<float> *Phi_K_2_eta, *Phi_K_2_phi, *Phi_K_2_pt;
    vector<float> *Phi_K_1_fromPV, *Phi_K_2_fromPV;
    vector<float> *Phi_K_1_pvAssocQuality, *Phi_K_2_pvAssocQuality;
    vector<int>   *Phi_K_1_vertexId, *Phi_K_2_vertexId;
    vector<int>   *Phi_K_1_hasAssocPV, *Phi_K_2_hasAssocPV;
    vector<int>   *Phi_K_1_passDzPV, *Phi_K_2_passDzPV;
    vector<int>   *Phi_K_1_passDxyPV, *Phi_K_2_passDxyPV;
    vector<int>   *Phi_K_1_passTrackPV, *Phi_K_2_passTrackPV;
    vector<float> *Phi_K_1_dzPV, *Phi_K_1_dxyPV, *Phi_K_1_dzAssocPV, *Phi_K_1_dxyAssocPV;
    vector<float> *Phi_K_2_dzPV, *Phi_K_2_dxyPV, *Phi_K_2_dzAssocPV, *Phi_K_2_dxyAssocPV;
    vector<int>   *Phi_K_1_genMatchIdx, *Phi_K_1_genMatchSource;
    vector<int>   *Phi_K_2_genMatchIdx, *Phi_K_2_genMatchSource;
    vector<float> *Phi_K_1_genMatchChi2, *Phi_K_2_genMatchChi2;

    // -- Upsilon candidate (for JpsiJpsiUps and JpsiUpsPhi modes) --
    vector<float> *Ups_mu_1_Idx, *Ups_mu_2_Idx;
    vector<float> *Ups_mass, *Ups_massErr, *Ups_massDiff;
    vector<float> *Ups_ctau, *Ups_ctauErr, *Ups_Chi2, *Ups_ndof, *Ups_VtxProb;
    vector<float> *Ups_px, *Ups_py, *Ups_pz;
    vector<float> *Ups_phi, *Ups_eta, *Ups_pt;
    vector<float> *Ups_pxErr, *Ups_pyErr, *Ups_pzErr, *Ups_ptErr;

    // ======================== MC gen-level branches ========================
    // Flat gen-particle lists (all J/psi, Upsilon, phi, mu, K in event)
    vector<int>   *MC_GenPart_pdgId;
    vector<int>   *MC_GenPart_status;
    vector<int>   *MC_GenPart_motherPdgId;
    vector<int>   *MC_GenPart_motherGenIdx;
    vector<int>   *MC_GenPart_handleIndex;
    vector<float> *MC_GenPart_px, *MC_GenPart_py, *MC_GenPart_pz, *MC_GenPart_mass;
    vector<float> *MC_GenPart_pt, *MC_GenPart_eta, *MC_GenPart_phi;
    
    // Legacy MC branches (backward-compatible with old ntuples)
    vector<float> *MC_X_px, *MC_X_py, *MC_X_pz, *MC_X_mass;
    vector<float> *MC_Dau_Jpsipx, *MC_Dau_Jpsipy, *MC_Dau_Jpsipz, *MC_Dau_Jpsimass;
    vector<float> *MC_Dau_psi2spx, *MC_Dau_psi2spy, *MC_Dau_psi2spz, *MC_Dau_psi2smass;
    vector<float> *MC_Granddau_mu1px, *MC_Granddau_mu1py, *MC_Granddau_mu1pz;
    vector<float> *MC_Granddau_mu2px, *MC_Granddau_mu2py, *MC_Granddau_mu2pz;
    vector<float> *MC_Granddau_Jpsipx, *MC_Granddau_Jpsipy, *MC_Granddau_Jpsipz, *MC_Granddau_Jpsimass;
    vector<float> *MC_Granddau_pi1px, *MC_Granddau_pi1py, *MC_Granddau_pi1pz;
    vector<float> *MC_Granddau_pi2px, *MC_Granddau_pi2py, *MC_Granddau_pi2pz;
    vector<float> *MC_Grandgranddau_mu3px, *MC_Grandgranddau_mu3py, *MC_Grandgranddau_mu3pz;
    vector<float> *MC_Grandgranddau_mu4px, *MC_Grandgranddau_mu4py, *MC_Grandgranddau_mu4pz;
    
    vector<int> *MC_X_chg;
    vector<int> *MC_Dau_JpsipdgId, *MC_Dau_psi2spdgId;
    vector<int> *MC_Granddau_mu1pdgId, *MC_Granddau_mu2pdgId;
    vector<int> *MC_Granddau_JpsipdgId;
    vector<int> *MC_Granddau_pi1pdgId, *MC_Granddau_pi2pdgId;
    vector<int> *MC_Grandgranddau_mu3pdgId, *MC_Grandgranddau_mu4pdgId;
    
    // RECO-to-GEN matching
    vector<float> *Match_mu1px, *Match_mu1py, *Match_mu1pz;
    vector<float> *Match_mu2px, *Match_mu2py, *Match_mu2pz;
    vector<float> *Match_mu3px, *Match_mu3py, *Match_mu3pz;
    vector<float> *Match_mu4px, *Match_mu4py, *Match_mu4pz;
    vector<float> *Match_pi1px, *Match_pi1py, *Match_pi1pz;
    vector<float> *Match_pi2px, *Match_pi2py, *Match_pi2pz;
};

#endif
