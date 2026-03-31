###############################################################################
# runMultiLepPAT_MCRun3_miniAOD_Run2022.py
# MC-specific test config for MultiLepPAT refactored analyzer
#
# Runs on DPS J/psi+J/psi+phi MINIAOD MC samples.
# Includes PAT sequence + muon trigger matching.
#
# Usage:
#   cmsRun runMultiLepPAT_MCRun3_miniAOD_Run2022.py \
#       inputFiles=file:/path/to/MINIAOD.root \
#       outputFile=output.root \
#       maxEvents=10
###############################################################################

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('analysis')
ivars.inputFiles = ()
ivars.outputFile = 'mymultilep_MC_DPS1.root'
ivars.parseArguments()

### Configuration switches
AddCaloMuon = False
runOnMC = True
HIFormat = False
UseGenPlusSim = False
UsepatMuonsWithTrigger = False

process = cms.Process("mkcands")

# --- Message logger ---
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.suppressInfo = cms.untracked.vstring("mkcands")
process.MessageLogger.suppressWarning = cms.untracked.vstring("mkcands")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# --- Geometry / B-field / TransientTrack ---
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

# --- Output module (not used but kept for compatibility) ---
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring('drop *'),
)

# --- MaxEvents ---
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(ivars.maxEvents))

# --- Global Tag ---
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnMC:
    process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v5', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_PromptAnalysis_v1', '')

# --- Source ---
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(ivars.inputFiles),
)

# --- Vertex filter ---
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2),
)
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25),
)
process.filter = cms.Sequence(process.primaryVertexFilter + process.noscraping)

# --- PAT + muon trigger matching ---
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)

# Gen-level SIM particles
process.genParticlePlusGEANT = cms.EDProducer("GenPlusSimParticleProducer",
    src = cms.InputTag("g4SimHits"),
    setStatus = cms.int32(8),
    filter = cms.vstring("pt > 0.0"),
    genParticles = cms.InputTag("genParticles"),
)

# PAT track candidates
process.load("PhysicsTools.PatAlgos.patSequences_cff")
if HIFormat:
    process.muonMatch.matched = cms.InputTag("hiGenParticles")
    process.genParticlePlusGEANT.genParticles = cms.InputTag("hiGenParticles")
if UseGenPlusSim:
    process.muonMatch.matched = cms.InputTag("genParticlePlusGEANT")

from PhysicsTools.PatAlgos.tools.trackTools import *
if runOnMC:
    makeTrackCandidates(process,
        label='TrackCands',
        tracks=cms.InputTag('generalTracks'),
        particleType='pi+',
        preselection='pt > 0.3',
        selection='pt > 0.3',
        isolation={},
        isoDeposits=[],
        mcAs='muon',
    )
    process.patTrackCandsMCMatch.resolveByMatchQuality = cms.bool(True)
    process.patTrackCandsMCMatch.resolveAmbiguities = cms.bool(True)
    process.patTrackCandsMCMatch.checkCharge = cms.bool(True)
    process.patTrackCandsMCMatch.maxDPtRel = cms.double(0.5)
    process.patTrackCandsMCMatch.maxDeltaR = cms.double(0.7)
    process.patTrackCandsMCMatch.mcPdgId = cms.vint32(111, 211, 311, 321)
    process.patTrackCandsMCMatch.mcStatus = cms.vint32(1)
    l1cands = getattr(process, 'patTrackCands')
    l1cands.addGenMatch = True
else:
    makeTrackCandidates(process,
        label='TrackCands',
        tracks=cms.InputTag('generalTracks'),
        particleType='pi+',
        preselection='pt > 0.3',
        selection='pt > 0.3',
        isolation={},
        isoDeposits=[],
        mcAs=None,
    )
    l1cands = getattr(process, 'patTrackCands')
    l1cands.addGenMatch = False

from PhysicsTools.PatAlgos.tools.coreTools import *

# Muon trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
if runOnMC:
    addMCinfo(process)
    process.muonMatch.resolveByMatchQuality = True
changeTriggerProcessName(process, "HLT")
switchOffAmbiguityResolution(process)
process.muonL1Info.maxDeltaR = 0.3
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR = 0.3
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0

# --- MultiLepPAT analyzer (refactored) ---
process.mkcands = cms.EDAnalyzer('MultiLepPAT',
    HLTriggerResults = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    inputGEN = cms.untracked.InputTag("prunedGenParticles"),
    MuonLabel = cms.untracked.InputTag("slimmedMuons"),
    TrackLabel = cms.untracked.InputTag("packedPFCandidates"),

    # ====== Analysis mode ======
    AnalysisMode = cms.untracked.string("JpsiJpsiPhi"),

    # ====== StringCutObjectSelector-based cuts ======
    MuonSelection  = cms.untracked.string("pt > 2.5 && abs(eta) < 2.4"),
    TrackSelection = cms.untracked.string("pt > 2.0 && abs(eta) < 2.5 && numberOfHits > 4"),

    # ====== Primary vertex cuts ======
    PVNdofMin = cms.untracked.int32(5),
    PVMaxAbsZ = cms.untracked.double(24.0),
    PVMaxRho  = cms.untracked.double(2.0),

    # ====== Onia mass windows (GeV) ======
    JpsiMassMin = cms.untracked.double(1.0),
    JpsiMassMax = cms.untracked.double(4.0),
    UpsMassMin  = cms.untracked.double(8.0),
    UpsMassMax  = cms.untracked.double(12.0),

    # ====== Meson mass window (GeV) ======
    PhiMassMin = cms.untracked.double(0.8),
    PhiMassMax = cms.untracked.double(1.2),

    # ====== Track kinematics ======
    TrackPtMin = cms.untracked.double(2.0),
    TrackDRMax = cms.untracked.double(0.7),

    # ====== Vertex probability cuts ======
    OniaDecayVtxProbCut = cms.untracked.double(0.01),
    PriVtxProbCut       = cms.untracked.double(0.0),
    PriRequireCommonAssocPV = cms.untracked.bool(True),
    PriRequireTrackPVCompatibility = cms.untracked.bool(True),
    PriTrackDzPVMax = cms.untracked.double(2.0),
    PriTrackDxyPVMax = cms.untracked.double(0.1),

    # ====== Per-resonance candidate pT/eta pre-cuts ======
    JpsiCandPtMin  = cms.untracked.double(0.0),
    JpsiCandEtaMax = cms.untracked.double(999.0),
    UpsCandPtMin   = cms.untracked.double(0.0),
    UpsCandEtaMax  = cms.untracked.double(999.0),
    PhiCandPtMin   = cms.untracked.double(0.0),
    PhiCandEtaMax  = cms.untracked.double(999.0),

    # ====== Minimum muon count ======
    MinMuonCount = cms.untracked.uint32(4),

    # ====== MC-specific ======
    # DoMonteCarloTree enables MC branches; the retention switch below controls
    # whether MC events without kept candidates are still written to the tree.
    DoJPsiMassConstraint = cms.untracked.bool(True),
    DoMonteCarloTree = cms.untracked.bool(True),
    RequireAcceptedCandidatesForMonteCarloTree = cms.untracked.bool(False),
    Debug_Output = cms.untracked.bool(False),

    # ====== Muon matching ======
    MuTrkMatchMethod = cms.untracked.string("sourceCandidatePtr"),
    MuonPackedMatchVectorRelPMax = cms.untracked.double(0.01),
    MuonPackedMatchChi2Max = cms.untracked.double(25.0),
    MuonPackedMatchDzPvChi2Max = cms.untracked.double(25.0),
    MuonPackedMatchDzAssocChi2Max = cms.untracked.double(25.0),
    RecoGenMuonMatchChi2Max = cms.untracked.double(25.0),
    RecoGenKaonMatchChi2Max = cms.untracked.double(25.0),

    # ====== Vertex ambiguity ======
    resolvePileUpAmbiguity   = cms.untracked.bool(True),
    addXlessPrimaryVertex    = cms.untracked.bool(True),

    # ====== Trigger paths ======
    TriggersForJpsi = cms.untracked.vstring(
        "HLT_Dimuon0_Jpsi3p5_Muon2_v",
        "HLT_DoubleMu4_3_LowMass_v",
    ),
    FiltersForJpsi = cms.untracked.vstring(
        "hltVertexmumuFilterJpsiMuon3p5",
        "hltDisplacedmumuFilterDoubleMu43LowMass",
    ),
    TriggersForUpsilon = cms.untracked.vstring(
        "HLT_Trimuon5_3p5_2_Upsilon_Muon_v",
    ),
    FiltersForUpsilon = cms.untracked.vstring(
        "hltVertexmumuFilterUpsilonMuon",
    ),
)

if HIFormat:
    process.mkcands.inputGEN = cms.untracked.InputTag('hiGenParticles')
if UseGenPlusSim:
    process.mkcands.inputGEN = cms.untracked.InputTag('genParticlePlusGEANT')
if UsepatMuonsWithTrigger:
    process.mkcands.MuonLabel = cms.untracked.InputTag('patMuonsWithTrigger')

# --- Output ---
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ivars.outputFile),
)

if runOnMC and UseGenPlusSim:
    process.patDefaultSequence *= process.genParticlePlusGEANT
if UsepatMuonsWithTrigger:
    process.patDefaultSequence *= process.patMuonsWithTriggerSequence

process.p = cms.Path(process.mkcands)
process.schedule = cms.Schedule(process.p)
