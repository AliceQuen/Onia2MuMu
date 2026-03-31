###############################################################################
# ConfFile_cfg.py - Default CMSSW config for MultiLepPAT EDAnalyzer
#
# Analysis Mode:
#   "JpsiJpsiPhi" - J/psi + J/psi + phi (default)
#   "JpsiJpsiUps" - J/psi + J/psi + Upsilon
#   "JpsiUpsPhi"  - J/psi + Upsilon + phi
#
# All selection cuts externalized via StringCutObjectSelector syntax:
#   MuonSelection  - String cut on pat::Muon
#   TrackSelection - String cut on pat::PackedCandidate
#
# Usage:
#   cmsRun ConfFile_cfg.py
#   cmsRun  ConfFile_cfg.py \
#           inputFiles_load=myInputs.list \
#           outputFile=output.root \
#           maxEvents=10000 \
#           analysisMode=JpsiJpsiUps \
#           runOnMC=False \
#           era=Run2023C \
#           duplicateCheckMode=noDuplicateCheck
###############################################################################

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# --- Include input arguments ---

ivars = VarParsing.VarParsing('analysis')

# --- Register additional arguments ---

ivars.register('analysisMode',
    default='JpsiJpsiPhi',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info='Analysis mode: JpsiJpsiPhi (default), JpsiJpsiUps, JpsiUpsPhi'
)
ivars.register('runOnMC',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info='Whether to run on MC (default: False)'
)
ivars.register('era',
    default='Run2023C',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info='Data-taking era (e.g. Run2023C, Run2023BPix). Required for GT selection.'
)
ivars.register('eventRange',
    default='',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info='Event range to process, in the format "run:event" or "run:event-run:event" (e.g. "369943:103642411")'
)
ivars.register('reportEvery',
    default=1,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info='Framework report frequency for processed events (default: 1)'
)
ivars.register('duplicateCheckMode',
    default='',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info='PoolSource duplicate check mode override, e.g. noDuplicateCheck'
)
ivars.register('requireAcceptedCandidatesForMonteCarloTree',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info='For MC, keep tree entries only when at least one candidate is stored (default: False)'
)

# --- Default values ---

ivars.inputFiles = ()
ivars.outputFile = 'mymultilep.root'

ivars.analysisMode = 'JpsiJpsiPhi'
ivars.runOnMC = False
ivars.era = 'Run2023C'

# ------ Dictionary to Handle Global Tag selection based on era and MC/data ------
globalTagDict = {
    # Source: PdmV twiki 'CMS/PdmVRun3Analysis'
    'data': {
        # Data from 2022, eras C-G, E-G corresponding to the 'postEE' period.
        "Run2022C": '124X_dataRun3_PromptAnalysis_v1',
        "Run2022D": '124X_dataRun3_PromptAnalysis_v1',
        "Run2022E": '124X_dataRun3_Prompt_v10',
        "Run2022F": '124X_dataRun3_PromptAnalysis_v2',
        "Run2022G": '124X_dataRun3_PromptAnalysis_v2',
        # Data from 2023, eras C-D, D corresponding to the 'postBPix' period.
        'Run2023C': '130X_dataRun3_PromptAnalysis_v1',
        'Run2023D': '130X_dataRun3_PromptAnalysis_v1',
        # Data from 2024, eras C-I. (Consistent run throughout the year. One GT.)
        'Run2024C': '150X_dataRun3_v2',
        'Run2024D': '150X_dataRun3_v2',
        'Run2024E': '150X_dataRun3_v2',
        'Run2024F': '150X_dataRun3_v2',
        'Run2024G': '150X_dataRun3_v2',
        'Run2024H': '150X_dataRun3_v2',
        'Run2024I': '150X_dataRun3_v2',
        # Data from 2025, eras C-G. (Consistent run throughout the year. One GT.)
        'Run2025C': '150X_dataRun3_v2',
        'Run2025D': '150X_dataRun3_v2',
        'Run2025E': '150X_dataRun3_v2',
        'Run2025F': '150X_dataRun3_v2',
        'Run2025G': '150X_dataRun3_v2',
    },
    'MC': {
        # MC for 2022, with 'postEE' period.
        "Run2022": '130X_mcRun3_2022_realistic_v5',
        "Run2022EE": '130X_mcRun3_2022_realistic_postEE_v6',
        # MC for 2023, with 'postBPix' period.
        'Run2023': '130X_mcRun3_2023_realistic_v14',
        'Run2023BPix': '130X_mcRun3_2023_realistic_postBPix_v2',
        # MC for 2024, consistent run throughout the year. One GT.
        'Run2024': '150X_mcRun3_2024_realistic_v2',
        # MC for 2025, consistent run throughout the year. One GT.
        'Run2025': '150X_mcRun3_2024_realistic_v2',
    }
}

# --- Parse and check input arguments ---

ivars.parseArguments()

if ivars.runOnMC and ivars.era == 'Run2023C':
    ivars.era = 'Run2022'

if ivars.analysisMode not in ['JpsiJpsiPhi', 'JpsiJpsiUps', 'JpsiUpsPhi']:
    raise ValueError(f"Invalid analysis mode: {ivars.analysisMode}. "
                     f"Valid options are: JpsiJpsiPhi, JpsiJpsiUps, JpsiUpsPhi")

if ivars.runOnMC and ivars.era not in globalTagDict['MC']:
    raise ValueError(f"Invalid era '{ivars.era}' for MC. Available options: {list(globalTagDict['MC'].keys())}")

if not ivars.runOnMC and ivars.era not in globalTagDict['data']:
    raise ValueError(f"Invalid era '{ivars.era}' for data. Available options: {list(globalTagDict['data'].keys())}")

# --- Begin process configuration ---

process = cms.Process("mkcands")

# --- Message logger ---
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.suppressInfo = cms.untracked.vstring("mkcands")
process.MessageLogger.suppressWarning = cms.untracked.vstring("mkcands")
process.MessageLogger.cerr.FwkReport.reportEvery = ivars.reportEvery

# --- Geometry / B-field / TransientTrack ---
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

# --- Global Tag ---
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag


# Configure Global Tag based on era and MC/data.
if ivars.runOnMC:
    myGlobalTag = globalTagDict['MC'][ivars.era]
else:
    myGlobalTag = globalTagDict['data'][ivars.era]
# Pass the selected Global Tag to the process configuration.
process.GlobalTag = GlobalTag(process.GlobalTag, myGlobalTag, '')

# --- Input ---
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(ivars.maxEvents))
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(ivars.inputFiles),
)
if ivars.eventRange != '':
    process.source.eventRange = cms.untracked.VEventRange(ivars.eventRange)
if ivars.duplicateCheckMode != '':
    process.source.duplicateCheckMode = cms.untracked.string(ivars.duplicateCheckMode)

# --- Primary vertex filter ---
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

# --- MultiLepPAT analyzer ---
process.mkcands = cms.EDAnalyzer('MultiLepPAT',
    HLTriggerResults = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    inputGEN = cms.untracked.InputTag("prunedGenParticles"),
    MuonLabel = cms.untracked.InputTag("slimmedMuons"),
    TrackLabel = cms.untracked.InputTag("packedPFCandidates"),

    # ====== Analysis mode ======
    # Options: "JpsiJpsiPhi", "JpsiJpsiUps", "JpsiUpsPhi"
    AnalysisMode = cms.untracked.string(ivars.analysisMode),

    # ====== StringCutObjectSelector cuts (runtime-configurable) ======
    MuonSelection  = cms.untracked.string("pt > 2.5 && abs(eta) < 2.4"),
    TrackSelection = cms.untracked.string("pt > 1.0 && abs(eta) < 2.5 && numberOfHits > 4"),

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
    PhiMassMin = cms.untracked.double(0.95),
    PhiMassMax = cms.untracked.double(1.07),

    # ====== Track kinematics ======
    TrackPtMin = cms.untracked.double(1.0),
    TrackDRMax = cms.untracked.double(0.7),

    # ====== Vertex probability cuts ======
    OniaDecayVtxProbCut = cms.untracked.double(0.0),
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

    # ====== PV selection mode ======
    # Options: "firstVertex" (default), "mostTracks", "highestSumPt2"
    PVSelectionMode = cms.untracked.string("firstVertex"),

    # ====== Track fromPV requirement ======
    # 0=none, 1=PVLoose, 2=PVTight, 3=PVUsedInFit
    MinTrackFromPV = cms.untracked.int32(1),

    # ====== Minimum muon count ======
    MinMuonCount = cms.untracked.uint32(4),

    # ====== MC toggle ======
    # DoMonteCarloTree enables MC branches; the retention switch below controls
    # whether MC events without kept candidates are still written to the tree.
    DoMonteCarloTree = cms.untracked.bool(ivars.runOnMC),
    RequireAcceptedCandidatesForMonteCarloTree = cms.untracked.bool(ivars.requireAcceptedCandidatesForMonteCarloTree),
    DoJPsiMassConstraint = cms.untracked.bool(True),
    Debug_Output = cms.untracked.bool(False),

    # ====== Muon matching ======
    MuonPackedMatchVectorRelPMax = cms.untracked.double(0.01),
    MuonPackedMatchChi2Max = cms.untracked.double(25.0),
    MuonPackedMatchDzPvChi2Max = cms.untracked.double(25.0),
    MuonPackedMatchDzAssocChi2Max = cms.untracked.double(25.0),
    RecoGenMuonMatchChi2Max = cms.untracked.double(25.0),
    RecoGenKaonMatchChi2Max = cms.untracked.double(25.0),

    # ====== Muon-track matching method ======
    # Methods: "sourceCandidatePtr" (default), "vector", "chi2", "dzAssoc", "dzPv"
    MuTrkMatchMethod = cms.untracked.string("sourceCandidatePtr"),
    MuTrkMatchDebug = cms.untracked.bool(True),

    # ====== Store all primary vertices and muon quantities ======
    StoreAllPVs = cms.untracked.bool(True),
    StoreMuonMomentumErrors = cms.untracked.bool(True),
    StoreMuonPVAssoc = cms.untracked.bool(True),

    # ====== Final fitted mass window check ======
    CheckFinalMass = cms.untracked.bool(True),

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

# --- Output ---
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ivars.outputFile),
)

process.p = cms.Path(process.mkcands)
process.schedule = cms.Schedule(process.p)
