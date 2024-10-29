import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')

ivars.inputFiles=(
'/store/data/Run2023B/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/367/079/00000/e20ba4ff-963d-4d1e-8933-a5e9ca8f9559.root'
)

ivars.outputFile='mymultilep.root'
# get and parse the command line arguments
ivars.parseArguments()

### Add Calo muons
#AddCaloMuon = True
AddCaloMuon = False

### Run on MC?
#runOnMC = True
runOnMC = False

### Switching between "hiGenParticles"(pPb MC) and "genParticles" (pp MC)
HIFormat = False

### Include SIM tracks for matching?
UseGenPlusSim = False

### Using pat muon with trigger or not
UsepatMuonsWithTrigger = False

process = cms.Process("mkcands")
process.load("FWCore.MessageService.MessageLogger_cfi")
#added by yik
process.MessageLogger.suppressInfo = cms.untracked.vstring( "mkcands" )
process.MessageLogger.suppressWarning = cms.untracked.vstring( "mkcands" )
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#added by yik

### Set TransientTrackBuilder 
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
### Set Geometry/GlobalTag/BField
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")


### output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *',
    )
)


### Set maxEvents  -1 means to run all events in root file
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

### Set global tag
if runOnMC:
    process.GlobalTag.globaltag = cms.string( 'MCRUN2_74_V9::All' ) 
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '') 
### PoolSource will be ignored when running crab
process.source = cms.Source("PoolSource",
    skipEvents=cms.untracked.uint32(0),
	fileNames = cms.untracked.vstring(ivars.inputFiles),
        eventsToProcess = cms.untracked.VEventRange("367079:199863171")
)

### Set basic filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
	vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
	minimumNDOF = cms.uint32(4) ,
	maxAbsZ = cms.double(24),	
	maxd0 = cms.double(2)	
)

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
	debugOn = cms.untracked.bool(False),
	numtrack = cms.untracked.uint32(10),
	thresh = cms.untracked.double(0.25)
)

#added by yik
# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)
#added by yik


# Common offline event selection

process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)

##Producing Gen list with SIM particles
process.genParticlePlusGEANT = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(8),             # set status = 8 for GEANT GPs
        filter        = cms.vstring("pt > 0.0"),  # just for testing (optional)
	genParticles   = cms.InputTag("genParticles") # original genParticle list
)

### Setup Pat
### Ref: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatching
process.load("PhysicsTools.PatAlgos.patSequences_cff")
### keep only Pat:: part 
if HIFormat:
	process.muonMatch.matched = cms.InputTag("hiGenParticles")
	process.genParticlePlusGEANT.genParticles = cms.InputTag("hiGenParticles")

##Using GEN plus SIM list for matching
if UseGenPlusSim:
	process.muonMatch.matched = cms.InputTag("genParticlePlusGEANT")

from PhysicsTools.PatAlgos.tools.trackTools import *
if runOnMC:
    makeTrackCandidates(process,              # patAODTrackCands
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        tracks=cms.InputTag('generalTracks'), # input track collection
    	particleType='pi+',                   # particle type (for assigning a mass)
        preselection='pt > 0.3',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
        selection='pt > 0.3',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
    	isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
       	isoDeposits=[],
        mcAs='muon'                           # Replicate MC match as the one used for Muons
    );                                        # you can specify more than one collection for this
    ### MC+mcAs+Match/pat_label options
    process.patTrackCandsMCMatch.resolveByMatchQuality = cms.bool(True)
    process.patTrackCandsMCMatch.resolveAmbiguities = cms.bool(True)
    process.patTrackCandsMCMatch.checkCharge = cms.bool(True)
    process.patTrackCandsMCMatch.maxDPtRel = cms.double(0.5)
    process.patTrackCandsMCMatch.maxDeltaR = cms.double(0.7)
    process.patTrackCandsMCMatch.mcPdgId = cms.vint32(111, 211, 311, 321)
    process.patTrackCandsMCMatch.mcStatus = cms.vint32(1)
    l1cands = getattr(process, 'patTrackCands')
    l1cands.addGenMatch = True

else :
    makeTrackCandidates(process,              # patAODTrackCands
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        tracks=cms.InputTag('generalTracks'), # input track collection
        particleType='pi+',                   # particle type (for assigning a mass)
        preselection='pt > 0.3',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
        selection='pt > 0.3',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
        isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
        isoDeposits=[],
        mcAs=None                             # Replicate MC match as the one used for Muons
    );                                        # you can specify more than one collection for this
    l1cands = getattr(process, 'patTrackCands')
    l1cands.addGenMatch = False
from PhysicsTools.PatAlgos.tools.coreTools import *


process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
if runOnMC:
	addMCinfo(process)
	process.muonMatch.resolveByMatchQuality = True
changeTriggerProcessName(process, "HLT")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
###Criterias from Hyunchul's 
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

# Merge muons, calomuons in a single collection for T&P
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    muons     = cms.InputTag("muons"),
    mergeCaloMuons = cms.bool(True),  ### NEEDED TO RUN ON AOD
    caloMuons = cms.InputTag("calomuons"),
    minCaloCompatibility = cms.double(0.6),
    mergeTracks = cms.bool(False),
    tracks = cms.InputTag("generalTracks"),
)
if AddCaloMuon:
    process.patMuonsWithoutTrigger.muonSource = cms.InputTag("mergedMuons")
    process.patMuonsWithoutTriggerMatch = PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi.muonMatch.clone( src = cms.InputTag("mergedMuons"))
    if runOnMC:
        process.patMuonsWithTriggerSequence.replace(process.patMuonsWithoutTrigger, process.patMuonsWithoutTriggerMatch + process.patMuonsWithoutTrigger)
        process.patMuonsWithoutTrigger.genParticleMatch = 'patMuonsWithoutTriggerMatch'
    process.patDefaultSequence = cms.Sequence(process.mergedMuons*process.patDefaultSequence)

### Set MultiLepPAT option
process.mkcands = cms.EDAnalyzer('MultiLepPAT',
        HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
        inputGEN  = cms.untracked.InputTag("genParticles"),
        VtxSample   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
        DoJPsiMassConstraint = cms.untracked.bool(True),  #False;
        DoMonteCarloTree = cms.untracked.bool(False),
        MonteCarloParticleId = cms.untracked.int32(20443),
        trackQualities = cms.untracked.vstring('loose','tight','highPurity'),
        MinNumMuPixHits = cms.untracked.int32(1),
        MinNumMuSiHits = cms.untracked.int32(3),
        MaxMuNormChi2 = cms.untracked.double(99999),  #bascially remove this cut by change 15 to 99999
        MaxMuD0 = cms.untracked.double(10.0),
        MaxJPsiMass = cms.untracked.double(3.4),
        MinJPsiMass = cms.untracked.double(2.7),
        MinNumTrSiHits = cms.untracked.int32(4),
        MinMuPt = cms.untracked.double(1.95),  # changed to 0.500 from 0.300 by yik
        JPsiKKKMaxDR = cms.untracked.double(1.5),
        XCandPiPiMaxDR = cms.untracked.double(1.5),
        UseXDr = cms.untracked.bool(False),
        JPsiKKKMaxMass = cms.untracked.double(5.6),
        JPsiKKKMinMass = cms.untracked.double(5.0),
        resolvePileUpAmbiguity = cms.untracked.bool(True),
        addXlessPrimaryVertex = cms.untracked.bool(True),
        Debug_Output = cms.untracked.bool(False),

        TriggersForJpsi = cms.untracked.vstring("HLT_Dimuon0_Jpsi3p5_Muon2_v"),
#        FiltersForJpsi = cms.untracked.vstring("hltTriggerType", "hltL1TripleMu5SQ3SQ0OQDoubleMu53SQOSMassMax9"),
        FiltersForJpsi = cms.untracked.vstring("hltVertexmumuFilterJpsiMuon3p5"),
 
        Chi2NDF_Track =  cms.untracked.double(15.0)
)

if HIFormat:
	process.mkcands.GenLabel = cms.InputTag('hiGenParticles')
if UseGenPlusSim:
	process.mkcands.GenLabel = cms.InputTag('genParticlePlusGEANT')
if UsepatMuonsWithTrigger:
	process.mkcands.MuonLabel = cms.InputTag('patMuonsWithTrigger')	
### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)

if runOnMC and UseGenPlusSim:
	process.patDefaultSequence *= process.genParticlePlusGEANT
if UsepatMuonsWithTrigger:
	process.patDefaultSequence *= process.patMuonsWithTriggerSequence

process.p = cms.Path(	
    process.mkcands
)
process.schedule = cms.Schedule(
	process.p
)
