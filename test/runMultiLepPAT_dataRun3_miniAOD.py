import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')

ivars.inputFiles=(
#'file:BPHSkim_UL.root',
#'/store/user/zhenhu/MuOnia/BPHSkim-v3-Run2018D-12Nov2019_UL2018-v1/210321_010747/0000/BPHSkim_UL_556.root',
#'/store/data/Run2023B/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/366/406/00000/7a0c3c42-f92e-4d1d-952c-e13f90b0a8e3.root',
#'/store/data/Run2023B/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/366/403/00000/5be75e01-f0e6-4690-ba57-da041d5437a0.root'
'/store/data/Run2023B/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/367/079/00000/e20ba4ff-963d-4d1e-8933-a5e9ca8f9559.root'
#'/store/data/Run2023D/ParkingDoubleMuonLowMass1/MINIAOD/PromptReco-v1/000/370/093/00000/9b491ed3-c0e4-4870-be7e-f93d79a5ae01.root'
#'/store/data/Run2023D/ParkingDoubleMuonLowMass3/MINIAOD/22Sep2023_v2-v1/60000/02862e37-eb36-4b7e-82c8-648274e03648.root'
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
#HIFormat = True
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
#process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
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
    #process.GlobalTag.globaltag = cms.string( 'START53_V7F::All' )  #Summer12_DR53X
    #process.GlobalTag.globaltag = cms.string( 'STARTHI53_V26::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START52_V5::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START52_V7::All' )
    #process.GlobalTag.globaltag = cms.string( 'START53_V17::All' )
#    process.GlobalTag.globaltag = cms.string( 'START53_V27::All' ) ##pPb
    process.GlobalTag.globaltag = cms.string( 'MCRUN2_74_V9::All' ) 
else:
#    process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v12', '')
    #process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v4', '')
    process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '') 
### PoolSource will be ignored when running crab
process.source = cms.Source("PoolSource",
    skipEvents=cms.untracked.uint32(0),
	fileNames = cms.untracked.vstring(ivars.inputFiles),
	#eventsToProcess = cms.untracked.VEventRange("367079:791559619-367079:MAX")
	#eventsToProcess = cms.untracked.VEventRange("367079:970546777")
)

### Set basic filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
	vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
	#vertexCollection = cms.InputTag('offlineSlimmedPrimaryVerticesWithBS'),
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
#process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)
#process.filter = cms.Sequence(process.noscraping)
#process.filter = cms.Sequence(process.PAcollisionEventSelection)

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
#from PhysicsTools.PatAlgos.patEventContent_cff import *
if HIFormat:
	process.muonMatch.matched = cms.InputTag("hiGenParticles")
	process.genParticlePlusGEANT.genParticles = cms.InputTag("hiGenParticles")

##Using GEN plus SIM list for matching
if UseGenPlusSim:
	process.muonMatch.matched = cms.InputTag("genParticlePlusGEANT")

#process.allLayer1Jets.addJetCorrFactors = False
from PhysicsTools.PatAlgos.tools.trackTools import *######MODIFIED
#from Bfinder.tempTools.trackTools import *######MODIFIED
#process.load( 'PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cff' )
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
    #process.patTrackCandsMCMatch.matched = cms.InputTag("hiGenParticles")
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
#from Bfinder.tempTools.tempTools import *######MODIFIED
#from UserCode.MultiLepPAT.tempTools import *######MODIFIED
#removeSpecificPATObjects(process, ['Photons', 'Electrons', 'Taus', 'Jets', 'METs', 'Muons'])######MODIFIED
#removeAllPATObjectsBut(process, ['Muons'])
#removeSpecificPATObjects(process, ['Jets'])



#if not runOnMC :
#	removeMCMatching(process, ['All'] )
	#removeMCMatching(process, ['Photons', 'Electrons', 'Muons', 'Taus', 'Jets', 'METs', 'PFAll', 'PFElectrons', 'PFTaus', 'PFMuons','OOTPhotons'] )


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
    #changeRecoMuonInput(process, "mergedMuons")#Add calo muon to the collection
    #process.patMuons.muonSource = cms.InputTag("mergedMuons")#Need to use the same collection as they are internally entengled
    #process.patMuons.embedCaloMETMuonCorrs = cms.bool(False)
    #process.patMuons.embedTcMETMuonCorrs   = cms.bool(False)

    #Or we change the muonMatch source of our patMuonsWithoutTrigger
    process.patMuonsWithoutTrigger.muonSource = cms.InputTag("mergedMuons")
    process.patMuonsWithoutTriggerMatch = PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi.muonMatch.clone( src = cms.InputTag("mergedMuons"))
    if runOnMC:
        process.patMuonsWithTriggerSequence.replace(process.patMuonsWithoutTrigger, process.patMuonsWithoutTriggerMatch + process.patMuonsWithoutTrigger)
        process.patMuonsWithoutTrigger.genParticleMatch = 'patMuonsWithoutTriggerMatch'
    process.patDefaultSequence = cms.Sequence(process.mergedMuons*process.patDefaultSequence)

### Set MultiLepPAT option
process.mkcands = cms.EDAnalyzer('MultiLepPAT',
        HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
        # HLTTriggerSummary = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
        inputGEN  = cms.untracked.InputTag("genParticles"),
        #VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
        VtxSample   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
	#PFCandidates = cms.untracked.InputTag('packedPFCandidates'),
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

TriggersForJpsi = cms.untracked.vstring("HLT_Dimuon0_Jpsi_Muon_v18","HLT_Dimuon0_Jpsi_Muon_v17","HLT_Dimuon0_Jpsi_Muon_v16","HLT_Dimuon0_Jpsi_Muon_v15","HLT_Dimuon0_Jpsi_Muon_v14","HLT_Dimuon0_Jpsi_Muon_v13","HLT_Dimuon0_Jpsi_Muon_v12","HLT_Dimuon0_Jpsi_Muon_v11","HLT_Dimuon0_Jpsi_Muon_v10","HLT_Dimuon0_Jpsi_Muon_v9","HLT_Dimuon0_Jpsi_Muon_v8","HLT_Dimuon0_Jpsi_Muon_v7","HLT_Dimuon0_Jpsi_Muon_v6","HLT_Dimuon0_Jpsi_Muon_v5","HLT_Dimuon0_Jpsi_Muon_v4","HLT_Dimuon0_Jpsi_Muon_v3","HLT_Dimuon0_Jpsi_Muon_v2","HLT_Dimuon0_Jpsi_Muon_v1"),
 FiltersForJpsi = cms.untracked.vstring("hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon"),

#TriggersForUpsilon = cms.untracked.vstring("HLT_Dimuon0_Upsilon_Muon_v18","HLT_Dimuon0_Upsilon_Muon_v17","HLT_Dimuon0_Upsilon_Muon_v16","HLT_Dimuon0_Upsilon_Muon_v15","HLT_Dimuon0_Upsilon_Muon_v14","HLT_Dimuon0_Upsilon_Muon_v13","HLT_Dimuon0_Upsilon_Muon_v12","HLT_Dimuon0_Upsilon_Muon_v11","HLT_Dimuon0_Upsilon_Muon_v10","HLT_Dimuon0_Upsilon_Muon_v9","HLT_Dimuon0_Upsilon_Muon_v8","HLT_Dimuon0_Upsilon_Muon_v7","HLT_Dimuon0_Upsilon_Muon_v6","HLT_Dimuon0_Upsilon_Muon_v5","HLT_Dimuon0_Upsilon_Muon_v4","HLT_Dimuon0_Upsilon_Muon_v3","HLT_Dimuon0_Upsilon_Muon_v2","HLT_Dimuon0_Upsilon_Muon_v1"),
#FiltersForUpsilon = cms.untracked.vstring("hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon"),

TriggersForUpsilon = cms.untracked.vstring("HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v1","HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v2","HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v3","HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v3"),
FiltersForUpsilon = cms.untracked.vstring("hltVertexmumuFilterUpsilonMuon53p52OpenMuon","hltVertexmumuFilterUpsilonMuon53p52OpenMuon","hltVertexmumuFilterUpsilonMuon53p52OpenMuon","hltVertexmumuFilterUpsilonMuon53p52OpenMuon"),


#        TriggersForJpsi = cms.untracked.vstring("HLT_Dimuon0_Jpsi_Muon_v18","HLT_Dimuon0_Jpsi_Muon_v17","HLT_Dimuon0_Jpsi_Muon_v16","HLT_Dimuon0_J$
#        FiltersForJpsi = cms.untracked.vstring("hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVertexmumuFilterJpsiMuon","hltVerte$

#        TriggersForUpsilon = cms.untracked.vstring("HLT_Dimuon0_Upsilon_Muon_v18","HLT_Dimuon0_Upsilon_Muon_v17","HLT_Dimuon0_Upsilon_Muon_v16","H$
#        FiltersForUpsilon = cms.untracked.vstring("hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuon","hltVertexmumuFilterUpsilonMuo$
 
        Chi2NDF_Track =  cms.untracked.double(15.0)
)

# process.mkcands = cms.EDAnalyzer('Bfinder',
# 	Bchannel 		= cms.vint32(
# 		1,#RECONSTRUCTION: J/psi + K
# 		1,#RECONSTRUCTION: J/psi + Pi
# 		1,#RECONSTRUCTION: J/psi + Ks 
# 		1,#RECONSTRUCTION: J/psi + K* (K+, Pi-)
# 		1,#RECONSTRUCTION: J/psi + K* (K-, Pi+)
# 		1,#RECONSTRUCTION: J/psi + phi
# 		1,#RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
# 		1,#RECONSTRUCTION: J/psi + lambda (p+, pi-) 
# 		1,),#RECONSTRUCTION: J/psi + lambda (p-, pi+) 
#     MuonTriggerMatchingPath = cms.vstring("HLT_Dimuon*", "HLT_DoubleMu*"),
# 	HLTLabel        = cms.InputTag('TriggerResults::HLT'),
#     GenLabel        = cms.InputTag('genParticles'),
# 	MuonLabel       = cms.InputTag('selectedPatMuons'),         #selectedPatMuons
# 	TrackLabel      = cms.InputTag('selectedPatTrackCands'),    #selectedPat
#     PUInfoLabel     = cms.InputTag("addPileupInfo"),
#     BSLabel     = cms.InputTag("offlineBeamSpot"),
#     PVLabel     = cms.InputTag("offlinePrimaryVerticesWithBS"),
#     tkPtCut = cms.double(0.4),
#     jpsiPtCut = cms.double(0.0),
#     bPtCut = cms.double(0.0),
#     RunOnMC = cms.bool(False),
#     doTkPreCut = cms.bool(True),
#     doMuPreCut = cms.bool(True)
# )


if HIFormat:
	process.mkcands.GenLabel = cms.InputTag('hiGenParticles')
if UseGenPlusSim:
	process.mkcands.GenLabel = cms.InputTag('genParticlePlusGEANT')
if UsepatMuonsWithTrigger:
	process.mkcands.MuonLabel = cms.InputTag('patMuonsWithTrigger')	

#### SetUp HLT info
##process.load('Bfinder.HiHLTAlgos.hltanalysis_cff')
#process.load('Bfinder.EventAnalysis.hltanalysis_cff')
#process.hltanalysis.dummyBranches = cms.untracked.vstring()
##if HIFormat:
#	#process.hltanalysis.mctruth = cms.InputTag("hiGenParticles")
#	#process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
#	#process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
#	#process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
#process.hltAna = cms.Path(process.filter*process.hltanalysis)


### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)

if runOnMC and UseGenPlusSim:
	process.patDefaultSequence *= process.genParticlePlusGEANT
if UsepatMuonsWithTrigger:
	process.patDefaultSequence *= process.patMuonsWithTriggerSequence

process.p = cms.Path(	
#    process.filter*
#    process.patDefaultSequence*
    process.mkcands
)
#process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule(
	process.p
#	,process.hltAna
#    ,process.e
)
