import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')

ivars.inputFiles=('/store/data/Run2018D/Charmonium/MINIAOD/15Feb2022_UL2018-v1/2820000/0086F269-BEC0-0549-9790-AF83C021A86D.root'
)

ivars.outputFile='mymultilep.root'
# get and parse the command line arguments
ivars.parseArguments()

process = cms.Process("mkcands")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.suppressInfo = cms.untracked.vstring("mkcands")
process.MessageLogger.suppressWarning = cms.untracked.vstring("mkcands")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

### Set global tag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37', '') 
### PoolSource will be ignored when running crab
process.source = cms.Source("PoolSource",
    skipEvents=cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(ivars.inputFiles),
)

### Set MultiLepPAT option
process.mkcands = cms.EDAnalyzer('MultiLepPAT',
        TriggersForJpsi = cms.untracked.vstring("HLT_Dimuon0_Jpsi3p5_Muon2"),
        FiltersForJpsi = cms.untracked.vstring("hltVertexmumuFilterJpsiMuon3p5"),
)

### Set output
process.TFileService = cms.Service("TFileService", fileName = cms.string(ivars.outputFile))

process.p = cms.Path(process.mkcands)
process.schedule = cms.Schedule(process.p)
