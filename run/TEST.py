import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles','file:./input.root')
options.setDefault('maxEvents',-1)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.parseArguments()

process = cms.Process('TEST')

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(options.skipEvents),
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '') # 106X_mc2017_realistic_v6

process.load('LowPtElectrons.LowPtElectrons.PATLowPtElectrons_UnitTest_cfi')
process.testPATLowPtElectrons.verbose = True
process.testPATLowPtElectrons_step = cms.Path(process.testPATLowPtElectrons)
process.schedule = cms.Schedule(process.testPATLowPtElectrons_step)
