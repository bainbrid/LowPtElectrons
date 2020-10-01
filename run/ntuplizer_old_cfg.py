import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles','file:input.root')
options.setDefault('outputFile','output_old.root')
options.setDefault('maxEvents',-1)
options.parseArguments()

process = cms.Process('TEST')

from files import files
process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring(options.inputFiles),
    fileNames = cms.untracked.vstring(files[:2]),
    secondaryFileNames = cms.untracked.vstring()
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
    )

#process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.SimL1EmulatorDM_cff')
#process.load('Configuration.StandardSequences.DigiToRawDM_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.RecoSim_cff')
#process.load('CommonTools.ParticleFlow.EITopPAG_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    )

process.load('LowPtElectrons.LowPtElectrons.TrackerElectronsFeatures_cfi')
checkFromB = cms.bool(True),
drMax = cms.double(0.02),
fakesMultiplier = cms.double(6.),
rho = cms.InputTag('fixedGridRhoFastjetAll'),
beamspot = cms.InputTag("offlineBeamSpot"),
prunedGenParticles = cms.InputTag("prunedGenParticles"),
packedGenParticles = cms.InputTag("packedGenParticles"),
gsfTracks = cms.InputTag("lowPtGsfEleGsfTracks"),
electrons = cms.InputTag("slimmedLowPtElectrons"),
MVASeedUnbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
MVASeedPtbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased"),
MVAIDLowPt = cms.InputTag('lowPtGsfElectronID'),

process.ntuplizer_step = cms.Path(process.features)
process.schedule = cms.Schedule(process.ntuplizer_step)
