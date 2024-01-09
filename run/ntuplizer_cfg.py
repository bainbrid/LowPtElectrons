import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles','file:input.root')
options.setDefault('outputFile','output.root')
options.setDefault('maxEvents',-1)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.register('useAOD',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.register('addSkim',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.parseArguments()

process = cms.Process('TEST') 

# /store/data/Run2018B/ParkingBPH5/RAW/v1/000/317/626/00001/E0032C29-4C6D-E811-91CA-FA163E8EFCB8.root

# AOD: 9558 + 6230 + 1661 = 17449 events
# MINOAOD: 17449 events
default_files = [
    #'/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/792E843F-F82F-4642-97D7-A7625ABA83D2.root',
    #'/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/DB36F8A5-E77B-4547-9320-6DD1FA5BA58F.root',
    #'/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/F617C5AF-E3AB-754D-BFFC-C004AAB2887E.root',
    'file:/afs/cern.ch/user/b/bainbrid/eos/WORK/dev/revisit-seed/CMSSW_10_2_14/src/AODSIM_mc.root'
] if options.useAOD else [
    '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/5A729845-ECB4-8E47-AF73-52C79F53B0B2.root'
]

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(default_files),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
    )

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('output_filtered.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('ntuplizer_path'))
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
if options.useAOD is True : switchOnVIDElectronIdProducer(process,DataFormat.AOD)
else :                      switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
for idmod in [
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
    ] : 
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.ntuplizer_seq = cms.Sequence()

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff')
if options.useAOD is False : 
    process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
process.ntuplizer_seq *= process.lowPtGsfElectronID

process.load('LowPtElectrons.LowPtElectrons.IDNtuplizer_cfi')
process.ntuplizer_seq *= process.ntuplizer

process.ntuplizer_path = cms.Path(process.egmGsfElectronIDSequence*
                                  process.ntuplizer_seq)
process.output_path = cms.EndPath(process.output)
if options.addSkim is False : process.schedule = cms.Schedule(process.ntuplizer_path)
else : process.schedule = cms.Schedule(process.ntuplizer_path,process.output_path)
