import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles','file:input.root')
options.setDefault('outputFile','output_data.root')
options.setDefault('maxEvents',-1)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.register('useAOD',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.register('addSkim',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.register('simple',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.parseArguments()

if options.simple: print("Using IDNtuplizerSimple...!!!")

process = cms.Process('TEST') 

################################################################################
# Various sets of input files:
################################################################################

# Command line 
fileNames_cli = options.inputFiles

# MINIAOD:
fileNames_miniaod = [
    #'/store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/50002/A16E8E15-1DBF-5044-A16D-88FC78887F59.root',
    #'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root',
    '/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/130000/00141A1E-02F3-0B4F-9128-2CE371D069B3.root',
]
fileNames_miniaod_UL = [
    '/store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2530000/521F0FA1-D91E-ED4D-8920-5BB1E7CDDE38.root'
]

# AOD:
fileNames_aod = []

################################################################################
# Input source, number of events, GT, ...
################################################################################

file_maod = ["file:./root_files/example_files/miniaod.root"]
file_aod = ["file:./root_files/example_files/aod.root"]

fileNames_to_use = None
if options.useAOD: fileNames_to_use = fileNames_aod
else:              fileNames_to_use = fileNames_miniaod

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(fileNames_to_use),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v16', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37', '') # UL

# BParkingNano chain (electron stuff overridden below??????)
process.load('PhysicsTools.BParkingNano.nanoBPark_cff')

outputFile = options.outputFile
if options.simple: outputFile = outputFile.replace(".root","_simple.root")
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(outputFile)
)

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('output_filtered.root'),
    outputCommands = cms.untracked.vstring('keep *','drop recoTransientTracks_*_*_*'),
    #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('ntuplizer_path')),
)

################################################################################
# Energy regression
################################################################################

# !!! Requires soft link to local DB file, e.g.
# lowPtEleReg_2018_02062020_nv.db -> $CMSSW_BASE/src/PhysicsTools/BParkingNano/test/lowPtEleReg_2018_02062020_nv.0Cdb

## Edit the GlobalTag for the low-pT energy regression
process.GlobalTag.toGet = cms.VPSet(
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")))

from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XUL
from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XULLP

process.regressionForEle = cms.EDProducer(
    'ElectronRegresser',
    lowptSrc = cms.InputTag('slimmedLowPtElectrons'),
    pfSrc    = cms.InputTag('slimmedElectrons'),
    lowPtRegressionConfig = cms.PSet(
        modifierName = cms.string('EGRegressionModifierLPV1'),
        rhoTag = cms.string('fixedGridRhoFastjetAll'),
        useClosestToCentreSeedCrysDef = cms.bool(False),
        maxRawEnergyForLowPtEBSigma = cms.double(-1),
        maxRawEnergyForLowPtEESigma = cms.double(1200.),
        eleRegs = cms.PSet(
            ecalOnlyMean = cms.PSet(
                rangeMinLowEt = cms.double(0.2),
                rangeMaxLowEt = cms.double(2.0),
                rangeMinHighEt = cms.double(-1.),
                rangeMaxHighEt = cms.double(3.0),
                forceHighEnergyTrainingIfSaturated = cms.bool(True),
                lowEtHighEtBoundary = cms.double(20.),
                ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
                ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
                eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
                eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
            ),
            ecalOnlySigma = cms.PSet(
                rangeMinLowEt = cms.double(0.0002),
                rangeMaxLowEt = cms.double(0.5),
                rangeMinHighEt = cms.double(0.0002),
                rangeMaxHighEt = cms.double(0.5),
                forceHighEnergyTrainingIfSaturated = cms.bool(True),
                lowEtHighEtBoundary = cms.double(20.),
                ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
                ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
                eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
                eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
            ),
            epComb = cms.PSet(
                ecalTrkRegressionConfig = cms.PSet(
                    rangeMinLowEt = cms.double(0.2),
                    rangeMaxLowEt = cms.double(2.0),
                    rangeMinHighEt = cms.double(0.2),
                    rangeMaxHighEt = cms.double(2.0),
                    lowEtHighEtBoundary = cms.double(20.),
                    forceHighEnergyTrainingIfSaturated = cms.bool(False),
                    ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To20_mean'),
                    ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_20To50_mean'),
                    eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To20_mean'),
                    eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_20To50_mean'),
                ),
                ecalTrkRegressionUncertConfig = cms.PSet(
                    rangeMinLowEt = cms.double(0.0002),
                    rangeMaxLowEt = cms.double(0.5),
                    rangeMinHighEt = cms.double(0.0002),
                    rangeMaxHighEt = cms.double(0.5),
                    lowEtHighEtBoundary = cms.double(20.),
                    forceHighEnergyTrainingIfSaturated = cms.bool(False),
                    ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To20_sigma'),
                    ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_20To50_sigma'),
                    eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To20_sigma'),
                    eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_20To50_sigma'),
                ),
                maxEcalEnergyForComb=cms.double(200.),
                minEOverPForComb=cms.double(0.025),
                maxEPDiffInSigmaForComb=cms.double(15.),
                maxRelTrkMomErrForComb=cms.double(10.),
            )
        ),
        phoRegs = regressionModifier106XUL.phoRegs.clone()
    ),
    gsfRegressionConfig = cms.PSet(
        modifierName = cms.string('EGRegressionModifierV3'),
        rhoTag = cms.string('fixedGridRhoFastjetAll'),
        useClosestToCentreSeedCrysDef = cms.bool(False),
        maxRawEnergyForLowPtEBSigma = cms.double(-1),
        maxRawEnergyForLowPtEESigma = cms.double(1200.),
        eleRegs = cms.PSet(
            ecalOnlyMean = cms.PSet(
                rangeMinLowEt = cms.double(0.2),
                rangeMaxLowEt = cms.double(2.0),
                rangeMinHighEt = cms.double(-1.),
                rangeMaxHighEt = cms.double(3.0),
                forceHighEnergyTrainingIfSaturated = cms.bool(True),
                lowEtHighEtBoundary = cms.double(999999.),
                ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
                ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
                eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
                eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
            ),
            ecalOnlySigma = cms.PSet(
                rangeMinLowEt = cms.double(0.0002),
                rangeMaxLowEt = cms.double(0.5),
                rangeMinHighEt = cms.double(0.0002),
                rangeMaxHighEt = cms.double(0.5),
                forceHighEnergyTrainingIfSaturated = cms.bool(True),
                lowEtHighEtBoundary = cms.double(999999.),
                ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
                ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
                eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
                eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
            ),
            epComb = cms.PSet(
                ecalTrkRegressionConfig = cms.PSet(
                    rangeMinLowEt = cms.double(0.2),
                    rangeMaxLowEt = cms.double(2.0),
                    rangeMinHighEt = cms.double(0.2),
                    rangeMaxHighEt = cms.double(2.0),
                    lowEtHighEtBoundary = cms.double(999999.),
                    forceHighEnergyTrainingIfSaturated = cms.bool(False),
                    ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
                    ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
                    eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
                eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
                ),
                ecalTrkRegressionUncertConfig = cms.PSet(
                    rangeMinLowEt = cms.double(0.0002),
                    rangeMaxLowEt = cms.double(0.5),
                    rangeMinHighEt = cms.double(0.0002),
                    rangeMaxHighEt = cms.double(0.5),
                    lowEtHighEtBoundary = cms.double(999999.),
                    forceHighEnergyTrainingIfSaturated = cms.bool(False),
                    ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
                    ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
                    eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
                    eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
                ),
                maxEcalEnergyForComb=cms.double(200.),
                minEOverPForComb=cms.double(0.025),
                maxEPDiffInSigmaForComb=cms.double(15.),
                maxRelTrkMomErrForComb=cms.double(10.),
            )
        ),
        phoRegs = regressionModifier106XUL.phoRegs.clone()
    )
)

################################################################################
# PF MVA ID (incl. Otto's retraining)
################################################################################

# Otto's retraining:
# https://indico.cern.ch/event/732971/contributions/3022864/attachments/1658765/2656595/180530_egamma.pdf

# vvv BEGIN: (was in: RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_BParkRetrain_cff.py)
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_tools import *
from os import path
mvaTag = "BParkRetrain"
weightFileDir = "RecoEgamma/ElectronIdentification/data/LowPtElectrons"
mvaWeightFiles = cms.vstring(
    path.join(weightFileDir, "BParkRetrain_LowPt_unbiased.xml.gz"),
    path.join(weightFileDir, "BParkRetrain_HighPt_unbiased.xml.gz"),
)
categoryCuts = cms.vstring(
    "pt < 5.",
    "pt >= 5.",
)
mvaEleID_BParkRetrain_producer_config = cms.PSet(
    mvaName             = cms.string(mvaClassName),
    mvaTag              = cms.string(mvaTag),
    nCategories         = cms.int32(2),
    categoryCuts        = categoryCuts,
    weightFileNames     = mvaWeightFiles,
    variableDefinition  = cms.string(mvaVariablesFile)
)
# ^^^ END

mvaConfigsForEleProducer = cms.VPSet( )

#from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_BParkRetrain_cff \
#    import mvaEleID_BParkRetrain_producer_config
mvaConfigsForEleProducer.append( mvaEleID_BParkRetrain_producer_config )

from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff \
    import mvaEleID_Fall17_noIso_V2_producer_config
mvaConfigsForEleProducer.append( mvaEleID_Fall17_noIso_V2_producer_config )

process.electronMVAVariableHelper = cms.EDProducer(
    'GsfElectronMVAVariableHelper',
    src = cms.InputTag('gedGsfElectrons'),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    beamSpot         = cms.InputTag("offlineBeamSpot"),
    conversions      = cms.InputTag("allConversions"),
    srcMiniAOD              = cms.InputTag('slimmedElectrons'),#processName=cms.InputTag.skipCurrentProcess()),
    vertexCollectionMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpotMiniAOD         = cms.InputTag("offlineBeamSpot"),
    conversionsMiniAOD      = cms.InputTag("reducedEgamma:reducedConversions"),
)

process.electronMVAValueMapProducer = cms.EDProducer(
    'ElectronMVAValueMapProducer',
    src = cms.InputTag('gedGsfElectrons'),
    srcMiniAOD = cms.InputTag('slimmedElectrons'),#processName=cms.InputTag.skipCurrentProcess()),
    mvaConfigurations = mvaConfigsForEleProducer
)

process.egmGsfElectronIDs = cms.EDProducer(
    "VersionedGsfElectronIdProducer",
    physicsObjectSrc = cms.InputTag('gedGsfElectrons'),
    physicsObjectIDs = cms.VPSet( )
)

process.egmGsfElectronIDTask = cms.Task(
    process.electronMVAVariableHelper,
    process.electronMVAValueMapProducer,
    #process.egmGsfElectronIDs,
)

################################################################################
# Default model, 2019Aug07
################################################################################

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff')
if options.useAOD == False: 
    process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
else:
    process.lowPtGsfElectronID.electrons = 'lowPtGsfElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAllTmp'
process.lowPtGsfElectronID.ModelNames = [
    '2019Aug07',
]
process.lowPtGsfElectronID.ModelWeights = [
    'RecoEgamma/ElectronIdentification/data/LowPtElectrons/RunII_Autumn18_LowPtElectrons_mva_id.root',
]
process.lowPtGsfElectronID.ModelThresholds = cms.vdouble([
    -99.,
])

################################################################################
# ROME models, 2020Sep15: depth13 and depth15 (default)
################################################################################

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronIDExtra_cff')
if options.useAOD == True : 
    process.lowPtGsfElectronIDExtra.electrons = 'lowPtGsfElectrons'
    process.lowPtGsfElectronIDExtra.rho = 'fixedGridRhoFastjetAllTmp'
#else:
    # use defaults?
    #process.lowPtGsfElectronIDExtra.electrons = 'regressionForEle:regressedLowPtElectrons'
    #process.lowPtGsfElectronIDExtra.rho = 'fixedGridRhoFastjetAll'

process.lowPtGsfElectronIDExtra.ModelNames = [
    '2020Sept15',
    '2020Nov28',
    '2021May17',
]
process.lowPtGsfElectronIDExtra.ModelWeights = [
    'RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root', # ele_mva_value_depth10
    'RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Nov28.root',  # ele_mva_value_depth11
    'RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2021May17.root',  # ele_mva_value_depth13
]
process.lowPtGsfElectronIDExtra.ModelThresholds = cms.vdouble([
    -99.,
    -99.,
    -99.,
])

################################################################################
# Ntuplizer code
################################################################################

if options.simple:
    process.load('LowPtElectrons.LowPtElectrons.IDNtuplizerSimple_cfi')
else:
    process.load('LowPtElectrons.LowPtElectrons.IDNtuplizer_cfi')

# Running on data !!!
process.ntuplizer.isMC = cms.bool(False)
process.ntuplizer.verbose = 0

from_tracks = False
if from_tracks : # Evaluating models for BParking studies
    process.ntuplizer.tagMuonPtThreshold  = 5.
    process.ntuplizer.tagMuonEtaThreshold = 2.5
    process.ntuplizer.filterNtupleContent = False
    process.ntuplizer.prescale = 0.#-2.94 # Poisson mean number of fakes/event
else : # Skim to keep just low-pT electrons (no seeds, no PF, etc, ...) for Max Hart training
    process.ntuplizer.tagMuonPtThreshold  = 7.
    process.ntuplizer.tagMuonEtaThreshold = 1.5
    process.ntuplizer.filterNtupleContent = True
    process.ntuplizer.prescale = 100. # Poisson mean number of fakes/event

################################################################################
# BParking Analysis sequences to define data control regions
################################################################################

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#
#from PhysicsTools.BParkingNano.muonsBPark_cff import muonTrgSelector
#process.muonTrgSelector = muonTrgSelector.clone()
#
#from PhysicsTools.BParkingNano.electronsBPark_cff import electronsForAnalysis
#process.electronsForAnalysis = electronsForAnalysis.clone(
#    lowptSrc = cms.InputTag('regressionForEle:regressedLowPtElectrons'),
#    pfSrc = cms.InputTag('slimmedElectrons'), #@@ NOT regressionForEle:regressedElectrons for PF!
#    mvaId = cms.InputTag("lowPtGsfElectronIDExtra:depth15"),
#    pfmvaId = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues"),
#    drForCleaning_wrtTrgMuon = cms.double(-1.),
#    dzForCleaning_wrtTrgMuon = cms.double(-1.),
#    pf_ptMin = cms.double(0.5),
#    bdtMin = cms.double(-99.),
#    )
#
#from PhysicsTools.BParkingNano.BToKLL_cff import electronPairsForKee
#process.electronPairsForKee = electronPairsForKee.clone(
#    lep1Selection = cms.string('pt > 0.5'),
#    lep2Selection = cms.string('pt > 0.5'),
#    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1.'\
#                                 ' && mass() > 0.'\
#                                 ' && mass() < 5.'\
#                                 ' && charge() == 0'\
#                                 ' && userFloat("lep_deltaR") > 0.03'\
#                                 ' && userInt("nlowpt")<3'
#                             ),
#    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
#)

from PhysicsTools.BParkingNano.nanoBPark_cff import *
process = nanoAOD_customizeMuonTriggerBPark(process)
process = nanoAOD_customizeElectronFilteredBPark(process)
process = nanoAOD_customizeTrackFilteredBPark(process)
process = nanoAOD_customizeBToKLL(process)
#process = nanoAOD_customizeBToKstarEE(process)
#process = nanoAOD_customizeBToKstarMuMu(process)
process = nanoAOD_customizeTriggerBitsBPark(process)

process.nanoAOD_KMuMu_step = cms.Path(process.nanoSequence + process.nanoBKMuMuSequence + CountBToKmumu )

################################################################################
# Paths, Sequences, etc
################################################################################

process.egamma_path = cms.Path(process.egmGsfElectronIDTask)
process.path = cms.Path(process.lowPtGsfElectronID
                        #+process.regressionForEle
                        +process.lowPtGsfElectronIDExtra
                        #+process.muonTrgSelector
                        #+process.electronsForAnalysis
                        #+process.electronPairsForKee
                    )
process.ntuplizer_path = cms.Path(process.ntuplizer)

process.path_tmp = cms.Path(process.nanoSequence
                            +process.nanoBKMuMuSequence
                            +process.CountBToKmumu 
                            +process.electronMVAVariableHelper
                            +process.electronMVAValueMapProducer
                            #+process.egmGsfElectronIDs
                            +process.lowPtGsfElectronID
                            +process.lowPtGsfElectronIDExtra
                            +process.ntuplizer
                        )

if options.addSkim is False : 
    process.schedule = cms.Schedule(
        #process.path_tmp,
        #process.nanoAOD_KMuMu_step,
        process.egamma_path,process.path,process.ntuplizer_path
    )
else : 
    process.output_path = cms.EndPath(process.output)
    process.schedule = cms.Schedule(process.egamma_path,process.path,process.ntuplizer_path,process.output_path)
