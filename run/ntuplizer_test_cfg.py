import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles','file:input.root')
options.setDefault('outputFile','output.root')
options.setDefault('maxEvents',-1)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.register('useAOD',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.register('addSkim',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.parseArguments()

TEST = not options.useAOD
AOD = options.useAOD
DERIVED = False

process = cms.Process('TEST') 

#from files import files
base1='file:/afs/cern.ch/user/b/bainbrid/work/public/6-ntuplizer/CMSSW_10_2_14/src/1-miniaod-from-crab/step3_inAODSIM_{:.0f}_numEvent200.root'
base2='file:/afs/cern.ch/user/b/bainbrid/work/public/6-ntuplizer/CMSSW_10_2_14/src/1-miniaod-from-crab/lowpteleid/190502_224210/crab_lowpteleid/results/step3_inAODSIM_{:.0f}.root'
process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring(options.inputFiles), # Command line
    #fileNames = cms.untracked.vstring(files[:1]), # Test of recent production
    fileNames = cms.untracked.vstring(
        [
            #'file:./5955A84F-6A21-9C49-9948-0697A60262A0.root' # earlier MINIAOD test file
            #'file:./MINIAOD_from_AOD_TEST.root' # 9558 events
            #'file:./output_filtered_MINIAOD.root'

            '/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/280000/C842E5BC-A265-424C-921B-8BF1913D73C4.root'

#            '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/5A729845-ECB4-8E47-AF73-52C79F53B0B2.root' # 17449 events

#            '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/5EA6D799-A04A-1D42-A76B-6A564B28B7AF.root',
#            '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/5A729845-ECB4-8E47-AF73-52C79F53B0B2.root',
#            '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/5270FDD3-DA34-C444-B7DD-B8A84A257BCE.root',
#            '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/4D1D7852-44FB-0E4F-882D-1F42B6B81722.root',
#            '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/3200A640-1F5C-B44E-89C1-F6CFD4FC3BD9.root',
#            '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/26FF8646-13A8-6648-9717-340AE968A408.root',
#            '/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/0E7ECF47-95EE-6449-9CB2-E06BAF618B3C.root',
            
        ]
        if TEST else \
        [
            #'file:./output_filtered_AOD.root',
            #'/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/4EFA4099-1C34-EB45-ABFF-39AA91656DB3.root' # 9558 events

            # 17449 events from AOD:
            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/792E843F-F82F-4642-97D7-A7625ABA83D2.root', # 9558 +
            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/DB36F8A5-E77B-4547-9320-6DD1FA5BA58F.root', # 6230 +
            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/F617C5AF-E3AB-754D-BFFC-C004AAB2887E.root', #1661 = 17449 events

#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/F52F5093-2728-B040-89B2-14D17E748D33.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/CB392D4D-FCE6-7A4A-8FE2-A7F4703A6F30.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/41DEF1B3-643F-164B-9AEE-8E8AB320297F.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/0D29D10E-3D38-BF4C-9CDA-98F19271C810.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/405295DF-F22D-D24B-8A89-0400F0CAA195.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/4EFA4099-1C34-EB45-ABFF-39AA91656DB3.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/098E956C-DF60-CD47-9A46-D61FC51A09ED.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/47AA1872-4C04-DE47-9B34-4B9A5E17A4E7.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/DB36F8A5-E77B-4547-9320-6DD1FA5BA58F.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/AAD087E2-399E-414A-A223-55227DF09B16.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/792E843F-F82F-4642-97D7-A7625ABA83D2.root',
#            '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/F617C5AF-E3AB-754D-BFFC-C004AAB2887E.root',
#           '/store/mc/RunIIAutumn18RECOBParking/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/E08C5049-036B-D34C-8F90-E71EA84C4D47.root',

        ]
        if AOD else \
        ['file:/afs/cern.ch/user/b/bainbrid/work/public/7-slc7/CMSSW_10_2_15/src/2-ntuples-from-crab/MINIAOD_from_AOD_test.root',] \
        if DERIVED else \
        ['/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/5955A84F-6A21-9C49-9948-0697A60262A0.root']
    ),

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

#####################################
# PF MVA ID (incl. Otto's retraining)
#####################################

mvaConfigsForEleProducer = cms.VPSet( )
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_BParkRetrain_cff \
    import mvaEleID_BParkRetrain_producer_config
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
    srcMiniAOD              = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),
    vertexCollectionMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpotMiniAOD         = cms.InputTag("offlineBeamSpot"),
    conversionsMiniAOD      = cms.InputTag("reducedEgamma:reducedConversions"),
)

process.electronMVAValueMapProducer = cms.EDProducer(
    'ElectronMVAValueMapProducer',
    src = cms.InputTag('gedGsfElectrons'),
    srcMiniAOD = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),
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
    process.egmGsfElectronIDs,
)

process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDTask)

#################
# 2019Aug07 model
#################

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff')
if AOD is False : 
    process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'

#################
# 2019Aug07 model
#################

# !!! Requires soft link to local DB file, e.g.
# lowPtEleReg_2018_02062020_nv.db -> $CMSSW_BASE/src/PhysicsTools/BParkingNano/test/lowPtEleReg_2018_02062020_nv.db

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

#############
# ROME models
#############

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronIDExtra_cff')
if AOD is False : 
    #process.lowPtGsfElectronIDExtra.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronIDExtra.electrons = 'regressionForEle:regressedLowPtElectrons'
    process.lowPtGsfElectronIDExtra.rho = 'fixedGridRhoFastjetAll'
    process.lowPtGsfElectronIDExtra.ModelNames = [
        'depth13',
        'depth15',
    ]
    process.lowPtGsfElectronIDExtra.ModelWeights = [
        'RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15_depth13.root',
        'RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root',
    ]
    process.lowPtGsfElectronIDExtra.ModelThresholds = cms.vdouble([
        -99.,
        -99.
    ])

########################
# .xml.gz or .root files
########################

use_root_files = True
if AOD is False and use_root_files :
    process.lowPtGsfElectronID.ModelWeights = [x.replace(".xml.gz",".root") 
                                               for x in process.lowPtGsfElectronID.ModelWeights]
    process.lowPtGsfElectronIDExtra.ModelWeights = [x.replace(".xml.gz",".root") 
                                                    for x in process.lowPtGsfElectronIDExtra.ModelWeights]
print(process.lowPtGsfElectronID.ModelWeights) 
print(process.lowPtGsfElectronIDExtra.ModelWeights) 

################
# Ntuplizer code
################

process.load('LowPtElectrons.LowPtElectrons.IDNtuplizer_cfi')

################
# Sequences, etc
################

process.ntuplizer_seq = cms.Sequence(process.lowPtGsfElectronID *
                                     process.regressionForEle *
                                     process.lowPtGsfElectronIDExtra *
                                     process.ntuplizer)
process.ntuplizer_path = cms.Path(process.egmGsfElectronIDSequence*
                                  process.ntuplizer_seq)
process.output_path = cms.EndPath(process.output)
if options.addSkim is False : process.schedule = cms.Schedule(process.ntuplizer_path)
else : process.schedule = cms.Schedule(process.ntuplizer_path,process.output_path)
