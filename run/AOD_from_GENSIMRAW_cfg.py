# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --filein file:input.root --fileout file:AOD_from_GEN-SIM-RAW.root --mc --eventcontent AODSIM --runUnscheduled --datatier AODSIM --conditions 102X_upgrade2018_realistic_v15 --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI --nThreads 1 --geometry DB:Extended --era Run2_2018,bParking --python_filename AOD_from_GEN-SIM-RAW_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 10
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles',[
    '/store/mc/RunIIAutumn18DR/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/GEN-SIM-RAW/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/AC0EFF4C-C75F-274E-A561-A19D20DB2668.root'

    '/store/mc/RunIIAutumn18DR/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/GEN-SIM-RAW/PUPoissonAve20_102X_upgrade2018_realistic_v15_ext1-v1/90000/FF6D06A9-FA1A-5A49-B448-7888E1B98299.root', # B --> K,Jpsi(-->ee)
    '/store/mc/RunIIAutumn18DR/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/GEN-SIM-RAW/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/FF53BD1A-5F6D-764D-ADCB-2909D3C0657D.root', # B --> K,ee
    '/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/GEN-SIM-RAW/PUPoissonAve20_102X_upgrade2018_realistic_v15_ext1-v1/10000/8737CFD9-0E66-554C-BCEA-2D48B10E6C23.root', # B --> K*,Jpsi(-->ee)
    '/store/mc/RunIIAutumn18DR/BdToKstar_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/GEN-SIM-RAW/PUPoissonAve20_102X_upgrade2018_realistic_v15_ext1-v1/90000/FF9D2AF2-88B6-CF4A-BFAA-88EBE6231394.root', # B --> K*,ee
][0]
)
options.setDefault('outputFile','AOD_from_GEN-SIM-RAW.root')
options.setDefault('maxEvents',10)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.parseArguments()

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2018,eras.bParking)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.AODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.eventinterpretaion_step,process.endjob_step,process.AODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(1)
process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
