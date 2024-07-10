from WMCore.Configuration import Configuration
import CRABClient
from CRABClient.UserUtilities import config
from datetime import date
import time

dataset = ["mc","data"][0]

conf = {
    "mc":{
        "psetName":"ntuplizer_test_cfg.py",
        "inputDataset":"/BuToKee_MufilterPt6_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "inputBlocks":[""],
    },
    "data":{
        "psetName":"ntuplizer_data_cfg.py",
        "inputDataset":"/ParkingBPH1/Run2018D-05May2019promptD-v1/MINIAOD",
        "inputBlocks":["/ParkingBPH1/Run2018D-05May2019promptD-v1/MINIAOD#06bb9204-7074-4799-bcb5-0d31f495eaea"],
    },
}.get(dataset)

config = Configuration()

folder = 'crab'
timestamp = date.today().strftime('_%d%m%y')+time.strftime('_%H%M%S')

config.section_('General')
config.General.requestName = folder+timestamp
config.General.workArea = folder
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName=conf['psetName']
config.JobType.maxMemoryMB=2500
config.JobType.numCores = 1

config.section_('Data')
config.Data.inputDBS='global'
config.Data.inputDataset=conf['inputDataset']
if dataset=="data": config.Data.inputBlocks=conf['inputBlocks']
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/bainbrid/work/ntuples/'+folder # maps onto /eos/user/b/bainbrid/
config.Data.publication = False

config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ["../../../PhysicsTools/BParkingNano/test/lowPtEleReg_2018_02062020_nv.db"]

config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX'
