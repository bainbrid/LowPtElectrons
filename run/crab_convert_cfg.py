from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()
folder = 'lowpteleid'

config.General.requestName = folder
config.General.workArea = folder
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'MINIAOD_from_AOD_cfg.py'
config.JobType.maxMemoryMB=3000
config.JobType.numCores = 1

config.Data.inputDBS = 'global'
config.Data.inputDataset  ='/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18RECOBParking-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/AODSIM'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/bainbrid/'+folder
config.Data.publication = True
config.JobType.allowUndistributedCMSSW = True # e.g. use CMSSW_10_2_15_patch2 on SLC7

config.Site.storageSite = 'T2_UK_London_IC'
