from CRABClient.UserUtilities import config, getUsernameFromSiteDB

index = 3
folder = 'lowpteleid_{:.0f}'.format(index)

config = config()

config.General.requestName = folder
config.General.workArea = folder
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step3_RECOSIM_cfg.py'
config.JobType.maxMemoryMB=3500
config.JobType.numCores=1

datasets = [
    '/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18DR-PUPoissonAve20_Bparking_102X_upgrade2018_realistic_v15-v1/GEN-SIM-DIGI-RAW',
    '/BsToPhiJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18DR-PUPoissonAve20_Bparking_102X_upgrade2018_realistic_v15-v1/GEN-SIM-DIGI-RAW',
    '/BdToKstar_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18DR-PUPoissonAve20_Bparking_102X_upgrade2018_realistic_v15-v1/GEN-SIM-DIGI-RAW',
    '/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18DR-PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/GEN-SIM-RAW',
    ]

config.Data.inputDataset = datasets[index]
config.Data.inputDBS = 'global'

if (False) : # test
    config.Data.splitting = 'EventAwareLumiBased'
    config.Data.unitsPerJob = 1
    config.Data.totalUnits = 1
else :
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1
    config.Data.totalUnits = -1

config.Data.outLFNDirBase = '/store/user/bainbrid/'+folder # '/eos/user/b/bainbrid/'+folder
config.Data.publication = False

#config.Site.whitelist =['T2_CH_CERNBOX']
config.Site.storageSite = 'T2_CH_CERNBOX'
#config.Site.storageSite = 'T2_UK_London_IC'
