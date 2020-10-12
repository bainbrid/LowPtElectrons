from CRABClient.UserUtilities import config#, getUsernameFromSiteDB

#import CRABClient
#from WMCore.Configuration import Configuration
#config = Configuration()

####################
# DEFINE DATA SET

data_sets = {
    "aod":{
        # Non-resonant
        'nonres_vsmall':'/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18RECOBParking-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/AODSIM', # 26.0GB, 91109 events, 13 files (T2_CH_CERN)
        'nonres_med':'/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18RECOBParking-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/AODSIM', # 2.6TB, 9072646 events, 959 files (requested...)
        'nonres_large':'/BuToKee_MufilterPt6_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18RECOBParking-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v1/AODSIM', # 9.1TB, 31772522 events, 3384 files (requested...)
        # Resonant (J/psi)
        'res_small':'/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18RECOBParking-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/AODSIM', # 602.7GB, 2122456 events, 255 files (requested @ T2_US_Vanderbilt)
        'res_med':'/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18RECOBParking-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/AODSIM', # 3.0TB, 10578394 events, 1118 files (requested...)
    },
    "miniaod":{
        # Non-resonant
        'nonres_vsmall':'/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM', # 6.3GB, 91109 events, 7 files (on multiple sites)
        'nonres_med':'/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', #  619.9GB, 9072646 events, 222 files (on multiple sites)
        'nonres_large':'/BuToKee_MufilterPt6_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v1/MINIAODSIM', # : 2.2TB, 31585558 events, 786 files (on multiple sites)
        # Resonant (J/psi)
        'res_small':'/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM', # 145.3GB, 2122456 events, 74 files (on multiple sites)
        'res_med':'/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM', # 696.6GB, 10193124 events, 253 files (on multiple sites)
    },
    "derived":'/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/bainbrid-crab_lowpteleid-0bd58594e6ade05f64e0c3a8301c3139/USER'
}

data_tier = ['aod','miniaod','derived'][1]
data_set = {'aod':    data_sets.get(data_tier,'').get('',''), # usuall 'res_small'
            'miniaod':data_sets.get(data_tier,'').get('nonres_large',''),
            'derived':data_sets.get(data_tier,''),
}.get(data_tier,'')

print('data_tier:',data_tier)
print('data_set:',data_set)
if data_set == '' : 
    print "DATASET NOT KNOWN!"
    quit()

# END
####################

config = config()

from datetime import date
import time
folder = 'lowpteleid'
timestamp = date.today().strftime('_%d%m%y')+time.strftime('_%H%M%S')

config.General.requestName = folder+timestamp
config.General.workArea = folder
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuplizer_test_cfg.py'
#config.JobType.psetName = 'MINIAOD_from_AOD_cfg.py'
config.JobType.maxMemoryMB=3000
config.JobType.numCores = 1

config.Data.inputDBS = {'aod':'global','miniaod':'global','derived':'phys03'}.get(data_tier,'')
config.Data.inputDataset = data_set

#from files import files
#config.Data.userInputFiles=files

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/bainbrid/'+folder
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True # e.g. use CMSSW_10_2_15_patch2 on SLC7
config.JobType.inputFiles = ["../PhysicsTools/BParkingNano/test/lowPtEleReg_2018_02062020_nv.db"]

#config.Site.whitelist =['T2_UK_London_IC']
config.Site.storageSite = 'T2_CH_CERNBOX'
#config.Site.storageSite = 'T2_UK_London_IC'
