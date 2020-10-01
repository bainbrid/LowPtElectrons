# Unknown production of AODSIM (and MINIAODSIM)
# config.Site.whitelist =['T2_UK_London_IC']
# ~400 files, ~1M events (?)
path_190328_152903='root://cms-xrd-global.cern.ch//store/user/bainbrid/lowpteleid/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_lowpteleid/190328_152903/0000/step3_inAODSIM_XXX.root'
jobids = [ id+1 for id in range(400) ]
for id in [1,9,152,155,157,162,164,177,196,202,204,212,281,321,372,389] : jobids.remove(id) 
files_190328_152903 = [ path_190328_152903.replace('XXX',str(id)) for id in jobids ]

# (Buggy) Test production of AODSIM (and MINIAODSIM) using BuToKJpsi_Toee based on CMSSW_10_2_14
# Store Ref to seed track in GsfElectronCore for low pT but NOT PF (due to bug!)
# Keep recoElectronSeeds and recoTrackExtras in AODSIM
# config.Site.whitelist =['T2_UK_London_IC']
# ~100 files, ~25k events (?)
path_190503_223344='root://cms-xrd-global.cern.ch//store/user/bainbrid/lowpteleid/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_lowpteleid/190503_223344/0000/step3_inAODSIM_XXX.root'
jobids = [ id+1 for id in range(100) ]
for id in [8,20,27,28,53,55,86] : jobids.remove(id) 
files_190503_223344 = [ path_190503_223344.replace('XXX',str(id)) for id in jobids ]

# Test production of AODSIM (and MINIAODSIM) using BuToKJpsi_Toee based on CMSSW_10_2_14
# Store Ref to seed track in GsfElectronCore for low pT and PF
# Keep recoElectronSeeds and recoTrackExtras in AODSIM
# config.Site.whitelist =['T2_UK_London_IC']
# ~1000 files, ~4M events (?)
path_190511_200059='root://cms-xrd-global.cern.ch//store/user/bainbrid/lowpteleid/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_lowpteleid/190511_200059/0000/step3_inAODSIM_XXX.root'
jobids = [ id+1 for id in range(6,999) ]
for id in [18,21,88,165,246,247,283,296,369,445,474,506,564,573,611,637,715,782,854,929] : jobids.remove(id) 
files_190511_200059 = [ path_190511_200059.replace('XXX',str(id)) for id in jobids ]

# Test production of MINIAODSIM (only) using BuToKJpsi_Toee based on CMSSW_10_2_14, at CERNBOX !!!
# config.Site.whitelist =['T2_CH_CERNBOX']
# ~400 files, ~1M events
# Code as above
path_190624_111044='root://cms-xrd-global.cern.ch//store/user/bainbrid/lowpteleid/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_lowpteleid/190624_111044/0000/step3_inMINIAODSIM_XXX.root'
jobids = [ id+1 for id in range(0,400) ]
for id in [63,70,77,265] : jobids.remove(id) 
files_190624_111044 = [ path_190624_111044.replace('XXX',str(id)) for id in jobids ]

# New production of AODSIM (and MINIAOD) using BuToKJpsi_Toee based on CMSSW_10_2_14, at CERNBOX !!!
# config.Site.whitelist =['T2_CH_CERNBOX']
# ~400 files, ~1M events
# Code as above
path_190625_113922='root://eosuser.cern.ch//eos/user/b/bainbrid/lowpteleid/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_lowpteleid/190625_113922/0000/step3_inAODSIM_XXX.root'
jobids = [ id+1 for id in range(0,399) ]
for id in [308,314] : jobids.remove(id) 
files_190625_113922 = [ path_190625_113922.replace('XXX',str(id)) for id in jobids ]

# Choose production here!
files = files_190625_113922
#for f in files : print f
