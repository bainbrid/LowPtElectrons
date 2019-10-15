# Recipe to produce ntuples used for BDT training 

## Release area
```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

## Get CMSSW package and model (2019Jun28) for electron ID
```
git cms-addpkg RecoEgamma/EgammaElectronProducers
git cms-merge-topic CMSBParking:from-CMSSW_10_2_15_LowPtElectronsID
git clone --branch 102X_LowPtElectrons_2019Jun28 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

## This is required if running on CRAB!
```
git cms-addpkg RecoEgamma/ElectronIdentification
mv $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data # this is required if running on CRAB
 ```

## Install ntuplizer code to produce (2019Jun28) ntuples
```
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons 
git checkout 2019Jul22_ntuples -b ntuplizer_dev
scram b
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun ntuplizer.py 
```
