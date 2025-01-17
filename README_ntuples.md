# Recipe to produce ntuples used for BDT training 

## Release area
```
cmssw-el6
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

## Add packages and ntuplizer 
```
git cms-addpkg RecoEgamma/EgammaElectronProducers
git cms-addpkg RecoEgamma/EgammaTools
git cms-merge-topic -u CMSBParking:from-CMSSW_10_2_15_2020Sept15
# old branches^^^: from-CMSSW_10_2_15_2019Aug07, from-CMSSW_10_2_15_2019Jul22, from-CMSSW_10_2_15_2019Jun28

git cms-addpkg CommonTools/MVAUtils
git cms-merge-topic -u CMSBParking:convertXMLToGBRForestROOT

git cms-addpkg RecoEgamma/ElectronIdentification
git clone --branch 102X_LowPtElectrons git@github.com:bainbrid/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
# old branches^^^: 102X_LowPtElectrons_2019Aug07, 102X_LowPtElectrons_2019Jul22, 102X_LowPtElectrons_2019Jun28
```

## This is required if running on CRAB!
```
mv $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data
 ```

## Build and run
```
scram b -j8
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun ntuplizer_test_cfg.py maxEvents=10
```
