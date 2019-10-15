# Recipe to run macros used to train/evaluate the BDTs 

## Release area

```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

## Install latest macros to train/evaluate using 2019Jul22 ntuples

```
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/macros
git checkout tags/2019Jul22_model -b macros_dev
```

## Evaluating with the 2019Jul22 ntuples

First, create a soft link to the existing XML model file.

```
ln -s /eos/cms/store/cmst3/group/bpark/electron_training/2019Jul22/models models/2019Jul22
```

Then:

```
python train_bdt.py cmssw_mva_id --load_model --config models/2019Jul22/mauro.json
```

ROC curves are found in ```plots/2019Jul22```. 

Acceptance times efficiency and fake rate values:

```
EGamma GSF trk:              AxE = 0.232, FR = 0.001
EGamma PF electron:          AxE = 0.231, FR = 0.001
Low pT GSF trk (PreId):      AxE = 0.535, FR = 0.171
Low pT electron (CMSSW):     AxE = 0.496, FR = 0.121
Low pT electron (2019Jun22): AxE = 0.496, FR = 0.121
```

Further plots:

```
python eval_bdt.py models/2019Jul22/bdt_cmssw_mva_id --what cmssw_mva_id --dataset test --noxml --plot plots/2019Jul22
```
