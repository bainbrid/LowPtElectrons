import FWCore.ParameterSet.Config as cms

ntuplizer = cms.EDFilter( # cms.EDAnalyzer
    "IDNtuplizerSimple",
    verbose = cms.int32(0),
    checkFromB = cms.bool(True),
    drMax = cms.double(0.1),
    drThreshold = cms.double(0.02),
    prescale = cms.double(0), # zero: no prescale, +ve: use 1/prescale, -ve: (poisson) mean number of fakes/event
    minPt = cms.double(0.5),
    tagMuonPtThreshold = cms.double(7.),
    tagMuonEtaThreshold = cms.double(1.5),
    filterNtupleContent = cms.bool(False), #@@ not needed, just keep for now to run cfg
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    #beamspot = cms.InputTag("offlineBeamSpot"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"), # MINIAOD
    patElectronsEGamma = cms.InputTag("slimmedElectrons"), # MINIAOD
    mvaValueEGamma = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2RawValues'),
    mvaValueEGammaRetrained = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'),
    )
