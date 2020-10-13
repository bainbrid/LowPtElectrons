import FWCore.ParameterSet.Config as cms

testPATLowPtElectrons = cms.EDAnalyzer("PATLowPtElectrons_UnitTest",
                                       patElectrons = cms.InputTag("slimmedLowPtElectrons"),
                                       gsfElectrons = cms.InputTag("lowPtGsfElectrons"),
                                       mvaId = cms.InputTag("lowPtGsfElectronID"),
                                       mvaUnbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
                                       mvaPtbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased"),
                                       lowPtGsfLinks = cms.InputTag("lowPtGsfLinks"),
                                       verbose = cms.bool(True),
                                   )
