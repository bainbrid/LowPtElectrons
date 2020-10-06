import FWCore.ParameterSet.Config as cms

testPATLowPtElectrons = cms.EDAnalyzer("PATLowPtElectrons_UnitTest",
                                       patElectrons = cms.InputTag("slimmedLowPtElectrons"),
                                       gsfElectrons = cms.InputTag("lowPtGsfElectrons"),
                                       mvaId = cms.InputTag("lowPtGsfElectronID"),
                                       verbose = cms.bool(True),
                                   )
