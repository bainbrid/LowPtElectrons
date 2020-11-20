import FWCore.ParameterSet.Config as cms

testPATLowPtElectrons = cms.EDAnalyzer("PATLowPtElectrons_UnitTest",
                                       patElectrons = cms.InputTag("slimmedLowPtElectrons"),
                                       gsfElectrons = cms.InputTag("lowPtGsfElectrons"),
                                       gsfTracks = cms.InputTag("lowPtGsfEleGsfTracks"),
                                       mvaId = cms.InputTag("lowPtGsfElectronID"),
                                       mvaUnbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
                                       mvaPtbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased"),
                                       lowPtGsfLinks1 = cms.InputTag("lowPtGsfLinks:gsf2packed"),
                                       lowPtGsfLinks2 = cms.InputTag("lowPtGsfLinks:gsf2lost"),
                                       verbose = cms.bool(True),
                                       PFCandidates =  cms.InputTag("particleFlow"),
                                       packedCandidates = cms.InputTag("packedPFCandidates"),
                                       patPfElectrons = cms.InputTag("slimmedElectrons"),
                                   )
