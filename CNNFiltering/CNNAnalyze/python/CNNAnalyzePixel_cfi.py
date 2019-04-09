import FWCore.ParameterSet.Config as cms

# ++IntermediateHitDoublets "detachedTripletStepHitDoublets" "" "RECO" (productId = 3:311)
# ++IntermediateHitDoublets "initialStepHitDoublets" "" "RECO" (productId = 3:312)
# ++IntermediateHitDoublets "initialStepHitDoubletsPreSplitting" "" "RECO" (productId = 3:313)
# ++IntermediateHitDoublets "lowPtTripletStepHitDoublets" "" "RECO" (productId = 3:314)
# ++IntermediateHitDoublets "mixedTripletStepHitDoubletsA" "" "RECO" (productId = 3:315)
# ++IntermediateHitDoublets "mixedTripletStepHitDoubletsB" "" "RECO" (productId = 3:316)
# ++IntermediateHitDoublets "pixelLessStepHitDoublets" "" "RECO" (productId = 3:317)
# ++IntermediateHitDoublets "tobTecStepHitDoubletsTripl" "" "RECO" (productId = 3:318)
# ++IntermediateHitDoublets "tripletElectronHitDoublets" "" "RECO" (productId = 3:319)

# 2018
# ++IntermediateHitDoublets "detachedQuadStepHitDoublets" "" "RECO" (productId = 3:309)
# ++IntermediateHitDoublets "detachedTripletStepHitDoublets" "" "RECO" (productId = 3:310)
# ++IntermediateHitDoublets "highPtTripletStepHitDoublets" "" "RECO" (productId = 3:311)
# ++IntermediateHitDoublets "initialStepHitDoublets" "" "RECO" (productId = 3:312)
# ++IntermediateHitDoublets "initialStepHitDoubletsPreSplitting" "" "RECO" (productId = 3:313)
# ++IntermediateHitDoublets "lowPtQuadStepHitDoublets" "" "RECO" (productId = 3:314)
# ++IntermediateHitDoublets "lowPtTripletStepHitDoublets" "" "RECO" (productId = 3:315)
# ++IntermediateHitDoublets "mixedTripletStepHitDoubletsA" "" "RECO" (productId = 3:316)
# ++IntermediateHitDoublets "mixedTripletStepHitDoubletsB" "" "RECO" (productId = 3:317)
# ++IntermediateHitDoublets "pixelLessStepHitDoublets" "" "RECO" (productId = 3:318)
# ++IntermediateHitDoublets "tobTecStepHitDoubletsTripl" "" "RECO" (productId = 3:319)
# ++IntermediateHitDoublets "tripletElectronHitDoublets" "" "RECO" (productId = 3:320)
#pixelPairStepHitDoublets
#jetCoreRegionalStepHitDoublets

pixelPairStepCNN = cms.EDAnalyzer('CNNAnalyze',
        processName = cms.string( "pixelTracksHitDoublets"),
        doublets    = cms.InputTag("pixelTracksHitDoublets" ),
        tpMap       = cms.InputTag("tpClusterProducerPixelTrackingOnly" ),
        beamSpot    = cms.InputTag("offlineBeamSpot"),
        infoPileUp  = cms.InputTag("addPileupInfo")
)


CNNDoubletsPixelSequence = cms.Sequence(pixelPairStepCNN)
                                   #pixelPairStepCNN *
                                   #jetCoreRegionalStepCNN)
