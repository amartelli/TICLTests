import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'file:/tmp/amartell/TICL_test_ttB_dR_photon.root'
'file:/tmp/amartell/TICL_test_ttB_dR_pre6.root'

        ),
                            secondaryFileNames = cms.untracked.vstring(),
                             noEventSort = cms.untracked.bool(True),
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )



process.ticlAnalyzer = cms.EDAnalyzer("AnalyzeTracksters")


process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string("file:ticl_photon.root")
                                   fileName = cms.string("file:ticl_ttbar.root")
                                   )


process.p = cms.Path(process.ticlAnalyzer)


