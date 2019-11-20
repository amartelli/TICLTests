import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023D28Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)



process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'file:../../TICL/test/stepTICL_211_Pt10_dummy_testTime.root'
#'file:../../TICL/test/stepTICL_211_Pt10_dummy_testTime_v2.root'
'file:/afs/cern.ch/work/a/amartell/HGCAL/2019_ticl/CMSSW_11_0_X_2019-06-23-2300/src/RecoHGCal/TICL/test/stepTICL_22_Pt150_LE_timeOK_PU.root'

        ),
                            secondaryFileNames = cms.untracked.vstring(),
                             noEventSort = cms.untracked.bool(True),
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )



process.ticlAnalyzer = cms.EDAnalyzer("CheckResolutionV2")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:ticl_singlePhoton_test_TIME_PU.root")
                                   )


process.p = cms.Path(process.ticlAnalyzer)


