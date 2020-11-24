import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)



process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'file:../../TICL/test/stepTICL_211_Pt10_dummy_testTime.root'
#'file:../../TICL/test/stepTICL_211_Pt10_dummy_testTime_v2.root'
#'file:/afs/cern.ch/work/a/amartell/HGCAL/2019_ticl/CMSSW_11_0_X_2019-06-23-2300/src/RecoHGCal/TICL/test/stepTICL_22_Pt150_LE_timeOK_PU.root'
#'file:/tmp/amartell/step3_singlepi_pT0to100GeV_scanEta_nopu_5k.root'
'file:/tmp/amartell/TTbar_0PU_306eaecf-87a4-4d5a-9d3f-16e6f98346f7.root'
        ),
                            secondaryFileNames = cms.untracked.vstring(),
                             noEventSort = cms.untracked.bool(True),
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )



process.ticlAnalyzer = cms.EDAnalyzer("CheckResolution")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:ticl_ttbar.root")
                                   )


process.p = cms.Path(process.ticlAnalyzer)


