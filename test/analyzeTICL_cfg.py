import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.source = cms.Source ("PoolSource",    
                             fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/RECO_211Pt10VtxNoSmear.root'
        ),
                             secondaryFileNames = cms.untracked.vstring(),
                             noEventSort = cms.untracked.bool(True),
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                             )




process.ticlAnalyzer = cms.EDAnalyzer("TICLAnalyzer")

process.p = cms.Path(process.ticlAnalyzer)


