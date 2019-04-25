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
    input = cms.untracked.int32(-1)
)



process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/muchLooserTune/stepTICL_211_Pt10_HighEta.root',
        'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/muchLooserTune/stepTICL_211_Pt10_LowEta.root'
        #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/looserTune/stepTICL_211_HighEta_Pt10.root',
        #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/looserTune/stepTICL_211_LowEta_Pt10.root'
        #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/stepTICL_211Pt10_HighEta.root', 
        #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/stepTICL_211Pt10_LowEta.root'   

        #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/looserTune/stepTICL_211_211_Pt10.root',
        #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/looserTune/stepTICL_211_130_Pt10.root'
        #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/stepTICL_211_211_Pt10.root',
        #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/looserTune/stepTICL_211_130_Pt10.root'

        ),
                            secondaryFileNames = cms.untracked.vstring(),
                             noEventSort = cms.untracked.bool(True),
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )


#process.source = cms.Source ("PoolSource",    
#                             fileNames = cms.untracked.vstring(options.input),
#                             secondaryFileNames = cms.untracked.vstring(),
#                             noEventSort = cms.untracked.bool(True),
#                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
#                             )




process.ticlAnalyzer = cms.EDAnalyzer("BuildParticles")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:ticl_analysis_BuildP_211_test.root")
                                   )


process.p = cms.Path(process.ticlAnalyzer)


