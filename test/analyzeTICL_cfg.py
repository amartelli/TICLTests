import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('input',None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file"
                 )
options.parseArguments()


process = cms.Process("ANA")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.source = cms.Source ("PoolSource",    
                             fileNames = cms.untracked.vstring(options.input),
                             secondaryFileNames = cms.untracked.vstring(),
                             noEventSort = cms.untracked.bool(True),
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                             )




process.ticlAnalyzer = cms.EDAnalyzer("TICLAnalyzer")

process.p = cms.Path(process.ticlAnalyzer)


