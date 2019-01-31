# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_14TeV_TuneCUETP8M1_cfi --conditions auto:phase2_realistic -n 10 --era Phase2C4 --eventcontent FEVTDEBUG --relval 9000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC14TeV --geometry Extended2023D28 --no_exec --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Phase2C4)

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('deltaR', 5,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 "deltaR between two candidates"
                 )
options.register('pdgId', 22,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "pdg id"
                 )
options.parseArguments()



# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D28Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D28_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from IOMC.EventVertexGenerators.VtxSmearedParameters_cfi import *
process.VtxSmeared = cms.EDProducer("FlatEvtVtxGenerator",
                                    MaxT   = cms.double(11.+0.177),
                                    MaxX   = cms.double(1e-08),
                                    MaxY   = cms.double(89),
                                    MaxZ   = cms.double(320.7),
                                    MinT   = cms.double(11.-0.177),
                                    MinX   = cms.double(-1e-08),
                                    MinY   = cms.double(87.),
                                    MinZ   = cms.double(320.6),
                                    readDB = cms.bool(False),
                                    src    = cms.InputTag("generator","unsmeared") 
                                    )


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('single particle gun_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.generator = cms.EDProducer("CloseByParticleGunProducer",
                                   PGunParameters = cms.PSet(PartID = cms.vint32(options.pdgId,options.pdgId),
                                                             En = cms.double(50.),
                                                             R = cms.double(options.deltaR),
                                                             MinEta = cms.double(1.5), #dummy values
                                                             MaxEta = cms.double(2.9),
                                                             MinPhi = cms.double(-3.1415),
                                                             MaxPhi = cms.double(+3.1415)
                                                             ),
                                   Verbosity = cms.untracked.int32(0),
                                   psethack = cms.string('two close by particles'),
                                   AddAntiParticle = cms.bool(False),
                                   firstRun = cms.untracked.uint32(1)
                                   )

#process.generator = cms.EDProducer("FlatRandomEGunProducer",
#                                   AddAntiParticle = cms.bool(False),
#                                   PGunParameters = cms.PSet(MaxE   = cms.double(50.01),
#                                                             MaxEta = cms.double(2.9),
#                                                             MaxPhi = cms.double(3.14159265359/6.),
#                                                             MinE   = cms.double(49.99),
#                                                             MinEta = cms.double(1.5),
#                                                             MinPhi = cms.double(-3.14159265359/6.),
#                                                             PartID = cms.vint32(22)
#                                                             ),
#                                   Verbosity = cms.untracked.int32(0),
#                                   firstRun = cms.untracked.uint32(1),
#                                   psethack = cms.string('single particle E 50')
#                                   )


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
