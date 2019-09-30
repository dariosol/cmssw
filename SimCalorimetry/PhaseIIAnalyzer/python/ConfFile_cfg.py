import FWCore.ParameterSet.Config as cms

process = cms.Process("PhaseII")

#Print out:
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PhaseIIAnalyzer')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )



#Max event
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#input:
process.source = cms.Source("PoolSource",

                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/j/jobereng/CMSSW_10_3_1/src/SimCalorimetry/EcalSimProducers/test/SingleElectronPt10_pythia8_cfi_py_GEN_SIM_DIGI_Pt10000.root'),

                            )

process.phaseII = cms.EDAnalyzer('PhaseIIAnalyzer')
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('Digitizer2.root')
                                   )
process.load("SimCalorimetry.PhaseIIAnalyzer.CfiFile_cfi")

#Running process:
process.p = cms.Path(process.phaseII)
