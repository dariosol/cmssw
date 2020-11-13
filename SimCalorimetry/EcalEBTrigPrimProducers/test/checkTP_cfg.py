import FWCore.ParameterSet.Config as cms

process = cms.Process("PROdTPA")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("CalibCalorimetry.EcalTPGTools.ecalTPGScale_cff")

process.source = cms.Source("PoolSource",

fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/n/nancy/private/EcalL1/workTPHPhaseII/workWithNewDigis/Oct21/CMSSW_11_2_0_pre7/src/step1_UpToDigi.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
)

from SimCalorimetry.EcalTrigPrimProducers.ecalTrigPrimESProducer_cff import *

process.tpAnalyzer = cms.EDAnalyzer("EcalEBTrigPrimAnalyzer",
                                    inputRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                                    barrelEcalDigis = cms.InputTag("simEcalDigis","ebDigis"),
                                    AnalyzeRecHits = cms.bool(False),
                                    Debug = cms.bool(False),
                                    inputTP = cms.InputTag("simEcalEBTriggerPrimitiveDigis")
    
)



process.p = cms.Path(process.tpAnalyzer)
