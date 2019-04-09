import FWCore.ParameterSet.Config as cms

phaseII = cms.EDAnalyzer('PhaseIIAnalyzer',
                            BarrelDigis=cms.InputTag('simEcalDigis','ebDigis','DIGI'),
                            EndcapDigis=cms.InputTag('simEcalDigis','eeDigis','DIGI')
)
