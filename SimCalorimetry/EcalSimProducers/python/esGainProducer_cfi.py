import FWCore.ParameterSet.Config as cms

EcalGainRatiosRcd =  cms.ESSource("EmptyESSource",
                                recordName = cms.string("EcalGainRatiosRcd"),
                                firstValid = cms.vuint32(1),
                                iovIsRunNotTime = cms.bool(True)
                                )

EcalGainRatios = cms.ESProducer("EcalGainRatiosESProducer"),
 
timeThresh=cms.double(2.0),
