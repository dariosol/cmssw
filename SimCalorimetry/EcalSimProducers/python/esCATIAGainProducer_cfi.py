import FWCore.ParameterSet.Config as cms

EcalCATIAGainRatiosRcd =  cms.ESSource("EmptyESSource",
                                recordName = cms.string("EcalCATIAGainRatiosRcd"),
                                firstValid = cms.vuint32(1),
                                iovIsRunNotTime = cms.bool(True)
                                )

EcalCATIAGainRatios = cms.ESProducer("EcalCATIAGainRatiosESProducer"),
 
timeThresh=cms.double(2.0),
