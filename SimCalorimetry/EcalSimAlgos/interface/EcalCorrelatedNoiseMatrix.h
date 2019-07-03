#ifndef EcalSimAlgos_EcalCorrelatedNoiseMatrix_h
#define EcalSimAlgos_EcalCorrelatedNoiseMatrix_h

#include "DataFormats/Math/interface/Error.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "CondFormats/EcalObjects/interface/EcalConstants.h"

typedef math::ErrorD<ecalPh2::sampleSize>::type EcalCorrMatrix;

#endif
