#ifndef EcalPulseCovariances_h
#define EcalPulseCovariances_h

#include "CondFormats/Serialization/interface/Serializable.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"

template <size_t nsamples> struct EcalPulseCovarianceT {
  static const size_t TEMPLATESAMPLES = nsamples;
public:
  EcalPulseCovarianceT() {
    for (size_t i = 0; i < TEMPLATESAMPLES; ++i) {
      for (size_t j = 0; j < TEMPLATESAMPLES; ++j) {
        covval[i][j] = 0.;
      }
    }
  };
  float covval[TEMPLATESAMPLES][TEMPLATESAMPLES];
  float val(size_t i, size_t j) const { return covval[i][j]; }

  COND_SERIALIZABLE;
};

typedef EcalPulseCovarianceT<12> EcalPulseCovariance;
typedef EcalPulseCovarianceT<16> EcalPhase2PulseCovariance;

typedef EcalCondObjectContainer<EcalPulseCovariance > EcalPulseCovariancesMap;
typedef EcalPulseCovariancesMap::const_iterator EcalPulseCovariancesMapIterator;
typedef EcalPulseCovariancesMap EcalPulseCovariances;


typedef EcalCondObjectContainer<EcalPhase2PulseCovariance > EcalPhase2PulseCovariancesMap;
typedef EcalPhase2PulseCovariancesMap::const_iterator EcalPhase2PulseCovariancesMapIterator;
typedef EcalPhase2PulseCovariancesMap EcalPhase2PulseCovariances;


#endif
