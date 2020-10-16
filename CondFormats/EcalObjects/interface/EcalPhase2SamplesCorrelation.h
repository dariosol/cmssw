#ifndef CondFormats_EcalObjects_Ecal2SamplesCorrelation_HH
#define CondFormats_EcalObjects_Ecal2SamplesCorrelation_HH


#include "CondFormats/Serialization/interface/Serializable.h"

#include "DataFormats/Math/interface/Matrix.h"
#include <iostream>
#include <vector>

class EcalPhase2SamplesCorrelation {
public:

  std::vector<double> EBG10SamplesCorrelation;
  std::vector<double> EBG1SamplesCorrelation;

  void print(std::ostream& o) const;

  COND_SERIALIZABLE;
};

#endif
