#include "CondFormats/EcalObjects/interface/EcalPhase2SamplesCorrelation.h"


template <typename T>
static inline void print_vector(std::ostream& o, const std::vector<T>& vect) {
  o << "[";
  for (std::vector<double>::const_iterator i = vect.begin(); i != vect.end(); ++i) {
    std::cout << *i << ", ";
  }
  o << "]";
}

void EcalPhase2SamplesCorrelation::print(std::ostream& o) const {
  o << "EB Gain 10 correlation:";
  print_vector<double>(o, this->EBG10SamplesCorrelation);
  o << std::endl;
 
  o << "EB Gain 1 correlation:";
  print_vector<double>(o, this->EBG1SamplesCorrelation);
  o << std::endl;

 }
