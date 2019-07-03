//Namespaces for Phase1 and Phase2
#ifndef EcalObject_EcalConstants_h
#define EcalObject_EcalConstants_h

namespace ecalPh2 
{ 
  constexpr double BUNCHSPACE = 6.25; 
  constexpr int NGAINS = 2;
  constexpr float gains[2] = {10.,1.};
  constexpr int sampleSize = 16;
} 


namespace ecalPh1 
{ 
  constexpr double BUNCHSPACE = 25.; 
  constexpr int NGAINS = 4; 
  constexpr float gains[4] = {0.,12.,6.,1.};
  constexpr int sampleSize = 10;
} 
#endif
