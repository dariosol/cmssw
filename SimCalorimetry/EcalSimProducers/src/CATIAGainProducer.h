// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

//
// class decleration
//
class CATIAGainProducer: public edm::ESProducer 
{
public:
  using ReturnType = std::map<int,float>;
  CATIAGainProducer(const edm::ParameterSet&); 
  ~CATIAGainProducer() override {} ;

  ReturnType produceGain() ;

private:

};
