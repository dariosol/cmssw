#include <memory>
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESProductHost.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/DataRecord/interface/EcalCATIAGainRatiosRcd.h"
#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "CondFormats/EcalObjects/src/classes.h"
//
// class declaration
//

class EcalCATIAGainRatioESProducer : public edm::ESProducer {

public:

  EcalCATIAGainRatioESProducer(const edm::ParameterSet& iConfig);

  typedef std::shared_ptr<EcalCATIAGainRatios> ReturnType;

  ReturnType produce(const EcalCATIAGainRatiosRcd& iRecord);

  //~EcalCATIAGainRatioESProducer();

private:
  edm::ParameterSet pset_;
  void setupChannelGain(const EcalCATIAGainRatiosRcd&);
  EcalCATIAGainRatio m_params;
};

EcalCATIAGainRatioESProducer::EcalCATIAGainRatioESProducer(const edm::ParameterSet& iConfig) : 
  pset_(iConfig) {
  //the following line is needed to tell the framework what
  // data is being produced
  setWhatProduced(this);
}
////
EcalCATIAGainRatioESProducer::ReturnType
EcalCATIAGainRatioESProducer::produce(const EcalCATIAGainRatiosRcd& iRecord){
  m_params=10; 
  std::shared_ptr<EcalCATIAGainRatios> prod;
  for(int i =0; i<61200;++i) {
    prod->setValue(i,10);
  }
  return prod;
}

void EcalCATIAGainRatioESProducer::setupChannelGain(const EcalCATIAGainRatiosRcd& rhs){
 
}


//Define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(EcalCATIAGainRatioESProducer);
