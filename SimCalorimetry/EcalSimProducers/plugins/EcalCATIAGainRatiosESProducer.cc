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
const int kEBChannels = 61200, kEEChannels = 14648;

class EcalCATIAGainRatiosESProducer : public edm::ESProducer {

public:

  EcalCATIAGainRatiosESProducer(const edm::ParameterSet& iConfig);

  typedef std::shared_ptr<EcalCATIAGainRatios> ReturnType;

  ReturnType produce(const EcalCATIAGainRatiosRcd& iRecord);


private:
  edm::ParameterSet pset_;
};

EcalCATIAGainRatiosESProducer::EcalCATIAGainRatiosESProducer(const edm::ParameterSet& iConfig) : 
  pset_(iConfig) {
  //the following line is needed to tell the framework what
  // data is being produced
  //std::cout<<"*********Creating EcalCATIAGainRatiosESProducer"<<std::endl;
  setWhatProduced(this);
}
////
EcalCATIAGainRatiosESProducer::ReturnType
EcalCATIAGainRatiosESProducer::produce(const EcalCATIAGainRatiosRcd& iRecord){
  //std::cout<<"********Starting Production"<<std::endl;
  auto prod = std::make_shared<EcalCATIAGainRatios>();

  //std::cout<<"**********Set EB Values "<<std::endl;

  for (int iChannel = 0; iChannel < kEBChannels; iChannel++) {
    EBDetId myEBDetId = EBDetId::unhashIndex(iChannel);
    float val = 10.;
    prod->setValue(myEBDetId.rawId(), val);
  }

  //std::cout<<"**********Set EE Values "<<std::endl;

  for (int iChannel = 0; iChannel < kEEChannels; iChannel++) {
    EEDetId myEEDetId = EEDetId::unhashIndex(iChannel);     
    float val = 10.;
    prod->setValue(myEEDetId.rawId(), val);
  }
   
  
  
  
  //std::cout<<prod->size()<<std::endl;
 
  //std::cout<<"***********Returning"<<std::endl;
  return prod;
}



//Define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(EcalCATIAGainRatiosESProducer);
