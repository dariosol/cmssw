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
  std::cout<<"*********Creating EcalCATIAGainRatiosESProducer"<<std::endl;
  setWhatProduced(this);
}
////
EcalCATIAGainRatiosESProducer::ReturnType
EcalCATIAGainRatiosESProducer::produce(const EcalCATIAGainRatiosRcd& iRecord){
  std::cout<<"********Starting Production"<<std::endl;
  std::shared_ptr<EcalCATIAGainRatios> prod;
  EcalCATIAGainRatios *test = new EcalCATIAGainRatios();

   for (int iChannel = 0; iChannel < kEBChannels; iChannel++) {
     std::cout<<"**********Set EB Values "<<iChannel<<std::endl;
     EBDetId myEBDetId = EBDetId::unhashIndex(iChannel);
     std::cout<<"**********EB ID "<<myEBDetId.rawId()<<std::endl;
    float val = 10;
    test->setValue(myEBDetId.rawId(), val);
  }
   
   for (int iChannel = 0; iChannel < kEEChannels; iChannel++) {
     std::cout<<"**********Set EE Values "<<iChannel<<std::endl;
     EEDetId myEEDetId = EEDetId::unhashIndex(iChannel);     
     float val = 10;
     prod->setValue(myEEDetId.rawId(), val);
   }
   
  
  
  
  std::cout<<prod->size()<<std::endl;
 
  std::cout<<"***********Returning"<<std::endl;
  return prod;
}



//Define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(EcalCATIAGainRatiosESProducer);
