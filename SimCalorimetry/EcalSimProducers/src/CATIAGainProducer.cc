// user include files
#include "CATIAGainProducer.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HcalDetId/interface/HcalCastorDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/CaloTowerGeometry.h"
#include "Geometry/ForwardGeometry/interface/CastorGeometry.h"
#include "Geometry/ForwardGeometry/interface/ZdcGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
// member functions
//
CATIAGainProducer::CATIAGainProducer(const edm::ParameterSet&)
{
  //the following line is needed to tell the framework what
  // data is being produced
  setWhatProduced( this, &CATIAGainProducer::produceGain ); 
 }

// ------------ method called to produce the data  ------------

CATIAGainProducer::ReturnType
CATIAGainProducer::produceGain( )
{

  ReturnType pGain;
  
  for (int ichannel=0;ichannel<61200;++ichannel ) 
    {
      // loop on ECAL channels
      edm::LogInfo("CATIAGainProducer") << "Building CATIA gain values";
      pGain[EBDetId::unhashIndex(ichannel)] = 10;
    }
  return pGain ;
}
