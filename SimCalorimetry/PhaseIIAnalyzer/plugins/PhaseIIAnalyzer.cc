// -*- C++ -*-
//
// Package:    PhaseII/PhaseIIAnalyzer
// Class:      PhaseIIAnalyzer
// 
/**\class PhaseIIAnalyzer PhaseIIAnalyzer.cc PhaseII/PhaseIIAnalyzer/plugins/PhaseIIAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dario Soldi
//         Created:  Tue, 16 Jan 2018 11:56:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//My includes
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "PhaseIIAnalyzer.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/Track/interface/SimTrack.h"


using namespace std;
using namespace edm;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PhaseIIAnalyzer::PhaseIIAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
usesResource("TFileService");

 digiTagEB_= iConfig.getParameter<edm::InputTag>("BarrelDigis");
 digiTagEE_= iConfig.getParameter<edm::InputTag>("EndcapDigis");

 digiTokenEB_ = consumes<EBDigiCollection>(digiTagEB_);
 digiTokenEE_ = consumes<EEDigiCollection>(digiTagEE_);
 
 //Files:
 edm::Service<TFileService> fs;
 //Histograms

 for(int isample=0;isample<10;isample++){
   EBEnergyHisto[isample] = fs->make<TH1I>(Form("EnergyEB_%d", isample), Form("Energy sample %d  Barrel;ADC",isample), 950 , 100 , 2000 );
   EEEnergyHisto[isample] = fs->make<TH1I>(Form("EnergyEE_%d", isample), Form("Energy sample %d  Endcap;ADC",isample), 950 , 100 , 2000 );
   EBGainHisto[isample]   = fs->make<TH1I>(Form("GainEB_%d", isample), Form("Gain Barrel sample %d;Gain",isample) , 5 , 0 , 4 );
   EEGainHisto[isample]   = fs->make<TH1I>(Form("GainEE_%d", isample), Form("Gain Endcap sample %d;Gain",isample) , 5 , 0 , 4 );


   Char_t histo[200];


  
   sprintf (histo, "EcalDigiTask Barrel occupancy" ) ;
   meEBDigiOccupancy_ = fs->make<TH2D>(histo, histo, 360, 0., 360., 170, -85., 85.);

   sprintf (histo, "EcalDigiTask Barrel digis multiplicity" ) ;
   meEBDigiMultiplicity_ = fs->make<TH1D>(histo, histo, 612, 0., 61200);
  
    
   for (int i = 0; i < 10 ; i++ ) {

     sprintf (histo, "EcalDigiTask Barrel analog pulse %02d", i+1) ;
     meEBDigiADCAnalog_[i] = fs->make<TH1D>(histo, histo, 4000, 0., 400.);

     sprintf (histo, "EcalDigiTask Barrel ADC pulse %02d Gain 0 - Saturated", i+1) ;
     meEBDigiADCgS_[i] = fs->make<TH1D>(histo, histo, 4096, -0.5, 4095.5);

     sprintf (histo, "EcalDigiTask Barrel ADC pulse %02d Gain 1", i+1) ;
     meEBDigiADCg1_[i] = fs->make<TH1D>(histo, histo, 4096, -0.5, 4095.5);

     sprintf (histo, "EcalDigiTask Barrel ADC pulse %02d Gain 6", i+1) ;
     meEBDigiADCg6_[i] = fs->make<TH1D>(histo, histo, 4096, -0.5, 4095.5);

     sprintf (histo, "EcalDigiTask Barrel ADC pulse %02d Gain 12", i+1) ;
     meEBDigiADCg12_[i] = fs->make<TH1D>(histo, histo, 4096, -0.5, 4095.5);

     sprintf (histo, "EcalDigiTask Barrel gain pulse %02d", i+1) ;
     meEBDigiGain_[i] = fs->make<TH1D>(histo, histo, 4, 0, 4);

   }
    
   sprintf (histo, "EcalDigiTask Barrel pedestal for pre-sample" ) ;
   meEBPedestal_ = fs->make<TH1D>(histo, histo, 4096, -0.5, 4095.5) ;

   sprintf (histo, "EcalDigiTask Barrel maximum position gt 100 ADC" ) ;
   meEBMaximumgt100ADC_ = fs->make<TH1D>(histo, histo, 10, 0., 10.) ;

   sprintf (histo, "EcalDigiTask Barrel maximum position gt 10 ADC" ) ;
   meEBMaximumgt10ADC_ = fs->make<TH1D>(histo, histo, 10, 0., 10.) ;

   //   sprintf (histo, "EcalDigiTask Barrel ADC counts after gain switch" ) ;
   //meEBnADCafterSwitch_ = fs->make<TH1D>histo, histo, 10, 0., 10.) ;

 }
 
 
}
  

PhaseIIAnalyzer::~PhaseIIAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhaseIIAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  //LogInfo("PhaseII") << "new event ";

   Handle<EBDigiCollection> pDigiEB;
   iEvent.getByToken(digiTokenEB_,pDigiEB);
   
   Handle<EEDigiCollection> pDigiEE;
   iEvent.getByToken(digiTokenEE_,pDigiEE);
   
   const int MAXSAMPLES=10; 
   std::vector<double> ebAnalogSignal ;
   std::vector<double> ebADCCounts ;
   std::vector<double> ebADCGains ;
   ebAnalogSignal.reserve(EBDataFrame::MAXSAMPLES);
   ebADCCounts.reserve(EBDataFrame::MAXSAMPLES);
   ebADCGains.reserve(EBDataFrame::MAXSAMPLES);
   int nDigis=0;
   for (EBDigiCollection::const_iterator pDigi=pDigiEB->begin(); pDigi!=pDigiEB->end(); ++pDigi) {
     EBDataFrame digi( *pDigi );
     int nrSamples = digi.size();
     EBDetId ebid = digi.id () ;

     nDigis++;
     if (meEBDigiOccupancy_) meEBDigiOccupancy_->Fill( ebid.iphi(), ebid.ieta() ); 

     double Emax = 0. ;
     int Pmax = 0 ;
     double pedestalPreSample = 0.;
     double pedestalPreSampleAnalog = 0.;
     int countsAfterGainSwitch = -1;
     double higherGain = 1.;
     int higherGainSample = 0;
     
     for (int sample = 0 ; sample < nrSamples; ++sample) {
       ebAnalogSignal[sample] = 0.;
       ebADCCounts[sample] = 0.;
       ebADCGains[sample] = 0.;
     }

   double  gainConv_[2]={10,1};
   // saturated channels
   double barrelADCtoGeV_ = 0.048; //GeV


     
for (int sample = 0 ; sample < nrSamples; ++sample) 
        {
	  int thisSample = digi[sample];
	  
          ebADCCounts[sample] = (thisSample&0xFFF);
          ebADCGains[sample]  = (thisSample&(0x3<<12));
          ebAnalogSignal[sample] = (ebADCCounts[sample]*gainConv_[(int)ebADCGains[sample]]*barrelADCtoGeV_);
	  if((thisSample&0xfff)>250) {
	  cout<<"Full data "<<thisSample<<endl;
	  cout<<"Sample "<<sample<<" E: "<<(thisSample&0xfff)<<" gain: "<< (thisSample&(0x3<<12))<<" Analog "<<ebAnalogSignal[sample]<<endl;
	    }
          if (Emax < ebAnalogSignal[sample] ) {
            Emax = ebAnalogSignal[sample] ;
            Pmax = sample ;
          }
	}
 if(0==1)cout<<"P max "<<Pmax<<endl;
   }//end digi
   
   

   
}


// ------------ method called once each job just before starting event loop  ------------
void 
PhaseIIAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhaseIIAnalyzer::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhaseIIAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhaseIIAnalyzer);
