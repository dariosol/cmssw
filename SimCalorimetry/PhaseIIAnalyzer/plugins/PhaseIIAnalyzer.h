#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "TH1.h"
#include "TH2.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include <fstream>

using namespace std;
//
// class declaration
//

class PhaseIIAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
 public:
  explicit PhaseIIAnalyzer(const edm::ParameterSet&);
  ~PhaseIIAnalyzer();     

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


 private:
  virtual void beginJob() ;
  void analyze(const edm::Event&, const edm::EventSetup&) override ;
  virtual void endJob() ;

  edm::InputTag digiTagEB_;
  edm::InputTag digiTagEE_;
  
  edm::EDGetTokenT<EBDigiCollection> digiTokenEB_; 
  edm::EDGetTokenT<EEDigiCollection> digiTokenEE_;
  

  edm::InputTag hitTagEB_;
  edm::InputTag hitTagEE_;
  
  edm::EDGetTokenT<vector<PCaloHit>> hitTokenEB_; 
  edm::EDGetTokenT<vector<PCaloHit>> hitTokenEE_;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<vector<SimTrack>> trackToken_; 
  
  //Histograms
   TH1I *EBEnergyHisto[10];
   TH1I *EEEnergyHisto[10];
   TH1I *EBGainHisto[10];
   TH1I *EEGainHisto[10];

   TH2D* meEBDigiOccupancy_;

   TH1D* meEBDigiMultiplicity_;

  TH1D*  meEBDigiADCGlobal_;


  TH1D*  meEBDigiADCAnalog_[10];
  TH1D*  meEBDigiADCgS_[10];
  TH1D*  meEBDigiADCg1_[10]; 
  TH1D*  meEBDigiADCg6_[10]; 
  TH1D*  meEBDigiADCg12_[10];
  TH1D*  meEBDigiGain_[10]; 


 TH1D*  meEBPedestal_;
 
 TH1D*   meEBMaximumgt100ADC_; 
 
 TH1D* meEBMaximumgt10ADC_; 
 
 TH1D*  meEBnADCafterSwitch_;
 








  
};
