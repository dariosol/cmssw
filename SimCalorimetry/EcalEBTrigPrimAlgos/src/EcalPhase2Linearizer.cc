#include <SimCalorimetry/EcalEBTrigPrimAlgos/interface/EcalPhase2Linearizer.h>

#include <CondFormats/EcalObjects/interface/EcalTPGLinearizationConst.h>
#include <CondFormats/EcalObjects/interface/EcalTPGPedestals.h>
#include <CondFormats/EcalObjects/interface/EcalTPGCrystalStatus.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

EcalPhase2Linearizer::EcalPhase2Linearizer(bool famos)
  : famos_(famos), init_(false),
    linConsts_(nullptr),
    peds_(nullptr),
    badXStatus_(nullptr)
{
}

EcalPhase2Linearizer::~EcalPhase2Linearizer(){
  if (init_) {
    for (int i=0;i<(int)vectorbadXStatus_.size();i++){
      delete vectorbadXStatus_[i];
    }
  }    
}

void EcalPhase2Linearizer::setParameters(uint32_t raw)
{

  std::cout << " EcalPhase2Linearizer::setParameters() " << std::endl;
  std::cout << " Raw data  " << raw << std::endl; 

}

int EcalPhase2Linearizer::doIt()
{
  
  int output = 0;
  
  std::cout << " EcalPhase2Linearizer::doIt() output non bit shifted " << output << std::endl; 
  if(famos_ || output<0) return 0;
  
  return output;
}
 
int EcalPhase2Linearizer::setInput(const EcalLiteDTUSample &RawSam)
//int EcalPhase2Linearizer::setInput(const EcalMGPASample  &RawSam)
{

  std::cout << " EcalPhase2Linearizer::setInput() RawSam.raw() " << RawSam.raw() <<  std::endl;
  if (famos_)
    base_ = 200;  //FIXME by preparing a correct TPG.txt for Famos

  return 1;

}

void EcalPhase2Linearizer::process(const EBDigiCollectionPh2::Digi &df, std::vector<int> & output_percry)
{

  //We know a tower numbering is:                                                                                                                               
  // S1 S2 S3 S4 S5                                                                                                                                               
  //                                                                                                                                                              
  // 4  5  14 15 24                                                                                                                                               
  // 3  6  13 16 23                                                                                                                                               
  // 2  7  12 17 22                                                                                                                                               
  // 1  8  11 18 21                                                                                                                                               
  // 0  9  10 19 20                                                                                                                                               
                                                                                                                                                                
  std::cout << " EcalPhase2Linearizer::process(const  .. DataFrame size  " << df.size() <<  std::endl;                                                          
  for (int i=0;i<df.size();i++) {                                                                                                                               
    std::cout <<  df[i] << " ";                                                                                                                                 
  }                                                                                                                                                             
                                                                                                                                                                
  for (int i=0;i<df.size();i++) {                                                                                                                      
    EcalLiteDTUSample thisSample = df[i];         
    setInput(thisSample);                                                                                                                                      
    output_percry[i]=doIt();                                                                                                                                    
  }                                                                                                                                                             
                                                                                                                                                                
  std::cout << " EcalPhase2Linearizer::process(const  .. Final output " << std::endl;                                                                           
  for (int i=0;i<df.size();i++) {                                                                                                                               
    std::cout << " output_percry " << output_percry[i]<< " ";                                                                                                    
  }                                                                                                                                                             
                                                                                                                                                                
  return;           

}
