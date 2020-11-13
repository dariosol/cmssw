#ifndef ECAL_PHASE2_LINEARIZER_H
#define ECAL_PHASE2_LINEARIZER_H


#include "DataFormats/EcalDigi/interface/EcalLiteDTUSample.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include <CondFormats/EcalObjects/interface/EcalTPGPedestals.h>
#include <CondFormats/EcalObjects/interface/EcalTPGLinearizationConst.h>
#include <CondFormats/EcalObjects/interface/EcalTPGCrystalStatus.h>

#include <vector> 

  /** 
   \class EcalPhase2Linearizer
   \brief Linearisation for Phase2 
   *  input: ??  bits  corresponding to input EBDataFrame
   *  output: ?? bits 
   *  
   */

  class EcalPhase2Linearizer  {


  private:
    bool famos_;
    int uncorrectedSample_;
    int gainID_;
    int base_;
    int mult_;
    int shift_;
    int strip_;
    bool init_;
    
    const EcalTPGLinearizationConstant  *linConsts_;
    const EcalTPGPedestal *peds_;
    const EcalTPGCrystalStatusCode *badXStatus_;
    
    std::vector<const EcalTPGCrystalStatusCode *> vectorbadXStatus_;
     	
    int setInput(const EcalLiteDTUSample  &RawSam) ;
    //int setInput(const EcalMGPASample  &RawSam) ;
 

    int doIt() ;


  public:
    EcalPhase2Linearizer(bool famos);
    virtual ~EcalPhase2Linearizer();

    void process(const EBDigiCollectionPh2::Digi &, std::vector<int>&); 
    void setParameters(uint32_t raw); /// to be filled 

};




#endif
