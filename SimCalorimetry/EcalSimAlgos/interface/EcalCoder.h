
#ifndef EcalSimAlgos_EcalCoder_h
#define EcalSimAlgos_EcalCoder_h 1

#include "CalibFormats/CaloObjects/interface/CaloTSamples.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
//#include "CondFormats/EcalObjects/interface/EcalLiteDTUPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalCorrelatedNoiseMatrix.h"

template<typename M> class CorrelatedNoisifier ;
class EcalMGPASample;
class EcalDataFrame;
class DetId;
class EcalLiteDTUPed;

#include<vector>

namespace CLHEP {
  class HepRandomEngine;
}


class EcalCoder
{
   public:
#warning hard coded sample size
      typedef CaloTSamples<float,10> EcalSamples ;
      
      typedef CorrelatedNoisifier<EcalCorrMatrix> Noisifier ;

      enum { NBITS         =   12 , // number of available bits
             MAXADC        = 4095,  // 2^12 -1,  adc max range
	     NGAINS        =    2   // number of electronic gains
      };

      /// ctor
      EcalCoder( bool        addNoise        ,
		 bool        isPhase2        , //////////////////////////////////////////////////////NEW
                 bool        PreMix1         ,
                 Noisifier* ebCorrNoise0     ,
                 Noisifier* ebCorrNoise1 = nullptr ) ; 

      /// dtor
      virtual ~EcalCoder() ;

      /// can be fetched every event from the EventSetup
      void setPedestals( const EcalPedestals* pedestals ) ;

      void setGainRatios( const EcalGainRatios* gainRatios ) ;

      void setGainRatiosPhase2( const EcalCATIAGainRatios* gainRatios ) ;

      void setFullScaleEnergy( double EBscale ,
			       double EEscale   ) ;

      void setIntercalibConstants( const EcalIntercalibConstantsMC* ical ) ; 
 

      /// from EcalSamples to EcalDataFrame
      virtual void analogToDigital( CLHEP::HepRandomEngine*,
                                    const EcalSamples& clf ,
                                    EcalDataFrame&     df    ) const;
 
   private:

      /// limit on the energy scale due to the electronics range
      double fullScaleEnergy( const DetId & did ) const ;

      /// produce the pulse-shape
      void encode( const EcalSamples& ecalSamples , 
                   EcalDataFrame&     df,
                   CLHEP::HepRandomEngine* ) const ;

      

      void findPedestal( const DetId& detId    , 
			 int          gainId   , 
			 double&      pedestal ,
			 double&      width      ) const ;
    
      void findGains( const DetId& detId, float theGains[] ) const ;

      void findIntercalibConstant( const DetId& detId ,
				   double&      icalconst ) const ;
   
      const EcalPedestals* m_peds ;
      
      const EcalGainRatios* m_gainRatios ; // the electronics gains //////////////NEW
      const EcalCATIAGainRatios* m_gainRatiosPhase2 ; // the electronics gains for Phase2

      const EcalIntercalibConstantsMC* m_intercals ; //record specific for simulation of gain variation in MC

      double m_maxEneEB ; // max attainable energy in the ecal barrel
      
      bool m_addNoise ;   // whether add noise to the pedestals and the gains
      bool m_isPhase2 ; /////////////////////////////////////////////////////////////////////////////////////NEW
      bool m_PreMix1 ;   // Follow necessary steps for PreMixing input

      const Noisifier* m_ebCorrNoise[NGAINS] ;


};

#endif
