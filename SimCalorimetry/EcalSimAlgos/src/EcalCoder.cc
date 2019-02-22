#include "SimCalorimetry/EcalSimAlgos/interface/EcalCoder.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimGeneral/NoiseGenerators/interface/CorrelatedNoisifier.h"
#include "DataFormats/EcalDigi/interface/EcalMGPASample.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"

#include <iostream>

//#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"


EcalCoder::EcalCoder( bool                  addNoise     , 
		      bool                  PreMix1      ,
		      EcalCoder::Noisifier* ebCorrNoise0 ,
		      EcalCoder::Noisifier* ebCorrNoise1 ) :
   m_peds        (           nullptr ) ,
   m_gainRatios  (           nullptr ) ,
   m_intercals   (           nullptr ) ,
   m_maxEneEB    (      1668.3 ) , // 4095(MAXADC)*12(gain 2)*0.035(GeVtoADC)*0.97
   m_addNoise    ( addNoise    ) ,
   m_PreMix1     ( PreMix1     ) 
   {
   m_ebCorrNoise[0] = ebCorrNoise0 ;
   assert( nullptr != m_ebCorrNoise[0] ) ;
   m_ebCorrNoise[1] = ebCorrNoise1 ;

}  

EcalCoder::~EcalCoder()
{
}

void 
EcalCoder::setFullScaleEnergy( double EBscale ,
                               double EEscale   )
{
   m_maxEneEB = EBscale ;
   
}


void  
EcalCoder::setPedestals( const EcalPedestals* pedestals ) 
{
   m_peds = pedestals ;
}

void  
EcalCoder::setGainRatios( const EcalCATIAGainRatios* gainRatios ) 
{
   m_gainRatios = gainRatios ; 
}

void 
EcalCoder::setIntercalibConstants( const EcalIntercalibConstantsMC* ical ) 
{
   m_intercals = ical ;
}

double 
EcalCoder::fullScaleEnergy( const DetId & detId ) const 
{
    //return detId.subdetId() == EcalBarrel ? m_maxEneEB : m_maxEneEE ;
    return m_maxEneEB ;
}

void 
EcalCoder::analogToDigital( CLHEP::HepRandomEngine* engine,
                            const EcalSamples& clf ,
                            EcalDataFrame&     df    ) const 
{
   df.setSize( clf.size() ) ;
   encode( clf, df, engine );

}

void 
EcalCoder::encode( const EcalSamples& ecalSamples , 
                   EcalDataFrame&     df,
                   CLHEP::HepRandomEngine* engine ) const
{
   assert( nullptr != m_peds ) ;

   const unsigned int csize ( ecalSamples.size() ) ;

  
   DetId detId = ecalSamples.id();             
   double Emax = fullScaleEnergy(detId);       

   //....initialisation
   if ( ecalSamples[5] > 0. ) LogDebug("EcalCoder") << "Input caloSample" << "\n" << ecalSamples;
  
   double LSB[NGAINS];
   double pedestals[NGAINS];
   double widths[NGAINS];
   float  gains[NGAINS];
   double trueRMS[NGAINS];

   double icalconst = 1. ;
   findIntercalibConstant( detId, icalconst );

   for( unsigned int igain ( 0 ); igain < NGAINS ; ++igain ) 
   {
      // fill in the pedestal and width
      findPedestal( detId ,
                    igain , 
                    pedestals[igain] , 
                    widths[igain]      ) ;

     
      trueRMS[igain] = std::sqrt( widths[igain]*widths[igain] - 1./12. ) ;

      // set nominal value first
      findGains( detId , gains  );               

      LSB[igain]= Emax/(MAXADC*gains[igain]);
      
   }

   CaloSamples noiseframe[] = { CaloSamples( detId , csize ) ,
                                CaloSamples( detId , csize ) ,
   } ;

   const Noisifier* noisy[NGAINS] = { m_ebCorrNoise[0],m_ebCorrNoise[1]} ;

   if( m_addNoise ){

       #warning noise generation to be checked
       noisy[0]->noisify( noiseframe[0], engine ) ; // high gain
       //if( nullptr == noisy[1] ) noisy[0]->noisify( noiseframe[1] ,
       //                                             engine,
       //                                              &noisy[0]->vecgau() ) ; // lwo gain
       noisy[1]->noisify( noiseframe[1], engine ) ; // low gain
   }


   bool isSaturated[]={false,false};
   std::vector<int> adctrace(csize);
   // fill ADC trace in gain 0 (x10) and gain 1 (x1)
   for (unsigned int igain=0; igain<NGAINS; ++igain){
      for( unsigned int i ( 0 ) ; i != csize ; ++i ){
       
	      double asignal =0;
          if (!m_PreMix1){
       
             asignal = pedestals[igain] +
                       ecalSamples[i] /( LSB[igain]*icalconst ) +
                       trueRMS[igain]*noiseframe[igain][i]    ;
             
          } else {
             //  no noise nor pedestal when premixing
             asignal = ecalSamples[i] /( LSB[igain]*icalconst ) ;
            
          }
          int isignal =  asignal ;
          int adc =  asignal - (double) isignal < 0.5 ? isignal : isignal + 1;
          if (adc>MAXADC) {
              adc = MAXADC;
              isSaturated[i] = true;
          }
          if (isSaturated[0] && igain==0) break; // gain 0 (x10) channel is saturated, readout will use gain 1 (x1)
          else adctrace[i] = adc;
      } // for adc
   
      if (!isSaturated[0]) break; //  gain 0 (x10) is not saturated, so don't bother with gain 1
   } // for igain

   int igain =0;
   if (isSaturated[0] ) igain =1;

   // Note: we assume that Pileup generates small signals, and we will not saturate when adding pedestals

   for ( unsigned int j =0; j < ecalSamples.size(); ++j ) {
      df.setSample(j, EcalMGPASample(adctrace[j], igain));   
   }
}
                                      
                                      


void 
EcalCoder::findPedestal( const DetId & detId  , 
			 int           gainId , 
			 double&       ped    , 
			 double&       width     ) const
{


   EcalPedestalsMap::const_iterator itped = m_peds->getMap().find( detId );
   ped   = (*itped).mean(gainId);
   width = (*itped).rms(gainId);
  
   if ( (detId.subdetId() != EcalBarrel) && (detId.subdetId() != EcalEndcap) ) 
   { 
      edm::LogError("EcalCoder") << "Could not find pedestal for " << detId.rawId() << " among the " << m_peds->getMap().size();
   } 


   LogDebug("EcalCoder") << "Pedestals for " << detId.rawId() << " gain range " << gainId << " : \n" << "Mean = " << ped << " rms = " << width;
}


void 
EcalCoder::findGains( const DetId & detId , float Gains[]        ) const
{

   EcalCATIAGainRatioMap::const_iterator grit = m_gainRatios->getMap().find( detId );

   Gains[1] = 1.;
   Gains[0] = Gains[1]*(*grit);
   
  
   if ( (detId.subdetId() != EcalBarrel) && (detId.subdetId() != EcalEndcap) ) 
   { 
      edm::LogError("EcalCoder") << "Could not find gain ratios for " << detId.rawId() << " among the " << m_gainRatios->getMap().size();
   }   
  
}

void 
EcalCoder::findIntercalibConstant( const DetId& detId, 
				   double&      icalconst ) const
{
   EcalIntercalibConstantMC thisconst = 1.;
   // find intercalib constant for this xtal
   const EcalIntercalibConstantMCMap &icalMap = m_intercals->getMap();
   EcalIntercalibConstantMCMap::const_iterator icalit = icalMap.find(detId);
   if( icalit!=icalMap.end() )
   {
      thisconst = (*icalit);
      if ( icalconst == 0. ) { thisconst = 1.; }
   } 
   else
   {
      edm::LogError("EcalCoder") << "No intercalib const found for xtal " << detId.rawId() << "! something wrong with EcalIntercalibConstants in your DB? ";
   }
   icalconst = thisconst;
}
