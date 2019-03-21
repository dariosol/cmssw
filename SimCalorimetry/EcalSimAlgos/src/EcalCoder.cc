#include "SimCalorimetry/EcalSimAlgos/interface/EcalCoder.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimGeneral/NoiseGenerators/interface/CorrelatedNoisifier.h"
#include "DataFormats/EcalDigi/interface/EcalMGPASample.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"

#include <iostream>

//#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"


EcalCoder::EcalCoder( bool                  addNoise     , 
		      bool                  isPhase2     , ////////////////////////////////////////////////////////////NEW
		      bool                  PreMix1      ,
		      EcalCoder::Noisifier* ebCorrNoise0 ,
		      EcalCoder::Noisifier* ebCorrNoise1 ) :
   m_peds        (           nullptr ) ,
   m_gainRatios  (           nullptr ) ,
   m_gainRatiosPhase2  (     nullptr ) ,
   m_intercals   (           nullptr ) ,
   m_maxEneEB    (      1668.3 ) , // 4095(MAXADC)*12(gain 2)*0.035(GeVtoADC)*0.97
  
   m_addNoise    ( addNoise    ) ,
   m_isPhase2    ( isPhase2    ) ,
   //   m_PreMix1     ( PreMix1     ) 
   //m_PreMix1     ( 1     ) 
   //{
   //m_maxEneEB    (      2000. ) , // Maximum for CATIA: LSB gain 10: 0.048 MeV
   //m_addNoise    ( addNoise    ) ,
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
  //   m_maxEneEB = EBscale ;
  m_maxEneEB = 2000. ; //I don 't know where is setFullScaleEnergy first call
   
}


void  
EcalCoder::setPedestals( const EcalPedestals* pedestals ) 
{
   m_peds = pedestals ;
}

void  
EcalCoder::setGainRatios( const EcalGainRatios* gainRatios ) 
{
   m_gainRatios = gainRatios ; 
}

void  
EcalCoder::setGainRatiosPhase2( const EcalCATIAGainRatios* gainRatiosPhase2 ) 
{
   m_gainRatiosPhase2 = gainRatiosPhase2 ; 
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
  

   //N Gains set to 2 in the .h
   double pedestals[NGAINS];
   double widths[NGAINS];
   float  gains[NGAINS];
   double LSB[NGAINS];
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
      //I insert an absolute value in the trueRMS
      trueRMS[igain] = std::sqrt( std::fabs(widths[igain]*widths[igain] - 1./12.) ) ;

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

   for (unsigned int igain=0; igain<NGAINS; ++igain) {

     for( unsigned int i ( 0 ) ; i != csize ; ++i ) {
       
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
       if(ecalSamples[i]>0){
	 std::cout<<"NGAIN "<<igain<<std::endl;
	 std::cout<<"icalconst "<<icalconst<<std::endl;
	 std::cout<<"Emax "<<Emax<<std::endl;
	 for(int j=0;j<2;++j) {
	   std::cout<<"index "<<j<<"pedestals "<<pedestals[j]<<std::endl;
	   std::cout<<"index "<<j<<"widths "<<widths[j]<<std::endl;
	 }
	 

	 std::cout<<"trueRMS "<<trueRMS[igain]<<std::endl;
	 std::cout<<"LSB Gain"<<igain<<" "<<LSB[igain]<<"\n";
	 std::cout<<"Sample " <<i<< " NoiseFrame: " << noiseframe[igain][i] <<"\n";
	 std::cout<<"Sample " <<i<< " Analog Samples: " << ecalSamples[i] <<"\n";
	 std::cout<<"Sample " <<i<< " ADC: " << (unsigned int)adc <<"\n";
	    

       }
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
EcalCoder::findGains( const DetId & detId , float Gains[]       ) const
{
   
   if(m_isPhase2){
     EcalCATIAGainRatioMap::const_iterator grit = m_gainRatiosPhase2->getMap().find( detId );
     Gains[1] = 1.;
     Gains[0] = Gains[1]*(*grit);
   }

   else{
     EcalGainRatioMap::const_iterator grit = m_gainRatios->getMap().find( detId );
     Gains[0] = 0.;
     Gains[3] = 1.;
     Gains[2] = (*grit).gain6Over1();
     Gains[1] = Gains[2]*( (*grit).gain12Over6() );  
   }
  
   if ( (detId.subdetId() != EcalBarrel) && (detId.subdetId() != EcalEndcap) ) 
   { 
      edm::LogError("EcalCoder") << "Could not find gain ratios for " << detId.rawId() << " among the " << m_gainRatios->getMap().size();
   }   
  
}

// void 
// EcalCoder::findGains( const DetId & detId , float Gains[]   ) const
// {

//    EcalCATIAGainRatioMap::const_iterator grit = m_gainRatios->getMap().find( detId );

//    Gains[1] = 1.;
//    Gains[0] = Gains[1]*(*grit);
   
  
//    if ( (detId.subdetId() != EcalBarrel) && (detId.subdetId() != EcalEndcap) ) 
//    { 
//       edm::LogError("EcalCoder") << "Could not find gain ratios for " << detId.rawId() << " among the " << m_gainRatios->getMap().size();
//    }   
  
// }

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
