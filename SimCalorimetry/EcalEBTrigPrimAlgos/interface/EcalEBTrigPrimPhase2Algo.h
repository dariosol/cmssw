#ifndef EcalEBTrigPrimPhase2Algo_h
#define EcalEBTrigPrimPhase2Algo_h
/** \class EcalEBTrigPrimPhase2Algo
\author N. Marinelli - Univ. of Notre Dame
 * forPhase II 
 * While the new digitization is not yet implemented, we use the old Digis to make TP per crystal
 *
 ************************************************************/
#include <sys/time.h>
#include <iostream>
#include <vector>

#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "DataFormats/EcalDetId/interface/EcalTriggerElectronicsId.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame_Ph2.h"

#include "SimCalorimetry/EcalEBTrigPrimAlgos/interface/EcalPhase2Linearizer.h"


#include <map>
#include <utility>


class EcalTrigTowerDetId;
class ETPCoherenceTest;
class EcalEBTriggerPrimitiveSample;
class CaloSubdetectorGeometry;
class EBDataFrame_Ph2;


 
class EcalEBTrigPrimPhase2Algo
{  
 public:
  
  explicit EcalEBTrigPrimPhase2Algo(const edm::EventSetup & setup, int nSamples, int binofmax, bool tcpFormat, bool barrelOnly, bool debug, bool famos);

  virtual ~EcalEBTrigPrimPhase2Algo();

  void run(const edm::EventSetup &, const EBDigiCollectionPh2 *col, EcalEBTrigPrimDigiCollection & result);

  /*  
  void setPointers(const EcalTPGLinearizationConst *ecaltpLin,
		   const EcalTPGPedestals *ecaltpPed,
		   const EcalTPGCrystalStatus * ecaltpgBadX,
		   const EcalTPGWeightIdMap *ecaltpgWeightMap,
		   const EcalTPGWeightGroup *ecaltpgWeightGroup,
		   const EcalTPGSlidingWindow* ecaltpgSlidW,
		   const EcalTPGLutGroup *ecaltpgLutGroup,
		   const EcalTPGLutIdMap *ecaltpgLut,
		   const EcalTPGTowerStatus *ecaltpgBadTT,
		   const EcalTPGSpike * ecaltpgSpike )  
{
    ecaltpPed_=ecaltpPed;
    ecaltpLin_=ecaltpLin;
    ecaltpgBadX_=ecaltpgBadX;
    ecaltpgWeightMap_=ecaltpgWeightMap;
    ecaltpgWeightGroup_=ecaltpgWeightGroup;
    ecaltpgSlidW_=ecaltpgSlidW;
    ecaltpgLutGroup_=ecaltpgLutGroup;
    ecaltpgLut_=ecaltpgLut;
    ecaltpgBadTT_=ecaltpgBadTT;
    ecaltpgSpike_= ecaltpgSpike;
  }
  */

  void setPointers() {;}

 private:

  void init(const edm::EventSetup&);
  template <class T>  
    void initStructures(std::vector<std::vector<std::pair<int,std::vector<T> > > > & towMap);
  template <class T> 
    void clean(std::vector<std::vector<std::pair<int,std::vector<T> > > > &towerMap);


  void fillMap(EBDigiCollectionPh2 const * col, std::vector<std::vector<std::pair<int,std::vector<EBDigiCollectionPh2::Digi> > >  > &towerMap);
  //  void fillMap(EBDigiCollectionPh2 const * col, std::vector<std::vector<std::pair<int,std::vector<EBDataFrame_Ph2> > >  > &towerMap);
 


  int findStripNr(const EBDetId &id);
  
  int findStripNr(const EEDetId &id);

  // FIXME: temporary until hashedIndex works alsom for endcap
  int getIndex(const  EBDigiCollectionPh2 *, EcalTrigTowerDetId& id) {return id.hashedIndex();}
  // mind that eta is continuous between barrel+endcap
  //  int getIndex(const  EEDigiCollectionPh2 *, EcalTrigTowerDetId& id) {
  // int ind=(id.ietaAbs()-18)*72 + id.iphi();
  // if (id.zside()<0) ind+=792;
  // return ind;
  // }

  edm::ESHandle<EcalTrigTowerConstituentsMap> eTTmap_;
  //  const CaloSubdetectorGeometry *theEndcapGeometry;
  edm::ESHandle<CaloGeometry> theGeometry;


  float threshold;
  int nSamples_;
  int binOfMaximum_;
  int maxNrSamples_;


  bool tcpFormat_;
  bool barrelOnly_;
  bool debug_;  
  bool famos_; 


  int nrTowers_;   // nr of towers found by fillmap method
  static const unsigned int maxNrTowers_;
  static const unsigned int maxNrSamplesOut_; 
  static const unsigned int nrSamples_; 

  // data structures kept during the whole run
  std::vector<std::vector<int> > striptp_;
  //std::vector<std::vector<std::pair<int,std::vector<EBDataFrame> > > > towerMapEB_;
  std::vector<std::vector<std::pair<int,std::vector<EBDataFrame_Ph2> > > > towerMapEB_; 
  std::vector<std::vector<std::pair<int,std::vector<EEDataFrame> > > > towerMapEE_;
  std::vector<std::pair<int,EcalTrigTowerDetId> > hitTowers_;
  std::vector<EcalEBTriggerPrimitiveSample> towtp_;
  std::vector<EcalEBTriggerPrimitiveSample> towtp2_;

  enum {nbMaxStrips_=5};
  enum {nbMaxXtals_=5};

  const EcalElectronicsMapping* theMapping_;


  std::vector <EcalPhase2Linearizer *> linearizer_;
  //  EcalPhase2AmplitudeReconstruction *amplitude_reconstructor_; 
  //EcalPhase2TimeReconstruction *time_reconstructor_;
  //EcalPhase2SpikeFinder *spike_finder_;
  
  //
  /*
  const EcalTPGPedestals * ecaltpPed_;
  const EcalTPGLinearizationConst *ecaltpLin_;
  const EcalTPGCrystalStatus *ecaltpgBadX_;
  const EcalTPGWeightIdMap *ecaltpgWeightMap_;
  const EcalTPGWeightGroup *ecaltpgWeightGroup_;
  const EcalTPGSlidingWindow *ecaltpgSlidW_;
  const EcalTPGLutGroup *ecaltpgLutGroup_;
  const EcalTPGLutIdMap *ecaltpgLut_;
  const EcalTPGTowerStatus *ecaltpgBadTT_;
  const EcalTPGSpike * ecaltpgSpike_;   
  */

  EcalPhase2Linearizer *getLinearizer (int i) const { return linearizer_[i];}
  std::vector<std::vector<int> > lin_out_;
  //
  //EcalPhase2AmplitudeReconstruction *getAmplitudeFinder() const { return amplitude_reconstructor_;}
  std::vector<int> filt_out_;
  std::vector<int> peak_out_;
  std::vector<int> format_out_;
  // these two are dummy
  std::vector<int> fgvb_out_;
  std::vector<int> fgvb_out_temp_;

  //EcalPhase2TimeReconstruction *getTimeFinder() const {return time_reconstructor_};
  //EcalPhase2SpikeFinder *getSpikeFinder() const {return spike_finder};
  

  //

  std::vector<int> tcpformat_out_;



};


template <class T> 
void EcalEBTrigPrimPhase2Algo::clean( std::vector<std::vector<std::pair<int,std::vector<T> > > > & towMap) {  
  // clean internal data structures
  for (unsigned int i=0;i<maxNrTowers_;++i) 
    for (int j=0;j<nbMaxStrips_ ;++j) (towMap[i])[j].first=0;
  return;
}



void EcalEBTrigPrimPhase2Algo::fillMap(EBDigiCollectionPh2 const * col, std::vector<std::vector<std::pair<int,std::vector<EBDigiCollectionPh2::Digi> >  > >  &towerMap)
//void EcalEBTrigPrimPhase2Algo::fillMap(EBDigiCollectionPh2 const * col, std::vector<std::vector<std::pair<int,std::vector<EBDataFrame_Ph2> >  > >  &towerMap)



{
  // typedef typename Coll::Digi Digi;

  // implementation for Barrel 
  if (col) {
    nrTowers_=0;
    if ( debug_) std::cout  <<"Fill mapping, Collection size = "<< col->size() << std::endl;;
    for(unsigned int i = 0; i < col->size() ; ++i) {
      EBDigiCollectionPh2::Digi samples((*col)[i]);
      //EcalDataFrame_Ph2  samples((*col)[i]); 
      EcalTrigTowerDetId coarser=(*eTTmap_).towerOf(samples.id());
      int index=getIndex(col,coarser);
      EBDetId id=samples.id();
      int stripnr=findStripNr(id);

      int filled=0;
      for (unsigned int ij=0;ij<towerMap[index].size();++ij) filled+=towerMap[index][ij].first;
      if (!filled) {
	hitTowers_[nrTowers_++]=std::pair <int,EcalTrigTowerDetId>(index,coarser);
      }

      //FIXME: temporary protection
      int ncryst=towerMap[index][stripnr-1].first;
      if (ncryst>=nbMaxXtals_ ) {
        edm::LogError("EcalTrigPrimFunctionAlgo")<<"! Too many xtals for TT "<<coarser<<" stripnr "<<stripnr<<" xtalid "<<samples.id() ;
	continue;
      }
      ((towerMap[index])[stripnr-1].second)[ncryst]=samples;
      (towerMap[index])[stripnr-1].first++;
    }
  
    if (debug_) std::cout << "fillMap"<<"[EcalTrigPrimFunctionalAlgo] (found " 
			  << col->size() << " frames in "<< towerMap.size() << " towers) " << std::endl;
  }
  else {
    if (debug_) std::cout <<"FillMap - FillMap Collection size=0 !!!!" << std::endl;;
  }
}

template <class T> 
void EcalEBTrigPrimPhase2Algo::initStructures( std::vector<std::vector<std::pair<int,std::vector<T> > > > & towMap) {  
  //initialise internal data structures

  std::vector <T> vec0(nbMaxXtals_ );
  std::vector<std::pair<int,std::vector<T> > > vec1(nbMaxStrips_);
  for (int i=0;i<nbMaxStrips_ ;++i) vec1[i]=std::pair<int,std::vector<T> >(0,vec0);
  towMap.resize(maxNrTowers_); 
  for (unsigned int i=0;i<maxNrTowers_;++i) towMap[i]=vec1;
  
  std::vector<int> vecint(maxNrSamples_);
  striptp_.resize(nbMaxStrips_);
  for (int i=0;i<nbMaxStrips_;++i) striptp_[i]=vecint;
  
}



#endif
