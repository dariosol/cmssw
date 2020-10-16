#ifndef RecoLocalCalo_EcalRecProducers_EcalUncalibRecHitRecWorkerGlobal_hh
#define RecoLocalCalo_EcalRecProducers_EcalUncalibRecHitRecWorkerGlobal_hh

/** \class EcalUncalibRecHitRecGlobalAlgo                                                                                                                                           
 *  Template used to compute amplitude, pedestal using a weights method                                                                                                            
 *                           time using a ratio method                                                                                                                             
 *                           chi2 using express method  
 *
 *
 */

#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerBaseClass.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitTimeWeightsAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecChi2Algo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRatioMethodAlgo.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalTimeOffsetConstant.h"
#include "CondFormats/EcalObjects/interface/EcalLiteDTUPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalWeightXtalGroups.h"
#include "CondFormats/EcalObjects/interface/EcalTBWeights.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"
#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"
#include "CondFormats/EcalObjects/interface/EcalPhase2SamplesCorrelation.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm

class EcalPhase2UncalibRecHitWorkerMultiFit final : public EcalUncalibRecHitWorkerBaseClass {
public:
  EcalPhase2UncalibRecHitWorkerMultiFit(const edm::ParameterSet&, edm::ConsumesCollector& c);

private:
  void set(const edm::EventSetup& es) override;
  void set(const edm::Event& evt) override;
  void run(const edm::Event& evt, const EcalDigiCollectionPh2& digis, EcalUncalibratedRecHitCollection& result) ;

public:
  edm::ParameterSetDescription getAlgoDescription() override;

private:
  edm::ESHandle<EcalLiteDTUPedestals> peds;
  edm::ESHandle<EcalCATIAGainRatios> gains;
  edm::ESHandle<EcalPhase2SamplesCorrelation> noisecovariances;
  edm::ESHandle<EcalPhase2PulseShapes> pulseshapes;
  edm::ESHandle<EcalPhase2PulseCovariances> pulsecovariances;

  double timeCorrection(float ampli, const std::vector<float>& amplitudeBins, const std::vector<float>& shiftBins);

  const ecalph2::SampleMatrix& noisecor(int gain) const { return noisecors_[gain]; }
  const ecalph2::SampleMatrixGainArray& noisecor() const { return noisecors_; }

  // multifit method
  ecalph2::SampleMatrixGainArray noisecors_;

  BXVector activeBX;
  bool ampErrorCalculation_;
  bool useLumiInfoRunHeader_;
  EcalUncalibRecHitMultiFitAlgo multiFitMethod_;

  int bunchSpacingManual_;
  edm::EDGetTokenT<unsigned int> bunchSpacing_;

  // determine which of the samples must actually be used by ECAL local reco
  edm::ESHandle<EcalSampleMask> sampleMaskHand_;

  // time algorithm to be used to set the jitter and its uncertainty
  enum TimeAlgo { noMethod, ratioMethod, weightsMethod };
  TimeAlgo timealgo_ = noMethod;

  // time weights method
  edm::ESHandle<EcalWeightXtalGroups> grps;
  edm::ESHandle<EcalTBWeights> wgts;
  const EcalWeightSet::EcalWeightMatrix* weights[2];
  EcalUncalibRecHitTimeWeightsAlgo<EBDataFrame> weightsMethod_barrel_;
  
  bool doPrefitEB_;
  double prefitMaxChiSqEB_;
  bool dynamicPedestalsEB_;
  bool mitigateBadSamplesEB_;
  bool gainSwitchUseMaxSampleEB_;
  bool selectiveBadSampleCriteriaEB_;
  double addPedestalUncertaintyEB_;
  bool simplifiedNoiseModelForGainSwitch_;

  // ratio method
  std::vector<double> EBtimeFitParameters_;
  std::vector<double> EBamplitudeFitParameters_;
  std::pair<double, double> EBtimeFitLimits_;
  
  EcalUncalibRecHitRatioMethodAlgo<EcalDataFrame_Ph2> ratioMethod_barrel_;
  

  double EBtimeConstantTerm_;
  double EBtimeNconst_;
  double outOfTimeThreshG10pEB_;
  double outOfTimeThreshG1pEB_;
  double outOfTimeThreshG10mEB_;
  double outOfTimeThreshG1mEB_;
  double amplitudeThreshEB_;
  double ebSpikeThresh_;

  edm::ESHandle<EcalTimeBiasCorrections> timeCorrBias_;

  edm::ESHandle<EcalTimeCalibConstants> itime;
  edm::ESHandle<EcalTimeOffsetConstant> offtime;
  std::vector<double> ebPulseShape_;


  // chi2 thresholds for flags settings
  bool kPoorRecoFlagEB_;
  double chi2ThreshEB_;

};

#endif
