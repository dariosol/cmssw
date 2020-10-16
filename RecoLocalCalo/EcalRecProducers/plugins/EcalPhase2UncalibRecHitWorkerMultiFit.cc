#include "RecoLocalCalo/EcalRecProducers/plugins/EcalPhase2UncalibRecHitWorkerMultiFit.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/DataRecord/interface/EcalCATIAGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalLiteDTUPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalWeightXtalGroupsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTBWeightsRcd.h"
#include "CondFormats/DataRecord/interface/EcalSampleMaskRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeOffsetConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeBiasCorrectionsRcd.h"
#include "CondFormats/DataRecord/interface/EcalPhase2SamplesCorrelationRcd.h"
#include "CondFormats/DataRecord/interface/EcalPhase2PulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPhase2PulseCovariancesRcd.h"

#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include <FWCore/ParameterSet/interface/EmptyGroupDescription.h>

EcalPhase2UncalibRecHitWorkerMultiFit::EcalPhase2UncalibRecHitWorkerMultiFit(const edm::ParameterSet& ps, edm::ConsumesCollector& c)
    : EcalUncalibRecHitWorkerBaseClass(ps, c) {
  // get the BX for the pulses to be activated
  std::vector<int32_t> activeBXs = ps.getParameter<std::vector<int32_t>>("activeBXs");
  activeBX.resize(activeBXs.size());
  for (unsigned int ibx = 0; ibx < activeBXs.size(); ++ibx) {
    activeBX.coeffRef(ibx) = activeBXs[ibx];
  }

  // uncertainty calculation (CPU intensive)
  ampErrorCalculation_ = ps.getParameter<bool>("ampErrorCalculation");
  useLumiInfoRunHeader_ = ps.getParameter<bool>("useLumiInfoRunHeader");

  if (useLumiInfoRunHeader_) {
    bunchSpacing_ = c.consumes<unsigned int>(edm::InputTag("bunchSpacingProducer"));
    bunchSpacingManual_ = 0;
  } else {
    bunchSpacingManual_ = ps.getParameter<int>("bunchSpacing");
  }

  doPrefitEB_ = ps.getParameter<bool>("doPrefitEB");
  
  prefitMaxChiSqEB_ = ps.getParameter<double>("prefitMaxChiSqEB");
 
  dynamicPedestalsEB_ = ps.getParameter<bool>("dynamicPedestalsEB");
  mitigateBadSamplesEB_ = ps.getParameter<bool>("mitigateBadSamplesEB");
  gainSwitchUseMaxSampleEB_ = ps.getParameter<bool>("gainSwitchUseMaxSampleEB");
  selectiveBadSampleCriteriaEB_ = ps.getParameter<bool>("selectiveBadSampleCriteriaEB");
  addPedestalUncertaintyEB_ = ps.getParameter<double>("addPedestalUncertaintyEB");
  
  simplifiedNoiseModelForGainSwitch_ = ps.getParameter<bool>("simplifiedNoiseModelForGainSwitch");

  // algorithm to be used for timing
  auto const& timeAlgoName = ps.getParameter<std::string>("timealgo");
  if (timeAlgoName == "RatioMethod")
    timealgo_ = ratioMethod;
  else if (timeAlgoName == "WeightsMethod")
    timealgo_ = weightsMethod;
  else if (timeAlgoName != "None")
    edm::LogError("EcalUncalibRecHitError") << "No time estimation algorithm defined";

  // ratio method parameters
  EBtimeFitParameters_ = ps.getParameter<std::vector<double>>("EBtimeFitParameters");
  
  EBamplitudeFitParameters_ = ps.getParameter<std::vector<double>>("EBamplitudeFitParameters");
  
  EBtimeFitLimits_.first = ps.getParameter<double>("EBtimeFitLimits_Lower");
  EBtimeFitLimits_.second = ps.getParameter<double>("EBtimeFitLimits_Upper");
  
  
  EBtimeConstantTerm_ = ps.getParameter<double>("EBtimeConstantTerm");
  
  EBtimeNconst_ = ps.getParameter<double>("EBtimeNconst");
  
  outOfTimeThreshG10pEB_ = ps.getParameter<double>("outOfTimeThresholdGain10pEB");
  outOfTimeThreshG10mEB_ = ps.getParameter<double>("outOfTimeThresholdGain10mEB");
  outOfTimeThreshG1pEB_ = ps.getParameter<double>("outOfTimeThresholdGain1pEB");
  outOfTimeThreshG1mEB_ = ps.getParameter<double>("outOfTimeThresholdGain1mEB");
  
  amplitudeThreshEB_ = ps.getParameter<double>("amplitudeThresholdEB");
  

  // spike threshold
  ebSpikeThresh_ = ps.getParameter<double>("ebSpikeThreshold");

  ebPulseShape_ = ps.getParameter<std::vector<double>>("ebPulseShape");
  

  // chi2 parameters for flags determination
  kPoorRecoFlagEB_ = ps.getParameter<bool>("kPoorRecoFlagEB");
  
  chi2ThreshEB_ = ps.getParameter<double>("chi2ThreshEB_");
  
}

void EcalPhase2UncalibRecHitWorkerMultiFit::set(const edm::EventSetup& es) {
  // common setup
  es.get<EcalCATIAGainRatiosRcd>().get(gains);
  es.get<EcalLiteDTUPedestalsRcd>().get(peds);

  // for the multifit method
  if (!ampErrorCalculation_)
    multiFitMethod_.disableErrorCalculation();
  es.get<EcalPhase2SamplesCorrelationRcd>().get(noisecovariances);
  es.get<EcalPhase2PulseShapesRcd>().get(pulseshapes);
  es.get<EcalPhase2PulseCovariancesRcd>().get(pulsecovariances);

  // weights parameters for the time
  es.get<EcalWeightXtalGroupsRcd>().get(grps);
  es.get<EcalTBWeightsRcd>().get(wgts);

  // which of the samples need be used
  es.get<EcalSampleMaskRcd>().get(sampleMaskHand_);

  // for the ratio method
  es.get<EcalTimeCalibConstantsRcd>().get(itime);
  es.get<EcalTimeOffsetConstantRcd>().get(offtime);

  // for the time correction methods
  es.get<EcalTimeBiasCorrectionsRcd>().get(timeCorrBias_);

  int nnoise = SampleVector::RowsAtCompileTime;
  ecalph2::SampleMatrix& noisecorEBg10 = noisecors_[0];
  ecalph2::SampleMatrix& noisecorEBg1 = noisecors_[1];

  

  for (int i = 0; i < nnoise; ++i) {
    for (int j = 0; j < nnoise; ++j) {
      int vidx = std::abs(j - i);

      noisecorEBg10(i, j) = noisecovariances->EBG10SamplesCorrelation[vidx];
      noisecorEBg1(i, j) = noisecovariances->EBG1SamplesCorrelation[vidx];
  
    }
  }
}

void EcalPhase2UncalibRecHitWorkerMultiFit::set(const edm::Event& evt) {
  unsigned int bunchspacing = 450;

  if (useLumiInfoRunHeader_) {
    edm::Handle<unsigned int> bunchSpacingH;
    evt.getByToken(bunchSpacing_, bunchSpacingH);
    bunchspacing = *bunchSpacingH;
  } else {
    bunchspacing = bunchSpacingManual_;
  }

  if (useLumiInfoRunHeader_ || bunchSpacingManual_ > 0) {
    if (bunchspacing == 25) {
      activeBX.resize(10);
      activeBX << -5, -4, -3, -2, -1, 0, 1, 2, 3, 4;
    } else {
      //50ns configuration otherwise (also for no pileup)
      activeBX.resize(5);
      activeBX << -4, -2, 0, 2, 4;
    }
  }
}

/**
 * Amplitude-dependent time corrections; 
 * EXtimeCorrAmplitudes (ADC) and EXtimeCorrShifts (ns) need to have the same number of elements
 * Bins must be ordered in amplitude. First-last bins take care of under-overflows.
 *
 *
 * @return Jitter (in clock cycles) which will be added to UncalibRechit.setJitter(), 0 if no correction is applied.
 */
double EcalPhase2UncalibRecHitWorkerMultiFit::timeCorrection(float ampli,
                                                       const std::vector<float>& amplitudeBins,
                                                       const std::vector<float>& shiftBins) {
  // computed initially in ns. Than turned in the BX's, as
  // EcalUncalibratedRecHit need be.
  double theCorrection = 0;

  // sanity check for arrays
  if (amplitudeBins.empty()) {
    edm::LogError("EcalRecHitError") << "timeCorrAmplitudeBins is empty, forcing no time bias corrections.";

    return 0;
  }

  if (amplitudeBins.size() != shiftBins.size()) {
    edm::LogError("EcalRecHitError") << "Size of timeCorrAmplitudeBins different from "
                                        "timeCorrShiftBins. Forcing no time bias corrections. ";

    return 0;
  }

  // FIXME? what about a binary search?
  int myBin = -1;
  for (int bin = 0; bin < (int)amplitudeBins.size(); bin++) {
    if (ampli > amplitudeBins[bin]) {
      myBin = bin;
    } else {
      break;
    }
  }

  if (myBin == -1) {
    theCorrection = shiftBins[0];
  } else if (myBin == ((int)(amplitudeBins.size() - 1))) {
    theCorrection = shiftBins[myBin];
  } else {
    // interpolate linearly between two assingned points
    theCorrection = (shiftBins[myBin + 1] - shiftBins[myBin]);
    theCorrection *= (((double)ampli) - amplitudeBins[myBin]) / (amplitudeBins[myBin + 1] - amplitudeBins[myBin]);
    theCorrection += shiftBins[myBin];
  }

  // convert ns into clocks
  constexpr double inv25 = 1. / 25.;
  return theCorrection * inv25;
}

void EcalPhase2UncalibRecHitWorkerMultiFit::run(const edm::Event& evt,
                                          const EcalDigiCollectionPh2& digis,
                                          EcalUncalibratedRecHitCollection& result) {
  if (digis.empty())
    return;

  // assume all digis come from the same subdetector (either barrel or endcap)
  DetId detid(digis.begin()->id());
  

  multiFitMethod_.setSimplifiedNoiseModelForGainSwitch(simplifiedNoiseModelForGainSwitch_);
  
  multiFitMethod_.setDoPrefit(doPrefitEB_);
  multiFitMethod_.setPrefitMaxChiSq(prefitMaxChiSqEB_);
  multiFitMethod_.setDynamicPedestals(dynamicPedestalsEB_);
  multiFitMethod_.setMitigateBadSamples(mitigateBadSamplesEB_);
  multiFitMethod_.setGainSwitchUseMaxSample(gainSwitchUseMaxSampleEB_);
  multiFitMethod_.setSelectiveBadSampleCriteria(selectiveBadSampleCriteriaEB_);
  multiFitMethod_.setAddPedestalUncertainty(addPedestalUncertaintyEB_);
  

  FullSampleVector fullpulse(FullSampleVector::Zero());
  FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());

  result.reserve(result.size() + digis.size());
  for (auto itdg = digis.begin(); itdg != digis.end(); ++itdg) {
    DetId detid(itdg->id());

    const EcalSampleMask* sampleMask_ = sampleMaskHand_.product();

    // intelligence for recHit computation
    float offsetTime = 0;

    //    const EcalLiteDTUPedestals* aped = nullptr;
    const EcalCATIAGainRatio* aGain = nullptr;
    const EcalXtalGroupId* gid = nullptr;
    const EcalPhase2PulseShapes::Item* aPulse = nullptr;
    const EcalPhase2PulseCovariances::Item* aPulseCov = nullptr;

    const EcalLiteDTUPedestalsMap* aped;
    unsigned int hashedIndex = EBDetId(detid).hashedIndex();
    EcalLiteDTUPedestalsMap::const_iterator itped = aped->getMap().find(hashedIndex);
    aGain = &gains->barrel(hashedIndex);
    gid = &grps->barrel(hashedIndex);
    aPulse = &pulseshapes->barrel(hashedIndex);
    aPulseCov = &pulsecovariances->barrel(hashedIndex);
    offsetTime = offtime->getEBValue();
    
    double pedVec[ecalPh2::NGAINS] = {(*itped).mean(ecalPh2::gainId10), (*itped).mean(ecalPh2::gainId1)};
    double pedRMSVec[ecalPh2::NGAINS] = {(*itped).rms(ecalPh2::gainId10),  (*itped).rms(ecalPh2::gainId1)};
    double gainRatios[ecalPh2::NGAINS] = {1., *aGain};

    for (unsigned int i = 0; i < EcalPhase2PulseShape::TEMPLATESAMPLES; ++i)
      fullpulse(i + 7) = aPulse->pdfval[i];

    for (unsigned int i = 0; i < EcalPhase2PulseShape::TEMPLATESAMPLES; i++)
      for (unsigned int j = 0; j < EcalPhase2PulseShape::TEMPLATESAMPLES; j++)
        fullpulsecov(i + 7, j + 7) = aPulseCov->covval[i][j];

    // compute the right bin of the pulse shape using time calibration constants
    EcalTimeCalibConstantMap::const_iterator it = itime->find(detid);
    EcalTimeCalibConstant itimeconst = 0;
    if (it != itime->end()) {
      itimeconst = (*it);
    } else {
      edm::LogError("EcalRecHitError") << "No time intercalib const found for xtal " << detid.rawId()
                                       << "! something wrong with EcalTimeCalibConstants in your DB? ";
    }

    int lastSampleBeforeSaturation = -2;
    for (unsigned int iSample = 0; iSample < EcalDataFrame_Ph2::MAXSAMPLES; iSample++) {
      if (((EcalDataFrame_Ph2)(*itdg)).sample(iSample).gainId() == 0) {
        lastSampleBeforeSaturation = iSample - 1;
        break;
      }
    }

    // === amplitude computation ===

    if (lastSampleBeforeSaturation == 4) {  // saturation on the expected max sample
      result.emplace_back((*itdg).id(), 4095 * 12, 0, 0, 0);
      auto& uncalibRecHit = result.back();
      uncalibRecHit.setFlagBit(EcalUncalibratedRecHit::kSaturated);
      // do not propagate the default chi2 = -1 value to the calib rechit (mapped to 64), set it to 0 when saturation
      uncalibRecHit.setChi2(0);
    } else if (lastSampleBeforeSaturation >=
               -1) {  // saturation on other samples: cannot extrapolate from the fourth one
      int gainId = ((EcalDataFrame)(*itdg)).sample(5).gainId();
      if (gainId == 0)
        gainId = 3;
      auto pedestal = pedVec[gainId - 1];
      auto gainratio = gainRatios[gainId - 1];
      double amplitude = ((double)(((EcalDataFrame)(*itdg)).sample(5).adc()) - pedestal) * gainratio;
      result.emplace_back((*itdg).id(), amplitude, 0, 0, 0);
      auto& uncalibRecHit = result.back();
      uncalibRecHit.setFlagBit(EcalUncalibratedRecHit::kSaturated);
      // do not propagate the default chi2 = -1 value to the calib rechit (mapped to 64), set it to 0 when saturation
      uncalibRecHit.setChi2(0);
    } else {
      // multifit
        const ecalph2::SampleMatrixGainArray& noisecors = noisecor();

      result.push_back(multiFitMethod_.makeRecHit(*itdg, aped, aGain, noisecors, fullpulse, fullpulsecov, activeBX));
      auto& uncalibRecHit = result.back();

      // === time computation ===
      if (timealgo_ == ratioMethod) {
        // ratio method
        constexpr float clockToNsConstant = 25.;
        constexpr float invClockToNs = 1. / clockToNsConstant;
       
        ratioMethod_barrel_.init(*itdg, *sampleMask_, pedVec, pedRMSVec, gainRatios);
        ratioMethod_barrel_.fixMGPAslew(*itdg);
        ratioMethod_barrel_.computeTime(EBtimeFitParameters_, EBtimeFitLimits_, EBamplitudeFitParameters_);
        ratioMethod_barrel_.computeAmplitude(EBamplitudeFitParameters_);
        EcalUncalibRecHitRatioMethodAlgo<EBDataFrame>::CalculatedRecHit crh =
            ratioMethod_barrel_.getCalculatedRecHit();

        double theTimeCorrectionEB = timeCorrection(
                      uncalibRecHit.amplitude(), timeCorrBias_->EBTimeCorrAmplitudeBins, timeCorrBias_->EBTimeCorrShiftBins);

        uncalibRecHit.setJitter(crh.timeMax - 5 + theTimeCorrectionEB);
        uncalibRecHit.setJitterError(std::hypot(crh.timeError, EBtimeConstantTerm_ / clockToNsConstant));

        // consider flagging as kOutOfTime only if above noise
        if (uncalibRecHit.amplitude() > pedRMSVec[0] * amplitudeThreshEB_) {
            float outOfTimeThreshP = outOfTimeThreshG12pEB_;
            float outOfTimeThreshM = outOfTimeThreshG12mEB_;
            // determine if gain has switched away from gainId==1 (x12 gain)
            // and determine cuts (number of 'sigmas') to ose for kOutOfTime
            // >3k ADC is necessasry condition for gain switch to occur
            if (uncalibRecHit.amplitude() > 3000.) {
              for (int iSample = 0; iSample < EBDataFrame::MAXSAMPLES; iSample++) {
                int GainId = ((EcalDataFrame)(*itdg)).sample(iSample).gainId();
                if (GainId != 1) {
                  outOfTimeThreshP = outOfTimeThreshG61pEB_;
                  outOfTimeThreshM = outOfTimeThreshG61mEB_;
                  break;
                }
              }
            }
            float correctedTime = (crh.timeMax - 5) * clockToNsConstant + itimeconst + offsetTime;
            float cterm = EBtimeConstantTerm_;
            float sigmaped = pedRMSVec[0];  // approx for lower gains
            float nterm = EBtimeNconst_ * sigmaped / uncalibRecHit.amplitude();
            float sigmat = std::sqrt(nterm * nterm + cterm * cterm);
            if ((correctedTime > sigmat * outOfTimeThreshP) || (correctedTime < -sigmat * outOfTimeThreshM)) {
              uncalibRecHit.setFlagBit(EcalUncalibratedRecHit::kOutOfTime);
            }
          }
        }
      } else if (timealgo_ == weightsMethod) {
        //  weights method on the PU subtracted pulse shape
        std::vector<double> amplitudes;
        for (unsigned int ibx = 0; ibx < activeBX.size(); ++ibx)
          amplitudes.push_back(uncalibRecHit.outOfTimeAmplitude(ibx));

        EcalTBWeights::EcalTDCId tdcid(1);
        EcalTBWeights::EcalTBWeightMap const& wgtsMap = wgts->getMap();
        EcalTBWeights::EcalTBWeightMap::const_iterator wit;
        wit = wgtsMap.find(std::make_pair(*gid, tdcid));
        if (wit == wgtsMap.end()) {
          edm::LogError("EcalUncalibRecHitError")
              << "No weights found for EcalGroupId: " << gid->id() << " and  EcalTDCId: " << tdcid
              << "\n  skipping digi with id: " << detid.rawId();
          result.pop_back();
          continue;
        }
        const EcalWeightSet& wset = wit->second;  // this is the EcalWeightSet

        const EcalWeightSet::EcalWeightMatrix& mat1 = wset.getWeightsBeforeGainSwitch();
        const EcalWeightSet::EcalWeightMatrix& mat2 = wset.getWeightsAfterGainSwitch();

        weights[0] = &mat1;
        weights[1] = &mat2;

        double timerh;
        if (detid.subdetId() == EcalEndcap) {
          timerh = weightsMethod_endcap_.time(*itdg, amplitudes, aped, aGain, fullpulse, weights);
        } else {
          timerh = weightsMethod_barrel_.time(*itdg, amplitudes, aped, aGain, fullpulse, weights);
        }
        uncalibRecHit.setJitter(timerh);
        uncalibRecHit.setJitterError(0.);  // not computed with weights
      } else {                             // no time method;
        uncalibRecHit.setJitter(0.);
        uncalibRecHit.setJitterError(0.);
      }
    }

    // set flags if gain switch has occurred
    auto& uncalibRecHit = result.back();
    if (((EcalDataFrame)(*itdg)).hasSwitchToGain6())
      uncalibRecHit.setFlagBit(EcalUncalibratedRecHit::kHasSwitchToGain6);
    if (((EcalDataFrame)(*itdg)).hasSwitchToGain1())
      uncalibRecHit.setFlagBit(EcalUncalibratedRecHit::kHasSwitchToGain1);
  }
}

edm::ParameterSetDescription EcalPhase2UncalibRecHitWorkerMultiFit::getAlgoDescription() {
  edm::ParameterSetDescription psd0;
  psd0.addNode((edm::ParameterDescription<std::vector<double>>("EBPulseShapeTemplate",
                                                               {1.13979e-02,
                                                                7.58151e-01,
                                                                1.00000e+00,
                                                                8.87744e-01,
                                                                6.73548e-01,
                                                                4.74332e-01,
                                                                3.19561e-01,
                                                                2.15144e-01,
                                                                1.47464e-01,
                                                                1.01087e-01,
                                                                6.93181e-02,
                                                                4.75044e-02},
                                                               true) and
                edm::ParameterDescription<std::vector<double>>("EEPulseShapeTemplate",
                                                               {1.16442e-01,
                                                                7.56246e-01,
                                                                1.00000e+00,
                                                                8.97182e-01,
                                                                6.86831e-01,
                                                                4.91506e-01,
                                                                3.44111e-01,
                                                                2.45731e-01,
                                                                1.74115e-01,
                                                                1.23361e-01,
                                                                8.74288e-02,
                                                                6.19570e-02},
                                                               true)));

  psd0.addNode((edm::ParameterDescription<std::string>("EEdigiCollection", "", true) and
                edm::ParameterDescription<std::string>("EBdigiCollection", "", true) and
                edm::ParameterDescription<std::string>("ESdigiCollection", "", true) and
                edm::ParameterDescription<bool>("UseLCcorrection", true, false) and
                edm::ParameterDescription<std::vector<double>>(
                    "EBCorrNoiseMatrixG12",
                    {1.00000, 0.71073, 0.55721, 0.46089, 0.40449, 0.35931, 0.33924, 0.32439, 0.31581, 0.30481},
                    true) and
                edm::ParameterDescription<std::vector<double>>(
                    "EECorrNoiseMatrixG12",
                    {1.00000, 0.71373, 0.44825, 0.30152, 0.21609, 0.14786, 0.11772, 0.10165, 0.09465, 0.08098},
                    true) and
                edm::ParameterDescription<std::vector<double>>(
                    "EBCorrNoiseMatrixG06",
                    {1.00000, 0.70946, 0.58021, 0.49846, 0.45006, 0.41366, 0.39699, 0.38478, 0.37847, 0.37055},
                    true) and
                edm::ParameterDescription<std::vector<double>>(
                    "EECorrNoiseMatrixG06",
                    {1.00000, 0.71217, 0.47464, 0.34056, 0.26282, 0.20287, 0.17734, 0.16256, 0.15618, 0.14443},
                    true) and
                edm::ParameterDescription<std::vector<double>>(
                    "EBCorrNoiseMatrixG01",
                    {1.00000, 0.73354, 0.64442, 0.58851, 0.55425, 0.53082, 0.51916, 0.51097, 0.50732, 0.50409},
                    true) and
                edm::ParameterDescription<std::vector<double>>(
                    "EECorrNoiseMatrixG01",
                    {1.00000, 0.72698, 0.62048, 0.55691, 0.51848, 0.49147, 0.47813, 0.47007, 0.46621, 0.46265},
                    true) and
                edm::ParameterDescription<bool>("EcalPreMixStage1", false, true) and
                edm::ParameterDescription<bool>("EcalPreMixStage2", false, true)));

  psd0.addOptionalNode(
      (edm::ParameterDescription<std::vector<double>>(
           "EBPulseShapeCovariance",
           {3.001e-06,  1.233e-05,  0.000e+00,  -4.416e-06, -4.571e-06, -3.614e-06, -2.636e-06, -1.286e-06, -8.410e-07,
            -5.296e-07, 0.000e+00,  0.000e+00,  1.233e-05,  6.154e-05,  0.000e+00,  -2.200e-05, -2.309e-05, -1.838e-05,
            -1.373e-05, -7.334e-06, -5.088e-06, -3.745e-06, -2.428e-06, 0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
            0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
            -4.416e-06, -2.200e-05, 0.000e+00,  8.319e-06,  8.545e-06,  6.792e-06,  5.059e-06,  2.678e-06,  1.816e-06,
            1.223e-06,  8.245e-07,  5.589e-07,  -4.571e-06, -2.309e-05, 0.000e+00,  8.545e-06,  9.182e-06,  7.219e-06,
            5.388e-06,  2.853e-06,  1.944e-06,  1.324e-06,  9.083e-07,  6.335e-07,  -3.614e-06, -1.838e-05, 0.000e+00,
            6.792e-06,  7.219e-06,  6.016e-06,  4.437e-06,  2.385e-06,  1.636e-06,  1.118e-06,  7.754e-07,  5.556e-07,
            -2.636e-06, -1.373e-05, 0.000e+00,  5.059e-06,  5.388e-06,  4.437e-06,  3.602e-06,  1.917e-06,  1.322e-06,
            9.079e-07,  6.529e-07,  4.752e-07,  -1.286e-06, -7.334e-06, 0.000e+00,  2.678e-06,  2.853e-06,  2.385e-06,
            1.917e-06,  1.375e-06,  9.100e-07,  6.455e-07,  4.693e-07,  3.657e-07,  -8.410e-07, -5.088e-06, 0.000e+00,
            1.816e-06,  1.944e-06,  1.636e-06,  1.322e-06,  9.100e-07,  9.115e-07,  6.062e-07,  4.436e-07,  3.422e-07,
            -5.296e-07, -3.745e-06, 0.000e+00,  1.223e-06,  1.324e-06,  1.118e-06,  9.079e-07,  6.455e-07,  6.062e-07,
            7.217e-07,  4.862e-07,  3.768e-07,  0.000e+00,  -2.428e-06, 0.000e+00,  8.245e-07,  9.083e-07,  7.754e-07,
            6.529e-07,  4.693e-07,  4.436e-07,  4.862e-07,  6.509e-07,  4.418e-07,  0.000e+00,  0.000e+00,  0.000e+00,
            5.589e-07,  6.335e-07,  5.556e-07,  4.752e-07,  3.657e-07,  3.422e-07,  3.768e-07,  4.418e-07,  6.142e-07},
           true) and
       edm::ParameterDescription<std::vector<double>>(
           "EEPulseShapeCovariance",
           {3.941e-05,  3.333e-05,  0.000e+00,  -1.449e-05, -1.661e-05, -1.424e-05, -1.183e-05, -6.842e-06, -4.915e-06,
            -3.411e-06, 0.000e+00,  0.000e+00,  3.333e-05,  2.862e-05,  0.000e+00,  -1.244e-05, -1.431e-05, -1.233e-05,
            -1.032e-05, -5.883e-06, -4.154e-06, -2.902e-06, -2.128e-06, 0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
            0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
            -1.449e-05, -1.244e-05, 0.000e+00,  5.840e-06,  6.649e-06,  5.720e-06,  4.812e-06,  2.708e-06,  1.869e-06,
            1.330e-06,  9.186e-07,  6.446e-07,  -1.661e-05, -1.431e-05, 0.000e+00,  6.649e-06,  7.966e-06,  6.898e-06,
            5.794e-06,  3.157e-06,  2.184e-06,  1.567e-06,  1.084e-06,  7.575e-07,  -1.424e-05, -1.233e-05, 0.000e+00,
            5.720e-06,  6.898e-06,  6.341e-06,  5.347e-06,  2.859e-06,  1.991e-06,  1.431e-06,  9.839e-07,  6.886e-07,
            -1.183e-05, -1.032e-05, 0.000e+00,  4.812e-06,  5.794e-06,  5.347e-06,  4.854e-06,  2.628e-06,  1.809e-06,
            1.289e-06,  9.020e-07,  6.146e-07,  -6.842e-06, -5.883e-06, 0.000e+00,  2.708e-06,  3.157e-06,  2.859e-06,
            2.628e-06,  1.863e-06,  1.296e-06,  8.882e-07,  6.108e-07,  4.283e-07,  -4.915e-06, -4.154e-06, 0.000e+00,
            1.869e-06,  2.184e-06,  1.991e-06,  1.809e-06,  1.296e-06,  1.217e-06,  8.669e-07,  5.751e-07,  3.882e-07,
            -3.411e-06, -2.902e-06, 0.000e+00,  1.330e-06,  1.567e-06,  1.431e-06,  1.289e-06,  8.882e-07,  8.669e-07,
            9.522e-07,  6.717e-07,  4.293e-07,  0.000e+00,  -2.128e-06, 0.000e+00,  9.186e-07,  1.084e-06,  9.839e-07,
            9.020e-07,  6.108e-07,  5.751e-07,  6.717e-07,  7.911e-07,  5.493e-07,  0.000e+00,  0.000e+00,  0.000e+00,
            6.446e-07,  7.575e-07,  6.886e-07,  6.146e-07,  4.283e-07,  3.882e-07,  4.293e-07,  5.493e-07,  7.027e-07},
           true)),
      true);

  edm::ParameterSetDescription psd;
  psd.addNode(
      edm::ParameterDescription<std::vector<int>>("activeBXs", {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4}, true) and
      edm::ParameterDescription<bool>("ampErrorCalculation", true, true) and
      edm::ParameterDescription<bool>("useLumiInfoRunHeader", true, true) and
      edm::ParameterDescription<int>("bunchSpacing", 0, true) and
      edm::ParameterDescription<bool>("doPrefitEB", false, true) and
      edm::ParameterDescription<double>("prefitMaxChiSqEB", 25., true) and
      edm::ParameterDescription<bool>("dynamicPedestalsEB", false, true) and
      edm::ParameterDescription<bool>("mitigateBadSamplesEB", false, true) and
      edm::ParameterDescription<bool>("gainSwitchUseMaxSampleEB", false, true) and
      edm::ParameterDescription<bool>("selectiveBadSampleCriteriaEB", false, true) and
      edm::ParameterDescription<double>("addPedestalUncertaintyEB", 0., true) and
      edm::ParameterDescription<bool>("simplifiedNoiseModelForGainSwitch", true, true) and
      edm::ParameterDescription<std::string>("timealgo", "RatioMethod", true) and
      edm::ParameterDescription<std::vector<double>>("EBtimeFitParameters",
                                                     {-2.015452e+00,
                                                      3.130702e+00,
                                                      -1.234730e+01,
                                                      4.188921e+01,
                                                      -8.283944e+01,
                                                      9.101147e+01,
                                                      -5.035761e+01,
                                                      1.105621e+01},
                                                     true) and
      edm::ParameterDescription<double>("EBtimeFitLimits_Lower", 0.2, true) and
      edm::ParameterDescription<double>("EBtimeFitLimits_Upper", 1.4, true) and
      edm::ParameterDescription<double>("EBtimeConstantTerm", .6, true) and
      edm::ParameterDescription<double>("EBtimeNconst", 28.5, true) and
      edm::ParameterDescription<double>("outOfTimeThresholdGain10pEB", 5, true) and
      edm::ParameterDescription<double>("outOfTimeThresholdGain10mEB", 5, true) and
      edm::ParameterDescription<double>("outOfTimeThresholdGain1pEB", 5, true) and
      edm::ParameterDescription<double>("outOfTimeThresholdGain1mEB", 5, true) and
      edm::ParameterDescription<double>("amplitudeThresholdEB", 10, true) and
      edm::ParameterDescription<double>("ebSpikeThreshold", 1.042, true) and
      edm::ParameterDescription<std::vector<double>>(
          "ebPulseShape", {5.2e-05, -5.26e-05, 6.66e-05, 0.1168, 0.7575, 1., 0.8876, 0.6732, 0.4741, 0.3194}, true) and
      edm::ParameterDescription<bool>("kPoorRecoFlagEB", true, true) and
      edm::ParameterDescription<double>("chi2ThreshEB_", 65.0, true) and
      edm::ParameterDescription<edm::ParameterSetDescription>("EcalPulseShapeParameters", psd0, true));

  return psd;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN(EcalUncalibRecHitWorkerFactory, EcalPhase2UncalibRecHitWorkerMultiFit, "EcalPhase2UncalibRecHitWorkerMultiFit");
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitFillDescriptionWorkerFactory.h"
DEFINE_EDM_PLUGIN(EcalUncalibRecHitFillDescriptionWorkerFactory,
                  EcalPhase2UncalibRecHitWorkerMultiFit,
                  "EcalPhase2UncalibRecHitWorkerMultiFit");
