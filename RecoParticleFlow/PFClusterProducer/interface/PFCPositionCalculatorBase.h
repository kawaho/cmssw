#ifndef __PFCPositionCalculatorBase_H__
#define __PFCPositionCalculatorBase_H__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "CondFormats/DataRecord/interface/HcalPFCutsRcd.h"
#include "CondTools/Hcal/interface/HcalPFCutsHandler.h"

#include <string>

namespace edm {
  class EventSetup;
}

class PFCPositionCalculatorBase {
  typedef PFCPositionCalculatorBase PosCalc;

public:
  PFCPositionCalculatorBase(const edm::ParameterSet& conf, edm::ConsumesCollector& cc)
      : _minFractionInCalc(conf.getParameter<double>("minFractionInCalc")),
        _algoName(conf.getParameter<std::string>("algoName")) {}
  virtual ~PFCPositionCalculatorBase() = default;
  //get rid of things we should never use
  PFCPositionCalculatorBase(const PosCalc&) = delete;
  PosCalc& operator=(const PosCalc&) = delete;

  virtual void update(const edm::EventSetup&) {}

  // here we transform one PFCluster to use the new position calculation
  virtual void calculateAndSetPosition(reco::PFCluster&, const HcalPFCuts*) = 0;
  // here you call a loop inside to transform the whole vector
  virtual void calculateAndSetPositions(reco::PFClusterCollection&, const HcalPFCuts*) = 0;

  const std::string& name() const { return _algoName; }

protected:
  const float _minFractionInCalc;

private:
  const std::string _algoName;
};

// define the factory for this base class
#include "FWCore/PluginManager/interface/PluginFactory.h"
typedef edmplugin::PluginFactory<PFCPositionCalculatorBase*(const edm::ParameterSet&, edm::ConsumesCollector&)>
    PFCPositionCalculatorFactory;

#endif
