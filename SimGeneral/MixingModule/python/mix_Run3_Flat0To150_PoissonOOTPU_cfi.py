# A simple distribution for Run3 studies consisting of a flat distribution from 0
# to 150 for the average pileup.

import FWCore.ParameterSet.Config as cms
from SimGeneral.MixingModule.mix_probFunction_25ns_PoissonOOTPU_cfi import *
mix.input.nbPileupEvents.probFunctionVariable = cms.vint32(range(151))
mix.input.nbPileupEvents.probValue = cms.vdouble([1.0/150.0 for _ in range(151)])

