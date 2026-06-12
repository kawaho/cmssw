#ifndef PhysicsTools_IFNFlavour_IFNFlavourCalculator_h
#define PhysicsTools_IFNFlavour_IFNFlavourCalculator_h

// Fastjet-free interface to the Interleaved Flavour Neutralisation (IFN)
// algorithm (IFNPlugin / fjcontrib 1.0.4, requiring fastjet >= 3.4.1).
//
// This header deliberately exposes NO fastjet types. The vendored fastjet
// 3.4.1 + IFN code is compiled into this package's library with hidden symbol
// visibility (see BuildFile.xml) so it cannot collide with the cvmfs fastjet
// 3.3.0 that the rest of the cmsRun process links against. EDProducers should
// include only this header.

#include <vector>

// The whole library is built with -fvisibility=hidden so the vendored fastjet
// 3.4.1 symbols stay internal. The wrapper API below must remain callable from
// the plugin library, so mark it with default visibility explicitly.
#define IFNFLAVOUR_API __attribute__((visibility("default")))

namespace ifnflavour {

  // One input four-vector. If `isFlavourTag` is true the particle is treated as
  // a flavour carrier (its `pdgId` seeds the IFN FlavInfo); otherwise it is an
  // ordinary, flavourless clustering input (e.g. a light jet constituent).
  //
  // `resetToFlav` is used for hadron-level flavour tagging: when > 0 (a quark
  // flavour 1..6), after the pdgId is decoded into net quark content all flavours
  // EXCEPT `resetToFlav` are zeroed. This is the standard IFN hadron-level
  // prescription (FlavInfo::reset_all_but_flav): a B/C hadron PDG decodes to a
  // heavy quark PLUS a spurious light (anti)quark, and we keep only the heavy
  // flavour (5 for b-hadrons, 4 for c-hadrons) with its correct net sign.
  struct Particle {
    double px, py, pz, E;
    int pdgId;          // PDG id; sign matters (quark vs antiquark)
    bool isFlavourTag;  // true => seed IFN flavour from pdgId
    int resetToFlav = 0;  // if >0 (1..6), keep only this flavour after decoding pdgId
  };

  // Net flavour content of a clustered jet, indexed by |pdgId| 1..6
  // (d,u,s,c,b,t). Net = (#quark - #antiquark). Index 0 is a sentinel:
  //   netFlavour[0] == 21 -> the jet's IFN flavour is FLAVOURLESS (gluon-like).
  //   netFlavour[0] == 0  -> not flavourless; read flavour off netFlavour[1..6].
  struct JetFlavour {
    double px, py, pz, E;
    int netFlavour[7];  // netFlavour[5] = net b, etc.
  };

  // Cluster `inputs` with the IFN algorithm built on anti-kt(R), and return the
  // inclusive jets (pt > ptMin) with their net flavour content. Jets are
  // returned sorted by descending pt. alpha/omega are the IFN parameters
  // (recommended alpha=2, omega=1).
  IFNFLAVOUR_API std::vector<JetFlavour> clusterIFN(const std::vector<Particle>& inputs,
                                                    double R,
                                                    double ptMin,
                                                    double alpha = 2.0,
                                                    double omega = 1.0);

}  // namespace ifnflavour

#endif
