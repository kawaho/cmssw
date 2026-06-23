// Implementation of the fastjet-free IFN wrapper. Everything fastjet-related is
// confined to this translation unit and compiled with hidden visibility (see
// BuildFile.xml) so the vendored fastjet 3.4.1 symbols never leak out to clash
// with the cvmfs fastjet 3.3.0 used elsewhere in the process.

#include "PhysicsTools/IFNFlavour/interface/IFNFlavourCalculator.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/FlavInfo.hh"
#include "fastjet/contrib/IFNPlugin.hh"

#include <algorithm>

using namespace fastjet;
using namespace fastjet::contrib;

namespace ifnflavour {

  int partonFlavourFromNet(const int netFlavour[7], const bool strict_) {
    // require EXACTLY one non-zero quark-flavour entry (indices 1..6): return that
    // flavour, signed (positive net => quark, negative => antiquark). If more than
    // one is non-zero, return -2 (ambiguous / multi-flavour). If none, fall back to
    // the wrapper gluon sentinel (netFlavour[0]==21 => gluon, else 0).
    int found = 0;
    int foundHeavy = 0;
    int flavour = 0;
    for (int f = 0; f <= 6; ++f) {
      if (netFlavour[f] != 0) {
        ++found;
        if (f>=4) { ++foundHeavy; }
        if (f==0) { flavour = 21; }
        else { flavour = (netFlavour[f] > 0 ? f : -f); }
      }
    }
    if ((foundHeavy > 1) | ((found > 1) & strict_))
      return 11;
    else
      return flavour;
  }

  std::vector<JetFlavour> clusterIFN(const std::vector<Particle>& inputs,
                                     double R,
                                     double ptMin,
                                     double alpha,
                                     double omega) {
    // Build the IFN jet definition on top of anti-kt(R), tracking net flavour.
    JetDefinition base_jet_def(antikt_algorithm, R);
    FlavRecombiner flav_recombiner;
    base_jet_def.set_recombiner(&flav_recombiner);

    auto* ifn_plugin = new IFNPlugin(base_jet_def, alpha, omega, FlavRecombiner::net);
    JetDefinition ifn_jet_def(ifn_plugin);
    ifn_jet_def.delete_plugin_when_unused();

    // Build the input PseudoJets, seeding flavour on the tagged carriers.
    std::vector<PseudoJet> event;
    event.reserve(inputs.size());
    for (const auto& p : inputs) {
      PseudoJet pj(p.px, p.py, p.pz, p.E);
      // Every input needs flavour tracking; flavourless inputs get pdg 0.
      // Decode the pdgId into net flavour content via FlavInfo, then (for
      // hadron-level tagging) optionally keep only the requested heavy flavour.
      FlavInfo flav(p.pdgId, p.charge); //p.isFlavourTag ? p.pdgId : 0);
//      if (p.isFlavourTag && p.resetToFlav > 0)
//        flav.reset_all_but_flav(p.resetToFlav);
      pj.set_user_info(new FlavHistory(flav));
      event.push_back(pj);
    }

    ClusterSequence cs(event, ifn_jet_def);
    std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptMin));

    std::vector<JetFlavour> out;
    out.reserve(jets.size());
    for (const auto& jet : jets) {
      JetFlavour jf;
      jf.px = jet.px();
      jf.py = jet.py();
      jf.pz = jet.pz();
      jf.E = jet.E();
      const FlavInfo& flav = FlavHistory::current_flavour_of(jet);
      //std::cout << "Jet flavour: " << flav.description() << std::endl;
      //std::cout << "Jet flavour list : " << flav[0] << " " << flav[1] << " " << flav[2] << " " << flav[3] << " " << flav[4] << " " << flav[5] << " " << flav[6] << std::endl;
      for (int i = 0; i <= 6; ++i)
        jf.netFlavour[i] = flav[i];
      jf.constituents = jet.constituents();
      out.push_back(jf);
    }
    return out;
  }

}  // namespace ifnflavour
