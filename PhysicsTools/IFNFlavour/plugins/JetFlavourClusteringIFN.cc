// -*- C++ -*-
//
// Package:    PhysicsTools/IFNFlavour
// Class:      JetFlavourClusteringIFN
//
/**\class JetFlavourClusteringIFN JetFlavourClusteringIFN.cc PhysicsTools/IFNFlavour/plugins/JetFlavourClusteringIFN.cc

 Description: Determines GEN-jet PARTON flavour using the Interleaved Flavour
 Neutralisation (IFN) algorithm (IFNPlugin, fjcontrib 1.0.4, requiring
 fastjet >= 3.4.1). The producer supports TWO flavour definitions, selectable
 by the `useHadrons` boolean parameter. In both modes the input collection is
 wired via the SAME `genParticles` InputTag -- its CONTENT is what differs.

 PARTON-LEVEL (useHadrons=false, DEFAULT):
   - The clustering input IS the physics-parton collection
     (e.g. patJetPartons:physicsPartons) passed in via `genParticles`. Each
     parton enters with its REAL momentum and pdg-decoded NET flavour. No
     hadron-level particles or heavy hadrons enter the clustering -- this is the
     traditional parton-level flavour definition, but clustered with IFN instead
     of plain anti-kt.

 HADRON-LEVEL (useHadrons=true):
 Unlike the ghost-association approach used by PhysicsTools/JetMCAlgos
 JetFlavourClustering, NO ghost particles are inserted. The IFN flavour is
 defined on the HADRON-level final state, with the effect of heavy-hadron decays
 switched off:

   - The clustering input is the hadron-level final state (e.g.
     packedGenParticles, wired through `genParticles`; neutrinos removed to
     match the *NoNu gen jets), but EVERY final-state particle that descends
     from a b- or c-hadron is REMOVED and the parent heavy hadron is inserted
     in its place (with its REAL momentum). This "switches off" the
     heavy-hadron decay: a B meson is clustered as a single b-flavoured object
     instead of its many light decay products.
   - The heavy hadrons carry net b/c flavour, decoded from their PDG id (with the
     correct net sign) and reduced to the heavy flavour only (b-hadrons -> net b,
     c-hadrons -> net c). b takes priority over c, so a c-hadron from a b decay is
     NOT inserted separately (its b parent already represents that chain). All
     other (light) final-state particles enter the clustering flavourless.

 Common to both modes:
   - IFN (built on anti-kt(R)) tracks the NET flavour of each resulting jet
     through the recombination tree, applying flavour neutralisation when
     flavoured objects merge.
   - Each TARGET gen jet (e.g. slimmedGenJets) is matched to the nearest IFN jet
     within deltaR (default R/2) and inherits its signed parton flavour.

 The output is a reco::JetFlavourInfoMatchingCollection keyed by the target gen
 jets, in which only the parton flavour is set (hadron flavour left at 0). It is
 intended to run in parallel with the stock ghost-based flavour producers so the
 two flavour definitions can be compared per-jet.

 Reco-jet IFN flavour is obtained downstream by matching reco jets to the gen
 jets (and thus to this collection) by dR -- no reco-level clustering is done.

 The fastjet-free IFN core lives in PhysicsTools/IFNFlavour (vendored fastjet
 3.4.1 + IFNPlugin, symbol-isolated from the cvmfs fastjet 3.3.0); this producer
 talks to it only through the ifnflavour::clusterIFN wrapper.
*/

#include <memory>
#include <set>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "PhysicsTools/JetMCUtils/interface/CandMCTag.h"

#include "PhysicsTools/IFNFlavour/interface/IFNFlavourCalculator.h"

class JetFlavourClusteringIFN : public edm::stream::EDProducer<> {
public:
  explicit JetFlavourClusteringIFN(const edm::ParameterSet&);
  ~JetFlavourClusteringIFN() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // translate IFN net flavour content into a signed parton flavour code
  static int partonFlavourFromNet(const int netFlavour[7]);

  // true if `c` has an ancestor that is a b- (wantB=true) or c-hadron
  static bool hasHeavyAncestor(const reco::Candidate* c, bool wantB);

  const edm::EDGetTokenT<edm::View<reco::Jet> > jetsToken_;          // target gen jets
  // genParticles: hadron-level final state when useHadrons=true; physics-parton collection when useHadrons=false
  const edm::EDGetTokenT<edm::View<reco::Candidate> > genParticlesToken_;
  const edm::EDGetTokenT<reco::GenParticleRefVector> bHadronsToken_;
  const edm::EDGetTokenT<reco::GenParticleRefVector> cHadronsToken_;

  const bool useHadrons_;  // false (default): parton-level flavour; true: hadron-level
  const std::string jetAlgorithm_;
  const double rParam_;
  const double deltaR_;  // match IFN jets to target gen jets (default R/2)
  const double alpha_;
  const double omega_;
};

JetFlavourClusteringIFN::JetFlavourClusteringIFN(const edm::ParameterSet& iConfig)
    : jetsToken_(consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
      genParticlesToken_(consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
      bHadronsToken_(consumes<reco::GenParticleRefVector>(iConfig.getParameter<edm::InputTag>("bHadrons"))),
      cHadronsToken_(consumes<reco::GenParticleRefVector>(iConfig.getParameter<edm::InputTag>("cHadrons"))),
      useHadrons_(iConfig.getParameter<bool>("useHadrons")),
      jetAlgorithm_(iConfig.getParameter<std::string>("jetAlgorithm")),
      rParam_(iConfig.getParameter<double>("rParam")),
      deltaR_(iConfig.exists("deltaR") ? iConfig.getParameter<double>("deltaR") : 0.5 * rParam_),
      alpha_(iConfig.exists("alpha") ? iConfig.getParameter<double>("alpha") : 2.0),
      omega_(iConfig.exists("omega") ? iConfig.getParameter<double>("omega") : 1.0) {
  // The IFN wrapper currently builds the base algorithm as anti-kt.
  if (jetAlgorithm_ != "AntiKt")
    throw cms::Exception("InvalidJetAlgorithm")
        << "JetFlavourClusteringIFN currently supports only jetAlgorithm=\"AntiKt\", got: " << jetAlgorithm_;

  produces<reco::JetFlavourInfoMatchingCollection>();
}

int JetFlavourClusteringIFN::partonFlavourFromNet(const int netFlavour[7]) {
  // require EXACTLY one non-zero quark-flavour entry (indices 1..6): return that
  // flavour, signed (positive net => quark, negative => antiquark). If more than
  // one is non-zero, return -2 (ambiguous / multi-flavour). If none, fall back to
  // the wrapper gluon sentinel (netFlavour[0]==21 => gluon, else 0).
  int found = 0;
  int flavour = 0;
  for (int f = 0; f <= 6; ++f) {
    if (netFlavour[f] != 0) {
      ++found;
      if (f==0) { flavour = 21; }
      else { flavour = (netFlavour[f] > 0 ? f : -f); }
    }
  }
  if (found > 1)
    return -2;
  else
    return flavour;
}

bool JetFlavourClusteringIFN::hasHeavyAncestor(const reco::Candidate* c, bool wantB) {
  // breadth/depth-first climb of the mother chain; gen ancestry is shallow.
  std::set<const reco::Candidate*> seen;
  std::vector<const reco::Candidate*> stack(1, c);
  while (!stack.empty()) {
    const reco::Candidate* cur = stack.back();
    stack.pop_back();
    for (size_t i = 0; i < cur->numberOfMothers(); ++i) {
      const reco::Candidate* m = cur->mother(i);
      if (!m || !seen.insert(m).second)
        continue;
      if (wantB ? CandMCTagUtils::hasBottom(*m) : CandMCTagUtils::hasCharm(*m))
        return true;
      stack.push_back(m);
    }
  }
  return false;
}

void JetFlavourClusteringIFN::produce(edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<edm::View<reco::Jet> > jets;
  iEvent.getByToken(jetsToken_, jets);

  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  // Build the IFN clustering inputs.
  //   useHadrons_ = false (default): `genParticles` IS the physics-parton
  //     collection; each entry is inserted with its real momentum and
  //     pdg-decoded NET flavour. b/cHadrons are not used.
  //   useHadrons_ = true            : `genParticles` is the hadron-level final
  //     state; drop neutrinos and any descendant of a b/c hadron, then insert
  //     the parent b/c hadrons themselves as flavour carriers.
  std::vector<ifnflavour::Particle> inputs;

  auto jetFlavourInfos =
      std::make_unique<reco::JetFlavourInfoMatchingCollection>(reco::JetRefBaseProd(jets));

  if (!useHadrons_) {
    // PARTON mode: genParticles is the physics-parton collection.
    inputs.reserve(genParticles->size());
    for (edm::View<reco::Candidate>::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it) {
      if (it->pt() == 0)
        continue;
      inputs.push_back({it->px(), it->py(), it->pz(), it->energy(), it->pdgId(), true, 0});
    }
  } else {
    // HADRON mode: genParticles is the hadron-level final state.
    edm::Handle<reco::GenParticleRefVector> bHadrons;
    iEvent.getByToken(bHadronsToken_, bHadrons);

    edm::Handle<reco::GenParticleRefVector> cHadrons;
    iEvent.getByToken(cHadronsToken_, cHadrons);

    inputs.reserve(genParticles->size() + bHadrons->size() + cHadrons->size());

    // (1) hadron-level final state, minus neutrinos and minus everything that
    //     descends from a heavy hadron (those are represented by the parent
    //     hadron, inserted below with its real momentum and net b/c flavour).
    for (edm::View<reco::Candidate>::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it) {
      const int absId = std::abs(it->pdgId());
      if (absId == 12 || absId == 14 || absId == 16)
        continue;  // drop neutrinos to match the *NoNu gen jets
      if (it->pt() == 0)
        continue;
      if (hasHeavyAncestor(&*it, /*wantB=*/true) || hasHeavyAncestor(&*it, /*wantB=*/false))
        continue;  // decay product of a b/c hadron -> replaced by the hadron below
      inputs.push_back({it->px(), it->py(), it->pz(), it->energy(), 0, false, 0});
    }

    // (2) b-hadrons: inserted as net-b flavour carriers (real momenta).
    for (reco::GenParticleRefVector::const_iterator it = bHadrons->begin(); it != bHadrons->end(); ++it) {
      if ((*it)->pt() == 0)
        continue;
      inputs.push_back({(*it)->px(), (*it)->py(), (*it)->pz(), (*it)->energy(), (*it)->pdgId(), true, 5});
    }

    // (3) c-hadrons: inserted as net-c flavour carriers, but ONLY when they do not
    //     descend from a b-hadron (b has priority; the b chain already covers them).
    for (reco::GenParticleRefVector::const_iterator it = cHadrons->begin(); it != cHadrons->end(); ++it) {
      if ((*it)->pt() == 0)
        continue;
      if (hasHeavyAncestor(&**it, /*wantB=*/true))
        continue;
      inputs.push_back({(*it)->px(), (*it)->py(), (*it)->pz(), (*it)->energy(), (*it)->pdgId(), true, 4});
    }
  }

  // cluster the parton or hadron-level final state with IFN
  std::vector<ifnflavour::JetFlavour> ifnJets =
      ifnflavour::clusterIFN(inputs, rParam_, 0.0, alpha_, omega_);

  // match each target gen jet to the nearest IFN jet within deltaR
  for (size_t j = 0; j < jets->size(); ++j) {
    int partonFlavour = -1;
    double bestDR2 = deltaR_ * deltaR_;
    for (const auto& ij : ifnJets) {
      reco::Particle::LorentzVector p4(ij.px, ij.py, ij.pz, ij.E);
      double dr2 = reco::deltaR2(jets->at(j).rapidity(), jets->at(j).phi(), p4.Rapidity(), p4.phi());
      if (dr2 < bestDR2) {
        bestDR2 = dr2;
        partonFlavour = partonFlavourFromNet(ij.netFlavour);
      }
    }

    reco::GenParticleRefVector empty;
    (*jetFlavourInfos)[jets->refAt(j)] =
        reco::JetFlavourInfo(empty, empty, empty, empty, 0, partonFlavour);
  }

  iEvent.put(std::move(jetFlavourInfos));
}

void JetFlavourClusteringIFN::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets")->setComment("target gen jets to assign IFN parton flavour to");
  desc.add<edm::InputTag>("genParticles")
      ->setComment("useHadrons=false: physics-parton collection (e.g. patJetPartons:physicsPartons); "
                   "useHadrons=true : hadron-level final state (e.g. packedGenParticles)");
  desc.add<edm::InputTag>("bHadrons")
      ->setComment("selected b-hadrons (HadronAndPartonSelector:bHadrons); only used when useHadrons=true");
  desc.add<edm::InputTag>("cHadrons")
      ->setComment("selected c-hadrons (HadronAndPartonSelector:cHadrons); only used when useHadrons=true");
  desc.add<bool>("useHadrons", false)
      ->setComment("false (default): PARTON-level IFN flavour (genParticles = physics partons); "
                   "true            : HADRON-level IFN flavour (b/c hadrons as flavour carriers)");
  desc.add<std::string>("jetAlgorithm", "AntiKt");
  desc.add<double>("rParam");
  desc.add<double>("deltaR", 0.2)->setComment("dR to match IFN jets to target gen jets (default R/2)");
  desc.add<double>("alpha", 2.0);
  desc.add<double>("omega", 1.0);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(JetFlavourClusteringIFN);
