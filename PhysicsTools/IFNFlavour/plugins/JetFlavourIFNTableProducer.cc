// Configurable single-column FlatTable extension producer for the IFN parton
// flavour. The stock PhysicsTools/NanoAOD JetFlavourTableProducer hardcodes the
// column names "partonFlavour"/"hadronFlavour", so it cannot be reused to add a
// parallel IFN branch on the same Jet/GenJet table without a name clash. This
// producer writes ONE int column with a configurable name (default
// "partonFlavourIFN") from the parton flavour stored in the IFN
// JetFlavourInfoMatchingCollection.
//
// It works on edm::View<reco::Jet>, which covers both reco::Jet (PAT/reco jets)
// and reco::GenJet (which derives from reco::Jet), so the same plugin serves the
// reco and gen paths.

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class JetFlavourIFNTableProducer : public edm::stream::EDProducer<> {
public:
  explicit JetFlavourIFNTableProducer(const edm::ParameterSet& iConfig)
      : name_(iConfig.getParameter<std::string>("name")),
        branchName_(iConfig.getParameter<std::string>("branchName")),
        doc_(iConfig.getParameter<std::string>("doc")),
        src_(consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("src"))),
        cut_(iConfig.getParameter<std::string>("cut"), true),
        deltaR_(iConfig.getParameter<double>("deltaR")),
        jetFlavourInfosToken_(consumes<reco::JetFlavourInfoMatchingCollection>(
            iConfig.getParameter<edm::InputTag>("jetFlavourInfos"))) {
    produces<nanoaod::FlatTable>();
  }

  ~JetFlavourIFNTableProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("input Jet/GenJet collection (the table being extended)");
    desc.add<edm::InputTag>("jetFlavourInfos")->setComment("input IFN flavour info collection");
    desc.add<std::string>("name")->setComment("name of the (Gen)Jet FlatTable being extended");
    desc.add<std::string>("branchName", "partonFlavourIFN")->setComment("name of the added flavour column");
    desc.add<std::string>("doc", "IFN parton flavour")->setComment("documentation string for the column");
    desc.add<std::string>("cut", "")->setComment("cut on input collection (must match the base table)");
    desc.add<double>("deltaR", 0.1)->setComment("deltaR to match jets to the flavour info");
    descriptions.addDefault(desc);
  }

private:
  void produce(edm::Event& iEvent, edm::EventSetup const&) override {
    edm::Handle<edm::View<reco::Jet> > jets;
    iEvent.getByToken(src_, jets);

    edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavourInfos;
    iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfos);

    // If `src` and the flavour-info key product are the same collection (e.g.
    // both point at slimmedGenJets on the gen path), entry i of the flavour-info
    // AssociationVector IS the i-th jet of `src` -- no dR matching needed. The
    // dR fallback is used only when the two products differ (the reco path,
    // where `src` is reco PAT jets and the flavour info is keyed by gen jets).
    const bool sameProduct = (jets.id() == jetFlavourInfos->keyProduct().id());

    unsigned int ncand = 0;
    std::vector<int> partonFlavour;

    for (size_t i = 0; i < jets->size(); ++i) {
      const reco::Jet& jet = jets->at(i);
      if (!cut_(jet))
        continue;
      ++ncand;
      if (sameProduct) {
        partonFlavour.push_back((*jetFlavourInfos)[i].second.getPartonFlavour());
        continue;
      }
      // pick the SMALLEST-dR match within deltaR_ (caller sets this to half the
      // target jet cone size); jets outside that cone get flavour 0.
      int best = 0;
      double bestDR2 = deltaR_ * deltaR_;
      for (const reco::JetFlavourInfoMatching& m : *jetFlavourInfos) {
        const double dr2 = reco::deltaR2(jet.p4(), m.first->p4());
        if (dr2 < bestDR2) {
          bestDR2 = dr2;
          best = m.second.getPartonFlavour();
        }
      }
      partonFlavour.push_back(best);
    }

    auto tab = std::make_unique<nanoaod::FlatTable>(ncand, name_, false, true);
    tab->addColumn<int>(branchName_, partonFlavour, doc_, nanoaod::FlatTable::IntColumn);
    iEvent.put(std::move(tab));
  }

  const std::string name_;
  const std::string branchName_;
  const std::string doc_;
  const edm::EDGetTokenT<edm::View<reco::Jet> > src_;
  const StringCutObjectSelector<reco::Jet> cut_;
  const double deltaR_;
  const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
};

DEFINE_FWK_MODULE(JetFlavourIFNTableProducer);
