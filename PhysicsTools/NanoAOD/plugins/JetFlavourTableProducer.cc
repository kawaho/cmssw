// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

class JetFlavourTableProducer : public edm::stream::EDProducer<> {
    public:
        explicit JetFlavourTableProducer(const edm::ParameterSet &iConfig) :
            name_(iConfig.getParameter<std::string>("name")),
            src_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("src"))),
            cut_(iConfig.getParameter<std::string>("cut"), true),
            deltaR_(iConfig.getParameter<double>("deltaR")),
            jetFlavourInfosToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("jetFlavourInfos")))
        {
            produces<nanoaod::FlatTable>();
        }

        ~JetFlavourTableProducer() override {};

        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
            edm::ParameterSetDescription desc;
            desc.add<edm::InputTag>("src")->setComment("input Jet collection");
            desc.add<edm::InputTag>("jetFlavourInfos")->setComment("input flavour info collection");
            desc.add<std::string>("name")->setComment("name of the Jet FlatTable we are extending with flavour information");
            desc.add<std::string>("cut")->setComment("cut on input Jet collection");
            desc.add<double>("deltaR")->setComment("deltaR to match jets");
            descriptions.add("JetFlavourTable", desc);
        }

    private:
        void produce(edm::Event&, edm::EventSetup const&) override ;

        std::string name_;
        edm::EDGetTokenT<edm::View<reco::Jet>> src_;
        const StringCutObjectSelector<reco::Jet> cut_;
        const double deltaR_;
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

};

// ------------ method called to produce the data  ------------
void
JetFlavourTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::Jet>> jets;    
    iEvent.getByToken(src_, jets);
    
    edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavourInfos;
    iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfos);

    unsigned int ncand = 0;
    std::vector<int> partonFlavour;
    std::vector<uint8_t> hadronFlavour;

    for (edm::View<reco::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
//    for (const reco::Candidate & jet : *jets) {
      reco::Jet sel = dynamic_cast<const reco::Jet&>(*jet);
      if (!cut_(sel)) continue;
      ++ncand;
      bool matched = false;
      for (const reco::JetFlavourInfoMatching & jetFlavourInfoMatching : *jetFlavourInfos) {
        if (deltaR(jet->p4(), jetFlavourInfoMatching.first->p4()) < deltaR_) {
          partonFlavour.push_back(jetFlavourInfoMatching.second.getPartonFlavour());
          hadronFlavour.push_back(jetFlavourInfoMatching.second.getHadronFlavour());
          matched = true;
          break;
        }
      }
      if (!matched) {
        partonFlavour.push_back(0);
        hadronFlavour.push_back(0);
      }
    }

    auto tab  = std::make_unique<nanoaod::FlatTable>(ncand, name_, false, true);
    tab->addColumn<int>("partonFlavour", partonFlavour, "flavour from parton matching", nanoaod::FlatTable::IntColumn);
    tab->addColumn<uint8_t>("hadronFlavour", hadronFlavour, "flavour from hadron ghost clustering", nanoaod::FlatTable::UInt8Column);

    iEvent.put(std::move(tab));
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(JetFlavourTableProducer);
