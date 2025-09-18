#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

template<typename T>
class JetConstituentTableProducer : public edm::stream::EDProducer<> {
public:
  explicit JetConstituentTableProducer(const edm::ParameterSet &);
  ~JetConstituentTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  bool isBAncestor(const reco::Candidate &particle);

  const std::string name_;
  const std::string idx_name_;

  edm::EDGetTokenT<edm::View<T>> jet_token_;
};

//
// constructors and destructor
//
template< typename T>
JetConstituentTableProducer<T>::JetConstituentTableProducer(const edm::ParameterSet &iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      idx_name_(iConfig.getParameter<std::string>("idx_name")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))) {
  produces<nanoaod::FlatTable>(name_);
}

template< typename T>
JetConstituentTableProducer<T>::~JetConstituentTableProducer() {}

template< typename T>
bool JetConstituentTableProducer<T>::isBAncestor(const reco::Candidate &particle) {
   //particle is already the B
   if (((int)((abs(particle.pdgId()) / 100) % 10) == 5) || ((int)((abs(particle.pdgId()) / 1000) % 10) == 5)) return true;

   //otherwise loop on mothers, if any and return true if the B ancestor is found
   for (size_t i = 0; i < particle.numberOfMothers(); i++) {
      if (isBAncestor(*particle.mother(i))) return true;
                                        
   }

   return false;
}

template< typename T>
void JetConstituentTableProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  std::vector<int> jetIdx_pf, from_b;
  std::vector<float> cand_pt, cand_eta, cand_phi, cand_mass, cand_charge, cand_pdgId, jet_pt;

  auto jets = iEvent.getHandle(jet_token_);
  int size = 0;
  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet = jets->at(i_jet);
    //Bad practice !!! hard coded pt and assummed jet col is always sorted
    if (jet.pt() < 10) break;    
    // PF Cands    
    std::vector<reco::CandidatePtr> const & daughters = jet.daughterPtrVector();
    for (const auto &cand : daughters) {
      jetIdx_pf.push_back(i_jet);
      jet_pt.push_back(jet.pt());
      cand_pt.push_back(cand->pt());
      cand_eta.push_back(cand->eta());
      cand_phi.push_back(cand->phi());
      cand_mass.push_back(cand->mass());
      cand_charge.push_back(cand->charge());
      cand_pdgId.push_back(cand->pdgId());
      if (isBAncestor(*cand)) {
        from_b.push_back(1);
      } else {
        from_b.push_back(0);
      }
      size+=1;
    }  // end jet loop
  }

  auto candTable = std::make_unique<nanoaod::FlatTable>(size, name_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTableProducer
  candTable->addColumn<int>("jetIdx", jetIdx_pf, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
  candTable->addColumn<float>("pt", cand_pt, "pt", nanoaod::FlatTable::FloatColumn);
  candTable->addColumn<float>("eta", cand_eta, "eta", nanoaod::FlatTable::FloatColumn);
  candTable->addColumn<float>("phi", cand_phi, "phi", nanoaod::FlatTable::FloatColumn);
  candTable->addColumn<float>("mass", cand_mass, "mass", nanoaod::FlatTable::FloatColumn);
  candTable->addColumn<float>("charge", cand_charge, "charge", nanoaod::FlatTable::FloatColumn);
  candTable->addColumn<float>("pdgId", cand_pdgId, "pdgId", nanoaod::FlatTable::FloatColumn);
  candTable->addColumn<float>("jetPt", jet_pt, "jetPt", nanoaod::FlatTable::FloatColumn);
  candTable->addColumn<int>("isbProduct", from_b, "True if ancestor is a b hadron", nanoaod::FlatTable::IntColumn);
  iEvent.put(std::move(candTable), name_);
}

template< typename T>
void JetConstituentTableProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("name", "JetPFCands");
  desc.add<std::string>("idx_name", "candIdx");
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK8"));
  descriptions.addWithDefaultLabel(desc);
}

typedef JetConstituentTableProducer<pat::Jet> PatJetConstituentTableProducer;
typedef JetConstituentTableProducer<reco::GenJet> GenJetConstituentTableProducer;

DEFINE_FWK_MODULE(PatJetConstituentTableProducer);
DEFINE_FWK_MODULE(GenJetConstituentTableProducer);
