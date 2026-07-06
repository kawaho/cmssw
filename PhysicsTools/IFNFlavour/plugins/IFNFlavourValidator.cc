// -*- C++ -*-
//
// Package:    PhysicsTools/IFNFlavour
// Class:      IFNFlavourValidator
//
// Side-by-side comparison plots for the IFN parton flavour against the stock
// ghost-based parton/hadron flavour, plus eta-phi event displays. Designed to
// run on the fly inside a NanoAOD job (e.g. nano_mc_2018_UL_NANO.py) so the
// validation file is produced in the same cmsRun pass that produces the NanoAOD.
//
// Inputs (all already wired up by PhysicsTools/IFNFlavour/python/ifnFlavour_cff
// when addIFNFlavour() has been called):
//   genJets               -- slimmedGenJets (or the matching AK8 collection)
//   ghostFlavourInfos     -- JetFlavourInfoMatchingCollection from
//                            genJetFlavourAssociation (ghost-based)
//   ifnFlavourInfos       -- JetFlavourInfoMatchingCollection from
//                            genJetFlavourAssociationIFN (IFN, parton flavour
//                            stored in the partonFlavour() slot)
//   ifnJets               -- reco::BasicJetCollection published by the producer
//                            as <module>:ifnJets (the actual IFN jets)
//   ifnJetNetFlavour      -- flat std::vector<int> published as
//                            <module>:ifnJetNetFlavour; jet j's netFlavour[k]
//                            is at index 7*j + k. Parallel to ifnJets.
//   genJetIFNIndex        -- std::vector<int> published as
//                            <module>:genJetIFNIndex; parallel to genJets, gives
//                            the matched IFN jet index (-1 if no match). The
//                            matching is done in the producer; the validator
//                            just forwards it as a TTree branch.
//
// Output (via TFileService):
//   per-IFN-category TH2D composition_(parton|hadron)Flavour_vs_pt_IFN_<X>
//       (x = gen jet pt, y = ghost flavour bucket; user normalises columns in
//        ROOT to get a composition vs pt)
//   per-flavour TH2D etaphi_GenJet_<X> filled one entry per gen jet whose
//       ghost partonFlavour is in bucket X
//   per-flavour TH2D etaphi_IFNJet_<X> filled one entry per IFN jet whose
//       net flavour (reduced via the same partonFlavourFromNet rule used by
//       JetFlavourClusteringIFN) is in bucket X
//   TTree "jets" with one entry per event holding the GenJet collection
//       (pt, eta, phi, ghost parton/hadron flavour, matched IFN parton flavour)
//       and the IFN jet collection (pt, eta, phi, netFlavour[7], reduced
//       parton flavour) -- everything else can be derived from these in ROOT.
//
// IFN parton-flavour buckets used for the per-category histograms:
//   "gluon"   -> 21
//   "light"   -> |f| in {1,2,3}
//   "c"       -> |f| == 4
//   "b"       -> |f| == 5
//   "t"       -> |f| == 6
//   "ambig"   -> -2 (multi-flavour, per partonFlavourFromNet)
//   "none"    -> 0  (no IFN flavour assigned)
//   "unmatched" -> -1 (no IFN jet within deltaR of this gen jet)
// Ghost-based buckets follow the same coarse grouping for |partonFlavour|, with
// hadronFlavour always in {0,4,5}.

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "PhysicsTools/IFNFlavour/interface/IFNFlavourCalculator.h"

#include "TH2D.h"
#include "TTree.h"

class IFNFlavourValidator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit IFNFlavourValidator(const edm::ParameterSet&);
  ~IFNFlavourValidator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // Mirror of JetFlavourClusteringIFN::partonFlavourFromNet -- kept in sync so
  // the bucketing here matches what the producer assigns to each gen jet.
  static int partonFlavourFromNet(const int netFlavour[7]);

  // Coarse bucket key for hist names. Same scheme used for IFN and ghost.
  static std::string flavourBucket(int signedFlavour);
  static float hasBHadronAncestor(const reco::Candidate& p,
                        const reco::GenParticle& bHadron);
  // Returns true if the particle itself or any ancestor is a b-hadron
  // (identified via the PDG hundreds / thousands digit equalling 5).
  static bool isBAncestor(const reco::Candidate& particle);

  const edm::EDGetTokenT<GenEventInfoProduct> genTag_;
  const edm::EDGetTokenT<edm::View<reco::GenJet> > jetsToken_;
  const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> ghostInfosToken_;
  const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> ifnInfosToken_;
  const edm::EDGetTokenT<reco::BasicJetCollection> ifnJetsToken_;
  const edm::EDGetTokenT<std::vector<double>> ifnJetsConstToken_;
  const edm::EDGetTokenT<std::vector<double>> ifnJetsConstToken2_;
  const edm::EDGetTokenT<std::vector<double>> ifnJetsConstToken3_;
  const edm::EDGetTokenT<std::vector<double>> ifnJetsConstToken4_;
  const edm::EDGetTokenT<std::vector<int>> nifnJetsConstToken_;
  const edm::EDGetTokenT<std::vector<int> > ifnNetFlavourToken_;
  const edm::EDGetTokenT<std::vector<int> > genJetIFNIndexToken_;
  const edm::EDGetTokenT<reco::GenParticleRefVector> bHadronsToken_;
  const edm::EDGetTokenT<edm::View<reco::Candidate> > genParticlesToken_;

  const bool strict_;
  const double ptMin_;   // minimum pt for TTree/etaphi entries (composition uses all)
  const double jetR_;    // jet cone size used for bHadron <-> jet dR matching

  // Per-bucket histograms, lazily booked on first sight of a bucket key.
  // composition_(p|h)Flavour_vs_pt_IFN_<bucket>: x = pt, y = ghost flavour code
  std::map<std::string, TH2D*> hCompPartonVsPt_;
  std::map<std::string, TH2D*> hCompHadronVsPt_;
  std::map<std::string, TH2D*> hEtaPhiGen_;
  std::map<std::string, TH2D*> hEtaPhiIFN_;

  TTree* tree_ = nullptr;

  // Per-event scratch buffers for the TTree. Cleared at the start of analyze().
  unsigned int b_run_ = 0, b_lumi_ = 0;
  unsigned long long b_event_ = 0;
  double weight;
  std::vector<float> b_genJet_pt_, b_genJet_eta_, b_genJet_phi_;
  std::vector<int> b_genJet_partonFlavour_, b_genJet_hadronFlavour_, b_genJet_partonFlavourIFN_;
  std::vector<int> b_genJet_ifnJetIdx_;  // matched IFN jet index per gen jet, -1 if no match
  std::vector<float> b_ifnJet_pt_, b_ifnJet_eta_, b_ifnJet_phi_;
  std::vector<int> b_ifnJet_partonFlavour_;
  std::vector<int> b_ifnJet_netFlavour_;  // flat, length 7 * IFNJet_pt.size()

  // bHadron collection (per-event) and matched-bHadron indices per jet.
  std::vector<float> b_bHadron_e_, b_bHadron_pt_, b_bHadron_eta_, b_bHadron_phi_;
  std::vector<int> b_bHadron_pdgId_;
  std::vector<std::vector<int>> b_genJet_bHadronIdx_;
  std::vector<std::vector<float>> b_genJet_bHadron_e_;
  std::vector<std::vector<int>> b_ifnJet_bHadronIdx_;

  // Charged daughters of GenJets: one entry per charged daughter (all)
  std::vector<std::vector<float>> b_genJet_chDaughter_pt_;
  std::vector<std::vector<float>> b_genJet_chDaughter_eta_;
  std::vector<std::vector<float>> b_genJet_chDaughter_phi_;
  std::vector<std::vector<float>> b_genJet_chDaughter_mass_;
  std::vector<std::vector<int>> b_genJet_chDaughter_ch_;
  std::vector<std::vector<int>> b_genJet_chDaughter_B_;
  std::vector<std::vector<int>> b_genJet_chDaughter_pdgid_;
  // All daughters of IFN jets (one entry per daughter)
  std::vector<std::vector<float>> b_ifnJet_daughter_pt_;
  std::vector<std::vector<float>> b_ifnJet_daughter_eta_;
  std::vector<std::vector<float>> b_ifnJet_daughter_phi_;
  std::vector<std::vector<float>> b_ifnJet_daughter_mass_;
};

IFNFlavourValidator::IFNFlavourValidator(const edm::ParameterSet& iConfig)
    : genTag_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
      jetsToken_(consumes<edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJets"))),
      genParticlesToken_(consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
      ghostInfosToken_(
          consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("ghostFlavourInfos"))),
      ifnInfosToken_(
          consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("ifnFlavourInfos"))),
      ifnJetsToken_(consumes<reco::BasicJetCollection>(iConfig.getParameter<edm::InputTag>("ifnJets"))),
      ifnJetsConstToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("ifnJetsConst"))),
      ifnJetsConstToken2_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("ifnJetsConst2"))),
      ifnJetsConstToken3_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("ifnJetsConst3"))),
      ifnJetsConstToken4_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("ifnJetsConst4"))),
      nifnJetsConstToken_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("nifnJetsConst"))),
      ifnNetFlavourToken_(
          consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("ifnJetNetFlavour"))),
      genJetIFNIndexToken_(
          consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genJetIFNIndex"))),
      bHadronsToken_(consumes<reco::GenParticleRefVector>(iConfig.getParameter<edm::InputTag>("bHadrons"))),
      strict_(iConfig.getParameter<bool>("strict")),
      ptMin_(iConfig.getParameter<double>("ptMin")),
      jetR_(iConfig.getParameter<double>("jetR")) {
  usesResource("TFileService");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("jets", "GenJet + IFN jet info for IFN-flavour validation");
  tree_->Branch("run", &b_run_, "run/i");
  tree_->Branch("lumi", &b_lumi_, "lumi/i");
  tree_->Branch("weight", &weight);
  tree_->Branch("event", &b_event_, "event/l");
  tree_->Branch("GenJet_pt", &b_genJet_pt_);
  tree_->Branch("GenJet_eta", &b_genJet_eta_);
  tree_->Branch("GenJet_phi", &b_genJet_phi_);
  tree_->Branch("GenJet_partonFlavour", &b_genJet_partonFlavour_);
  tree_->Branch("GenJet_hadronFlavour", &b_genJet_hadronFlavour_);
  tree_->Branch("GenJet_partonFlavourIFN", &b_genJet_partonFlavourIFN_);
  tree_->Branch("GenJet_ifnJetIdx", &b_genJet_ifnJetIdx_);
  tree_->Branch("IFNJet_pt", &b_ifnJet_pt_);
  tree_->Branch("IFNJet_eta", &b_ifnJet_eta_);
  tree_->Branch("IFNJet_phi", &b_ifnJet_phi_);
  tree_->Branch("IFNJet_partonFlavour", &b_ifnJet_partonFlavour_);
  // Flat vector of length 7 * IFNJet_pt.size(): jet j's netFlavour[k] is at 7*j+k.
  tree_->Branch("IFNJet_netFlavour", &b_ifnJet_netFlavour_);
  tree_->Branch("bHadron_e", &b_bHadron_e_);
  tree_->Branch("bHadron_pt", &b_bHadron_pt_);
  tree_->Branch("bHadron_eta", &b_bHadron_eta_);
  tree_->Branch("bHadron_phi", &b_bHadron_phi_);
  tree_->Branch("bHadron_pdgId", &b_bHadron_pdgId_);
  tree_->Branch("GenJet_bHadronIdx", &b_genJet_bHadronIdx_);
  tree_->Branch("GenJet_bHadronE", &b_genJet_bHadron_e_);
  tree_->Branch("IFNJet_bHadronIdx", &b_ifnJet_bHadronIdx_);
  // Charged daughters of GenJets
  tree_->Branch("GenJet_chDaughter_pt",    &b_genJet_chDaughter_pt_);
  tree_->Branch("GenJet_chDaughter_eta",   &b_genJet_chDaughter_eta_);
  tree_->Branch("GenJet_chDaughter_phi",   &b_genJet_chDaughter_phi_);
  tree_->Branch("GenJet_chDaughter_mass",  &b_genJet_chDaughter_mass_);
  tree_->Branch("GenJet_chDaughter_ch",    &b_genJet_chDaughter_ch_);
  tree_->Branch("GenJet_chDaughter_B",   &b_genJet_chDaughter_B_);
  tree_->Branch("GenJet_chDaughter_pdgid",   &b_genJet_chDaughter_pdgid_);
  // All daughters of IFN jets
  tree_->Branch("IFNJet_daughter_pt",   &b_ifnJet_daughter_pt_);
  tree_->Branch("IFNJet_daughter_eta",  &b_ifnJet_daughter_eta_);
  tree_->Branch("IFNJet_daughter_phi",  &b_ifnJet_daughter_phi_);
  tree_->Branch("IFNJet_daughter_mass", &b_ifnJet_daughter_mass_);
}

float IFNFlavourValidator::hasBHadronAncestor(
    const reco::Candidate& p,
    const reco::GenParticle& bHadron)
{
    std::vector<const reco::Candidate*> stack;

    // start from input particle
    stack.push_back(&p);

    while (!stack.empty()) {
        const reco::Candidate* current = stack.back();
        stack.pop_back();

        if (!current) continue;

        for (size_t i = 0; i < current->numberOfMothers(); ++i) {
            const reco::Candidate* mom = current->mother(i);
            if (!mom) continue;

            if (mom == &bHadron) {
              return p.energy();
            }
            stack.push_back(mom);
        }
    }

    return 0;
}

// Returns true if the particle itself or any ancestor (walking the full mother
// chain) is a b-hadron. B-hadrons are identified by the standard PDG convention:
// the hundreds digit (mesons: |pdgId|/100 % 10) or the thousands digit
// (baryons: |pdgId|/1000 % 10) equals 5.
bool IFNFlavourValidator::isBAncestor(const reco::Candidate& particle) {
    const int absPdg = std::abs(particle.pdgId());
    if ((absPdg / 100) % 10 == 5 || (absPdg / 1000) % 10 == 5) return true;
    for (size_t i = 0; i < particle.numberOfMothers(); ++i) {
        if (particle.mother(i) && isBAncestor(*particle.mother(i))) return true;
    }
    return false;
}

std::string IFNFlavourValidator::flavourBucket(int signedFlavour) {
  if (signedFlavour == 21)
    return "gluon";
  if (signedFlavour == -1)
    return "unmatched";
  if (signedFlavour == -2)
    return "ambig";
  if (signedFlavour == 0)
    return "none";
  const int a = std::abs(signedFlavour);
  if (a >= 1 && a <= 3)
    return "light";
  if (a == 4)
    return "c";
  if (a == 5)
    return "b";
  if (a == 6)
    return "t";
  return "other";
}

void IFNFlavourValidator::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<edm::View<reco::GenJet> > genJets;
  iEvent.getByToken(jetsToken_, genJets);

  edm::Handle<reco::JetFlavourInfoMatchingCollection> ghostInfos;
  iEvent.getByToken(ghostInfosToken_, ghostInfos);

  edm::Handle<reco::JetFlavourInfoMatchingCollection> ifnInfos;
  iEvent.getByToken(ifnInfosToken_, ifnInfos);

  edm::Handle<reco::BasicJetCollection> ifnJets;
  iEvent.getByToken(ifnJetsToken_, ifnJets);

  edm::Handle<std::vector<double>> ifnJetsConst;
  iEvent.getByToken(ifnJetsConstToken_, ifnJetsConst);

  edm::Handle<std::vector<double>> ifnJetsConst2;
  iEvent.getByToken(ifnJetsConstToken2_, ifnJetsConst2);

  edm::Handle<std::vector<double>> ifnJetsConst3;
  iEvent.getByToken(ifnJetsConstToken3_, ifnJetsConst3);

  edm::Handle<std::vector<double>> ifnJetsConst4;
  iEvent.getByToken(ifnJetsConstToken4_, ifnJetsConst4);

  edm::Handle<std::vector<int>> nifnJetsConst;
  iEvent.getByToken(nifnJetsConstToken_, nifnJetsConst);

  edm::Handle<std::vector<int> > ifnNet;
  iEvent.getByToken(ifnNetFlavourToken_, ifnNet);

  edm::Handle<std::vector<int> > genJetIFNIdx;
  iEvent.getByToken(genJetIFNIndexToken_, genJetIFNIdx);

  edm::Handle<reco::GenParticleRefVector> bHadrons;
  iEvent.getByToken(bHadronsToken_, bHadrons);

  if (ifnNet->size() != 7 * ifnJets->size())
    throw cms::Exception("IFNFlavourValidator")
        << "ifnJetNetFlavour size (" << ifnNet->size() << ") != 7 * ifnJets size ("
        << ifnJets->size() << ")";

  if (genJetIFNIdx->size() != genJets->size())
    throw cms::Exception("IFNFlavourValidator")
        << "genJetIFNIndex size (" << genJetIFNIdx->size() << ") != genJets size ("
        << genJets->size() << ")";

  // helper: book histograms on first use of a bucket.
  edm::Service<TFileService> fs;
  auto getOrBook = [&fs](std::map<std::string, TH2D*>& m,
                         const std::string& base,
                         const std::string& key,
                         int nx, double xlo, double xhi,
                         int ny, double ylo, double yhi,
                         const std::string& xtitle,
                         const std::string& ytitle) -> TH2D* {
    auto it = m.find(key);
    if (it != m.end())
      return it->second;
    const std::string name = base + "_" + key;
    TH2D* h = fs->make<TH2D>(name.c_str(), name.c_str(), nx, xlo, xhi, ny, ylo, yhi);
    h->GetXaxis()->SetTitle(xtitle.c_str());
    h->GetYaxis()->SetTitle(ytitle.c_str());
    m[key] = h;
    return h;
  };

  // y-axis bins for the "ghost flavour" composition: discrete codes
  // {-21, -6,...,-1, 0, 1,...,6, 21}. Use ny=15 bins from -7.5 to 7.5 for
  // |flavour|<=6, plus put gluon (21) into the overflow bin via abs-bucketing.
  // Simpler: bin by absolute value 0..7 with 8 bins so 0=undef, 1-6=quark,
  // 7=gluon, and we collapse hadronFlavour {0,4,5} into the same bin set.
//  auto absCodeBin = [](int signedFlavour) -> int {
//    const int a = std::abs(signedFlavour);
//    if (signedFlavour == 21)
//      return 7;
//    if (a >= 1 && a <= 6)
//      return a;
//    return 0;
//  };

  // (2) Per-gen-jet loop: ghost + IFN flavour are AssociationVectors keyed by
  //     the same product (slimmedGenJets), so we just index them by i.
  b_run_ = iEvent.id().run();
  b_lumi_ = iEvent.id().luminosityBlock();
  b_event_ = iEvent.id().event();

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(genTag_, genInfo);

  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  weight = genInfo->weight();

  b_genJet_pt_.clear();
  b_genJet_eta_.clear();
  b_genJet_phi_.clear();
  b_genJet_partonFlavour_.clear();
  b_genJet_hadronFlavour_.clear();
  b_genJet_partonFlavourIFN_.clear();
  b_genJet_ifnJetIdx_.clear();
  b_ifnJet_pt_.clear();
  b_ifnJet_eta_.clear();
  b_ifnJet_phi_.clear();
  b_ifnJet_partonFlavour_.clear();
  b_ifnJet_netFlavour_.clear();
  b_bHadron_e_.clear();
  b_bHadron_pt_.clear();
  b_bHadron_eta_.clear();
  b_bHadron_phi_.clear();
  b_bHadron_pdgId_.clear();
  b_genJet_bHadronIdx_.clear();
  b_genJet_bHadron_e_.clear();
  b_ifnJet_bHadronIdx_.clear();
  b_genJet_chDaughter_pt_.clear();
  b_genJet_chDaughter_eta_.clear();
  b_genJet_chDaughter_phi_.clear();
  b_genJet_chDaughter_mass_.clear();
  b_genJet_chDaughter_ch_.clear();
  b_genJet_chDaughter_B_.clear();
  b_genJet_chDaughter_pdgid_.clear();
  b_ifnJet_daughter_pt_.clear();
  b_ifnJet_daughter_eta_.clear();
  b_ifnJet_daughter_phi_.clear();
  b_ifnJet_daughter_mass_.clear();

  // Fill the bHadron collection up front so jet<->bHadron matching can index it.
  b_bHadron_e_.reserve(bHadrons->size());
  b_bHadron_pt_.reserve(bHadrons->size());
  b_bHadron_eta_.reserve(bHadrons->size());
  b_bHadron_phi_.reserve(bHadrons->size());
  b_bHadron_pdgId_.reserve(bHadrons->size());
  for (const auto& bh : *bHadrons) {
    b_bHadron_e_.push_back(bh->energy());
    b_bHadron_pt_.push_back(bh->pt());
    b_bHadron_eta_.push_back(bh->eta());
    b_bHadron_phi_.push_back(bh->phi());
    b_bHadron_pdgId_.push_back(bh->pdgId());
  }
  const double dR2Max = jetR_ * jetR_;

  const bool ghostKeyed = (genJets.id() == ghostInfos->keyProduct().id());
  const bool ifnKeyed = (genJets.id() == ifnInfos->keyProduct().id());

  for (size_t i = 0; i < genJets->size(); ++i) {
    const reco::GenJet& j = genJets->at(i);

    const double pt = j.pt();
    const double eta = j.eta();
    const double phi = j.phi();

    int ghostParton = 0;
    int ghostHadron = 0;
    if (ghostKeyed && i < ghostInfos->size()) {
      const reco::JetFlavourInfo& gi = (*ghostInfos)[i].second;
      ghostParton = gi.getPartonFlavour();
      ghostHadron = gi.getHadronFlavour();
    }

    int ifnParton = -1;  // -1 == no IFN match assigned by the producer
    if (ifnKeyed && i < ifnInfos->size()) {
      const reco::JetFlavourInfo& ii = (*ifnInfos)[i].second;
      ifnParton = ii.getPartonFlavour();
    }

    const std::string ifnKey = flavourBucket(ifnParton);

//    // composition vs pt, per IFN bucket
//    TH2D* hp = getOrBook(hCompPartonVsPt_,
//                        "composition_partonFlavour_vs_pt_IFN",
//                        ifnKey,
//                        50, 0., 500.,
//                        8, -0.5, 7.5,
//                        "gen jet p_{T} [GeV]",
//                        "|ghost partonFlavour| (7=gluon, 0=undef)");
//    hp->Fill(pt, absCodeBin(ghostParton));
//
//    TH2D* hh = getOrBook(hCompHadronVsPt_,
//                        "composition_hadronFlavour_vs_pt_IFN",
//                        ifnKey,
//                        50, 0., 500.,
//                        8, -0.5, 7.5,
//                        "gen jet p_{T} [GeV]",
//                        "ghost hadronFlavour (0/4/5)");
//    hh->Fill(pt, absCodeBin(ghostHadron));
//
//    // eta-phi event display, gen jets bucketed by ghost partonFlavour
//    if (pt > ptMin_) {
//      TH2D* he = getOrBook(hEtaPhiGen_,
//                          "etaphi_GenJet",
//                          flavourBucket(ghostParton),
//                          80, -5., 5.,
//                          63, -M_PI, M_PI,
//                          "#eta",
//                          "#phi");
//      he->Fill(eta, phi);
//    }

    b_genJet_pt_.push_back(pt);
    b_genJet_eta_.push_back(eta);
    b_genJet_phi_.push_back(phi);
    b_genJet_partonFlavour_.push_back(ghostParton);
    b_genJet_hadronFlavour_.push_back(ghostHadron);
    b_genJet_partonFlavourIFN_.push_back(ifnParton);
    b_genJet_ifnJetIdx_.push_back((*genJetIFNIdx)[i]);

    std::vector<int> matched;
    std::vector<float> matched_e;

    std::vector<float> chPt, chEta, chPhi, chMass;
    std::vector<int> chCh, chPdgid, chB;

//      if (true) { //reco::deltaR2(eta, phi, b_bHadron_eta_[b], b_bHadron_phi_[b]) < dR2Max*2) {
    std::vector<reco::CandidatePtr> const & daughters = j.daughterPtrVector();
    for (size_t b = 0; b < bHadrons->size(); ++b) {
      float energyfromB = 0;
      const reco::GenParticle& bHadron = *((*bHadrons)[b]);
      for (size_t ic = 0; ic < daughters.size(); ++ic) {
          const auto &cand = daughters[ic];
          float energyfromBCand = hasBHadronAncestor(*cand, bHadron);
          energyfromB += energyfromBCand;
          if (b==0) {
            chPt.push_back((*cand).px());
            chEta.push_back((*cand).py());
            chPhi.push_back((*cand).pz());
            chMass.push_back((*cand).energy());
            chPdgid.push_back((*cand).pdgId());
            chCh.push_back((*cand).charge());
            chB.push_back(energyfromBCand > 0 ? 0 : -1);
          }
          else {
            if (energyfromBCand > 0) chB[ic] = b;
          }
       }

      if (energyfromB > 0) {
        matched.push_back(static_cast<int>(b));
        matched_e.push_back(energyfromB/bHadron.energy());
      }
    }

    if (matched.empty()) {
      matched.push_back(-1);
      matched_e.push_back(-1);
    }

    for (edm::View<reco::Candidate>::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it) {
      const int absId = std::abs(it->pdgId());
      if (!(absId==12 || absId==14 || absId==16)) continue;
      if (reco::deltaR2(eta, phi, it->eta(), it->phi()) < dR2Max) {
        float fromb = -1;
        for (size_t b = 0; b < bHadrons->size(); ++b) {
           if (hasBHadronAncestor(*it, *((*bHadrons)[b])) > 0) fromb = b;
        }
        chPt.push_back(it->px());
        chEta.push_back(it->py());
        chPhi.push_back(it->pz());
        chMass.push_back(it->energy());
        chPdgid.push_back(it->pdgId());
        chCh.push_back(it->charge());
        chB.push_back(fromb);
      }
    }

    b_genJet_bHadronIdx_.push_back(std::move(matched));
    b_genJet_bHadron_e_.push_back(std::move(matched_e));

    b_genJet_chDaughter_pt_.push_back(std::move(chPt));
    b_genJet_chDaughter_eta_.push_back(std::move(chEta));
    b_genJet_chDaughter_phi_.push_back(std::move(chPhi));
    b_genJet_chDaughter_mass_.push_back(std::move(chMass));
    b_genJet_chDaughter_ch_.push_back(std::move(chCh));
    b_genJet_chDaughter_B_.push_back(std::move(chB));
    b_genJet_chDaughter_pdgid_.push_back(std::move(chPdgid));
  }

  // (3) Per IFN-jet loop: eta-phi event display + TTree. The IFN jets and
  //     their netFlavour come straight from the producer's published products.
  int nconst = 0;
  for (size_t i = 0; i < ifnJets->size(); ++i) {
    const reco::BasicJet& ij = (*ifnJets)[i];
    const double pt = ij.pt();
    const double eta = ij.rapidity();  // matches what the producer uses for matching
    const double phi = ij.phi();
    int net[7];
    for (int k = 0; k < 7; ++k)
      net[k] = (*ifnNet)[7 * i + k];
    const int reduced = ifnflavour::partonFlavourFromNet(net, strict_);

    if (pt > ptMin_) {
      TH2D* he = getOrBook(hEtaPhiIFN_,
                          "etaphi_IFNJet",
                          flavourBucket(reduced),
                          80, -5., 5.,
                          63, -M_PI, M_PI,
                          "y",
                          "#phi");
      he->Fill(eta, phi);
    }

    b_ifnJet_pt_.push_back(pt);
    b_ifnJet_eta_.push_back(eta);
    b_ifnJet_phi_.push_back(phi);
    b_ifnJet_partonFlavour_.push_back(reduced);
    for (int k = 0; k < 7; ++k)
      b_ifnJet_netFlavour_.push_back(net[k]);

    std::vector<int> matched;
    for (size_t b = 0; b < bHadrons->size(); ++b) {
      if (reco::deltaR2(eta, phi, b_bHadron_eta_[b], b_bHadron_phi_[b]) < dR2Max)
        matched.push_back(static_cast<int>(b));
    }

    if (matched.empty()) matched.push_back(-1);
    b_ifnJet_bHadronIdx_.push_back(std::move(matched));

    // // Save all daughters of this IFN jet
    std::vector<float> dPt, dEta, dPhi, dMass;
     for (size_t d = nconst; d < nconst+nifnJetsConst->at(i); ++d) {
       dPt.push_back(ifnJetsConst->at(d));
       dEta.push_back(ifnJetsConst2->at(d));
       dPhi.push_back(ifnJetsConst3->at(d));
       dMass.push_back(ifnJetsConst4->at(d));
     }
     b_ifnJet_daughter_pt_.push_back(std::move(dPt));
     b_ifnJet_daughter_eta_.push_back(std::move(dEta));
     b_ifnJet_daughter_phi_.push_back(std::move(dPhi));
     b_ifnJet_daughter_mass_.push_back(std::move(dMass));

    nconst+=nifnJetsConst->at(i);
  }

  tree_->Fill();
}

void IFNFlavourValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("generator", edm::InputTag("generator"));
  desc.add<edm::InputTag>("genJets", edm::InputTag("slimmedGenJets"));
  desc.add<edm::InputTag>("ghostFlavourInfos", edm::InputTag("genJetFlavourAssociation"));
  desc.add<edm::InputTag>("ifnFlavourInfos", edm::InputTag("genJetFlavourAssociationIFN"));
  desc.add<edm::InputTag>("ifnJets", edm::InputTag("genJetFlavourAssociationIFN", "ifnJets"))
      ->setComment("BasicJetCollection published by JetFlavourClusteringIFN");
  desc.add<edm::InputTag>("ifnJetsConst", edm::InputTag("genJetFlavourAssociationIFN", "ifnJetsConst"));
  desc.add<edm::InputTag>("ifnJetsConst2", edm::InputTag("genJetFlavourAssociationIFN", "ifnJetsConst2"));
  desc.add<edm::InputTag>("ifnJetsConst3", edm::InputTag("genJetFlavourAssociationIFN", "ifnJetsConst3"));
  desc.add<edm::InputTag>("ifnJetsConst4", edm::InputTag("genJetFlavourAssociationIFN", "ifnJetsConst4"));
  desc.add<edm::InputTag>("nifnJetsConst", edm::InputTag("genJetFlavourAssociationIFN", "nifnJetsConst"));
  desc.add<edm::InputTag>("ifnJetNetFlavour", edm::InputTag("genJetFlavourAssociationIFN", "ifnJetNetFlavour"))
      ->setComment("flat vector<int> of length 7 * nIFNJets, parallel to ifnJets");
  desc.add<edm::InputTag>("genJetIFNIndex", edm::InputTag("genJetFlavourAssociationIFN", "genJetIFNIndex"))
      ->setComment("vector<int> parallel to genJets: index of the matched IFN jet, or -1 if no match");
  desc.add<edm::InputTag>("bHadrons", edm::InputTag("patJetPartons", "bHadrons"))
      ->setComment("selected b-hadrons (HadronAndPartonSelector:bHadrons)");
  desc.add<double>("ptMin", 10.0)->setComment("min pt for eta-phi map and TTree entries");
  desc.add<double>("jetR", 0.4)->setComment("jet cone size used for bHadron <-> jet dR matching");
  desc.add<bool>("strict", false);
  desc.add<edm::InputTag>("genParticles", edm::InputTag("packedGenParticles"));
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(IFNFlavourValidator);
