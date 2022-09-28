#ifndef RecoParticleFlow_PFProducer_interface_MLPFModel
#define RecoParticleFlow_PFProducer_interface_MLPFModel

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

namespace reco::mlpf {

  //The model takes the following number of features for each input PFElement
  static constexpr unsigned int NUM_ELEMENT_FEATURES = 41;

  struct ElementFeatures {
    float type;
    float pt;
    float eta;
    float phi;
    float energy;
    float layer;
    float depth;
    float charge;
    float trajpoint;
    float eta_ecal;
    float phi_ecal;
    float eta_hcal;
    float phi_hcal;
    float muon_dt_hits;
    float muon_csc_hits;
    float muon_type;
    float px;
    float py;
    float pz;
    float deltap;
    float sigmadeltap;
    float gsf_electronseed_trkorecal;
    float gsf_electronseed_dnn1;
    float gsf_electronseed_dnn2;
    float gsf_electronseed_dnn3;
    float gsf_electronseed_dnn4;
    float gsf_electronseed_dnn5;
    float num_hits;
    float cluster_flags;
    float corr_energy;
    float corr_energy_err;
    float vx;
    float vy;
    float vz;
    float pterror;
    float etaerror;
    float phierror;
    float lambda;
    float lambdaerror;
    float theta;
    float thetaerror;

    std::array<float, NUM_ELEMENT_FEATURES> as_array() {
      return {{type,
               pt,
               eta,
               phi,
               energy,
               layer,
               depth,
               charge,
               trajpoint,
               eta_ecal,
               phi_ecal,
               eta_hcal,
               phi_hcal,
               muon_dt_hits,
               muon_csc_hits,
               muon_type,
               px,
               py,
               pz,
               deltap,
               sigmadeltap,
               gsf_electronseed_trkorecal,
               gsf_electronseed_dnn1,
               gsf_electronseed_dnn2,
               gsf_electronseed_dnn3,
               gsf_electronseed_dnn4,
               gsf_electronseed_dnn5,
               num_hits,
               cluster_flags,
               corr_energy,
               corr_energy_err,
               vx,
               vy,
               vz,
               pterror,
               etaerror,
               phierror,
               lambda,
               lambdaerror,
               theta,
               thetaerror}};
    }
  };

  static constexpr unsigned int NUM_OUTPUT_FEATURES = 14;

  //these are defined at model creation time and set the random LSH codebook size
  static constexpr int LSH_BIN_SIZE = 100;
  static constexpr int NUM_MAX_ELEMENTS_BATCH = 200 * LSH_BIN_SIZE;

  //In CPU mode, we want to evaluate each event separately
  static constexpr int BATCH_SIZE = 1;

  //index [0, N_pdgids) -> PDGID
  //this maps the absolute values of the predicted PDGIDs to an array of ascending indices
  static constexpr std::array<int, 8> pdgid_encoding{{0, 211, 130, 1, 2, 22, 11, 13}};
  
  static constexpr unsigned int IDX_CLASS_LAST = pdgid_encoding.size();

  static constexpr unsigned int IDX_CHARGE = IDX_CLASS_LAST+1;
  static constexpr unsigned int IDX_PT = IDX_CLASS_LAST+2;
  static constexpr unsigned int IDX_ETA = IDX_CLASS_LAST+3;
  static constexpr unsigned int IDX_SIN_PHI = IDX_CLASS_LAST+4;
  static constexpr unsigned int IDX_COS_PHI = IDX_CLASS_LAST+5;
  static constexpr unsigned int IDX_ENERGY = IDX_CLASS_LAST+6;

  //for consistency with the baseline PFAlgo
  static constexpr float PI_MASS = 0.13957;


  //PFElement::type -> index [0, N_types)
  //this maps the type of the PFElement to an ascending index that is used by the model to distinguish between different elements
  static const std::map<int, int> elem_type_encoding = {
      {0, 0},
      {1, 1},
      {2, 2},
      {3, 3},
      {4, 4},
      {5, 5},
      {6, 6},
      {7, 7},
      {8, 8},
      {9, 9},
      {10, 10},
      {11, 11},
  };

  ElementFeatures getElementProperties(const reco::PFBlockElement& orig,
                                       const edm::View<reco::GsfElectron>& gsfElectrons);
  float normalize(float in);

  int argMax(std::vector<float> const& vec);

  reco::PFCandidate makeCandidate(int pred_pid,
                                  int pred_charge,
                                  float pred_pt,
                                  float pred_eta,
                                  float pred_sin_phi,
                                  float pred_cos_phi,
                                  float pred_e);

  const std::vector<const reco::PFBlockElement*> getPFElements(const reco::PFBlockCollection& blocks);

  void setCandidateRefs(reco::PFCandidate& cand,
                        const std::vector<const reco::PFBlockElement*> elems,
                        size_t ielem_originator);
};  // namespace reco::mlpf

#endif
