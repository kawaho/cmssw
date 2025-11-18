#ifndef RecoParticleFlow_PFProducer_interface_MLPFModel
#define RecoParticleFlow_PFProducer_interface_MLPFModel

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

namespace reco::mlpf {

  //The model takes the following number of features for each input PFElement
  static constexpr unsigned int NUM_ELEMENT_FEATURES = 72;
  static constexpr unsigned int NUM_VERTEX_FEATURES = 17;

  struct ElementFeatures {
    float type = 0.0;
    float pt = 0.0;
    float eta = 0.0;
    float phi = 0.0;
    float energy = 0.0;
    float layer = 0.0;
    float depth = 0.0;
    float charge = 0.0;
    float trajpoint = 0.0;
    float eta_ecal = 0.0;
    float phi_ecal = 0.0;
    float eta_hcal = 0.0;
    float phi_hcal = 0.0;
    float muon_dt_hits = 0.0;
    float muon_csc_hits = 0.0;
    float muon_type = 0.0;
    float px = 0.0;
    float py = 0.0;
    float pz = 0.0;
    float sigma_x = 0.0;
    float sigma_y = 0.0;
    float sigma_z = 0.0;
    float deltap = 0.0;
    float sigmadeltap = 0.0;
    float gsf_electronseed_trkorecal = 0.0;
    float gsf_electronseed_dnn1 = 0.0;
    float gsf_electronseed_dnn2 = 0.0;
    float gsf_electronseed_dnn3 = 0.0;
    float gsf_electronseed_dnn4 = 0.0;
    float gsf_electronseed_dnn5 = 0.0;
    float num_hits = 0.0;
    float cluster_flags = 0.0;
    float corr_energy = 0.0;
    float corr_energy_err = 0.0;
    float pca_x = 0.0;
    float pca_y = 0.0;
    float pca_z = 0.0;
    float pterror = 0.0;
    float etaerror = 0.0;
    float phierror = 0.0;
    float lambda = 0.0;
    float lambdaerror = 0.0;
    float theta = 0.0;
    float thetaerror = 0.0;
    float time = 0.0;
    float timeerror = 0.0;
    float etaerror1 = 0.0;
    float phierror1 = 0.0;
    float etaerror2 = 0.0;
    float phierror2 = 0.0;
    float etaerror3 = 0.0;
    float phierror3 = 0.0;
    float etaerror4 = 0.0;
    float phierror4 = 0.0;
    float vtx_x = 0.0;
    float vtx_y = 0.0;
    float vtx_z = 0.0;
    float ntracks = 0.0;
    float v_normalized_chi2 = 0.0;
    float vx = 0.0;
    float vy = 0.0;
    float vz = 0.0;
    float vt = 0.0;
    float vx_err = 0.0;
    float vy_err = 0.0;
    float vz_err = 0.0;
    float vt_err = 0.0;
    float vpx = 0.0;
    float vpy = 0.0;
    float vpz = 0.0;
    float ve = 0.0;

    // MLPF features in 2024
    // from particleflow/mlpf/heptfds/cms_pf/utils.py
    std::array<float, NUM_ELEMENT_FEATURES> as_array() {
      return {{type,
               pt,
               eta,
               std::sin(phi),
               std::cos(phi),
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
               pca_x,
               pca_y,
               pca_z,
               pterror,
               etaerror,
               phierror,
               lambda,
               lambdaerror,
               theta,
               thetaerror,
               time,
               timeerror,
               etaerror1,
               etaerror2,
               etaerror3,
               etaerror4,
               phierror1,
               phierror2,
               phierror3,
               phierror4,
	       sigma_x,
	       sigma_y,
	       sigma_z,
               vtx_x,
               vtx_y,
               vtx_z,
               ntracks,
               v_normalized_chi2,
               vx,
               vy,
               vz,
               vt,
               vx_err,
               vy_err,
               vz_err,
               vt_err,
               vpx,
               vpy,
               vpz,
               ve, 
        }};
    }
  };

  static constexpr unsigned int NUM_OUTPUT_FEATURES_CLS = 9;
  static constexpr unsigned int NUM_OUTPUT_FEATURES_P4 = 5;

  //In CPU mode, we want to evaluate each event separately
  static constexpr int BATCH_SIZE = 1;

  //index [0, N_pdgids) -> PDGID
  //this maps the absolute values of the predicted PDGIDs to an array of ascending indices
  static constexpr std::array<int, 9> pdgid_encoding{{0, 211, 130, 1, 2, 22, 11, 13, 15}};
  
  static constexpr unsigned int IDX_CLASS_LAST = pdgid_encoding.size()-1;

  static constexpr unsigned int IDX_PT = 0;
  static constexpr unsigned int IDX_ETA = 1;
  static constexpr unsigned int IDX_SIN_PHI = 2;
  static constexpr unsigned int IDX_COS_PHI = 3;
  static constexpr unsigned int IDX_ENERGY = 4;

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
                                       const edm::View<reco::GsfElectron>& gsfElectrons,
                                       const reco::VertexCollection& primaryVertices);

  ElementFeatures getVertexProperties(const reco::Vertex& primaryVertex);

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
