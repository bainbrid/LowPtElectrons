#ifndef LowPtElectrons_LowPtElectrons_IDNtuple
#define LowPtElectrons_LowPtElectrons_IDNtuple

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include <vector>

class TTree;

namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }

constexpr size_t NHITS_MAX = 30;
constexpr int NEG_INT = -10;
constexpr float NEG_FLOAT = -10.;

// Small class to provide fillers and hide tree I/O
class IDNtuple {

 public:

  IDNtuple() {}

  void reset() {
    IDNtuple dummy; // create a new object 
    *this = dummy; // use assignment to reset
  }
  
  void link_tree( TTree* tree );
  
  void set_weight( float w ) { weight_ = w; }
  void set_rho( float r ) { rho_ = r; }

  void is_e( bool t = true ) { is_e_ = t; }
  void is_e_not_matched( bool t = true ) { is_e_not_matched_ = t; }
  void is_other( bool t = true ) { is_other_ = t; }
  void is_egamma( bool t = true ) { is_egamma_ = t; }

  void has_trk( bool f = false ) { has_trk_ = f; }
  void has_seed( bool f = false ) { has_seed_ = f; }
  void has_gsf( bool f = false ) { has_gsf_ = f; }
  void has_ele( bool f = false ) { has_ele_ = f; }
  
  //void has_egamma_gsf( bool f = false ) { has_egamma_gsf_ = f; }
  //void has_egamma_ele( bool f = false ) { has_egamma_ele_ = f; }

  void fill_evt( const edm::EventID& id );

  void fill_gen( const pat::PackedGenParticleRef );
  void fill_gen( const reco::GenParticleRef ); //@@ AOD
  void fill_gen( const reco::CandidatePtr );

  void fill_trk( const reco::TrackRef& trk,
		 const reco::BeamSpot& spot );

  void fill_seed( double seed_unbiased,
		  double seed_ptbiased );
  
  void fill_preid( const reco::PreId& preid_ecal,
		   const reco::PreId& preid_hcal,
		   const reco::BeamSpot& spot, 
		   const double rho, 
		   noZS::EcalClusterLazyTools& ecalTools );

  void fill_gsf( const reco::GsfTrackRef trk,
		 const reco::BeamSpot& spot );

  void fill_ele( const reco::GsfElectronPtr ele,
		 float mva_value,
		 int mva_id,
		 float ele_conv_vtx_fit_prob,
		 const double rho );
  
 private:

  // Event
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long evt_ = 0;
  float weight_ = 1.;
  float rho_ = NEG_FLOAT;

  // Labels
  bool is_e_ = false;
  bool is_e_not_matched_ = false;
  bool is_other_ = false;
  bool is_egamma_ = false;

  // CMS PF reco
  bool has_trk_ = false;
  bool has_seed_ = false;
  bool has_gsf_ = false;
  bool has_ele_ = false;

  // GEN electrons
  float gen_pt_ = NEG_FLOAT;
  float gen_eta_ = NEG_FLOAT;
  float gen_phi_ = NEG_FLOAT;
  float gen_e_ = NEG_FLOAT;
  float gen_p_ = NEG_FLOAT;
  int gen_charge_ = NEG_INT;
  int gen_pdgid_ = 0;
  int gen_mom_pdgid_ = 0;
  int gen_gran_pdgid_ = 0;

  // KF tracks: kine
  float trk_pt_ = NEG_FLOAT;
  float trk_eta_ = NEG_FLOAT;
  float trk_phi_ = NEG_FLOAT;
  float trk_p_ = NEG_INT;
  int trk_charge_ = NEG_INT;
  float trk_inp_ = NEG_FLOAT;
  float trk_outp_ = NEG_FLOAT;
  float trk_dpt_ = NEG_FLOAT;

  // KF tracks: quality
  int trk_nhits_ = NEG_INT;
  int trk_missing_inner_hits_ = NEG_INT;
  int trk_high_purity_ = NEG_INT;
  float trk_chi2red_ = NEG_FLOAT;

  // KF tracks: displ
  float trk_dxy_ = NEG_FLOAT;
  float trk_dxy_err_ = NEG_FLOAT;
  float trk_dz_ = NEG_FLOAT;
  float trk_dz_err_ = NEG_FLOAT;
  
  // Seed BDT discriminator values
  float preid_unbiased_ = NEG_FLOAT;
  float preid_ptbiased_ = NEG_FLOAT;

  // Seed BDT discriminator values at GsfTrack level
  float seed_unbiased_ = NEG_FLOAT;
  float seed_ptbiased_ = NEG_FLOAT;

  // GSF tracks: kine
  float gsf_pt_ = NEG_FLOAT;
  float gsf_eta_ = NEG_FLOAT;
  float gsf_phi_ = NEG_FLOAT;
  float gsf_p_ = NEG_FLOAT;
  int gsf_charge_ = NEG_INT;
  float gsf_inp_ = NEG_FLOAT;
  float gsf_outp_ = NEG_FLOAT;
  float gsf_dpt_ = NEG_FLOAT;

  // GSF tracks: kine (mode)
  float gsf_mode_pt_ = NEG_FLOAT;
  float gsf_mode_eta_ = NEG_FLOAT;
  float gsf_mode_phi_ = NEG_FLOAT;
  float gsf_mode_p_ = NEG_FLOAT;

  // GSF tracks: quality
  int gsf_nhits_ = NEG_INT;
  int gsf_missing_inner_hits_ = NEG_INT;
  float gsf_chi2red_ = NEG_FLOAT;

  // GSF tracks: displacement
  float gsf_dxy_ = NEG_FLOAT;
  float gsf_dxy_err_ = NEG_FLOAT;
  float gsf_dz_ = NEG_FLOAT;
  float gsf_dz_err_ = NEG_FLOAT;

  // GSF tracks: tangents
  int gsf_ntangents_ = NEG_INT;
  float gsf_hit_dpt_[NHITS_MAX] = {NEG_FLOAT};
  float gsf_hit_dpt_unc_[NHITS_MAX] = {NEG_FLOAT};
  //std::vector<float> gsf_extapolated_eta_;
  //std::vector<float> gsf_extapolated_phi_;

  // GSF electrons: kinematics
  float ele_pt_ = NEG_FLOAT;
  float ele_eta_ = NEG_FLOAT;
  float ele_phi_ = NEG_FLOAT;
  float ele_p_ = NEG_FLOAT;

  // Electrons: IDs
  float ele_mva_value_ = NEG_FLOAT;
  int ele_mva_id_ = NEG_INT;
  float ele_conv_vtx_fit_prob_ = NEG_FLOAT;

  // Electrons: MVA variables
  float eid_rho_ = NEG_FLOAT;
  float eid_ele_pt_ = NEG_FLOAT;

  float eid_trk_p_ = NEG_FLOAT;
  float eid_trk_nhits_ = NEG_FLOAT;
  float eid_trk_chi2red_ = NEG_FLOAT;

  float eid_gsf_nhits_ = NEG_FLOAT;
  float eid_gsf_chi2red_ = NEG_FLOAT;

  float eid_sc_E_ = NEG_FLOAT;
  float eid_sc_eta_ = NEG_FLOAT;
  float eid_sc_etaWidth_ = NEG_FLOAT;
  float eid_sc_phiWidth_ = NEG_FLOAT;

  float eid_match_seed_dEta_ = NEG_FLOAT;
  float eid_match_eclu_EoverP_ = NEG_FLOAT;
  float eid_match_SC_EoverP_ = NEG_FLOAT;
  float eid_match_SC_dEta_ = NEG_FLOAT;
  float eid_match_SC_dPhi_ = NEG_FLOAT;

  float eid_shape_full5x5_sigmaIetaIeta_ = NEG_FLOAT;
  float eid_shape_full5x5_sigmaIphiIphi_ = NEG_FLOAT;
  float eid_shape_full5x5_HoverE_ = NEG_FLOAT;
  float eid_shape_full5x5_r9_ = NEG_FLOAT;
  float eid_shape_full5x5_circularity_ = NEG_FLOAT;

  float eid_brem_frac_ = NEG_FLOAT;
  
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace eleid {
  
  class Features {
  public:
    // KF track
    float trk_p_ = NEG_FLOAT;
    float trk_nhits_ = NEG_FLOAT;
    float trk_chi2red_ = NEG_FLOAT;
    // GSF track
    float gsf_nhits_ = NEG_FLOAT;
    float gsf_chi2red_ = NEG_FLOAT;
    // SC 
    float sc_E_ = NEG_FLOAT;
    float sc_eta_ = NEG_FLOAT;
    float sc_etaWidth_ = NEG_FLOAT;
    float sc_phiWidth_ = NEG_FLOAT;
    // Track-cluster matching
    float match_seed_dEta_ = NEG_FLOAT;
    float match_eclu_EoverP_ = NEG_FLOAT;
    float match_SC_EoverP_ = NEG_FLOAT;
    float match_SC_dEta_ = NEG_FLOAT;
    float match_SC_dPhi_ = NEG_FLOAT;
    // Shower shape vars
    float shape_full5x5_sigmaIetaIeta_ = NEG_FLOAT;
    float shape_full5x5_sigmaIphiIphi_ = NEG_FLOAT;
    float shape_full5x5_HoverE_ = NEG_FLOAT;
    float shape_full5x5_r9_ = NEG_FLOAT;
    float shape_full5x5_circularity_ = NEG_FLOAT;
    // Misc
    float rho_ = NEG_FLOAT;
    float brem_frac_ = NEG_FLOAT;
    float ele_pt_ = NEG_FLOAT;
  public:
    std::vector<float> get();
    void set( const pat::ElectronRef& ele, double rho );
    void set( const reco::GsfElectronPtr& ele, double rho );
  };

}

#endif // LowPtElectrons_LowPtElectrons_IDNtuple