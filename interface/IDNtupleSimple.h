#ifndef LowPtElectrons_LowPtElectrons_IDNtupleSimple
#define LowPtElectrons_LowPtElectrons_IDNtupleSimple

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Framework/interface/Event.h"
#include "LowPtElectrons/LowPtElectrons/interface/Common.h"
#include <vector>

class TTree;

constexpr size_t ARRAY_SIZE = 20;

// Small class to provide fillers and hide tree I/O
class IDNtupleSimple {

 public:
  
  IDNtupleSimple() {}

  void reset() {
    IDNtupleSimple dummy; // create a new object 
    *this = dummy; // use assignment to reset
  }
  
  void link_tree( TTree* tree );
  
  void set_weight( float w ) { weight_ = w; }
  void set_prescale( float p ) { prescale_ = p; }
  void set_rho( float r ) { rho_ = r; }

  void is_aod( int aod ) { is_aod_ = aod; }
  void is_mc( int mc ) { is_mc_ = mc; }

  void tag_pt( float x ) { tag_pt_ = x; }
  void tag_eta( float x ) { tag_eta_ = x; }

  void is_e( bool t = true ) { is_e_ = t; }
  void is_e_not_matched( bool t = true ) { is_e_not_matched_ = t; }
  void is_other( bool t = true ) { is_other_ = t; }
  void is_egamma( bool t = true ) { is_egamma_ = t; }

  void has_ele( bool f = false ) { has_ele_ = f; }
  void ele_dr( float dr ) { ele_dr_ = dr; }

  void fill_evt( const edm::EventID& id );

  //void fill_gen( const pat::PackedGenParticleRef );
  //void fill_gen( const reco::GenParticlePtr );
  void fill_gen( const reco::CandidatePtr );

  void fill_ele( const reco::GsfElectronPtr ele,
		 float mva_value,
		 float mva_value_retrained,
		 const double rho,
		 bool is_egamma = false );
  
 public:

  // Event
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long evt_ = 0;
  float prescale_ = 0.;
  float weight_ = 1.;
  float rho_ = id::NEG_FLOAT;

  // Data sample
  int is_aod_ = -1;
  int is_mc_ = -1;

  // Tag-side muon
  float tag_pt_ = id::NEG_FLOAT;
  float tag_eta_ = id::NEG_FLOAT;

  // Labels
  bool is_e_ = false;
  bool is_e_not_matched_ = false;
  bool is_other_ = false;
  bool is_egamma_ = false;

  // RECO steps
  bool has_ele_ = false;
  float ele_dr_ = id::NEG_FLOAT;

  // GEN electrons
  float gen_pt_ = id::NEG_FLOAT;
  float gen_eta_ = id::NEG_FLOAT;
  float gen_phi_ = id::NEG_FLOAT;
  float gen_e_ = id::NEG_FLOAT;
  float gen_p_ = id::NEG_FLOAT;
  int gen_charge_ = id::NEG_INT;
  int gen_pdgid_ = 0;
  int gen_mom_pdgid_ = 0;
  int gen_gran_pdgid_ = 0;

  // GSF electrons: kinematics
  float ele_pt_ = id::NEG_FLOAT;
  float ele_eta_ = id::NEG_FLOAT;
  float ele_phi_ = id::NEG_FLOAT;
  float ele_p_ = id::NEG_FLOAT;

  // Electrons: IDs
  float ele_mva_value_ = -999.; //@ id::NEG_FLOAT;
  float ele_mva_value_retrained_ = -999.; //@ id::NEG_FLOAT;

};

#endif // LowPtElectrons_LowPtElectrons_IDNtupleSimple
