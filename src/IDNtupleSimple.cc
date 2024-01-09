#include "LowPtElectrons/LowPtElectrons/interface/IDNtupleSimple.h"
#include "TTree.h"

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtupleSimple::link_tree( TTree *tree ) {

  tree->Branch("run",  &run_ , "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt",  &evt_ , "evt/i");
  tree->Branch("weight", &weight_, "weight/f");
  tree->Branch("prescale", &prescale_, "prescale/f");
  tree->Branch("rho", &rho_, "rho/f");

  tree->Branch("is_aod", &is_aod_, "is_aod/i");
  tree->Branch("is_mc", &is_mc_, "is_mc/i");

  tree->Branch("tag_pt", &tag_pt_, "tag_pt/f");
  tree->Branch("tag_eta", &tag_eta_, "tag_eta/f");

  tree->Branch("is_e", &is_e_, "is_e/O");
  tree->Branch("is_egamma", &is_egamma_, "is_egamma/O");

  tree->Branch("has_ele", &has_ele_, "has_ele/O");
  tree->Branch("ele_dr", &ele_dr_, "ele_dr/f");
  
  tree->Branch("gen_pt" , &gen_pt_ , "gen_pt/f" );
  tree->Branch("gen_eta", &gen_eta_, "gen_eta/f");
  tree->Branch("gen_phi", &gen_phi_, "gen_phi/f");
  tree->Branch("gen_e", &gen_e_, "gen_e/f");
  tree->Branch("gen_p", &gen_p_, "gen_p/f");
  tree->Branch("gen_charge", &gen_charge_, "gen_charge/I");
  tree->Branch("gen_pdgid", &gen_pdgid_, "gen_pdgid/I");
  tree->Branch("gen_mom_pdgid", &gen_mom_pdgid_, "gen_mom_pdgid/I");
  tree->Branch("gen_gran_pdgid", &gen_gran_pdgid_, "gen_gran_pdgid/I");
  
  tree->Branch("ele_pt", &ele_pt_, "ele_pt/f");
  tree->Branch("ele_p", &ele_p_, "ele_p/f");
  tree->Branch("ele_eta", &ele_eta_, "ele_eta/f");
  tree->Branch("ele_phi", &ele_phi_, "ele_phi/f");

  tree->Branch("ele_mva_value", &ele_mva_value_, "ele_mva_value/f");
  tree->Branch("ele_mva_value_retrained", &ele_mva_value_retrained_, "ele_mva_value_retrained/f");

}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtupleSimple::fill_evt( const edm::EventID& id ) {
  run_  = id.run();
  lumi_ = id.luminosityBlock();
  evt_  = id.event();
}

/////////////////////////////////////////////////////////////////////////////////
//
//void IDNtupleSimple::fill_gen( const reco::GenParticlePtr genp ) {
//  gen_pt_  = genp->pt();
//  gen_eta_ = genp->eta();
//  gen_phi_ = genp->phi();
//  gen_e_ = genp->energy();
//  gen_p_ = genp->p();
//  gen_charge_ = genp->charge();
//  gen_pdgid_ = 0;
//  gen_mom_pdgid_ = 0;
//  gen_gran_pdgid_ = 0;
//}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtupleSimple::fill_gen( const reco::CandidatePtr genp ) {
  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_e_ = genp->energy();
  gen_p_ = genp->p();
  gen_charge_ = genp->charge();
  gen_pdgid_ = 0;
  gen_mom_pdgid_ = 0;
  gen_gran_pdgid_ = 0;
}

/////////////////////////////////////////////////////////////////////////////////
//
//void IDNtupleSimple::fill_gen( const pat::PackedGenParticleRef genp ) {
//  gen_pt_  = genp->pt();
//  gen_eta_ = genp->eta();
//  gen_phi_ = genp->phi();
//  gen_e_ = genp->energy();
//  gen_p_ = genp->p();
//  gen_charge_ = genp->charge();
//  gen_pdgid_ = 0;
//  gen_mom_pdgid_ = 0;
//  gen_gran_pdgid_ = 0;
//}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtupleSimple::fill_ele( const reco::GsfElectronPtr ele,
			       float mva_value,
			       float mva_value_retrained,
			       const double rho,
			       bool is_egamma ) {
  
  // Kinematics
  if ( is_egamma ) {
    ele_p_ = ele->p();
    ele_pt_ = ele->pt();
    ele_eta_ = ele->eta();
    ele_phi_ = ele->phi();
  } else {
    reco::GsfTrackRef gsf = ele->gsfTrack();
    if ( gsf.isNonnull() && gsf.isAvailable() ) {
      ele_p_ = gsf->p();
      ele_pt_ = gsf->pt();
      ele_eta_ = gsf->eta();
      ele_phi_ = gsf->phi();
    } else {
      std::cout << "[IDNtupleSimple::fill_ele] ERROR: Null GsfTrackRef!" << std::endl;
    }
  }
  
  // MVA IDs: only filled if 'ValueMap->size() == electrons->size()' in IDFeatures::analyze()
  if ( mva_value > -666. ) { ele_mva_value_ = mva_value; }
  if ( mva_value_retrained > -666. ) { ele_mva_value_retrained_ = mva_value_retrained; }

}
