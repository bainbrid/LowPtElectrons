#include "LowPtElectrons/LowPtElectrons/interface/Common.h"

////////////////////////////////////////////////////////////////////////////////
//
template <typename T> 
void key_id( const edm::Ptr<T>& ptr, std::stringstream& ss ) {
  ss << "id/key: ";
  if ( ptr.isNull() || !ptr.isAvailable() ) { ss << "InvalidKey"; }
  else { 
    std::stringstream tmp; tmp << ptr.id();
    ss << std::setw(6) << tmp.str() << "/"
       << std::setw(3) << int( ptr.key() ); }
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T> 
void pt_eta_phi( const edm::Ptr<T>& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << ", pt/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->pt() << ", " 
     << std::setw(4) << ptr->eta() << ", " 
     << std::setw(4) << ptr->phi();
}

////////////////////////////////////////////////////////////////////////////////
//
void pt_eta_phi( const edm::Ptr<reco::GsfTrack>& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << ", pt/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->ptMode() << ", " 
     << std::setw(4) << ptr->etaMode() << ", " 
     << std::setw(4) << ptr->phiMode();
}

////////////////////////////////////////////////////////////////////////////////
//
void pt_eta_phi( const reco::CaloClusterPtr& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << ", Et/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->energy() / std::cosh(ptr->eta()) << ", " 
     << std::setw(4) << ptr->eta() << ", " 
     << std::setw(4) << ptr->phi();
}

////////////////////////////////////////////////////////////////////////////////
//
std::ostream& operator<< ( std::ostream& out, const ElectronChain& obj ) {
  std::stringstream ss;
  ss << "ElectronChain:"
     << " is_egamma: " << obj.is_egamma_
     << " is_e: " << obj.is_e_
     << " is_mc: " << obj.is_mc_
     << " is_aod: " << obj.is_aod_;

  ss << "\n TRG:   "
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << obj.tag_pt_ << ", " 
     << std::setw(4) << obj.tag_eta_ << ", " 
     << std::setw(4) << " n/a";
  //<< std::setprecision(2) << obj.tag_pt_ << "/"
  //<< std::setprecision(2) << obj.tag_eta_;

  ss << "\n SIG:   "; key_id(obj.sig_,ss); pt_eta_phi(obj.sig_,ss);
  ss << "\n TRK:   "; key_id(obj.trk_,ss); pt_eta_phi(obj.trk_,ss);
  ss << ", PdgId: " << obj.pdg_id_;
  if (obj.surrogate_) ss << " SURROGATE!";
  ss << "\n CALO:  "; key_id(obj.calo_,ss); pt_eta_phi(obj.calo_,ss);
  ss << "\n GSF:   "; key_id(obj.gsf_,ss); pt_eta_phi(obj.gsf_,ss);
  ss << "\n PFGSF: "; key_id(obj.pfgsf_,ss); pt_eta_phi(obj.pfgsf_,ss);
  ss << "\n ELE:   "; key_id(obj.ele_,ss); pt_eta_phi(obj.ele_,ss);

  ss << "\n SEED:  "; key_id(obj.seed_,ss);
  ss << ", TRK driven: " << obj.seed_tracker_driven_
     << ", ECAL driven: " << obj.seed_ecal_driven_;
  //ss << "\n PREID: "; key_id(obj.preid_ecal_,ss);
  ss << "\n BDTs:  " << "unbiased: " << std::setprecision(4) << obj.unbiased_ 
     << ", ptbiased: " << std::setprecision(4) << obj.ptbiased_
     << ", ID: " << std::setprecision(4) << obj.id_;
  ss << "\n MATCH: "
     << "trk: " << obj.trk_match_ << "/" << std::setprecision(4) << obj.trk_dr_ 
     << ", gsf: " << obj.gsf_match_ << "/" << std::setprecision(4) << obj.gsf_dr_ 
     << ", ele: " << obj.ele_match_ << "/" << std::setprecision(4) << obj.ele_dr_;

  out << ss.str();
  return out;
};
