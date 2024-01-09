#ifndef LowPtElectrons_LowPtElectrons_Common
#define LowPtElectrons_LowPtElectrons_Common

//#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
//#include "DataFormats/Common/interface/Association.h"
//#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
//#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "FWCore/Framework/interface/EDFilter.h" // EDAnalyzer.h
//#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "FWCore/Utilities/interface/EDGetToken.h"
//#include "LowPtElectrons/LowPtElectrons/interface/IDNtuple.h"
#include "TRandom3.h"
#include "TTree.h"
#include <set>
#include <vector>
#include <math.h>
#include <boost/core/demangle.hpp>
#include <algorithm>
#include <random>

////////////////////////////////////////////////////////////////////////////////
//
namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ref<CaloClusterCollection> CaloClusterRef; }
namespace reco { typedef edm::Ptr<PreId> PreIdPtr; }
namespace reco { typedef edm::Ptr<ElectronSeed> ElectronSeedPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }
namespace  pat { typedef edm::Ptr<PackedCandidate> PackedCandidatePtr; }

typedef std::map<unsigned long,int> PdgIds;
typedef std::pair<float,float> TagMuon;

namespace id {
  static constexpr size_t NHITS_MAX = 30;
  static constexpr int NEG_INT = -10;
  static constexpr float NEG_FLOAT = -10.;
  static constexpr float NEG_FLOATSQ = -1.*NEG_FLOAT*NEG_FLOAT;
}

////////////////////////////////////////////////////////////////////////////////
//
template <class T1, class T2> 
class DeltaR2 {
public:
  DeltaR2( const edm::Ptr<T1>& obj1, const edm::Ptr<T2>& obj2, const double dr2 ) {
    obj1_ = obj1; obj2_ = obj2; dr2_ = dr2;
  };
  edm::Ptr<T1> obj1_;
  edm::Ptr<T2> obj2_;
  double dr2_ = id::NEG_FLOATSQ; // Use ^2 because this is deltaR^2
  static bool compare_by_dr2( const DeltaR2<T1,T2>& a, const DeltaR2<T1,T2>& b ) {
    return a.dr2_ < b.dr2_;
  };
};

typedef DeltaR2<reco::Candidate,reco::Track> SigToTrkDR2;
typedef DeltaR2<reco::Candidate,reco::GsfTrack> SigToGsfDR2;
typedef DeltaR2<reco::Candidate,reco::GsfElectron> SigToEleDR2;
typedef DeltaR2<reco::Track,reco::GsfTrack> TrkToGsfDR2;
typedef DeltaR2<reco::GsfTrack,reco::GsfElectron> GsfToEleDR2;
typedef DeltaR2<reco::Track,reco::GsfElectron> TrkToEleDR2;
typedef DeltaR2<reco::GsfTrack,reco::GsfTrack> GsfToGsfDR2;

////////////////////////////////////////////////////////////////////////////////
//
class ElectronChain {

public :

  explicit ElectronChain() {;}
  ~ElectronChain() {;}

public:
  
  int is_mc_ = -1;
  int is_aod_ = -1;

  bool is_e_ = false;
  bool is_egamma_ = false;

  float tag_pt_ = id::NEG_FLOAT;
  float tag_eta_ = id::NEG_FLOAT;

  // "Signal electron" info 
  reco::CandidatePtr sig_;

  // Track info
  reco::TrackPtr trk_;
  float trk_dr_ = id::NEG_FLOAT;
  bool trk_match_ = false;
  int pdg_id_ = 0;
  bool surrogate_ = false;

  // CaloCluster info
  reco::CaloClusterPtr calo_;

  // ElectronSeed info
  reco::ElectronSeedPtr seed_;
  bool seed_tracker_driven_ = false;
  bool seed_ecal_driven_ = false;
  
  // PreId info
  reco::PreIdPtr preid_ecal_;
  reco::PreIdPtr preid_hcal_;

  // Seed BDTs
  float unbiased_ = id::NEG_FLOAT;
  float ptbiased_ = id::NEG_FLOAT;

  // GsfTrack info
  reco::GsfTrackPtr gsf_;
  float gsf_dr_ = id::NEG_FLOAT;
  bool gsf_match_ = false;

  // PF GsfTrack info
  reco::GsfTrackPtr pfgsf_;
  float pfgsf_dr_ = id::NEG_FLOAT;
  bool pfgsf_match_ = false;

  // GsfElectron info
  reco::GsfElectronPtr ele_;
  float ele_dr_ = id::NEG_FLOAT;
  bool ele_match_ = false;
  float id_ = id::NEG_FLOAT;

  // Inner GSF track P4: defines reference coords in eta-phi space 

  float gsf_ref_eta_ = id::NEG_FLOAT;
  float gsf_ref_phi_ = id::NEG_FLOAT;
  float gsf_ref_R_ = id::NEG_FLOAT;
  float gsf_ref_p_ = id::NEG_FLOAT;
  float gsf_ref_pt_ = id::NEG_FLOAT;

  // GEN channel (2 points: inner P4 + inner P4 extrapolated to ECAL)

  float gen_inner_eta_ = id::NEG_FLOAT;
  float gen_inner_phi_ = id::NEG_FLOAT;
  float gen_inner_R_ = id::NEG_FLOAT;
  float gen_inner_p_ = id::NEG_FLOAT;
  float gen_inner_pt_ = id::NEG_FLOAT;

  float gen_proj_eta_ = id::NEG_FLOAT;
  float gen_proj_phi_ = id::NEG_FLOAT;
  float gen_proj_R_ = id::NEG_FLOAT;

  // GSF channel (3 points: inner P4 + inner P4 extrapolated to ECAL + outer P4 at ECAL surface)

  float gsf_inner_eta_ = id::NEG_FLOAT;
  float gsf_inner_phi_ = id::NEG_FLOAT;
  float gsf_inner_R_ = id::NEG_FLOAT;
  float gsf_inner_p_ = id::NEG_FLOAT;
  float gsf_inner_pt_ = id::NEG_FLOAT;

  float gsf_proj_eta_ = id::NEG_FLOAT;
  float gsf_proj_phi_ = id::NEG_FLOAT;
  float gsf_proj_R_ = id::NEG_FLOAT;
  float gsf_proj_p_ = id::NEG_FLOAT;

  float gsf_atcalo_eta_ = id::NEG_FLOAT;
  float gsf_atcalo_phi_ = id::NEG_FLOAT;
  float gsf_atcalo_R_ = id::NEG_FLOAT;
  float gsf_atcalo_p_ = id::NEG_FLOAT;

  // Cluster constituents of SC
  std::vector<float> clu_eta_;
  std::vector<float> clu_phi_;
  std::vector<float> clu_e_;
  std::vector<int> clu_nhit_;

  // PF candidates
  std::vector<float> pf_eta_;
  std::vector<float> pf_phi_;
  std::vector<float> pf_p_;
  std::vector<int> pf_pdgid_;
  std::vector<int> pf_matched_;
  std::vector<int> pf_lost_;

};

////////////////////////////////////////////////////////////////////////////////
//
template <class T1, class T2> 
  std::ostream& operator<< ( std::ostream& out, const DeltaR2<T1,T2>& obj ) {
out << "Class type:        " << boost::core::demangle( typeid(obj).name() ) << "\n"
      << "  Obj1 type:       " << boost::core::demangle( typeid(obj.obj1_).name() ) << "\n"
      << "  Obj2 type:       " << boost::core::demangle( typeid(obj.obj2_).name() ) << "\n";
  out << "  Obj1 id/key:     "; 
  if ( obj.obj1_.isNull() || !obj.obj1_.isAvailable() ) { out << "InvalidKey"; }
  else { out << obj.obj1_.id() << "/" << obj.obj1_.key(); }
  out << "\n";
  out << "  Obj2 id/key:     ";
  if ( obj.obj2_.isNull() || !obj.obj2_.isAvailable() ) { out << "InvalidKey"; }
  else { out << obj.obj2_.id() << "/" << obj.obj2_.key(); }
  out << "\n";
  out << "  Obj1 pt/eta/phi: ";
  if ( obj.obj1_.isNull() || !obj.obj1_.isAvailable() ) { out << "InvalidKey"; }
  else { out << obj.obj1_->pt()  << "/" 
	     << obj.obj1_->eta() << "/"
	     << obj.obj1_->phi(); } 
  out  << "\n";
  out << "  Obj2 pt/eta/phi: ";
  if ( obj.obj2_.isNull() || !obj.obj2_.isAvailable() ) { out << "InvalidKey"; }
  else { out << obj.obj2_->pt()  << "/" 
	     << obj.obj2_->eta() << "/"
	     << obj.obj2_->phi(); } 
  out  << "\n";
  out << "  dR2(1,2):        " << obj.dr2_;
  return out;
};

////////////////////////////////////////////////////////////////////////////////
//
template <typename T> 
void key_id( const edm::Ptr<T>& ptr, std::stringstream& ss );

////////////////////////////////////////////////////////////////////////////////
//
template <typename T> 
void pt_eta_phi( const edm::Ptr<T>& ptr, std::stringstream& ss );

////////////////////////////////////////////////////////////////////////////////
//
void pt_eta_phi( const edm::Ptr<reco::GsfTrack>& ptr, std::stringstream& ss );

////////////////////////////////////////////////////////////////////////////////
//
void pt_eta_phi( const reco::CaloClusterPtr& ptr, std::stringstream& ss );

////////////////////////////////////////////////////////////////////////////////
//
std::ostream& operator<< ( std::ostream& out, const ElectronChain& obj );

#endif // LowPtElectrons_LowPtElectrons_Common
