#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
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
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDNtuple.h"
#include "TRandom3.h"
#include "TTree.h"
#include <set>
#include <vector>
#include <math.h>
#include <boost/core/demangle.hpp>
#include <algorithm>

namespace reco { typedef edm::Ref<CaloClusterCollection> CaloClusterRef; }
namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ptr<ElectronSeed> ElectronSeedPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
namespace reco { typedef edm::Ptr<PreId> PreIdPtr; }
namespace pat { typedef edm::Ptr<PackedCandidate> PackedCandidatePtr; }
typedef std::map<unsigned long,int> PdgIds;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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
  double dr2_ = IDNtuple::NEG_FLOAT;
  static bool compare_by_dr2( const DeltaR2<T1,T2>& a, const DeltaR2<T1,T2>& b ) {
    return a.dr2_ < b.dr2_;
  };
};

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
  out << "  dR(1,2):         " << ( obj.dr2_ >= 0. ? sqrt(obj.dr2_) : IDNtuple::NEG_FLOAT );
  return out;
};

typedef DeltaR2<reco::Candidate,reco::Track> SigToTrkDR2;
typedef DeltaR2<reco::Candidate,reco::GsfTrack> SigToGsfDR2;
typedef DeltaR2<reco::Candidate,reco::GsfElectron> SigToEleDR2;
typedef DeltaR2<reco::Track,reco::GsfTrack> TrkToGsfDR2;
typedef DeltaR2<reco::GsfTrack,reco::GsfElectron> GsfToEleDR2;
typedef DeltaR2<reco::Track,reco::GsfElectron> TrkToEleDR2;
typedef DeltaR2<reco::GsfTrack,reco::GsfTrack> GsfToGsfDR2;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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

  // "Signal electron" info 
  reco::CandidatePtr sig_;

  // Track info
  reco::TrackPtr trk_;
  float trk_dr_ = IDNtuple::NEG_FLOAT;
  bool trk_match_ = false;
  int pdg_id_ = 0;

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
  float unbiased_ = IDNtuple::NEG_FLOAT;
  float ptbiased_ = IDNtuple::NEG_FLOAT;

  // GsfTrack info
  reco::GsfTrackPtr gsf_;
  float gsf_dr_ = IDNtuple::NEG_FLOAT;
  bool gsf_match_ = false;

  // PF GsfTrack info
  reco::GsfTrackPtr pfgsf_;
  float pfgsf_dr_ = IDNtuple::NEG_FLOAT;
  bool pfgsf_match_ = false;

  // GsfElectron info
  reco::GsfElectronPtr ele_;
  float ele_dr_ = IDNtuple::NEG_FLOAT;
  bool ele_match_ = false;

};

template <typename T> 
void key_id( const edm::Ptr<T>& ptr, std::stringstream& ss ) {
  ss << "id/key: ";
  if ( ptr.isNull() || !ptr.isAvailable() ) { ss << "InvalidKey"; }
  else { 
    std::stringstream tmp; tmp << ptr.id();
    ss << std::setw(6) << tmp.str() << "/"
       << std::setw(3) << int( ptr.key() ); }
}

void pt_eta_phi( const edm::Ptr<reco::GsfTrack>& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << ", pt/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->ptMode() << ", " 
     << std::setw(4) << ptr->etaMode() << ", " 
     << std::setw(4) << ptr->phiMode();
}

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

void pt_eta_phi( const reco::CaloClusterPtr& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << ", Et/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->energy() / std::cosh(ptr->eta()) << ", " 
     << std::setw(4) << ptr->eta() << ", " 
     << std::setw(4) << ptr->phi();
}

std::ostream& operator<< ( std::ostream& out, const ElectronChain& obj ) {
  std::stringstream ss;
  ss << "ElectronChain:"
     << " is_egamma: " << obj.is_egamma_
     << " is_e: " << obj.is_e_
     << " is_mc: " << obj.is_mc_
     << " is_aod:  " << obj.is_aod_;
  ss << "\n SIG:   "; key_id(obj.sig_,ss); pt_eta_phi(obj.sig_,ss);
  ss << "\n TRK:   "; key_id(obj.trk_,ss); pt_eta_phi(obj.trk_,ss); 
  ss << ", PdgId: " << obj.pdg_id_;
  ss << "\n CALO:  "; key_id(obj.calo_,ss); pt_eta_phi(obj.calo_,ss);
  ss << "\n GSF:   "; key_id(obj.gsf_,ss); pt_eta_phi(obj.gsf_,ss);
  ss << "\n PFGSF: "; key_id(obj.pfgsf_,ss); pt_eta_phi(obj.pfgsf_,ss);
  ss << "\n ELE:   "; key_id(obj.ele_,ss); pt_eta_phi(obj.ele_,ss);
  ss << "\n SEED:  "; key_id(obj.seed_,ss);
  ss << ", TRK driven: " << obj.seed_tracker_driven_
     << ", ECAL driven: " << obj.seed_ecal_driven_;
  //ss << "\n PREID: "; key_id(obj.preid_ecal_,ss);
  ss << "\n BDTs:  " << "unbiased: " << std::setprecision(4) << obj.unbiased_ 
     << ", ptbiased: " << std::setprecision(4) << obj.ptbiased_;
  ss << "\n MATCH: "
     << "trk: " << obj.trk_match_ << "/" << std::setprecision(4) << obj.trk_dr_ 
     << ", gsf: " << obj.gsf_match_ << "/" << std::setprecision(4) << obj.gsf_dr_ 
     << ", ele: " << obj.ele_match_ << "/" << std::setprecision(4) << obj.ele_dr_;
  out << ss.str();
  return out;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
class IDNtuplizer : public edm::EDAnalyzer {
  
public:
  
  explicit IDNtuplizer( const edm::ParameterSet& );
  ~IDNtuplizer() {}
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
  
  ///////////////
  // Main methods
  ///////////////

  // Reads all collections from the Event
  void readCollections( const edm::Event&, const edm::EventSetup& );

  // Wraps other methods to provide a sample of "signal" electrons
  void signalElectrons( std::set<reco::CandidatePtr>& signal_electrons );
  
  // GEN-based method to provide a sample of "signal" electrons
  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B, 
			  float muon_pt = 5., float muon_eta = 2.5 );
  
  // Top-level method that creates ElectronChain objects for EGamma electrons
  void pfElectrons( std::set<reco::CandidatePtr>& signal_electrons,
		    std::vector<SigToTrkDR2>& sig2trk,
		    std::vector<SigToTrkDR2>& other_trk,
		    std::vector<GsfToGsfDR2>& gsf2pfgsf );

  // Method that creates ElectronChain objects for EGamma signal candidates
  void pfElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
			   std::vector<SigToTrkDR2>& sig2trk,
			   std::vector<SigToGsfDR2>& sig2gsf,
			   std::vector<SigToGsfDR2>& sig2pfgsf,
			   std::vector<GsfToGsfDR2>& gsf2pfgsf,
			   std::vector<GsfToEleDR2>& pfgsf2ele,
			   std::vector<SigToTrkDR2>& other_trk );

  // Method that creates ElectronChain objects for EGamma fake candidates
  void pfElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
			  std::vector<TrkToGsfDR2>& trk2gsf,
			  std::vector<GsfToGsfDR2>& gsf2pfgsf,
			  std::vector<GsfToEleDR2>& pfgsf2ele );
    
  // Top-level method that creates ElectronChain objects for low-pT electrons
  void lowPtElectrons( std::set<reco::CandidatePtr>& signal_electrons,
		       std::vector<SigToTrkDR2>& sig2trk,
		       std::vector<SigToTrkDR2>& other_trk,
		       std::vector<GsfToGsfDR2>& gsf2pfgsf );

  // Method that creates ElectronChain objects for low-pT electron signal candiates
  void lowPtElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
			      std::vector<SigToTrkDR2>& sig2trk,
			      std::vector<SigToGsfDR2>& sig2gsf,
			      std::vector<GsfToGsfDR2>& gsf2pfgsf,
			      std::vector<GsfToEleDR2>& gsf2ele,
			      std::vector<SigToTrkDR2>& other_trk );

  // Method that creates ElectronChain objects for low-pT electron fake candidates
  void lowPtElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
			     std::vector<TrkToGsfDR2>& trk2gsf,
			     std::vector<GsfToGsfDR2>& gsf2pfgsf,
			     std::vector<GsfToEleDR2>& gsf2ele );
    
  // Fills tree per ElectronChain object
  void fill( const edm::Event& event, const edm::EventSetup& setup );

  //////////////////
  // Utility methods
  //////////////////

  // Extracts TrackPtrs to common vector, obtained from AOD or mAOD collections via Handles
  void extractTrackPtrs();
  
  // Various navigation methods (high- to low-level objects)
  bool gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk, bool is_egamma );
  bool eleToGsf( reco::GsfElectronPtr& ele, reco::GsfTrackPtr& gsf );
  bool gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed );
  bool seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk );
  bool seedToCalo( reco::ElectronSeedPtr& seed, reco::CaloClusterPtr& calo );

  // Various navigation maps (low- to high-level objects)
  void trkToGsfLinks( std::vector<reco::TrackPtr>& ctfTracks,
		      edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::vector<TrkToGsfDR2>& trk2gsf );
  void gsfToEleLinks( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
		      std::vector<GsfToEleDR2>& gsf2ele );
  void gsfToPfGsfLinks( edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
			edm::Handle< std::vector<reco::GsfTrack> >& gsfTracksEGamma,
			std::vector<GsfToGsfDR2>& gsf2pfgsf );

  // Links "signal" electrons to reconstructed objects
  template <typename T> 
  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
		       std::vector< edm::Ptr<T> >& candidates,
		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand, 
		       bool append = false );

  // Wraps method above to allow use of Handle<View<T>>
  template <typename T> 
  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
		       edm::Handle< edm::View<T> >& candidates,
		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand, 
		       bool append = false );

  // Wraps method above to allow use of Handle<vector<T>>
  template <typename T> 
  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
		       edm::Handle< std::vector<T> >& candidates,
		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand, 
		       bool append = false );

  // Return by reference the 'cand' that is matched to 'sig' in 'sig2cand' map
  template <typename T>
  void match( reco::CandidatePtr& sig,
	      std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
	      edm::Ptr<T>& cand, float& dr, bool& match ); // pass by ref

  // Wrap reco::deltaR2 method for various types, specifically GsfTrack uses eta and phi 'Mode'
  template <typename T1, typename T2> 
  float deltaR2( edm::Ptr<T1>& cand1, edm::Ptr<T2>& cand2 ); // used by sigToCandLinks
  float deltaR2( reco::TrackPtr& trk, reco::GsfTrackPtr& gsf ); // overload
  float deltaR2( reco::GsfTrackPtr& pfgsf, reco::GsfTrackPtr& gsf ); // overload
  float deltaR2( reco::CandidatePtr& sig, reco::GsfTrackPtr& gsf ); // overload

  // Filter track candidates by quality flag, simple pass for all other cands
  template <typename T> 
  bool filterCand( edm::Ptr<T>& cand );
  bool filterCand( edm::Ptr<reco::Track>& trk );
  bool filterCand( edm::Ptr<reco::GsfTrack>& gsf );

  // Check is Ptr is valid and available
  template <typename T> bool validPtr( edm::Ptr<T>& ptr );
  
private:
  
  // Misc
  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDNtuple ntuple_;
  int verbose_;
  bool check_from_B_;
  double dr_max_; // Max DeltaR value considered
  double dr_threshold_; // Threshold for DeltaR matching
  double prescale_;
  int isAOD_;
  bool isMC_;
  double minTrackPt_;
  
  // Generic collections

  const edm::EDGetTokenT<double> rho_;
  edm::Handle<double> rhoH_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  edm::Handle<reco::BeamSpot> beamspotH_;

  const edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticles_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;

  const edm::EDGetTokenT< edm::View<reco::Track> > ctfTracks_; // AOD
  edm::Handle< edm::View<reco::Track> > ctfTracksH_;

  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > packedCands_; // MINIAOD
  edm::Handle< edm::View<pat::PackedCandidate> > packedCandsH_;
  
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > lostTracks_; // MINIAOD
  edm::Handle< edm::View<pat::PackedCandidate> > lostTracksH_;
  
  const edm::EDGetTokenT<EcalRecHitCollection> ebRecHits_;
  edm::Handle<EcalRecHitCollection> ebRecHitsH_;

  const edm::EDGetTokenT<EcalRecHitCollection> eeRecHits_;
  edm::Handle<EcalRecHitCollection> eeRecHitsH_;

  //noZS::EcalClusterLazyTools ecalTools_;
  
  const edm::EDGetTokenT<reco::SuperClusterCollection> barrelSCs_; // AOD
  edm::Handle<reco::SuperClusterCollection> barrelSCsH_;

  const edm::EDGetTokenT<reco::SuperClusterCollection> endcapSCs_; // AOD
  edm::Handle<reco::SuperClusterCollection> endcapSCsH_;

  // Low pT collections

  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeeds_; // AOD
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeedsH_;

  const edm::EDGetTokenT< std::vector<reco::PreId> > preIdsEcal_; // AOD
  edm::Handle< std::vector<reco::PreId> > preIdsEcalH_;

  const edm::EDGetTokenT< std::vector<reco::PreId> > preIdsHcal_; // AOD
  edm::Handle< std::vector<reco::PreId> > preIdsHcalH_;

  const edm::EDGetTokenT< edm::ValueMap<reco::PreIdRef> > preIdRefs_; // AOD
  edm::Handle< edm::ValueMap<reco::PreIdRef> > preIdRefsH_;

  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracks_; // AOD and MINIAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksH_;

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectrons_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsH_;

  const edm::EDGetTokenT< edm::Association<reco::TrackCollection> > gsfTrackLinks_; // AOD
  edm::Handle<edm::Association<reco::TrackCollection> > gsfTrackLinksH_;

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > packedCandLinks_; // MINIAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > packedCandLinksH_;

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > lostTrackLinks_; // MINIAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > lostTrackLinksH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaUnbiased_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaUnbiasedH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaPtbiased_; //  on the fly?
  edm::Handle< edm::ValueMap<float> > mvaPtbiasedH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPt_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtH_;

  // EGamma collections

  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeedsEGamma_; // AOD
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeedsEGammaH_; // AOD

  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_; // AOD
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_MAOD_; // MINIAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksEGammaH_;

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectronsEGamma_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectronsEGamma_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGamma_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<bool> > mvaIdEGamma_; // on the fly?
  edm::Handle< edm::ValueMap<bool> > mvaIdEGammaH_;

  // Conversions

  //@@ const edm::EDGetTokenT< edm::ValueMap<float> > convVtxFitProb_;

  std::vector<ElectronChain> chains_;

  std::vector<reco::TrackPtr> tracks_;
  PdgIds pdgids_;

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
IDNtuplizer::IDNtuplizer( const edm::ParameterSet& cfg ) 
  : tree_(nullptr),
    ntuple_(),
    verbose_(cfg.getParameter<int>("verbose")),
    check_from_B_(cfg.getParameter<bool>("checkFromB")),
    dr_max_(cfg.getParameter<double>("drMax")),
    dr_threshold_(cfg.getParameter<double>("drThreshold")),
    prescale_(cfg.getParameter<double>("prescale")),
    isAOD_(-1),
    isMC_(true),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    // Generic collections
    rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
    rhoH_(),
    beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
    beamspotH_(),
    genParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
    prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
    genParticlesH_(),
    ctfTracks_(consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
    ctfTracksH_(),
    packedCands_(consumes< edm::View<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedCands"))),
    packedCandsH_(),
    lostTracks_(consumes< edm::View<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("lostTracks"))),
    lostTracksH_(),
    ebRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ebRecHits"))),
    ebRecHitsH_(),
    eeRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("eeRecHits"))),
    eeRecHitsH_(),
    //ecalTools_(),
    barrelSCs_(consumes<reco::SuperClusterCollection>(cfg.getParameter<edm::InputTag>("barrelSuperClusters"))),
    barrelSCsH_(),
    endcapSCs_(consumes<reco::SuperClusterCollection>(cfg.getParameter<edm::InputTag>("endcapSuperClusters"))),
    endcapSCsH_(),
    // Low pT collections
    eleSeeds_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeeds"))),
    eleSeedsH_(),
    preIdsEcal_(consumes< std::vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preIdsEcal"))),
    preIdsEcalH_(),
    preIdsHcal_(consumes< std::vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preIdsHcal"))),
    preIdsHcalH_(),
    preIdRefs_(consumes< edm::ValueMap<reco::PreIdRef> >(cfg.getParameter<edm::InputTag>("preIdRefs"))),
    preIdRefsH_(),
    gsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
    gsfTracksH_(),
    gsfElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
    patElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("patElectrons"))),
    gsfElectronsH_(),
    gsfTrackLinks_(consumes< edm::Association<reco::TrackCollection> >(cfg.getParameter<edm::InputTag>("gsfTrackLinks"))),
    gsfTrackLinksH_(),
    packedCandLinks_(consumes< edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("packedCandLinks"))),
    packedCandLinksH_(),
    lostTrackLinks_(consumes< edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("lostTrackLinks"))),
    lostTrackLinksH_(),
    mvaUnbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaUnbiased"))),
    mvaUnbiasedH_(),
    mvaPtbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaPtbiased"))),
    mvaPtbiasedH_(),
    mvaValueLowPt_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueLowPt"))),
    mvaValueLowPtH_(),
    // EGamma collections
    eleSeedsEGamma_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeedsEGamma"))),
    eleSeedsEGammaH_(),
    gsfTracksEGamma_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma"))),
    gsfTracksEGamma_MAOD_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma_MAOD"))),
    gsfTracksEGammaH_(),
    gsfElectronsEGamma_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectronsEGamma"))),
    patElectronsEGamma_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("patElectronsEGamma"))),
    gsfElectronsEGammaH_(),
    mvaValueEGamma_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueEGamma"))),
    mvaValueEGammaH_(),
    mvaIdEGamma_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("mvaIdEGamma"))),
    mvaIdEGammaH_(),
    // Conversions
    //convVtxFitProb_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("convVtxFitProb")))
    chains_(),
    tracks_(),
    pdgids_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "[IDNtuplizer::IDNtuplizer] Verbosity level: "<< verbose_ << std::endl;
  }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Initialise the weights LUT to filter fake tracks
void IDNtuplizer::beginRun( const edm::Run& run, const edm::EventSetup& es ) {
  //@@ ?
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::analyze( const edm::Event& event, const edm::EventSetup& setup ) {

  // Update all handles - MUST be called every event!
  readCollections(event,setup);
  
  // Clear ElectronChain vector and populate
  chains_.clear();

  // Identify signal electrons, from (GEN) MC or data control regions
  std::set<reco::CandidatePtr> signal_electrons;
  signalElectrons(signal_electrons);

  // Populate std::vector<reco::TrackPtr> tracks_ (from reco::Tracks, PF candidates, lost tracks)
  // Populate std::map<unsigned long,int> pdgids_ (typedef'ed to PdgIds)
  extractTrackPtrs();

  // Match "signal" electrons to generalTracks
  std::vector<SigToTrkDR2> sig2trk;
  std::vector<SigToTrkDR2> other_trk;
  sigToCandLinks<reco::Track>( signal_electrons, tracks_, sig2trk, other_trk );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::analyze] sigToCandLinks<reco::Track>:" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " tracks_.size(): " << tracks_.size() << std::endl
	      << " sig2trk.size(): " << sig2trk.size() << std::endl
	      << " other_trk.size(): " << other_trk.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : sig2trk ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Hack: PF GsfTracks --> low-pT GsfTracks map
  std::vector<GsfToGsfDR2> gsf2pfgsf;
  gsfToPfGsfLinks( gsfTracksH_, //@@ low-pT GSF tracks!
		   gsfTracksEGammaH_, 
		   gsf2pfgsf );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::analyze] gsfToPfGsfLinks:" << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " gsf2pfgsf.size(): " << gsf2pfgsf.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : gsf2pfgsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Populate ElectronChain objects using low-pT electrons
  lowPtElectrons( signal_electrons, sig2trk, other_trk, gsf2pfgsf );

  // Populate ElectronChain objects using PF electrons
  pfElectrons( signal_electrons, sig2trk, other_trk, gsf2pfgsf );
  
  // Print ElectronChain objects
  if ( verbose_ > 0 ) {
    for ( auto iter : chains_ ) { std::cout << iter << std::endl; }
  }
  
  // Fill ntuple
  fill(event,setup);

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {

  // Low pT electrons (and identify if data or MC and RECO/AOD or MINIAOD)
  if ( isAOD_ == -1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_);
    if ( gsfElectronsH_.isValid() ) {
      isAOD_ = 1;
      std::cout << "[IDNtuplizer::readCollections] File contains AOD data tier!" << std::endl;
    } else {
      event.getByToken(patElectrons_,gsfElectronsH_);
      if ( gsfElectronsH_.isValid() ) { 
	isAOD_ = 0;
	std::cout << "[IDNtuplizer::readCollections] File contains MINIAOD data tier!" << std::endl;
      } else {
	throw cms::Exception(" Collection not found: ") 
	  << "[IDNtuplizer::readCollections] Failed to find a standard AOD or miniAOD particle collection " 
	  << std::endl;
      }
    }
  } else if ( isAOD_ == 1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_);
  } else if ( isAOD_ == 0 ) {
    event.getByToken(patElectrons_,gsfElectronsH_);
  } else {
    throw cms::Exception(" Invalid value for isAOD: ") 
      << isAOD_ 
      << std::endl;
  }

  // Generic collections 
  event.getByToken(rho_, rhoH_);
  event.getByToken(beamspot_, beamspotH_);
  
  // GEN particles
  if ( isMC_ ) {
    if ( isAOD_ == 1 ) { 
      event.getByToken(genParticles_, genParticlesH_);
      if ( !(genParticlesH_.isValid()) ) { 
	isMC_ = false;
	std::cout << "[IDNtuplizer::readCollections] No GEN info found in AOD data tier!" << std::endl;
      }
    } else if ( isAOD_ == 0 ) { 
      event.getByToken(prunedGenParticles_, genParticlesH_);
      if ( !(genParticlesH_.isValid()) ) { 
	isMC_ = false;
	std::cout << "[IDNtuplizer::readCollections] No GEN info found in MINIAOD data tier!" << std::endl;
      }
    }
  }

  // KF tracks
  if ( isAOD_ == 1 ) { 
    event.getByToken(ctfTracks_, ctfTracksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCands_,packedCandsH_);
    event.getByToken(lostTracks_,lostTracksH_);
  }

  // RecHits and SuperClusters
  if ( isAOD_ == 1 ) { 
    event.getByToken(ebRecHits_, ebRecHitsH_);
    event.getByToken(eeRecHits_, eeRecHitsH_);
    event.getByToken(barrelSCs_, barrelSCsH_);
    event.getByToken(endcapSCs_, endcapSCsH_);
    //ecalTools_ = noZS::EcalClusterLazyTools(event, setup, ebRecHitsH_, eeRecHitsH_);
  }

  // ElectronSeeds and PreIds
  if ( isAOD_ == 1 ) { 
    event.getByToken(eleSeeds_, eleSeedsH_); 
    event.getByToken(eleSeedsEGamma_, eleSeedsEGammaH_); 
    event.getByToken(preIdsEcal_, preIdsEcalH_); 
    event.getByToken(preIdsHcal_, preIdsHcalH_); 
    event.getByToken(preIdRefs_, preIdRefsH_); 
  }

  // GsfTracks
  event.getByToken(gsfTracks_, gsfTracksH_);

  // Links
  if      ( isAOD_ == 1 ) { 
    event.getByToken(gsfTrackLinks_, gsfTrackLinksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCandLinks_, packedCandLinksH_); 
    event.getByToken(lostTrackLinks_, lostTrackLinksH_); 
  }
  
  // EGamma collections 
  if      ( isAOD_ == 1 ) { event.getByToken(eleSeedsEGamma_, eleSeedsEGammaH_); }
  if      ( isAOD_ == 1 ) { event.getByToken(gsfTracksEGamma_, gsfTracksEGammaH_); }
  else if ( isAOD_ == 0 ) { event.getByToken(gsfTracksEGamma_MAOD_, gsfTracksEGammaH_); }
  if      ( isAOD_ == 1 ) { event.getByToken(gsfElectronsEGamma_, gsfElectronsEGammaH_); }
  else if ( isAOD_ == 0 ) { event.getByToken(patElectronsEGamma_, gsfElectronsEGammaH_); }

  // IDs
  event.getByToken(mvaUnbiased_, mvaUnbiasedH_);
  event.getByToken(mvaPtbiased_, mvaPtbiasedH_);
  event.getByToken(mvaValueLowPt_, mvaValueLowPtH_);
  event.getByToken(mvaValueEGamma_, mvaValueEGammaH_);
  event.getByToken(mvaIdEGamma_, mvaIdEGammaH_);

  // Conversions
  //event.getByToken(convVtxFitProb_, convVtxFitProbH_);

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::signalElectrons( std::set<reco::CandidatePtr>& signal_electrons ) {

  signal_electrons.clear();

  if ( isMC_ ) { // Identify "signal" (GEN) electrons from B decays
    std::set<reco::GenParticlePtr> electrons_from_B;
    genElectronsFromB(electrons_from_B);
    for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
  } else { // Identify "signal" electrons from data control regions
    //@@ FOR NOW, A HACK ...
    std::set<reco::GenParticlePtr> electrons_from_B;
    genElectronsFromB(electrons_from_B);
    for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// From AOD
void IDNtuplizer::genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,
				     float muon_pt, float muon_eta ) {
  
  electrons_from_B.clear();
  bool tag_side_muon = false;
  
  for ( size_t idx = 0; idx < genParticlesH_->size(); idx++ ) {
    
    reco::GenParticlePtr gen(genParticlesH_, idx);
    if ( !validPtr(gen) ) {
      std::cout << "[IDNtuplizer::genElectronsFromB] ERROR! GenParticlePtr:"
		<< " gen.isNull(): " << gen.isNull()
		<< " gen.isAvailable(): " << gen.isAvailable()
		<< std::endl;
      continue;
    }
    
    // Last copy of GEN electron 
    bool is_ele = std::abs(gen->pdgId()) == 11 && gen->isLastCopy(); //@@ not a method of Candidate
    
    // Does GEN ele comes from B decay?
    bool non_resonant = gen->numberOfMothers() >= 1 && gen->mother() &&   // has mother
      std::abs(gen->mother()->pdgId()) > 510 &&                           // mother is B
      std::abs(gen->mother()->pdgId()) < 546;                             // mother is B
    bool resonant = gen->numberOfMothers() >= 1 && gen->mother() &&       // has mother
      std::abs(gen->mother()->pdgId()) == 443 &&                          // mother is J/psi
      gen->mother()->numberOfMothers() >= 1 && gen->mother()->mother() && // has grandmother
      std::abs(gen->mother()->mother()->pdgId()) > 510 &&                 // grandmother is B
      std::abs(gen->mother()->mother()->pdgId()) < 546;                   // grandmother is B
    
    //  Check for tag side muon
    tag_side_muon |= ( std::abs(gen->pdgId()) == 13 && gen->isLastCopy() 
		       && 
		       ( ( gen->numberOfMothers() >= 1 &&
			   gen->mother() &&
			   std::abs(gen->mother()->pdgId()) > 510 &&
			   std::abs(gen->mother()->pdgId()) < 546 && 
			   gen->mother()->pt() > muon_pt && 
			   std::abs(gen->mother()->eta()) < muon_eta ) 
			 ||
			 ( gen->numberOfMothers() >= 1 && 
			   gen->mother() &&
			   gen->mother()->numberOfMothers() >= 1 && 
			   gen->mother()->mother() &&
			   std::abs(gen->mother()->mother()->pdgId()) > 510 &&
			   std::abs(gen->mother()->mother()->pdgId()) < 546 && 
			   gen->mother()->mother()->pt() > muon_pt &&
			   std::abs(gen->mother()->mother()->eta()) < muon_eta ) ) );
    
    // is coming from a B
    if ( is_ele && ( ( resonant || non_resonant ) || !check_from_B_ ) ) {
      electrons_from_B.insert(gen);
      if ( verbose_ > 0 ) {
	std::cout << "[IDNtuplizer::genElectronsFromB] "
		  << " #signal_electrons: " << electrons_from_B.size()
		  << " resonant? " << resonant
		  << " non resonant? " << non_resonant
		  << " tag-side muon? " << tag_side_muon
		  << std::endl;
      }
    }
    
  } // genParticles loop

  if ( !tag_side_muon ) { electrons_from_B.clear(); }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::lowPtElectrons( std::set<reco::CandidatePtr>& signal_electrons,
				  std::vector<SigToTrkDR2>& sig2trk,
				  std::vector<SigToTrkDR2>& other_trk,
				  std::vector<GsfToGsfDR2>& gsf2pfgsf ) {
  
  // Match "signal" electrons to low-pT GsfTracks
  std::vector<SigToGsfDR2> sig2gsf;
  std::vector<SigToGsfDR2> other_gsf;
  sigToCandLinks<reco::GsfTrack>( signal_electrons, 
				  gsfTracksH_, 
				  sig2gsf, 
				  other_gsf );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::lowPtElectrons] sigToCandLinks<reco::GsfTrack>:" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " sig2gsf.size(): " << sig2gsf.size() << std::endl
	      << " other_gsf.size(): " << other_gsf.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : sig2gsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Match Tracks (including surrogates) to low-pT GsfTracks
  std::vector<TrkToGsfDR2> trk2gsf;
  trkToGsfLinks( tracks_, 
		 gsfTracksH_, 
		 trk2gsf );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::lowPtElectrons] trkToGsfLinks:" << std::endl
	      << " tracks_.size(): " << tracks_.size() << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " trk2gsf.size(): " << trk2gsf.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : trk2gsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Match low-pT GsfTracks to low-pT GsfElectrons
  std::vector<GsfToEleDR2> gsf2ele;
  gsfToEleLinks( gsfTracksH_, 
		 gsfElectronsH_, 
		 gsf2ele );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::lowPtElectrons] gsfToEleLinks:" << std::endl 
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " gsfElectronsH_->size(): " << gsfElectronsH_->size() << std::endl
	      << " gsf2ele.size(): " << gsf2ele.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : gsf2ele ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  lowPtElectrons_signal( signal_electrons, sig2trk, sig2gsf, gsf2pfgsf, gsf2ele, other_trk );
  lowPtElectrons_fakes( other_trk, trk2gsf, gsf2pfgsf, gsf2ele );
    
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "signal electrons" 
void IDNtuplizer::lowPtElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
					 std::vector<SigToTrkDR2>& sig2trk,
					 std::vector<SigToGsfDR2>& sig2gsf,
					 std::vector<GsfToGsfDR2>& gsf2pfgsf,
					 std::vector<GsfToEleDR2>& gsf2ele,
					 std::vector<SigToTrkDR2>& other_trk ) {
  
  ////////////////////
  // Iterate through "signal electrons" 
  ////////////////////

  for ( auto sig : signal_electrons ) {

    // Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = true;
    chain.is_egamma_ = false;
    chain.sig_ = sig;

    // Store matches between "signal electron" and tracks
    match<reco::Track>(sig,
		       sig2trk,
		       chain.trk_, // by ref
		       chain.trk_dr_,
		       chain.trk_match_ );
    
    // Store matches between "signal electron" and GSF tracks
    match<reco::GsfTrack>(sig,
			  sig2gsf,
			  chain.gsf_, // by ref
			  chain.gsf_dr_,
			  chain.gsf_match_ );

    // If not matched to GsfTrack, then move onto next "signal electron"
    if ( !chain.gsf_match_ ) { continue; }
    // Otherwise ... 
    
    // Store ElectronSeed BDT discrimator outputs
    chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
    chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];

    // Update Track info (i.e. overwrite gen-trk match via deltaR with seed track)
    reco::TrackPtr trk; 
    if ( gsfToTrk(chain.gsf_,trk,false/*is_egamma*/) ) {
      chain.trk_ = trk; 
      chain.trk_match_ = true;
      chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
      PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
      if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
    } else {
      // If no seed track is found when GSF is matched to GEN, then reset (this shouldn't happen for low-pT)
      chain.trk_ = reco::TrackPtr(); 
      chain.trk_match_ = false;
      chain.trk_dr_ = IDNtuple::NEG_FLOAT;
      chain.pdg_id_ = 0;
    }

    // If "signal electron" is matched to a low-pT GSF track,
    // check if matched to a PF GSF track with dR < threshold ...
    auto match_pfgsf = std::find_if( gsf2pfgsf.begin(),
				     gsf2pfgsf.end(),
				     [chain](const GsfToGsfDR2& dr2) {
				       return dr2.obj1_ == chain.gsf_;
				     }
				     );
    if ( match_pfgsf != gsf2pfgsf.end() && 
	 validPtr(match_pfgsf->obj2_) && 
	 match_pfgsf->dr2_ < dr_threshold_*dr_threshold_ ) {
      // Update PF GSF track
      chain.pfgsf_ = match_pfgsf->obj2_;
      chain.pfgsf_match_ = true;
      chain.pfgsf_dr_ = match_pfgsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_pfgsf->dr2_);
    }
    
    // Store GsfElectron info
    auto match_ele = std::find_if( gsf2ele.begin(), 
				   gsf2ele.end(), 
				   [chain](const GsfToEleDR2& dr2) { 
				     return dr2.obj1_ == chain.gsf_; // find gsf
				   }
				   );
    if ( match_ele != gsf2ele.end() && validPtr(match_ele->obj2_) ) { 
      chain.ele_ = match_ele->obj2_; 
      chain.ele_match_ = true;;
      chain.ele_dr_ = sqrt(deltaR2(chain.sig_,chain.ele_));
    } else {
      chain.ele_ = reco::GsfElectronPtr();
      chain.ele_match_ = false;
      chain.ele_dr_ = IDNtuple::NEG_FLOAT;
    }

    if ( isAOD_ == 0 ) {
      // In the absence of AOD info, set as tracker-driven
      chain.seed_tracker_driven_ = true;
    } else if ( isAOD_ == 1 ) {
      // Store ElectronSeed info
      reco::ElectronSeedPtr seed;
      if ( gsfToSeed(chain.gsf_,seed) ) {
	chain.seed_ = seed;
	chain.seed_tracker_driven_ = seed->isTrackerDriven();
	chain.seed_ecal_driven_ = seed->isEcalDriven();
	// Store CaloCluster info
	reco::CaloClusterPtr calo;
	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; }
	// Store PreId info
	//chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
	//chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
      }
    }
    
  } // for ( auto sig : signal_electrons )

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "fake tracks" 
void IDNtuplizer::lowPtElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
					std::vector<TrkToGsfDR2>& trk2gsf,
					std::vector<GsfToGsfDR2>& gsf2pfgsf,
					std::vector<GsfToEleDR2>& gsf2ele ) {

  ////////////////////
  // Iterate through "fake tracks" 
  ////////////////////

  // Iterate through tracks
  for ( auto iter : other_trk ) {

    if ( !validPtr(iter.obj2_) ) { continue; } // Shouldn't happen

    // Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = false;
    chain.is_egamma_ = false;
    
    // Store Track info
    chain.trk_ = iter.obj2_;
    chain.trk_match_ = true;
    chain.trk_dr_ = IDNtuple::NEG_FLOAT;
    
    // Find matched GsfTrack match, either via ElectronSeed or just deltaR
    auto match_gsf = std::find_if( trk2gsf.begin(), 
				   trk2gsf.end(), 
				   [chain](const TrkToGsfDR2& dr2) { 
				     return dr2.obj1_ == chain.trk_; // find trk
				   }
				   );
    if ( match_gsf != trk2gsf.end() && validPtr(match_gsf->obj2_) ) {
      chain.gsf_ = match_gsf->obj2_;
      chain.gsf_match_ = true;
      chain.gsf_dr_ = match_gsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_gsf->dr2_); // deltaR(trk,gsf)
    } else {
      if ( verbose_ > 3 ) {
	std::cout << "[IDNtuplizer::lowPtElectrons] INFO! Cannot match TrackPtr to GsfTrackPtr:"
		  << " TrackPtr: " << chain.trk_.id() << "/" << chain.trk_.key()
		  << std::endl;
      }
    }

    // If not matched to GsfTrack, then move on to next "signal electron"
    if ( !chain.gsf_match_ ) { continue; }
    // Otherwise ... 

    // Store ElectronSeed BDT discrimator outputs
    chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
    chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];

    // If "signal electron" is matched to a low-pT GSF track,
    // check if matched to a PF GSF track with dR < threshold ...
    auto match_pfgsf = std::find_if( gsf2pfgsf.begin(),
				     gsf2pfgsf.end(),
				     [chain](const GsfToGsfDR2& dr2) {
				       return dr2.obj1_ == chain.gsf_;
				     }
				     );
    if ( match_pfgsf != gsf2pfgsf.end() && 
	 validPtr(match_pfgsf->obj2_) && 
	 match_pfgsf->dr2_ < dr_threshold_*dr_threshold_ ) {
      // Update PF GSF track
      chain.pfgsf_ = match_pfgsf->obj2_;
      chain.pfgsf_match_ = true;
      chain.pfgsf_dr_ = match_pfgsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_pfgsf->dr2_);
    }
    
    // Store GsfElectron info
    auto match_ele = std::find_if( gsf2ele.begin(), 
				   gsf2ele.end(), 
				   [chain](const GsfToEleDR2& dr2) { 
				     return dr2.obj1_ == chain.gsf_; // find gsf
				   }
				   );
    if ( match_ele != gsf2ele.end() && validPtr(match_ele->obj2_) ) { 
      chain.ele_ = match_ele->obj2_; 
      chain.ele_match_ = true;
      chain.ele_dr_ = match_ele->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_ele->dr2_);
    } else {
      chain.ele_ = reco::GsfElectronPtr();
      chain.ele_match_ = false;
      chain.ele_dr_ = IDNtuple::NEG_FLOAT;
    }

    if ( isAOD_ == 0 ) {
      // In the absence of AOD info, set as tracker-driven
      chain.seed_tracker_driven_ = true;
    } else if ( isAOD_ == 1 ) {
      // Store ElectronSeed info
      reco::ElectronSeedPtr seed;
      if ( gsfToSeed(chain.gsf_,seed) ) {
	chain.seed_ = seed;
	chain.seed_tracker_driven_ = seed->isTrackerDriven();
	chain.seed_ecal_driven_ = seed->isEcalDriven();
	// Store CaloCluster info
	reco::CaloClusterPtr calo;
	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; }
	// Store PreId info
	//chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
	//chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
      }
    }
    
  } // for ( auto iter : other_trk )

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::pfElectrons( std::set<reco::CandidatePtr>& signal_electrons,
			       std::vector<SigToTrkDR2>& sig2trk,
			       std::vector<SigToTrkDR2>& other_trk,
			       std::vector<GsfToGsfDR2>& gsf2pfgsf ) {
  
  // Match "signal" electrons to *low-pT* GsfTracks
  std::vector<SigToGsfDR2> sig2gsf;
  std::vector<SigToGsfDR2> other_gsf; // not used
  sigToCandLinks<reco::GsfTrack>( signal_electrons, 
				  gsfTracksH_, //@@ low-pT GSF tracks!
				  sig2gsf, 
				  other_gsf );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::pfElectrons] sigToCandLinks<reco::GsfTrack>:" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " sig2gsf.size(): " << sig2gsf.size() << std::endl
	      << " other_gsf.size(): " << other_gsf.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : sig2gsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Match Tracks (including surrogates) to *low-pT* GsfTracks
  std::vector<TrkToGsfDR2> trk2gsf;
  trkToGsfLinks( tracks_, 
		 gsfTracksH_, //@@ low-pT GSF tracks! 
		 trk2gsf );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::pfElectrons] trkToGsfLinks:" << std::endl
	      << " tracks_.size(): " << tracks_.size() << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " trk2gsf.size(): " << trk2gsf.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : trk2gsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }
  
  // Match "signal" electrons to PF GsfTracks
  std::vector<SigToGsfDR2> sig2pfgsf;
  std::vector<SigToGsfDR2> other_pfgsf;
  sigToCandLinks<reco::GsfTrack>( signal_electrons, 
				  gsfTracksEGammaH_,
				  sig2pfgsf, 
				  other_pfgsf );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::pfElectrons] sigToCandLinks<reco::GsfTrack>:" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " sig2pfgsf.size(): " << sig2pfgsf.size() << std::endl
	      << " other_pfgsf.size(): " << other_pfgsf.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : sig2pfgsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Match PF GsfTracks to PF GsfElectrons
  std::vector<GsfToEleDR2> pfgsf2ele;
  gsfToEleLinks( gsfTracksEGammaH_, 
		 gsfElectronsEGammaH_, 
		 pfgsf2ele );
  if ( verbose_ > 0 ) {
    std::cout << "[IDNtuplizer::pfElectrons] gsfToEleLinks:" << std::endl 
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " gsfElectronsEGammaH_->size(): " << gsfElectronsEGammaH_->size() << std::endl
	      << " pfgsf2ele.size(): " << pfgsf2ele.size() << std::endl;
    if ( verbose_ > 1 ) {
      for ( auto iter : pfgsf2ele ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }
  
  pfElectrons_signal( signal_electrons, sig2trk, sig2gsf, sig2pfgsf, gsf2pfgsf, pfgsf2ele, other_trk );
  pfElectrons_fakes( other_trk, trk2gsf, gsf2pfgsf, pfgsf2ele );
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "signal electrons" 
void IDNtuplizer::pfElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
				      std::vector<SigToTrkDR2>& sig2trk,
				      std::vector<SigToGsfDR2>& sig2gsf,
				      std::vector<SigToGsfDR2>& sig2pfgsf,
			       	      std::vector<GsfToGsfDR2>& gsf2pfgsf,
				      std::vector<GsfToEleDR2>& pfgsf2ele,
				      std::vector<SigToTrkDR2>& other_trk ) {

  for ( auto sig : signal_electrons ) {

    // Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = true;
    chain.is_egamma_ = true;
    chain.sig_ = sig;

    // Store matches between "signal electron" and tracks
    match<reco::Track>(sig,
		       sig2trk,
		       chain.trk_, // by ref
		       chain.trk_dr_,
		       chain.trk_match_ );
    
    // Store matches between "signal electron" and *low-pT* GSF tracks
    match<reco::GsfTrack>(sig,
			  sig2gsf,
			  chain.gsf_, // by ref
			  chain.gsf_dr_,
			  chain.gsf_match_ );
	
    // Store matches between "signal electron" and PF GSF tracks
    match<reco::GsfTrack>(sig, 
			  sig2pfgsf,
			  chain.pfgsf_, // by ref 
			  chain.pfgsf_dr_,
			  chain.pfgsf_match_ );

    // If "signal electron" is matched to a *low-pT* GSF track,
    // which is matched to a PF GSF track with dR < threshold ...
    if ( chain.gsf_match_ ) {

      // Store ElectronSeed BDT discrimator outputs
      chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
      chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
      
      auto match_gsf = std::find_if( gsf2pfgsf.begin(),
				     gsf2pfgsf.end(),
				     [chain](const GsfToGsfDR2& dr2) {
				       return dr2.obj1_ == chain.gsf_;
				     }
				     );
      if ( match_gsf != gsf2pfgsf.end() && 
	   validPtr(match_gsf->obj2_) && 
	   match_gsf->dr2_ < dr_threshold_*dr_threshold_ ) {
	// .. Then:

	// 1) Update PF GSF track
	chain.pfgsf_ = match_gsf->obj2_;
	chain.pfgsf_match_ = true;
	chain.pfgsf_dr_ = match_gsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_gsf->dr2_);

	// 2) Update Track info (i.e. overwrite deltaR match via with seed track)
	reco::TrackPtr trk; 
	if ( gsfToTrk(chain.gsf_,trk,false/*is_egamma*/) ) {
	  chain.trk_ = trk; 
	  chain.trk_match_ = true;
	  chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
	  PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
	  if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
	} else {
	  // If no seed track is found, then reset (shouldn't happen for low-pT)
	  chain.trk_ = reco::TrackPtr(); 
	  chain.trk_match_ = false;
	  chain.trk_dr_ = IDNtuple::NEG_FLOAT;
	  chain.pdg_id_ = 0;
	}
	
	// 3) Store ElectronSeed info
	chain.seed_tracker_driven_ = true;
	
      } 

    } else if  ( chain.pfgsf_match_ && !other_trk.empty() ) {
      
      // Identify "best" surrogate track to be the closest in pT to the GSF track
      reco::GsfTrackPtr gsf = chain.pfgsf_;
      auto best = std::min_element(tracks_.begin(),
				   tracks_.end(),
				   [gsf]( const reco::TrackPtr& trk1, 
					  const reco::TrackPtr& trk2 ) {
				     return				\
				     std::abs(trk1->pt()-gsf->ptMode())/gsf->ptMode() \
				     <					\
				     std::abs(trk2->pt()-gsf->ptMode())/gsf->ptMode();
				   }
				   );
      if ( best != tracks_.end() && validPtr(*best) ) {
	chain.trk_ = *best; // Store surrogate track
	chain.trk_match_ = true;
	chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
	PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
	if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
      } else {
	std::cout << "[IDNtuplizer::pfElectrons] " 
		  << "ERROR: Could find a valid 'best' surrogate track matched in pT to the PF GSF track!";
      }

      // Set as ECAL-driven (to record this unusual case)
      chain.seed_ecal_driven_ = true;

    } else { continue; } // if no matches, move onto next "signal electron"
    
    // Store GsfElectron info
    auto match_pfgsf = std::find_if( pfgsf2ele.begin(), 
				     pfgsf2ele.end(), 
				     [chain](const GsfToEleDR2& dr2) { 
				       return dr2.obj1_ == chain.pfgsf_; // find pfgsf
				     }
				     );
    if ( match_pfgsf != pfgsf2ele.end() && validPtr(match_pfgsf->obj2_) ) { 
      chain.ele_ = match_pfgsf->obj2_; 
      chain.ele_match_ = true;;
      chain.ele_dr_ = sqrt(deltaR2(chain.sig_,chain.ele_));
    } else {
      chain.ele_ = reco::GsfElectronPtr();
      chain.ele_match_ = false;
      chain.ele_dr_ = IDNtuple::NEG_FLOAT;
    }

//    if ( isAOD_ == 1 ) {
//      // Store ElectronSeed info
//      reco::ElectronSeedPtr seed;
//      if ( gsfToSeed(chain.gsf_,seed) ) {
//	chain.seed_ = seed;
//	chain.seed_tracker_driven_ = seed->isTrackerDriven();
//	chain.seed_ecal_driven_ = seed->isEcalDriven();
//	// Store CaloCluster info
//	reco::CaloClusterPtr calo;
//	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; }
//	// Store PreId info
//	//chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
//	//chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
//      }
//    }
    
  } // for ( auto sig : signal_electrons )

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "fake tracks" 
void IDNtuplizer::pfElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
				     std::vector<TrkToGsfDR2>& trk2gsf,
				     std::vector<GsfToGsfDR2>& gsf2pfgsf,
				     std::vector<GsfToEleDR2>& pfgsf2ele ) {
  
  // Iterate through tracks
  for ( auto iter : other_trk ) {
    
    if ( !validPtr(iter.obj2_) ) { continue; } // Shouldn't happen
    
    // Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = false;
    chain.is_egamma_ = true;
    
    // Store Track info
    chain.trk_ = iter.obj2_;
    chain.trk_match_ = true;
    chain.trk_dr_ = IDNtuple::NEG_FLOAT;
    
    // Find matched GsfTrack match, either via ElectronSeed or just deltaR
    auto match_gsf = std::find_if( trk2gsf.begin(), 
				   trk2gsf.end(), 
				   [chain](const TrkToGsfDR2& dr2) { 
				     return dr2.obj1_ == chain.trk_; // find trk
				   }
				   );
    if ( match_gsf != trk2gsf.end() && validPtr(match_gsf->obj2_) ) {
      
      // Store *low-pT* GSF track info
      chain.gsf_ = match_gsf->obj2_;
      chain.gsf_match_ = true;
      chain.gsf_dr_ = match_gsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(deltaR2(chain.gsf_,chain.trk_));
      
      // Store ElectronSeed BDT discrimator outputs
      chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
      chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];

      if ( isAOD_ == 0 ) {
	// In the absence of AOD info, set as tracker-driven
	chain.seed_tracker_driven_ = true;
      }
      
      // Attempt to match to PF GSF track
      auto match_gsf = std::find_if( gsf2pfgsf.begin(),
				     gsf2pfgsf.end(),
				     [chain](const GsfToGsfDR2& dr2) {
				       return dr2.obj1_ == chain.gsf_;
				     }
				     );
      if ( match_gsf != gsf2pfgsf.end() && 
	   validPtr(match_gsf->obj2_) && 
	   match_gsf->dr2_ < dr_threshold_*dr_threshold_ ) {
	
	// Store PF GSF track
	chain.pfgsf_ = match_gsf->obj2_;
	chain.pfgsf_match_ = true;
	chain.pfgsf_dr_ = match_gsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(deltaR2(chain.pfgsf_,chain.trk_));
	
      }

      // Store GsfElectron info
      auto match_pfgsf = std::find_if( pfgsf2ele.begin(), 
				       pfgsf2ele.end(), 
				       [chain](const GsfToEleDR2& dr2) { 
					 return dr2.obj1_ == chain.pfgsf_; // find pfgsf
				       }
				       );
      if ( match_pfgsf != pfgsf2ele.end() && validPtr(match_pfgsf->obj2_) ) { 
	chain.ele_ = match_pfgsf->obj2_; 
	chain.ele_match_ = true;;
	chain.ele_dr_ = match_pfgsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(deltaR2(chain.ele_,chain.trk_));
      }
      
    } else {
      if ( verbose_ > 3 ) {
	std::cout << "[IDNtuplizer::pfElectrons] INFO! Cannot match TrackPtr to GsfTrackPtr:"
		  << " TrackPtr: " << chain.trk_.id() << "/" << chain.trk_.key()
		  << std::endl;
      }
    }
    
  } // for ( auto iter : other_trk )

  // Iterate through unmatched PF GSF tracks and assign surrogate track
  for ( auto iter : gsf2pfgsf ) {
    
    // If set, apply prescale to select a random subset of PF GSF tracks
    if ( prescale_ > 0. && ( gRandom->Rndm() > prescale_ ) ) { continue; }

    if ( !validPtr(iter.obj2_) || 
	 iter.dr2_ < dr_threshold_*dr_threshold_ ) { continue; }
    
    // Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = false;
    chain.is_egamma_ = true;

    // Identify "best" surrogate track to be the closest in pT to the GSF track
    reco::GsfTrackPtr gsf = iter.obj2_;
    auto best = std::min_element(tracks_.begin(),
				 tracks_.end(),
				 [gsf]( const reco::TrackPtr& trk1, 
					const reco::TrackPtr& trk2 ) {
				   return				\
				   std::abs(trk1->pt()-gsf->ptMode())/gsf->ptMode() \
				   <					\
				   std::abs(trk2->pt()-gsf->ptMode())/gsf->ptMode();
				 }
				 );
    if ( best != tracks_.end() && validPtr(*best) ) {
      chain.trk_ = *best; // Store surrogate track
      chain.trk_match_ = true;
      chain.trk_dr_ = IDNtuple::NEG_FLOAT;
      PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
      if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
    } else {
      std::cout << "[IDNtuplizer::pfElectrons]"
		<< " ERROR: Could find a valid 'best' surrogate track matched in pT to the PF GSF track!";
    }

    // Store PF GSF track
    chain.pfgsf_ = iter.obj2_;
    chain.pfgsf_match_ = true;
    chain.pfgsf_dr_ = iter.dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(deltaR2(chain.pfgsf_,chain.trk_));
    
    // Set as ECAL-driven (to record this unusual case)
    chain.seed_ecal_driven_ = true;

    // Store GsfElectron info
    auto match_pfgsf = std::find_if( pfgsf2ele.begin(), 
				     pfgsf2ele.end(), 
				     [chain](const GsfToEleDR2& dr2) { 
				       return dr2.obj1_ == chain.pfgsf_; // find pfgsf
				     }
				     );
    if ( match_pfgsf != pfgsf2ele.end() && validPtr(match_pfgsf->obj2_) ) { 
      chain.ele_ = match_pfgsf->obj2_; 
      chain.ele_match_ = true;;
      chain.ele_dr_ = match_pfgsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(deltaR2(chain.ele_,chain.trk_));
    }
    
  } // for ( auto iter : gsf2pfgsf )

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::fill( const edm::Event& event,
			const edm::EventSetup& setup ) {
  
  for ( auto chain : chains_ ) {
    
    // Init tree here
    ntuple_.reset();

    // Data sample
    ntuple_.is_mc( chain.is_mc_ );
    ntuple_.is_aod( chain.is_aod_ );

    // Event stuff
    ntuple_.fill_evt( event.id() );
    ntuple_.set_rho( *rhoH_ );

    // EGamma of low-pT
    ntuple_.is_egamma( chain.is_egamma_ );

    // Truth label
    ntuple_.is_e( chain.is_e_ );
    ntuple_.is_other( !chain.is_e_ );
    
    // Set background weight
    if ( !chain.is_e_ ) { ntuple_.set_weight( prescale_ > 0. ? prescale_ : 1. ); }

    // "signal" info
    if ( validPtr(chain.sig_) ) {
      ntuple_.fill_gen( chain.sig_ ); // Add sig info
    }

    // Track info
    if ( validPtr(chain.trk_) ) {
      ntuple_.has_trk( chain.trk_match_ );
      ntuple_.fill_trk( chain.trk_, *beamspotH_ );
      ntuple_.trk_dr( chain.trk_dr_ );
      ntuple_.pdg_id( chain.pdg_id_ );
    }

    // ElectronSeed
    if ( validPtr(chain.seed_) ) { ntuple_.has_seed(true); }

    // PreId
//    noZS::EcalClusterLazyTools ecal_tools(event, setup, ebRecHits_, eeRecHits_);
//    ntuple_.fill_preid( *chain.preid_ecal_,
//			*chain.preid_hcal_,
//			*beamspotH_,
//			*rhoH_,
//			ecal_tools );

    // Tracker- or ECAL driven? 
    ntuple_.fill_seed( chain.seed_tracker_driven_, 
		       chain.seed_ecal_driven_ );
    
    // GsfTrack info
    if ( validPtr(chain.gsf_) ) {
      if ( !chain.is_egamma_ ) { ntuple_.has_trk( chain.gsf_match_ ); }
      ntuple_.has_gsf( chain.gsf_match_ );
      ntuple_.fill_gsf( chain.gsf_, *beamspotH_ );
      ntuple_.gsf_dr( chain.gsf_dr_ );
      ntuple_.fill_bdt( chain.unbiased_, chain.ptbiased_ );
    }
    
    // PF GsfTrack info
    if ( validPtr(chain.pfgsf_) ) {
      if ( chain.is_egamma_ ) { ntuple_.has_trk( chain.pfgsf_match_ ); }
      ntuple_.has_pfgsf( chain.pfgsf_match_ );
      ntuple_.fill_pfgsf( chain.pfgsf_, *beamspotH_ );
      ntuple_.pfgsf_dr( chain.pfgsf_dr_ );
    }

    // GsfElectron info
    if ( validPtr(chain.ele_) ) {

      ntuple_.has_trk( chain.ele_match_ );
      if ( chain.is_egamma_ ) { ntuple_.has_pfgsf( chain.ele_match_ ); }
      else { ntuple_.has_gsf( chain.ele_match_ ); }
      ntuple_.has_ele( chain.ele_match_ );
      ntuple_.ele_dr( chain.ele_dr_ );

      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
      float mva_value = -999.;
      int mva_id = -999;
      if ( !chain.is_egamma_ ) {
	if ( mvaValueLowPtH_.isValid() && 
	     mvaValueLowPtH_->size() == gsfElectronsH_->size() ) {
	  mva_value = mvaValueLowPtH_->get( chain.ele_.key() );
	} else {
	  std::cout << "[IDNtuplizer::fill] ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
	}
      } else {
	if ( mvaValueEGammaH_.isValid() && 
	     mvaValueEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
	  mva_value = mvaValueEGammaH_->get( chain.ele_.key() );
	} else {
	  std::cout << "[IDNtuplizer::fill] ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
	}
	if ( mvaIdEGammaH_.isValid() && 
	     mvaIdEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
	  mva_id = mvaIdEGammaH_->get( chain.ele_.key() ) ? 1 : 0;
	} else {
	  std::cout << "[IDNtuplizer::fill] ERROR! Issue matching MVA ID to GsfElectrons!" << std::endl;
	}
      }
      
      //@@ dirty hack as is not in Event nor embedded in pat::Electron
      float conv_vtx_fit_prob = -999.;
      //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == gsfElectrons->size() ) {
      //  conv_vtx_fit_prob = convVtxFitProb->get( chain.ele_.key() );
      //}
      
      ntuple_.fill_ele( chain.ele_, mva_value, mva_id, conv_vtx_fit_prob, *rhoH_, chain.is_egamma_ );

    }

    //ntuple_.fill_supercluster(chain.ele_);
    
    tree_->Fill(); 
    
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Utility methods ... /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::extractTrackPtrs() {
  tracks_.clear();
  pdgids_.clear();
  if ( isAOD_ == 1 ) {
    for ( size_t itrk = 0; itrk < ctfTracksH_->size(); ++itrk ) { 
      reco::TrackPtr ptr(ctfTracksH_,itrk);
      if ( !filterCand<reco::Track>(ptr) ) { continue; }
      tracks_.push_back(ptr); 
      pdgids_.insert(PdgIds::value_type(ptr.key(),0));
    }
  } else if ( isAOD_ == 0 ) {
    size_t iptr = 0;
    for ( const auto& ptr : *packedCandsH_ ) {
      if ( ptr.bestTrack() == nullptr ) { continue; }
      reco::TrackPtr trk(ptr.bestTrack(),iptr);
      if ( !filterCand<reco::Track>(trk) ) { continue; }
      tracks_.push_back(trk);
      pdgids_.insert(PdgIds::value_type(trk.key(),ptr.pdgId()));
      ++iptr;
    }
    for ( const auto& ptr : *lostTracksH_ ) { 
      if ( ptr.bestTrack() == nullptr ) { continue; }
      reco::TrackPtr trk(ptr.bestTrack(),iptr);
      if ( !filterCand<reco::Track>(trk) ) { continue; }
      tracks_.push_back(trk); 
      pdgids_.insert(PdgIds::value_type(trk.key(),ptr.pdgId()));
      ++iptr;
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::gsfToTrk( reco::GsfTrackPtr& gsf, 
			    reco::TrackPtr& trk, 
			    bool is_egamma ) {
  
  // Shouldn't happen (but "link maps" may contain a invalid ptr)
  if ( !validPtr(gsf) ) {
    if ( verbose_ > 0 ) {
      std::cout << "[IDNtuplizer::gsfToTrk] ERROR! GsfTrackPtr:"
		<< " gsf.isNull(): " << gsf.isNull()
		<< " gsf.isAvailable(): " << gsf.isAvailable()
		<< std::endl;
    }
    return false;
  }
  
  // Attempt to navigate via Seed (and TrackExtra) to Track
  reco::ElectronSeedPtr seed;
  if ( gsfToSeed(gsf,seed) && seedToTrk(seed,trk) ) { return true; }
  //@@ else { return false; }

  // In the case of low-pT electrons ...
  if ( !is_egamma ) {
    // ... if above fails (e.g. TrackExtra missing), attempt to use following maps:
    //@@ Warning: linked tracks can be missing and/or duplicated!
    
    reco::GsfTrackRef gsf_ref(gsf.id(),gsf.get(),gsf.key());
    if ( gsf_ref.isNull() || !gsf_ref.isAvailable() ) { return false; }

    if ( isAOD_ == 1 ) {
      // ... gsf->track maps in AOD
      reco::TrackRef trk_ref = (*gsfTrackLinksH_)[gsf_ref];
      trk = edm::refToPtr(trk_ref);
      if ( validPtr(trk) ) { return true; }
    } else if ( isAOD_ == 0 ) {
      // ... gsf->packed maps in MINIAOD
      pat::PackedCandidateRef packed_ref = (*packedCandLinksH_)[gsf_ref];
      if ( packed_ref.isAvailable() && 
	   packed_ref.isNonnull() && 
	   packed_ref->bestTrack() != nullptr ) { 
	trk = reco::TrackPtr(packed_ref->bestTrack(),packed_ref.key());
	if ( validPtr(trk) ) { return true; }
      }
      // ... gsf->lost maps in MINIAOD
      pat::PackedCandidateRef lost_ref = (*lostTrackLinksH_)[gsf_ref];
      if ( lost_ref.isAvailable() && 
	   lost_ref.isNonnull() && 
	   lost_ref->bestTrack() != nullptr ) { 
	trk = reco::TrackPtr(lost_ref->bestTrack(),lost_ref.key());
	if ( validPtr(trk) ) { return true; }
      }
    }
  }

  return false;
}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::eleToGsf( reco::GsfElectronPtr& ele, reco::GsfTrackPtr& gsf ) {
  if ( !validPtr(ele) ) {
    std::cout << "[IDNtuplizer::eleToGsf] ERROR! GsfElectronPtr:"
	      << " ele.isNull(): " << ele.isNull()
	      << " ele.isAvailable(): " << ele.isAvailable()
	      << std::endl;
    return false;
  }
  gsf = edm::refToPtr(ele->gsfTrack());
  if ( !validPtr(gsf) ) {
    std::cout << "[IDNtuplizer::eleToGsf] ERROR! GsfTrackPtr:"
	      << " gsf.isNull(): " << gsf.isNull()
	      << " gsf.isAvailable(): " << gsf.isAvailable()
	      << std::endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed ) {
  if ( !validPtr(gsf) ) {
    if ( verbose_ > 0 ) {
      std::cout << "[IDNtuplizer::gsfToSeed] ERROR! GsfTrackPtr:"
		<< " gsf.isNull(): " << gsf.isNull()
		<< " gsf.isAvailable(): " << gsf.isAvailable()
		<< std::endl;
    }
    return false;
  }
  edm::RefToBase<TrajectorySeed> traj;
  if ( gsf->extra().isNonnull() && gsf->extra().isAvailable() ) { 
    traj = gsf->seedRef(); 
  } else {
    if ( verbose_ > 3 ) { // TrackExtra are not stored by default in MINIAOD
      std::cout << "[IDNtuplizer::gsfToSeed] ERROR: TrackExtra:" 
		<< " gsf->extra().isNull(): " << gsf->extra().isNull()
		<< " gsf->extra().isAvailable(): " << gsf->extra().isAvailable()
		<< std::endl; 
    }
    return false;
  }
  if ( traj.isNull() || !traj.isAvailable() ) { 
    if ( verbose_ > 0 ) {
      std::cout << "[IDNtuplizer::gsfToSeed] ERROR: TrajectorySeedRef:" 
		<< " traj.isNull(): " << traj.isNull()
		<< " traj.isAvailable(): " << traj.isAvailable()
		<< std::endl; 
    }
    return false;
  }
  seed = edm::refToPtr(traj.castTo<reco::ElectronSeedRef>());
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "[IDNtuplizer::gsfToSeed] ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk ) {
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "[IDNtuplizer::seedToTrk] ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  if ( !seed->isTrackerDriven() ) {
    if ( verbose_ > 3 ) {
      std::cout << "[IDNtuplizer::seedToTrk] INFO! ElectronSeedPtr:"
		<< " seed->isTrackerDriven(): " << seed->isTrackerDriven()
		<< std::endl;
    }
  }
  trk = edm::refToPtr(seed->ctfTrack());
  if ( !validPtr(trk) ) { 
    if ( verbose_ > 3 ) {
      std::cout << "[IDNtuplizer::seedToTrk] INFO! TrackPtr:"
		<< " trk.isNull(): " << trk.isNull()
		<< " trk.isAvailable(): " << trk.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::seedToCalo( reco::ElectronSeedPtr& seed, reco::CaloClusterPtr& calo ) {
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "[IDNtuplizer::seedToCalo] ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }

  edm::RefToBase<reco::CaloCluster> base = seed->caloCluster();
  if ( base.isNull() || !base.isAvailable() ) { 
    if ( verbose_ > 3 ) {
      std::cout << "[IDNtuplizer::seedToCalo] INFO! edm::RefToBase<reco::CaloCluster>:"
		<< " base.isNull(): " << base.isNull()
		<< " base.isAvailable(): " << base.isAvailable()
		<< std::endl;
    }
    return false;
  }
  calo = edm::refToPtr(seed->caloCluster().castTo<reco::CaloClusterRef>());
  if ( !validPtr(calo) ) { 
    std::cout << "[IDNtuplizer::seedToCalo] ERROR! CaloClusterPtr:"
	      << " calo.isNull(): " << calo.isNull()
	      << " calo.isAvailable(): " << calo.isAvailable()
	      << std::endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::trkToGsfLinks( std::vector<reco::TrackPtr>& ctfTracks,
				 edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				 std::vector<TrkToGsfDR2>& trk2gsf ) {

  trk2gsf.clear();

  // Record matched tracks via their keys
  std::vector<size_t> keys(ctfTracks.size());
  std::generate( keys.begin(), keys.end(), [n=0] () mutable { return n++; } );

  // Shuffled to allow unbiased selection of surrogate tracks
  std::random_shuffle( keys.begin(), keys.end() );

  //////////
  // 1) Add GsfTracks to trk2gsf map if link (via provenence information) can be made
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackPtr gsf(gsfTracks, idx); 

    // Check validity of GsfTrackPtr (shouldn't ever fail?!)
    if ( !validPtr(gsf) ) { continue; }
    
    // If gsf not linked to trk (via ElectronSeed or "link maps"), then continue
    reco::TrackPtr trk;
    if ( !gsfToTrk(gsf,trk,false) ) { continue; }

    // Check if gsf is already stored in map and, if not, add entry (warning: occasionally points to duplicate trks!)
    auto match_gsf = std::find_if( trk2gsf.begin(), 
				   trk2gsf.end(), 
				   [gsf](const TrkToGsfDR2& dr2) { 
				     return dr2.obj2_ == gsf; 
				   }
				   );
    if ( match_gsf == trk2gsf.end() ) {
      trk2gsf.emplace_back( trk, gsf, deltaR2(trk,gsf) );
      keys.erase( std::remove( keys.begin(), 
			       keys.end(), 
			       trk.key() ), 
		  keys.end() ); // move key to end and erase
    } else { 
      // This can happen on occasion with the "link maps"
      std::cout << "[IDNtuplizer::trkToGsfLinks] ERROR! GsfTrackPtr is already in the map!" << std::endl; 
    }
  }
  
  //////////
  // 2) Add entries with GsfTracks linked to (random) surrogate tracks if link cannot be made
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackPtr gsf(gsfTracks, idx); 

    if ( !validPtr(gsf) ) { continue; } //@@ shouldn't ever happen?!
    
    // If gsf is *linked* to trk (i.e. opposite to above), then continue
    reco::TrackPtr trk;
    if ( gsfToTrk(gsf,trk,false) ) { continue; }
    
    // Check if gsf already stored in map and, if not, add entry with surrogate track
    auto match_gsf = std::find_if( trk2gsf.begin(), 
				   trk2gsf.end(), 
				   [gsf](const TrkToGsfDR2& dr2) { 
				     return dr2.obj2_ == gsf; 
				   }
				   );
    if ( match_gsf == trk2gsf.end() ) {
      reco::TrackPtr trk(tracks_[keys.front()]); // use 1st track (in shuffled vector) as a surrogate ...
      keys.erase(keys.begin()); // ... and then remove
      trk2gsf.emplace_back( trk, gsf, IDNtuple::NEG_FLOAT ); //deltaR2(trk,gsf) );
    } else { std::cout << "[IDNtuplizer::trkToGsfLinks] ERROR! GsfTrackPtr is already in the map!" << std::endl; }
  }

  //////////
  // 3) Add (null) entries for all Tracks without a match to a GsfTrack
  for ( auto key : keys ) {
    reco::TrackPtr trk(tracks_[key]);

    // Check if trk already stored in map and, if not, add "empty" entry
    auto match_trk = std::find_if( trk2gsf.begin(), 
				   trk2gsf.end(), 
				   [trk](const TrkToGsfDR2& dr2) { 
				     return dr2.obj1_ == trk; 
				   }
				   );
    if ( match_trk == trk2gsf.end() ) {
      reco::GsfTrackPtr gsf;
      trk2gsf.emplace_back( trk, gsf, IDNtuple::NEG_FLOAT );
    }
  }    

}

////////////////////////////////////////////////////////////////////////////////
// 
void IDNtuplizer::gsfToEleLinks( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				 const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
				 std::vector<GsfToEleDR2>& gsf2ele ) {
  
  gsf2ele.clear();

  // Record matched GsfTracks via their keys
  std::vector<size_t> keys(gsfTracks->size());
  std::generate( keys.begin(), keys.end(), [n=0] () mutable { return n++; } );
  
  // 1) Add GsfElectrons to gsf2ele map if link can be made
  for ( size_t idx = 0; idx < gsfElectrons->size(); ++idx ) {
    reco::GsfElectronPtr ele(gsfElectrons, idx);
    if ( !validPtr(ele) ) { continue; } //@@ shouldn't ever happen?!
    
    // If ele not linked to trk, then continue
    reco::GsfTrackPtr gsf;
    if ( !eleToGsf(ele,gsf) ) { continue; }
    
    // Check if ele is already stored in map and, if not, add entry
    auto match_ele = std::find_if( gsf2ele.begin(), 
				   gsf2ele.end(), 
				   [ele](const GsfToEleDR2& dr2) { 
				     return dr2.obj2_ == ele; 
				   }
				   );
    if ( match_ele == gsf2ele.end() ) {
      gsf2ele.emplace_back( gsf, ele, deltaR2(gsf,ele) );
      keys.erase( std::remove( keys.begin(), 
			       keys.end(), 
			       gsf.key() ), 
		  keys.end() ); // move key to end and erase
    } else { std::cout << "[IDNtuplizer::trkToGsfLinks]" 
		       << " ERROR! GsfElectronPtr is already in the map!" << std::endl; }

  }

  // 2) Add (null) entries for all GsfTracks without a match to a GsfElectron
  for ( auto key : keys ) {
    reco::GsfTrackPtr gsf(gsfTracks, key);

    // Check if gsf already stored in map and, if not, add "empty" entry
    auto match_gsf = std::find_if( gsf2ele.begin(), 
				   gsf2ele.end(), 
				   [gsf](const GsfToEleDR2& dr2) { 
				     return dr2.obj1_ == gsf;
				   }
				   );
    if ( match_gsf == gsf2ele.end() ) {
      reco::GsfElectronPtr ele;
      gsf2ele.emplace_back( gsf, ele, IDNtuple::NEG_FLOAT );
    }
    
  } 
}

////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::gsfToPfGsfLinks( edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				   edm::Handle< std::vector<reco::GsfTrack> >& gsfTracksEGamma,
				   std::vector<GsfToGsfDR2>& gsf2pfgsf ) {

  if ( gsfTracks->empty() && gsfTracksEGamma->empty() ) { return; }

  // Record matched GsfTracks via their keys
  std::vector<size_t> keys(gsfTracks->size());
  std::generate( keys.begin(), keys.end(), [n=0] () mutable { return n++; } );
  
  //////////
  // 1) Determine deltaR(pfgsg,gsf) for all combinations  
  std::vector<GsfToGsfDR2> gsf2pfgsf_all;
  gsf2pfgsf_all.reserve( gsfTracks->size()*gsfTracksEGamma->size() );
  for ( size_t igsf = 0; igsf < gsfTracks->size(); ++igsf ) {
    reco::GsfTrackPtr gsf(gsfTracks, igsf);
    for ( size_t ipfgsf = 0; ipfgsf < gsfTracksEGamma->size(); ++ipfgsf ) {
      reco::GsfTrackPtr pfgsf(gsfTracksEGamma, ipfgsf);
      gsf2pfgsf_all.emplace_back( gsf, pfgsf, deltaR2(gsf,pfgsf) );
    }
  }

  //////////
  // 2) For each pfgsf, select best match according to DeltaR2 metric
  gsf2pfgsf.clear();
  std::sort( gsf2pfgsf_all.begin(), 
	     gsf2pfgsf_all.end(), 
	     GsfToGsfDR2::compare_by_dr2 );

  for ( size_t ipfgsf = 0; ipfgsf < gsfTracksEGamma->size(); ++ipfgsf ) {
    reco::GsfTrackPtr pfgsf(gsfTracksEGamma, ipfgsf);
    auto match_pfgsf = std::find_if( gsf2pfgsf_all.begin(),
				     gsf2pfgsf_all.end(),
				     [pfgsf](const GsfToGsfDR2& dr2) {
				       return dr2.obj2_ == pfgsf;
				     }
				     );
    if ( match_pfgsf != gsf2pfgsf.end() ) { 
      gsf2pfgsf.emplace_back(*match_pfgsf); 
      keys.erase( std::remove( keys.begin(), 
			       keys.end(), 
			       match_pfgsf->obj1_.key() ), 
		  keys.end() ); // move key to end and erase
    } else {
      std::cout << "[IDNtuplizer::trkToGsfLinks]"
		<< " ERROR! PF GsfTrackPtr is not found in the map!" << std::endl;
    }
    if ( gsf2pfgsf.size() >= gsfTracksEGamma->size() ) { break; }  // found unique match for all PF GsfTracks
  }

  //////////
  // 3) Add (null) entries for all low-pT GsfTracks without a match to a PF GsfTrack
  for ( auto key : keys ) {
    reco::GsfTrackPtr gsf(gsfTracks,key);
    // Check if trk already stored in map and, if not, add "empty" entry
    auto match_gsf = std::find_if( gsf2pfgsf.begin(), 
				   gsf2pfgsf.end(), 
				   [gsf](const GsfToGsfDR2& dr2) { 
				     return dr2.obj1_ == gsf; 
				   }
				   );
    if ( match_gsf == gsf2pfgsf.end() ) {
      reco::GsfTrackPtr pfgsf;
      gsf2pfgsf.emplace_back( gsf, pfgsf, IDNtuple::NEG_FLOAT );
    }
  }    

}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
				  edm::Handle< edm::View<T> >& candidates,
				  std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
				  std::vector< DeltaR2<reco::Candidate,T> >& other_cand,
				  bool append ) {
  std::vector< edm::Ptr<T> > cands;
  for ( size_t icand = 0; icand < candidates->size(); ++icand ) {
    edm::Ptr<T> cand(candidates,icand);
    if ( validPtr(cand) ) { cands.push_back(cand); }
    else {
      std::cout << "[IDNtuplizer::sigToCandLinks] ERROR! CandidatePtr:"
		<< " cand.isNull(): " << cand.isNull()
		<< " cand.isAvailable(): " << cand.isAvailable()
		<< std::endl;
    }
  }
  sigToCandLinks<T>( signal_electrons, cands, sig2cand, other_cand, append );
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
				  edm::Handle< std::vector<T> >& candidates,
				  std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
				  std::vector< DeltaR2<reco::Candidate,T> >& other_cand,
				  bool append ) {
  std::vector< edm::Ptr<T> > cands;
  for ( size_t icand = 0; icand < candidates->size(); ++icand ) {
    edm::Ptr<T> cand(candidates,icand);
    if ( validPtr(cand) ) { cands.push_back(cand); }
    else {
      std::cout << "[IDNtuplizer::sigToCandLinks] ERROR! CandidatePtr:"
		<< " cand.isNull(): " << cand.isNull()
		<< " cand.isAvailable(): " << cand.isAvailable()
		<< std::endl;
    }
  }
  sigToCandLinks<T>( signal_electrons, cands, sig2cand, other_cand, append );
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
				  std::vector< edm::Ptr<T> >& candidates,
				  std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
				  std::vector< DeltaR2<reco::Candidate,T> >& other_cand,
				  bool append ) {
  
  if ( !append ) { sig2cand.clear(); }
  if ( !append ) { other_cand.clear(); }
  
  // DeltaR2 matches for all combinations of signal electrons and reco::Candidate
  std::vector< DeltaR2<reco::Candidate,T> > sig2cand_all;
  sig2cand_all.reserve( sig2cand_all.size() + signal_electrons.size()*candidates.size() );
  for ( auto cand : candidates ) {
    if ( !validPtr(cand) ) { continue; }
    if ( !filterCand<T>(cand) ) { continue; }
    for ( auto sig : signal_electrons ) {
      sig2cand_all.emplace_back( sig,
				 cand,
				 deltaR2<reco::Candidate,T>(sig,cand) );
    }
    // Note: matched candidates are removed below
    if ( prescale_ < 0. || ( gRandom->Rndm() < prescale_ ) ) { 
      reco::CandidatePtr null;
      other_cand.emplace_back( null, cand, IDNtuple::NEG_FLOAT );
    }
  }
  
  // Sort by DeltaR2!!
  std::sort( sig2cand_all.begin(), 
	     sig2cand_all.end(), 
	     DeltaR2<reco::Candidate,T>::compare_by_dr2 );
  
  // Select best matches according to DeltaR2 metric
  sig2cand.clear();
  for ( auto iter : sig2cand_all ) {
    auto sig = std::find_if( sig2cand.begin(), 
			     sig2cand.end(), 
			     [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
			       return dr2.obj1_ == iter.obj1_; 
			     }
			     );
    auto cand = std::find_if( sig2cand.begin(), 
			      sig2cand.end(), 
			      [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
				return dr2.obj2_ == iter.obj2_; 
			      }
			      );
    if ( sig == sig2cand.end() && cand == sig2cand.end() ) {
      if ( iter.dr2_ < dr_max_*dr_max_ ) { sig2cand.push_back(iter); }
      else { 
	sig2cand.emplace_back( iter.obj1_, edm::Ptr<T>()/*null*/, IDNtuple::NEG_FLOAT );
      }
    }
    if ( sig2cand.size() >= signal_electrons.size() ) { break; }
  }
  
  // Remove matched candidates
  for ( auto iter : sig2cand ) {
    auto match = std::find_if( other_cand.begin(), 
			       other_cand.end(), 
			       [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
				 return dr2.obj2_ == iter.obj2_; 
			       }
			       );
    if ( match != other_cand.end() ) { other_cand.erase(match); }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::match( reco::CandidatePtr& sig,
			 std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
			 edm::Ptr<T>& cand_ptr, 
			 float& cand_dr, 
			 bool& cand_match ) { 
  
  // Identify GEN electron matches to candidate
  auto match = std::find_if( sig2cand.begin(), 
			     sig2cand.end(), 
			     [sig](const DeltaR2<reco::Candidate,T>& dr2) { 
			       return dr2.obj1_ == sig; 
			     }
			     );
  if ( match != sig2cand.end() ) {
    if ( validPtr(match->obj2_) ) {
      cand_ptr = match->obj2_; // pass by ref
      if ( deltaR2(sig,cand_ptr) >= 0. ) cand_dr = sqrt(deltaR2(sig,cand_ptr)); // pass by ref
      cand_match = ( cand_dr >= 0. ) && ( cand_dr < dr_threshold_ ); // pass by ref
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T1, typename T2>
float IDNtuplizer::deltaR2( edm::Ptr<T1>& cand1, edm::Ptr<T2>& cand2 ) {
  return reco::deltaR2( cand1->eta(), 
			cand1->phi(),
			cand2->eta(),
			cand2->phi() ); 
}

float IDNtuplizer::deltaR2( reco::TrackPtr& trk,
			    reco::GsfTrackPtr& gsf ) {
  return reco::deltaR2( trk->eta(), 
			trk->phi(),
			gsf->etaMode(),
			gsf->phiMode() ); 
}

float IDNtuplizer::deltaR2( reco::GsfTrackPtr& pfgsf,
			    reco::GsfTrackPtr& gsf ) {
  return reco::deltaR2( pfgsf->etaMode(), 
			pfgsf->phiMode(),
			gsf->etaMode(),
			gsf->phiMode() ); 
}

float IDNtuplizer::deltaR2( reco::CandidatePtr& sig,
			    reco::GsfTrackPtr& gsf ) {
  return reco::deltaR2( sig->eta(), 
			sig->phi(),
			gsf->etaMode(),
			gsf->phiMode() ); 
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
bool IDNtuplizer::filterCand( edm::Ptr<T>& cand ) { 
  if ( cand->pt() < minTrackPt_ ) { return false; }
  //if ( std::abs(cand->eta()) > 2.5 ) { return false; } //@@ Needed for PFCandidates?
  return true;
}

bool IDNtuplizer::filterCand( edm::Ptr<reco::Track>& trk ) { 
  if ( trk->pt() < minTrackPt_ ) { return false; }
  //if ( std::abs(trk->eta()) > 2.5 ) { return false; } //@@ Needed for PFCandidates?
  if ( !(trk->quality(reco::TrackBase::qualityByName("highPurity"))) ) { return false; }
  return true; 
}

bool IDNtuplizer::filterCand( edm::Ptr<reco::GsfTrack>& gsf ) { 
  if ( gsf->ptMode() < minTrackPt_ ) { return false; }
  //if ( std::abs(gsf->etaMode()) > 2.5 ) { return false; } //@@ Needed for PFCandidates?
  return true; 
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T> 
bool IDNtuplizer::validPtr( edm::Ptr<T>& ptr ) {
  return ( ptr.isNonnull() && ptr.isAvailable() );
}

////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDNtuplizer);




















//  // Lower-level method that creates ElectronChain objects for low-pT electrons
//  void signal( bool is_egamma,
//	       std::set<reco::CandidatePtr>& signal_electrons,
//	       std::vector<SigToTrkDR2>& sig2trk,
//	       std::vector<SigToGsfDR2>& sig2gsf,
//	       std::vector<SigToEleDR2>& sig2ele,
//	       std::vector<TrkToGsfDR2>& trk2gsf,
//	       std::vector<TrkToEleDR2>& trk2ele,
//	       std::vector<GsfToEleDR2>& gsf2ele);
//  
//  // Lower-level method that creates ElectronChain objects for low-pT electrons
//  void fakes( bool is_egamma,
//	      std::vector<SigToTrkDR2>& other_trk,
//	      std::vector<SigToGsfDR2>& other_gsf,
//	      std::vector<TrkToGsfDR2>& trk2gsf,
//	      std::vector<TrkToEleDR2>& trk2ele,
//	      std::vector<GsfToEleDR2>& gsf2ele );

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////
//void IDNtuplizer::signal( bool is_egamma,
//			  std::set<reco::CandidatePtr>& signal_electrons,
//			  std::vector<SigToTrkDR2>& sig2trk,
//			  std::vector<SigToGsfDR2>& sig2gsf,
//			  std::vector<SigToEleDR2>& sig2ele,
//			  std::vector<TrkToGsfDR2>& trk2gsf,
//			  std::vector<TrkToEleDR2>& trk2ele,
//			  std::vector<GsfToEleDR2>& gsf2ele ) {
//  
//  // Iterate through "signal electrons" and initialise ElectronChain objects
//  for ( auto sig : signal_electrons ) {
//    chains_.push_back(ElectronChain());
//    ElectronChain& chain = chains_.back();
//    chain.is_mc_ = isMC_;
//    chain.is_aod_ = isAOD_;
//    chain.is_e_ = true;
//    chain.is_egamma_ = is_egamma;
//    chain.sig_ = sig;
//
//    // Find match between "signal electron" and reconstructed objects (info passed by ref)
//    match<reco::GsfElectron>(sig,sig2ele,chain.ele_,chain.ele_dr_,chain.ele_match_);
//    match<reco::GsfTrack>(sig,sig2gsf,chain.gsf_,chain.gsf_dr_,chain.gsf_match_);
//    match<reco::Track>(sig,sig2trk,chain.trk_,chain.trk_dr_,chain.trk_match_);
//    
//    // Update ElectronChain if matched to GsfTrack
//    if ( chain.gsf_match_ ) {
//
//      // Update Seed BDTs
//      if ( !chain.is_egamma_ ) { 
//	chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
//	chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
//      }
//      
//      // Update (override?) GsfElectron info
//      auto match_ele = std::find_if( gsf2ele.begin(), 
//				     gsf2ele.end(), 
//				     [chain](const GsfToEleDR2& dr2) { 
//				       return dr2.obj1_ == chain.gsf_; // find gsf
//				     }
//				     );
//      if ( match_ele != gsf2ele.end() && validPtr(match_ele->obj2_) ) { 
//	chain.ele_ = match_ele->obj2_; 
//	chain.ele_match_ = true;;
//	chain.ele_dr_ = sqrt(deltaR2(chain.sig_,chain.ele_));
//      } else {
//	chain.ele_ = reco::GsfElectronPtr();
//	chain.ele_match_ = false;
//	chain.ele_dr_ = IDNtuple::NEG_FLOAT;
//      }
//
//      // Update Track info
//      reco::TrackPtr trk;
//      if ( gsfToTrk(chain.gsf_,trk,is_egamma) ) {  // Use "link maps" for low-pT if no TrackExtra (in mAOD)
//	//@@ what about PF and mAOD? how to navigate back to trk without eleSeed?
//	//@@ track pT? 
//	//@@ embedded track for mAOD???? 
//	chain.trk_ = trk; 
//	chain.trk_match_ = true;
//	chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
//	PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
//	if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
//      } else {
//	chain.trk_ = reco::TrackPtr(); 
//	chain.trk_match_ = false;
//	chain.trk_dr_ = IDNtuple::NEG_FLOAT;
//	chain.pdg_id_ = 0;
//      }
//
//      // Update ElectronSeed info
//      reco::ElectronSeedPtr seed;
//      if ( gsfToSeed(chain.gsf_,seed) ) {
//	chain.seed_ = seed;
//	chain.seed_tracker_driven_ = seed->isTrackerDriven();
//	chain.seed_ecal_driven_ = seed->isEcalDriven();
//	
////	// Update PreId info
////	if ( !chain.is_egamma_ ) { // && isAOD_ == 1
////	  chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
////	  chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
////	}
//
//	// Update CaloCluster info
//	reco::CaloClusterPtr calo;
//	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; }
//
//	//@@ ONLY ECAL DRIVEN????? WHAT ABOUT PF + MINIAOD?
//	// If ECAL-driven _only_, update "surrogate" Track info
//	if ( !seed->isTrackerDriven() && seed->isEcalDriven() ) {
//	  auto match_trk = std::find_if( trk2gsf.begin(), 
//					 trk2gsf.end(), 
//					 [chain](const TrkToGsfDR2& dr2) { 
//					   return dr2.obj2_ == chain.gsf_; // find gsf
//					 }
//					 );
//	  if ( match_trk != trk2gsf.end() && validPtr(match_trk->obj1_) ) { 
//	    chain.trk_ = match_trk->obj1_; 
//	    chain.trk_match_ = true;
//	    chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
//	    PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
//	    if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
//	  } else {
//	    chain.trk_ = reco::TrackPtr(); 
//	    chain.trk_match_ = false;
//	    chain.trk_dr_ = IDNtuple::NEG_FLOAT;
//	    chain.pdg_id_ = 0;
//	  }
//	}
//      }
//
//    } // if ( chain.gsf_match_ )
//
//  } // for ( auto sig : signal_electrons )
//    
//}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////
//void IDNtuplizer::fakes( bool is_egamma,
//			 std::vector<SigToTrkDR2>& other_trk,
//			 std::vector<SigToGsfDR2>& other_gsf,
//			 std::vector<TrkToGsfDR2>& trk2gsf,
//			 std::vector<TrkToEleDR2>& trk2ele,
//			 std::vector<GsfToEleDR2>& gsf2ele ) {
//  
//  // Iterate through generalTracks
//  for ( auto iter : other_trk ) {
//
//    if ( !validPtr(iter.obj2_) ) { continue; }
//
//    // Create ElectronChain object and start to populate
//    chains_.push_back(ElectronChain());
//    ElectronChain& chain = chains_.back();
//    chain.is_mc_ = isMC_;
//    chain.is_aod_ = isAOD_;
//    chain.is_e_ = false;
//    chain.is_egamma_ = is_egamma;
//    
//    // Store Track info
//    chain.trk_ = iter.obj2_;
//    chain.trk_match_ = true;
//    chain.trk_dr_ = IDNtuple::NEG_FLOAT;
//    
//    // Find matched GsfTrack match, either via ElectronSeed or just deltaR
//    auto match_gsf = std::find_if( trk2gsf.begin(), 
//				   trk2gsf.end(), 
//				   [chain](const TrkToGsfDR2& dr2) { 
//				     return dr2.obj1_ == chain.trk_; // find trk
//				   }
//				   );
//    if ( match_gsf != trk2gsf.end() && validPtr(match_gsf->obj2_) ) {
//      chain.gsf_ = match_gsf->obj2_;
//      chain.gsf_match_ = true;
//      chain.gsf_dr_ = match_gsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_gsf->dr2_);
//    } else {
//      chain.gsf_ = reco::GsfTrackPtr();
//      chain.gsf_match_ = false;
//      chain.gsf_dr_ = IDNtuple::NEG_FLOAT;
//      if ( verbose_ > 3 ) {
//	std::cout << "INFO! Cannot match TrackPtr to GsfTrackPtr:"
//		  << " TrackPtr: " << chain.trk_.id() << "/" << chain.trk_.key()
//		  << std::endl;
//      }
//    }
//
//    // Check GsfTrackPtr, then update ElectronChain info
//    if ( !validPtr(chain.gsf_) ) { continue; } 
//
//    // Update Seed BDTs
//    if ( !chain.is_egamma_ ) { 
//      chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
//      chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
//    }
//    
//    // Update (override?) GsfElectron info (should be identical to that above?!)
//    auto match_ele = std::find_if( gsf2ele.begin(), 
//				   gsf2ele.end(), 
//				   [chain](const GsfToEleDR2& dr2) { 
//				     return dr2.obj1_ == chain.gsf_; // find gsf
//				   }
//				   );
//    if ( match_ele != gsf2ele.end() && validPtr(match_ele->obj2_) ) {
//      if ( validPtr(chain.ele_) && chain.ele_ != match_ele->obj2_ ) {
//	std::cout << "ERROR! GsfElectrons are not consistent!"
//		  << " matched ele id/key: " << match_ele->obj2_.id() << "/" << match_ele->obj2_.key()
//		  << " chain.ele id/key: " << chain.ele_.id() << "/" << chain.ele_.key()
//		  << std::endl;
//      }
//      chain.ele_ = match_ele->obj2_;
//      chain.ele_match_ = true;
//      chain.ele_dr_ = match_ele->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_ele->dr2_);
//     } else {
//      chain.ele_ = reco::GsfElectronPtr();
//      chain.ele_match_ = false;
//      chain.ele_dr_ = IDNtuple::NEG_FLOAT;
//    }
//
//    // Check GsfElectron and GsfTrack are consistent
//    reco::GsfTrackPtr gsf;
//    if ( validPtr(chain.ele_) && 
//	 eleToGsf(chain.ele_,gsf) && 
//	 gsf != chain.gsf_ ) { 
//      std::cout << "ERROR! GsfElectron and GsfTrack are not consistent!"
//		<< " GsfElectron id/key: " << chain.ele_.id() << "/" << chain.ele_.key()
//		<< " GsfTrack id/key: " << chain.gsf_.id() << "/" << chain.gsf_.key()
//		<< " matched gsf id/key: " << gsf.id() << "/" << gsf.key() 
//		<< std::endl;
//    }
//    
//    // Add ElectronSeed, CaloCluster, check Track
//    if ( validPtr(chain.gsf_) ) { 
//      reco::ElectronSeedPtr seed;
//      if ( gsfToSeed(chain.gsf_,seed) ) {
//	chain.seed_ = seed;
//	chain.seed_tracker_driven_ = seed->isTrackerDriven();
//	chain.seed_ecal_driven_ = seed->isEcalDriven();
//      }
//      
//      reco::CaloClusterPtr calo;
//      if ( seedToCalo(seed,calo) ) { 
//	chain.calo_ = calo; 
//      }
//      if ( verbose_ > 3 ) {
//	reco::TrackPtr trk;
//	if ( seedToTrk(seed,trk) && ( trk != chain.trk_ ) && chain.gsf_dr_ < 0. ) { 
//	  std::cout << "WARNING! Track is not consistent!"
//		    << " Track id/key: " << chain.trk_.id() << "/" << chain.trk_.key()
//		    << " matched trk id/key: " << trk.id() << "/" << trk.key() 
//		    << ". Likely due to multiple GsfTracks originating from a"
//		    << " common ElectronSeed, coupled with the use of gsfToSeed()."
//		    << std::endl; //@@ gsfToTrk? Missing TrackExtra?
//	}
//      }
//      
//    }
//    
//  } // for ( auto trk : other_trk )
//  
//}











//  // Extract GsfTrackRef for use with gsf->track maps 
//  reco::GsfTrackRef gsf_ref;
//  if ( !is_egamma ) { 
//    gsf_ref = reco::GsfTrackRef(gsfTracksH_,(unsigned long)gsf.key());
//  } else {
//    pfgsf_ref = reco::GsfTrackRef(gsfTracksEGammaH_,(unsigned long)gsf.key());
//    gsf_ref = gsfToPfGsfLinks(pfgsf_ref); //@@ Hack: PF GSF track --> low-pT GSF track via deltaR matching
//  }
//  if ( !gsf_ref.isValid() || !gsf_ref.isAvailable() ) { return false; }
//
//  // If navigation via ElectronSeed fails (e.g. TrackExtra missing), attempt to use:
//  if ( isAOD_ == 1 ) {
//    // 1) low-pT "gsf->trk" map in AOD
//    reco::TrackRef trk_ref = (*gsfTrackLinksH_)[gsf_ref];
//    trk = edm::refToPtr(trk_ref);
//    if ( validPtr(trk) ) { return true; }
//  } else if ( isAOD_ == 0 ) {
//    // 2a) low-pT "gsf->packed" map in MINIAOD
//    pat::PackedCandidateRef packed_ref = (*packedCandLinksH_)[gsf_ref];
//    if ( packed_ref.isAvailable() && 
//	 packed_ref.isNonnull() && 
//	 packed_ref->bestTrack() != nullptr ) { 
//      trk = reco::TrackPtr(packed_ref->bestTrack(),packed_ref.key());
//      if ( validPtr(trk) ) { return true; }
//    }
//    // 2b) low-pT "gsf->lost" map in MINIAOD
//    pat::PackedCandidateRef lost_ref = (*lostTrackLinksH_)[gsf_ref];
//    if ( lost_ref.isAvailable() && 
//	 lost_ref.isNonnull() && 
//	 lost_ref->bestTrack() != nullptr ) { 
//      trk = reco::TrackPtr(lost_ref->bestTrack(),lost_ref.key());
//      if ( validPtr(trk) ) { return true; }
//    }
//  }

















//  void trkToGsfLinks( std::vector<reco::TrackPtr>& ctfTracks,
//		      edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
//		      std::vector<TrkToGsfDR2>& trk2gsf,
//		      bool is_egamma,
//		      std::vector<GsfToGsfDR2>& gsf2pfgsf );

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////
//void IDNtuplizer::trkToGsfLinks( std::vector<reco::TrackPtr>& ctfTracks,
//				 edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
//				 std::vector<TrkToGsfDR2>& trk2gsf,
//				 bool is_egamma,
//				 std::vector<GsfToGsfDR2>& gsf2pfgsf ) {
//  
//  // If not, use deltaR match between GsfTrack and "surrogate" Track
//
////  std::vector<TrkToGsfDR2> trk2gsf_all;
////  if ( ctfTracks.empty() || gsfTracks->empty() ) { return; }
////  trk2gsf_all.reserve( ctfTracks.size()*gsfTracks->size() );
////  for ( auto& trk : ctfTracks ) {
////    for ( size_t igsf = 0; igsf < gsfTracks->size(); ++igsf ) {
////      reco::GsfTrackPtr gsf(gsfTracks, igsf);
////      reco::TrackPtr ptr;
////      if ( gsfToTrk(gsf,ptr,is_egamma) && ptr == trk ) { 
////	trk2gsf_all.emplace_back( trk, gsf, IDNtuple::NEG_FLOAT ); // match via ElectronSeed
////      } else { 
////	trk2gsf_all.emplace_back( trk, gsf, deltaR2(trk,gsf) ); // match with a surrogate track
////      }
////    }
////  }
//  
//  std::vector<TrkToGsfDR2> trk2gsf_all;
//  if ( ctfTracks.empty() || gsfTracks->empty() ) { return; }
//  trk2gsf_all.reserve( ctfTracks.size()*gsfTracks->size() );
//  for ( auto& trk : ctfTracks ) {
//    for ( size_t igsf = 0; igsf < gsfTracks->size(); ++igsf ) {
//      reco::GsfTrackPtr gsf(gsfTracks, igsf);
//      if ( is_egamma ) { 
//	//@@ Hack! PF-->lowPt (via deltaR matching)
//	auto match_gsf = std::find_if( gsf2pfgsf.begin(), 
//				       gsf2pfgsf.end(), 
//				       [gsf](const GsfToGsfDR2& dr2) { 
//					 return dr2.obj1_ == gsf; // find pfgsf
//				       }
//				       );
//	if ( match_gsf == gsf2pfgsf.end() || !validPtr(match_gsf->obj2_) ) { continue; }
//	gsf = match_gsf->obj2_; 
//      }
//      reco::TrackPtr ptr;
//      if ( gsfToTrk(gsf,ptr,false/*is_egamma*/) && ptr == trk ) { //@@ Hack! Force is_gamma==false
//	trk2gsf_all.emplace_back( trk, gsf, IDNtuple::NEG_FLOAT ); // match via ElectronSeed
//      } else { 
//	//trk2gsf_all.emplace_back( trk, gsf, deltaR2(trk,gsf) ); // match with a surrogate track
//	trk2gsf_all.emplace_back( trk, gsf, gRandom->Rndm() ); // randomize surrogate order in sorted vector via dr2
//      }
//    }
//  }
//
//  // Sort by ascending deltaR2 (Note, matches via ElectronSeed have -ve deltaR2)
//  std::sort( trk2gsf_all.begin(), 
//	     trk2gsf_all.end(), 
//	     TrkToGsfDR2::compare_by_dr2 );
//  
//  // Select best matches according to DeltaR2 metric (i.e. keeps matches via ElectronSeed)
//  trk2gsf.clear();
//  for ( auto iter : trk2gsf_all ) {
//    auto trk = std::find_if( trk2gsf.begin(), 
//			     trk2gsf.end(), 
//			     [iter](const TrkToGsfDR2& dr2) { 
//			       return dr2.obj1_ == iter.obj1_; 
//			     }
//			     );
//    auto gsf = std::find_if( trk2gsf.begin(), 
//			     trk2gsf.end(), 
//			     [iter](const TrkToGsfDR2& dr2) { 
//			       return dr2.obj2_ == iter.obj2_; 
//			     }
//			     );
//    if ( trk == trk2gsf.end() && // Check trk not already used
//	 gsf == trk2gsf.end() ) { // Check gsf not already used
//      if ( iter.dr2_ < 0. ) {
//	// For matches via ElectronSeed, update to correct deltaR2
//	trk2gsf.emplace_back(iter.obj1_,
//			     iter.obj2_,
//			     deltaR2(iter.obj1_,iter.obj2_)); 
////      } else if ( iter.dr2_ == 0. ) {
////	// For matches via deltaR, update to correct deltaR2
////	trk2gsf.emplace_back(iter.obj1_,
////			     iter.obj2_,
////			     deltaR2(iter.obj1_,iter.obj2_)); 
//      } else {
//	// For matches using "surrogate" Tracks, update to IDNtuple::NEG_FLOAT
//	trk2gsf.emplace_back(iter.obj1_,
//			     iter.obj2_,
//			     IDNtuple::NEG_FLOAT); 
//      }
//    }
//    if ( trk2gsf.size() >= ctfTracks.size() ) { break; } // found unique match for all tracks
//  }
//  
//}
