#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
//#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
//#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
//#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
//#include "DataFormats/EgammaReco/interface/SuperCluster.h"
//#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
//#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
//#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
//#include "DataFormats/Math/interface/Vector.h"
//#include "DataFormats/Math/interface/Vector3D.h"
//#include "DataFormats/Math/interface/Point3D.h"
//#include "DataFormats/ParticleFlowReco/interface/PreId.h"
//#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "LowPtElectrons/LowPtElectrons/interface/Common.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDNtupleSimple.h"
#include "TRandom3.h"
#include "TTree.h"
#include <set>
#include <vector>
#include <math.h>
#include <memory>
#include <boost/core/demangle.hpp>
#include <algorithm>
#include <random>

////////////////////////////////////////////////////////////////////////////////
//
class IDNtuplizerSimple : public edm::stream::EDFilter<> {
  
public:

  explicit IDNtuplizerSimple( const edm::ParameterSet& );
  ~IDNtuplizerSimple() override;
  
  void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  bool filter( edm::Event&, const edm::EventSetup& ) override;

//  // Wraps other methods to provide a sample of "signal" electrons
//  void signalElectrons( std::set<reco::CandidatePtr>& signal_electrons,
//			std::set<reco::CandidatePtr>& tag_side_muons );
//  
//  // GEN-based method to provide a sample of "signal" electrons
//  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B, 
//			  std::set<reco::GenParticlePtr>& gen_muons, 
//			  float tag_muon_pt_threshold,
//			  float tag_muon_eta_threshold );
//
//  void extractGsfElectronPtrs();
//  
//  // Populate ElectronChain vector
//  void createChains( std::set<reco::CandidatePtr>& signal_electrons,
//		     std::set<reco::CandidatePtr>& tag_side_muons );
//
//  // Method that creates ElectronChain objects for signal candiates
//  void signal( std::set<reco::CandidatePtr>& signal_electrons,
//	       std::set<reco::CandidatePtr>& tag_side_muons,
//	       std::vector<SigToEleDR2>& sig2ele,
//	       std::vector<SigToEleDR2>& other_ele );
//  
//  // Method that creates ElectronChain objects for bkgd candidates
//  void bkgd( std::set<reco::CandidatePtr>& signal_electrons,
//	     std::set<reco::CandidatePtr>& tag_side_muons,
//	     std::vector<SigToEleDR2>& sig2ele,
//	     std::vector<SigToEleDR2>& other_ele );
//  
//  // Fills tree per ElectronChain object
//  void fill( const edm::Event& event, const edm::EventSetup& setup );
//
//  // Links "signal" electrons to reconstructed objects
//  template <typename T> 
//  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
//		       std::vector< edm::Ptr<T> >& candidates,
//		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
//		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand, 
//		       bool append = false );
//
//  // Wraps method above to allow use of Handle<View<T>>
//  template <typename T> 
//  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
//		       edm::Handle< edm::View<T> >& candidates,
//		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
//		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand, 
//		       bool append = false );
//
//  // Return by reference the 'cand' that is matched to 'sig' in 'sig2cand' map
//  template <typename T>
//  void match( reco::CandidatePtr& sig,
//	      std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
//	      edm::Ptr<T>& cand, float& dr, bool& match ); // pass by ref
//
//  // Wrap reco::deltaR2 method for various types, specifically GsfTrack uses eta and phi 'Mode'
//  template <typename T1, typename T2> 
//  float deltaR2( edm::Ptr<T1>& cand1, edm::Ptr<T2>& cand2 ); // used by sigToCandLinks
//
//  // Filter track candidates by quality flag, simple pass for all other cands
//  template <typename T> 
//  bool filterCand( edm::Ptr<T>& cand );
//  bool filterCand( edm::Ptr<reco::GsfElectron>& ele );
//
//  // Check is Ptr is valid and available
//  template <typename T> bool validPtr( edm::Ptr<T>& ptr );
  
private:
  
  // Misc
  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDNtupleSimple ntuple_;
  int verbose_;
  bool check_from_B_;
  double dr_max_; // Max DeltaR value considered
  double dr_threshold_; // Threshold for DeltaR matching
  double prescale_;
  int isAOD_;
  bool isMC_;
  double minPt_;
  float tagMuonPtThreshold_;
  float tagMuonEtaThreshold_;

  // Generic collections

  const edm::EDGetTokenT<double> rho_;
  edm::Handle<double> rhoH_;

//  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
//  edm::Handle<reco::BeamSpot> beamspotH_;

  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectronsEGamma_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGamma_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGammaRetrained_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaRetrainedH_;

  std::vector<reco::GsfElectronPtr> ptrs_;
  std::vector<ElectronChain> chains_;

};

////////////////////////////////////////////////////////////////////////////////
//
IDNtuplizerSimple::~IDNtuplizerSimple() {
}

////////////////////////////////////////////////////////////////////////////////
//
IDNtuplizerSimple::IDNtuplizerSimple( const edm::ParameterSet& cfg ) 
  : tree_(nullptr),
    ntuple_(),
    verbose_(cfg.getParameter<int>("verbose")),
    check_from_B_(cfg.getParameter<bool>("checkFromB")),
    dr_max_(cfg.getParameter<double>("drMax")),
    dr_threshold_(cfg.getParameter<double>("drThreshold")),
    prescale_(cfg.getParameter<double>("prescale")),
    isAOD_(0), // -1 = unknown, 0 = MINIAOD, 1 = AOD (force MINIAOD for now...)
    isMC_(true), // force this true for now ...
    minPt_(cfg.getParameter<double>("minPt")),
    tagMuonPtThreshold_(cfg.getParameter<double>("tagMuonPtThreshold")),
    tagMuonEtaThreshold_(cfg.getParameter<double>("tagMuonEtaThreshold")),
    // Generic collections
    rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
    rhoH_(),
//beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
//beamspotH_(),
    prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
    genParticlesH_(),
    patElectronsEGamma_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("patElectronsEGamma"))),
    gsfElectronsEGammaH_(),
    mvaValueEGamma_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueEGamma"))),
    mvaValueEGammaH_(),
    mvaValueEGammaRetrained_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueEGammaRetrained"))),
    mvaValueEGammaRetrainedH_(),
    ptrs_(),
    chains_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "[IDNtuplizerSimple::IDNtuplizerSimple] Verbosity level: "<< verbose_ << std::endl;
  }

////////////////////////////////////////////////////////////////////////////////
// Initialise the weights LUT to filter fake tracks
void IDNtuplizerSimple::beginRun( const edm::Run& run,
				  const edm::EventSetup& es ) {
  //@@ ?
}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizerSimple::filter( edm::Event& event, 
				const edm::EventSetup& setup ) {
  
//  // Read collections from Event
//  event.getByToken(prunedGenParticles_, genParticlesH_);
//  event.getByToken(patElectronsEGamma_, gsfElectronsEGammaH_);
//  event.getByToken(mvaValueEGamma_, mvaValueEGammaH_);
//  event.getByToken(mvaValueEGammaRetrained_, mvaValueEGammaRetrainedH_);
//  
//  // Identify signal electrons
//  std::set<reco::CandidatePtr> signal_electrons;
//  std::set<reco::CandidatePtr> tag_side_muons;
//  signalElectrons(signal_electrons,tag_side_muons);
//
//  // Trigger requirement (DO WE ALWAYS WANT THIS ???)
//  if ( tag_side_muons.empty() ) { return false; }
//
//  // Populate std::vector<reco::GsfElectronPtr> ptrs_
//  extractGsfElectronPtrs();
//  
//  // Populate ElectronChain vector
//  chains_.clear();
//  createChains(signal_electrons,tag_side_muons);
//
//  // Fill ntuple
//  fill(event,setup);
//  
//  // Print ElectronChain objects
//  if ( verbose_ > 0 ) {
//    std::cout << "[IDNtuplizerSimple::filter]"
//	      << " chains_.size() = " << chains_.size()
//	      << std::endl;
//    for ( auto iter : chains_ ) { 
//      //if ( !iter.is_e_ ) { continue; } // skip fakes
//      //if ( !iter.is_egamma_ || iter.is_e_ ) { continue; } // keep PF and fakes
//      if ( !iter.is_egamma_ ) { continue; } // keep PF and fakes
//      if ( iter.ele_match_ == false ) { continue; } // keep only electron candidates
//      std::cout << iter << std::endl; 
//    }
//  }
  
  // No filter logic used here ...
  return false;
  
}

//
//////////////////////////////////////////////////////////////////////////////////
////
//void IDNtuplizerSimple::signalElectrons( std::set<reco::CandidatePtr>& signal_electrons,
//					 std::set<reco::CandidatePtr>& tag_side_muons ) {
//  
//  signal_electrons.clear();
//  tag_side_muons.clear();
//
//  if ( isMC_ ) { // Identify "signal" (GEN) electrons from B decays
//    std::set<reco::GenParticlePtr> electrons_from_B;
//    std::set<reco::GenParticlePtr> gen_muons;
//    genElectronsFromB(electrons_from_B,gen_muons,tagMuonPtThreshold_,tagMuonEtaThreshold_);
//    for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
//    for ( auto gen : gen_muons ) { tag_side_muons.insert(gen); }
//  } else { // Identify "signal" electrons from data control regions
//    //@@ FOR NOW, A HACK ...
//    std::set<reco::GenParticlePtr> electrons_from_B;
//    std::set<reco::GenParticlePtr> gen_muons;
//    genElectronsFromB(electrons_from_B,gen_muons,tagMuonPtThreshold_,tagMuonEtaThreshold_);
//    for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
//    for ( auto gen : gen_muons ) { tag_side_muons.insert(gen); }
//  }
//  
//}
//
//////////////////////////////////////////////////////////////////////////////////
//// 
//void IDNtuplizerSimple::genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,
//					   std::set<reco::GenParticlePtr>& gen_muons,
//					   float tag_muon_pt_threshold, 
//					   float tag_muon_eta_threshold ) {
//  
//  electrons_from_B.clear();
//  gen_muons.clear();
//  
//  for ( size_t idx = 0; idx < genParticlesH_->size(); idx++ ) {
//    
//    reco::GenParticlePtr gen(genParticlesH_, idx);
//    if ( !validPtr(gen) ) {
//      std::cout << "[IDNtuplizerSimple::genElectronsFromB] ERROR! reco::GenParticlePtr:"
//		<< " gen.isNull(): " << gen.isNull()
//		<< " gen.isAvailable(): " << gen.isAvailable()
//		<< std::endl;
//      continue;
//    }
//    
//    // Last copy of GEN electron 
//    bool is_ele = std::abs(gen->pdgId()) == 11 && gen->isLastCopy(); //@@ not a method of Candidate
//    
//    // Does GEN ele comes from B decay?
//    bool non_resonant = gen->numberOfMothers() >= 1 && gen->mother() &&   // has mother
//      std::abs(gen->mother()->pdgId()) > 510 &&                           // mother is B
//      std::abs(gen->mother()->pdgId()) < 546;                             // mother is B
//    bool resonant = gen->numberOfMothers() >= 1 && gen->mother() &&       // has mother
//      std::abs(gen->mother()->pdgId()) == 443 &&                          // mother is J/psi
//      gen->mother()->numberOfMothers() >= 1 && gen->mother()->mother() && // has grandmother
//      std::abs(gen->mother()->mother()->pdgId()) > 510 &&                 // grandmother is B
//      std::abs(gen->mother()->mother()->pdgId()) < 546;                   // grandmother is B
//    
//    //  Check for tag side muon
//    bool is_muon = std::abs(gen->pdgId()) == 13 && gen->isLastCopy() && 
//      gen->pt() > tag_muon_pt_threshold && 
//      std::abs(gen->eta()) < tag_muon_eta_threshold;
//    
//    // Does GEN muon comes from B decay?
//    bool non_res_to_muons = gen->numberOfMothers() >= 1 && gen->mother() && // has mother
//      std::abs(gen->mother()->pdgId()) > 510 &&                             // mother is B
//      std::abs(gen->mother()->pdgId()) < 546;                               // mother is B
//    bool res_to_muons = gen->numberOfMothers() >= 1 && gen->mother() &&     // has mother
//      std::abs(gen->mother()->pdgId()) == 443 &&                            // mother is J/psi
//      gen->mother()->numberOfMothers() >= 1 && gen->mother()->mother() &&   // has grandmother
//      std::abs(gen->mother()->mother()->pdgId()) > 510 &&                   // grandmother is B
//      std::abs(gen->mother()->mother()->pdgId()) < 546;                     // grandmother is B
//    
//    bool tag_muon = is_muon && ( non_res_to_muons || res_to_muons );
//    if ( tag_muon ) { gen_muons.insert(gen); }
//    
//    // is coming from a B
//    if ( is_ele && ( ( resonant || non_resonant ) || !check_from_B_ ) ) {
//      electrons_from_B.insert(gen);
//      if ( verbose_ > 1 ) {
//	std::cout << "[IDNtuplizerSimple::genElectronsFromB] "
//		  << " #signal_electrons: " << electrons_from_B.size()
//		  << " resonant? " << resonant
//		  << " non resonant? " << non_resonant
//		  << " tag-side muon? " << !gen_muons.empty()
//		  << std::endl;
//      }
//    }
//    
//  } // genParticles loop
//
//  if ( gen_muons.empty() ) { electrons_from_B.clear(); }
//
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//void IDNtuplizerSimple::extractGsfElectronPtrs() {
//  ptrs_.clear();
//  for ( size_t iele = 0; iele < gsfElectronsEGammaH_->size(); ++iele ) { 
//    reco::GsfElectronPtr ptr(gsfElectronsEGammaH_,iele);
//    if ( !filterCand<reco::GsfElectron>(ptr) ) { continue; }
//    ptrs_.push_back(ptr); 
//  }
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//void IDNtuplizerSimple::createChains( std::set<reco::CandidatePtr>& signal_electrons,
//				      std::set<reco::CandidatePtr>& tag_side_muons ) {
//  
//  // Match "signal" electrons to reco::GsfElectrons
//  std::vector<SigToEleDR2> sig2ele;
//  std::vector<SigToEleDR2> other_ele;
//  sigToCandLinks<reco::GsfElectron>( signal_electrons, 
//				     ptrs_, 
//				     sig2ele, 
//				     other_ele );
//  if ( verbose_ > 1 ) {
//    std::cout << "[IDNtuplizerSimple::createChains] sigToCandLinks<reco::GsfElectron>:" << std::endl
//	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
//	      << " gsfElectronsEGammaH_->size(): " << gsfElectronsEGammaH_->size() << std::endl
//	      << " sig2ele.size(): " << sig2ele.size() << std::endl
//	      << " other_ele.size(): " << other_ele.size() << std::endl;
//    if ( verbose_ > 2 ) {
//      for ( auto iter : sig2ele ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
//    }
//    std::cout << std::endl;
//  }
//  
//  signal( signal_electrons, tag_side_muons, sig2ele, other_ele );
//  bkgd( signal_electrons, tag_side_muons, sig2ele, other_ele );
//  
//}
//
//////////////////////////////////////////////////////////////////////////////////
//// Iterate through "signal electrons" 
//void IDNtuplizerSimple::signal( std::set<reco::CandidatePtr>& signal_electrons,
//				std::set<reco::CandidatePtr>& tag_side_muons,
//				std::vector<SigToEleDR2>& sig2ele,
//				std::vector<SigToEleDR2>& other_ele ) {
//  
//  // Find highest-pT tag-side muon in event and store
//  float tag_pt = id::NEG_FLOAT;
//  float tag_eta = id::NEG_FLOAT;
//  for ( auto tag : tag_side_muons ) {
//    if ( tag->pt() > tag_pt ) { // Find highest-pT tag muon 
//      tag_pt = tag->pt();
//      tag_eta = tag->eta();
//    }
//  }
//  
//  // Iterate through signal electrons
//  for ( auto sig : signal_electrons ) {
//    
//    // Repeat for both low pT and EGamma
//    for ( auto is_egamma : std::vector<bool>{ true } ) { //@@ was {true, false} !!!!!!!!!!!!!!!!!
//      
//      // SIG: Initialise ElectronChain object
//      chains_.push_back(ElectronChain());
//      ElectronChain& chain = chains_.back();
//      chain.is_mc_ = isMC_;
//      chain.is_aod_ = isAOD_;
//      chain.is_e_ = true;
//      chain.is_egamma_ = is_egamma; // Set here if low pT or EGamma!
//      chain.sig_ = sig;
//      chain.tag_pt_ = tag_pt;
//      chain.tag_eta_ = tag_eta;
//
//      // ELE: Store matches between "signal electron" and reco'ed electrons
//      match<reco::GsfElectron>(sig,
//			       sig2ele,
//			       chain.ele_,
//			       chain.ele_dr_,
//			       chain.ele_match_ );
//      
//    }
//  }
//
//}
//
//////////////////////////////////////////////////////////////////////////////////
//// Iterate through fake candidates
//void IDNtuplizerSimple::bkgd( std::set<reco::CandidatePtr>& signal_electrons,
//			      std::set<reco::CandidatePtr>& tag_side_muons,
//			      std::vector<SigToEleDR2>& sig2ele,
//			      std::vector<SigToEleDR2>& other_ele ) {
//  
//  // Find highest-pT tag-side muon in event and store
//  float tag_pt = id::NEG_FLOAT;
//  float tag_eta = id::NEG_FLOAT;
//  for ( auto tag : tag_side_muons ) {
//    if ( tag->pt() > tag_pt ) { // Find highest-pT tag muon 
//      tag_pt = tag->pt();
//      tag_eta = tag->eta();
//    }
//  }
//  
//  // Iterate through tracks
//  for ( auto iter : other_ele ) {
//    
//    // Repeat for two options: low pT and PF EGamma reconstruction
//    for ( auto is_egamma : std::vector<bool>{ true } ) { //@@ was {true, false} !!!!!!!!!!!!!!!!!
//
//      // SIG: Initialise ElectronChain object
//      chains_.push_back(ElectronChain());
//      ElectronChain& chain = chains_.back();
//      chain.is_mc_ = isMC_;
//      chain.is_aod_ = isAOD_;
//      chain.is_e_ = false;
//      chain.is_egamma_ = is_egamma; // Set here if low pT or EGamma!
//      chain.tag_pt_ = tag_pt;
//      chain.tag_eta_ = tag_eta;
//      
//      // ELE: Store reco'ed electron info
//      chain.ele_ = iter.obj2_; 
//      chain.ele_match_ = true;
//      chain.ele_dr_ = 0.;
//      
//    }
//  }
//
//}
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////
//void IDNtuplizerSimple::fill( const edm::Event& event,
//			      const edm::EventSetup& setup ) {
//  
//  using value_type = unsigned int;
//
//  // Random seed generator
//  static value_type seed = 1; //std::chrono::system_clock::now().time_since_epoch().count(); //@@ fixed seed !!
//  static std::default_random_engine generator(seed);
//
//  // Randomise the order of the EletronChain objects
//  std::shuffle( chains_.begin(), chains_.end(), generator );
//
//   // Poisson with mean = prescale_ (used when prescale_ is -ve)
//  std::poisson_distribution<value_type> poisson(-1.*prescale_);
//
//  std::vector<value_type> cands;
//  std::vector<value_type> cands2;
//  cands.reserve(chains_.size()); 
//  cands2.reserve(chains_.size()); 
//
//  if ( prescale_ < -1.e-9 ) { 
//
//    // If prescale_ is -ve, select cands based on Poisson distr (prescale_ is mean)
//    // 1) Store indices of fakes (separately for low pT and EGamma)
//    for ( unsigned int idx = 0; idx < chains_.size(); ++idx ) { 
//      if ( chains_[idx].is_e_ == false ) { 
//	if ( chains_[idx].is_egamma_ == false ) { cands.push_back(idx); }
//	if ( chains_[idx].is_egamma_ == true ) { cands2.push_back(idx); }
//      }
//    }
//    // 2) Shuffle then truncate to give 'prescale' elements
//    unsigned int size = poisson(generator);
//    cands.resize(size);
//    cands2.resize(size);
//
//  } else if ( prescale_ > 1.e-9 ) { 
//
//    // If prescale_ is +ve, select cands based on flat probability prior (1./prescale_ is threshold)
//    for ( unsigned int idx = 0; idx < chains_.size(); ++idx ) {
//      if ( gRandom->Rndm() < (1./prescale_) ) {
//	if ( chains_[idx].is_e_ == false ) { 
//	  if ( chains_[idx].is_egamma_ == false ) { cands.push_back(idx); }
//	  if ( chains_[idx].is_egamma_ == true ) { cands2.push_back(idx); }
//	}
//      }
//    }
//    
//  } else {
//
//    // If prescale is 0, select all cands
//    for ( unsigned int idx = 0; idx < chains_.size(); ++idx ) {
//      if ( chains_[idx].is_e_ == false ) {
//	if ( chains_[idx].is_egamma_ == false ) { cands.push_back(idx); }
//	if ( chains_[idx].is_egamma_ == true ) { cands2.push_back(idx); }
//      }
//    }
//    
//  }
//
//  // Fill ntuple
//  for ( size_t idx = 0; idx < chains_.size(); ++idx ) {
//    ElectronChain chain = chains_[idx];
//
////    std::cout << "TEST"
////	      << " idx " << idx
////	      <<  " size " << chains_.size()
////	      << " is_e_ " << chain.is_e_
////	      << " is_egamma_ " << chain.is_egamma_
////	      << " is_e_ " << chain.is_e_
////	      << " validPtr(sig) " << validPtr(chain.sig_)
////	      << " validPtr(ele) " << validPtr(chain.ele_)
////	      << " pt " << ( validPtr(chain.sig_) ? chain.sig_->pt() : -1.)
////	      << std::endl;
//
//    // Init tree here
//    ntuple_.reset();
//
//    // Apply prescale to fake candidates
//    if ( chain.is_e_ == false && 
//	 std::find( cands.begin(), cands.end(), idx ) == cands.end() &&
//	 std::find( cands2.begin(), cands2.end(), idx ) == cands2.end() ) { continue; }
//    
//    // Data sample
//    ntuple_.is_mc( chain.is_mc_ );
//    ntuple_.is_aod( chain.is_aod_ );
//    
//    // Tag-side muon
//    ntuple_.tag_pt( chain.tag_pt_ );
//    ntuple_.tag_eta( chain.tag_eta_ );
//
//    // Truth label
//    ntuple_.is_e( chain.is_e_ );
//    ntuple_.is_other( !chain.is_e_ );
//    
//    // Event stuff
//    ntuple_.fill_evt( event.id() );
//    ntuple_.set_rho( -1. ); //*rhoH_ );
//
//    // EGamma of low-pT
//    ntuple_.is_egamma( chain.is_egamma_ );
//    
//    // Set prescale for BACKGROUND
//    if ( !chain.is_e_ ) { ntuple_.set_prescale( std::abs(prescale_) > 1.e-6 ? prescale_ : 1. ); }
//
//    //@@ Set weight??
//
//    // "signal" info
//    if ( validPtr(chain.sig_) ) {
//      ntuple_.fill_gen( chain.sig_ ); // Add sig info
//    }
//
//    // GsfElectron info
//    if ( validPtr(chain.ele_) ) {
//
//      float mva_value = -999.;//@@
//      float mva_value_retrained = -999.;//@@
//
//      ntuple_.has_ele( chain.ele_match_ );
//      ntuple_.ele_dr( chain.ele_dr_ );
//
//      if ( chain.is_egamma_ ) {
//	if ( mvaValueEGammaH_.isValid() && 
//	     mvaValueEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
//	  mva_value = mvaValueEGammaH_->get( chain.ele_.key() );
//	} else {
//	  std::cout << "[IDNtuplizerSimple::fill] ERROR! Issue matching MVA output to PF GsfElectrons!" << std::endl;
//	}
//	if ( mvaValueEGammaRetrainedH_.isValid() && 
//	     mvaValueEGammaRetrainedH_->size() == gsfElectronsEGammaH_->size() ) {
//	  mva_value_retrained = mvaValueEGammaRetrainedH_->get( chain.ele_.key() );
//	} else {
//	  std::cout << "[IDNtuplizerSimple::fill] ERROR! Issue matching retrained MVA to PF GsfElectrons! "
//	    	    << mvaValueEGammaRetrainedH_.isValid() << " " 
//		    << mvaValueEGammaRetrainedH_->size() << " " 
//		    << gsfElectronsEGammaH_->size() << " " 
//		    << std::endl;
//	}
//      }
//      
//      ntuple_.fill_ele( chain.ele_, 
//			mva_value,
//			mva_value_retrained,
//			-1.,//*rhoH_,
//			chain.is_egamma_ );
//      
//    }
//
//    tree_->Fill();
//    
//  }
//  
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//template <typename T>
//void IDNtuplizerSimple::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
//					edm::Handle< edm::View<T> >& candidates,
//					std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
//					std::vector< DeltaR2<reco::Candidate,T> >& other_cand,
//					bool append ) {
//  std::vector< edm::Ptr<T> > cands;
//  for ( size_t icand = 0; icand < candidates->size(); ++icand ) {
//    edm::Ptr<T> cand(candidates,icand);
//    if ( validPtr(cand) ) { cands.push_back(cand); }
//    else {
//      std::cout << "[IDNtuplizerSimple::sigToCandLinks] ERROR! CandidatePtr:"
//		<< " cand.isNull(): " << cand.isNull()
//		<< " cand.isAvailable(): " << cand.isAvailable()
//		<< std::endl;
//    }
//  }
//  sigToCandLinks<T>( signal_electrons, cands, sig2cand, other_cand, append );
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//template <typename T>
//void IDNtuplizerSimple::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
//					std::vector< edm::Ptr<T> >& candidates,
//					std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
//					std::vector< DeltaR2<reco::Candidate,T> >& other_cand,
//					bool append ) {
//  
//  if ( !append ) { sig2cand.clear(); }
//  if ( !append ) { other_cand.clear(); }
//  
//  // DeltaR2 matches for all combinations of signal electrons and reco::Candidate
//  std::vector< DeltaR2<reco::Candidate,T> > sig2cand_all;
//  sig2cand_all.reserve( sig2cand.size() + signal_electrons.size()*candidates.size() );
//  unsigned int icand = 0;
//  for ( auto cand : candidates ) {
//    ++icand;
//    if ( !validPtr(cand) ) { continue; }
//    if ( !filterCand<T>(cand) ) { continue; }
//    for ( auto sig : signal_electrons ) {
//      sig2cand_all.emplace_back( sig,
//				 cand,
//				 deltaR2<reco::Candidate,T>(sig,cand) );
//    }
//    other_cand.emplace_back( reco::CandidatePtr(), cand, id::NEG_FLOATSQ ); // Null signal
//  }
//
//  // Sort by DeltaR2!!
//  std::sort( sig2cand_all.begin(), 
//	     sig2cand_all.end(), 
//	     DeltaR2<reco::Candidate,T>::compare_by_dr2 );
//  
//  // Select best matches according to DeltaR2 metric
//  sig2cand.clear();
//  for ( auto iter : sig2cand_all ) {
//    auto sig = std::find_if( sig2cand.begin(), 
//			     sig2cand.end(), 
//			     [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
//			       return dr2.obj1_ == iter.obj1_; 
//			     }
//			     );
//    auto cand = std::find_if( sig2cand.begin(), 
//			      sig2cand.end(), 
//			      [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
//				return dr2.obj2_ == iter.obj2_; 
//			      }
//			      );
//    if ( sig == sig2cand.end() && cand == sig2cand.end() ) {
//      if ( iter.dr2_ < dr_max_*dr_max_ ) { sig2cand.push_back(iter); }
//      else { 
//	sig2cand.emplace_back( iter.obj1_, edm::Ptr<T>(), id::NEG_FLOATSQ ); // Null cand
//      }
//    }
//    if ( sig2cand.size() >= signal_electrons.size() ) { break; }
//  }
//  
//  // Remove matched candidates
//  for ( auto iter : sig2cand ) {
//    auto match = std::find_if( other_cand.begin(), 
//			       other_cand.end(), 
//			       [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
//				 return dr2.obj2_ == iter.obj2_; 
//			       }
//			       );
//    if ( match != other_cand.end() ) { other_cand.erase(match); }
//  }
//  
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//template <typename T>
//void IDNtuplizerSimple::match( reco::CandidatePtr& sig,
//			       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
//			       edm::Ptr<T>& cand_ptr, 
//			       float& cand_dr, 
//			       bool& cand_match ) { 
//  
//  // Identify GEN electron matches to candidate
//  auto match = std::find_if( sig2cand.begin(), 
//			     sig2cand.end(), 
//			     [sig](const DeltaR2<reco::Candidate,T>& dr2) { 
//			       return dr2.obj1_ == sig; 
//			     }
//			     );
//  if ( match != sig2cand.end() ) {
//    if ( validPtr(match->obj2_) ) {
//      cand_ptr = match->obj2_; // pass by ref
//      if ( deltaR2(sig,cand_ptr) >= 0. ) cand_dr = sqrt(deltaR2(sig,cand_ptr)); // pass by ref
//      cand_match = ( cand_dr >= 0. ) && ( cand_dr < dr_threshold_ ); // pass by ref
//    }
//  }
//  
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//template <typename T1, typename T2>
//float IDNtuplizerSimple::deltaR2( edm::Ptr<T1>& cand1, edm::Ptr<T2>& cand2 ) {
//  return reco::deltaR2( cand1->eta(), 
//			cand1->phi(),
//			cand2->eta(),
//			cand2->phi() ); 
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//template <typename T>
//bool IDNtuplizerSimple::filterCand( edm::Ptr<T>& cand ) { 
//  if ( cand->pt() < minPt_ ) { return false; }
//  //if ( std::abs(cand->eta()) > 2.5 ) { return false; }
//  return true;
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//bool IDNtuplizerSimple::filterCand( edm::Ptr<reco::GsfElectron>& ele ) { 
//  if ( ele->gsfTrack()->ptMode() < minPt_ ) { return false; }
//  //if ( std::abs(gsf->gsfTrack()->etaMode()) > 2.5 ) { return false; }
//  return true; 
//}
//
//////////////////////////////////////////////////////////////////////////////////
////
//template <typename T> 
//bool IDNtuplizerSimple::validPtr( edm::Ptr<T>& ptr ) {
//  return ( ptr.isNonnull() && ptr.isAvailable() );
//}
//
//

////////////////////////////////////////////////////////////////////////////////
//
DEFINE_FWK_MODULE(IDNtuplizerSimple);
