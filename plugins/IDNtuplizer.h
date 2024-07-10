#ifndef LowPtElectrons_LowPtElectrons_IDNtuplizer
#define LowPtElectrons_LowPtElectrons_IDNtuplizer

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
#include "FWCore/Framework/interface/EDFilter.h" // EDAnalyzer.h
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "LowPtElectrons/LowPtElectrons/interface/Common.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDNtuple.h"
#include "TRandom3.h"
#include "TTree.h"
#include <set>
#include <vector>
#include <math.h>
#include <boost/core/demangle.hpp>
#include <algorithm>
#include <random>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
class IDNtuplizer : public edm::EDFilter { // edm::EDAnalyzer
  
public:
  
  explicit IDNtuplizer( const edm::ParameterSet& );
  ~IDNtuplizer();
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual bool filter( edm::Event&, const edm::EventSetup& ) override; // analyze(const,const)
  
  ///////////////
  // Main methods
  ///////////////

  // Reads all collections from the Event
  void readCollections( const edm::Event&, const edm::EventSetup& );
  
  // Wraps other methods to provide a sample of "signal" electrons
  void signalElectrons( std::set<reco::CandidatePtr>& signal_electrons,
			std::set<reco::CandidatePtr>& tag_side_muons );
  
  // GEN-based method to provide a sample of "signal" electrons
  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B, 
			  std::set<reco::GenParticlePtr>& gen_muons, 
			  float tag_muon_pt_threshold,
			  float tag_muon_eta_threshold );
  
  // Data-based method to provide a sample of "signal" electrons
  void electronsFromB( std::set<reco::CandidatePtr>& electrons_from_B, 
		       std::set<reco::CandidatePtr>& muons, 
		       float tag_muon_pt_threshold,
		       float tag_muon_eta_threshold );
  
  // 
  void createChains( std::set<reco::CandidatePtr>& signal_electrons,
		     std::set<reco::CandidatePtr>& tag_side_muons );

  // Method that creates ElectronChain objects for both low-pT and PF electron signal candiates
  void signal( std::set<reco::CandidatePtr>& signal_electrons,
	       std::set<reco::CandidatePtr>& tag_side_muons,
	       std::vector<SigToTrkDR2>& sig2trk,
	       std::vector<SigToGsfDR2>& sig2gsf,
	       std::vector<SigToGsfDR2>& sig2pfgsf,
	       std::vector<SigToTrkDR2>& other_trk,
	       std::vector<TrkToGsfDR2>& trk2gsf,
	       std::vector<TrkToGsfDR2>& trk2pfgsf,
	       std::vector<GsfToGsfDR2>& gsf2pfgsf,
	       std::vector<GsfToEleDR2>& gsf2ele,
	       std::vector<GsfToEleDR2>& pfgsf2ele );

  // Method that creates ElectronChain objects for both low-pT and PF electron bkgd candidates
  void bkgd( std::set<reco::CandidatePtr>& signal_electrons,
	     std::set<reco::CandidatePtr>& tag_side_muons,
	     std::vector<SigToTrkDR2>& sig2trk,
	     std::vector<SigToGsfDR2>& sig2gsf,
	     std::vector<SigToGsfDR2>& sig2pfgsf,
	     std::vector<SigToTrkDR2>& other_trk,
	     std::vector<TrkToGsfDR2>& trk2gsf,
	     std::vector<TrkToGsfDR2>& trk2pfgsf,
	     std::vector<GsfToGsfDR2>& gsf2pfgsf,
	     std::vector<GsfToEleDR2>& gsf2ele,
	     std::vector<GsfToEleDR2>& pfgsf2ele );
  
  // Fills tree per ElectronChain object
  void fill( const edm::Event& event, const edm::EventSetup& setup );

  //////////////////
  // Utility methods
  //////////////////

  // Extracts TrackPtrs to common vector, obtained from AOD or mAOD collections via Handles
  void extractTrackPtrs();
  
  // Various navigation methods (high- to low-level objects)
  bool gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk, bool is_egamma);// = false );
  bool eleToGsf( reco::GsfElectronPtr& ele, reco::GsfTrackPtr& gsf );
  bool eleToTrk( reco::GsfElectronPtr& gsf, reco::TrackPtr& trk, bool is_egamma);// = false );
  bool gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed );
  bool seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk );
  bool seedToCalo( reco::ElectronSeedPtr& seed, reco::CaloClusterPtr& calo );

  // Various navigation maps (low- to high-level objects)
  void trkToGsfLinks( std::vector<reco::TrackPtr>& ctfTracks,
		      edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::vector<TrkToGsfDR2>& trk2gsf,
		      bool is_egamma);// = false );
  void trkToEleLinks( std::vector<reco::TrackPtr>& ctfTracks,
		      edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
		      std::vector<TrkToEleDR2>& trk2ele,
		      bool is_egamma);// = false );
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
  float deltaR2( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk ); // overload
  float deltaR2( reco::GsfElectronPtr& ele, reco::TrackPtr& trk ); // overload
  float deltaR2( reco::GsfTrackPtr& pfgsf, reco::GsfTrackPtr& gsf ); // overload
  float deltaR2( reco::CandidatePtr& sig, reco::GsfTrackPtr& gsf ); // overload

  // Filter track candidates by quality flag, simple pass for all other cands
  template <typename T> 
  bool filterCand( edm::Ptr<T>& cand );
  bool filterCand( edm::Ptr<reco::Track>& trk );
  bool filterCand( edm::Ptr<reco::GsfTrack>& gsf );

  // Check is Ptr is valid and available
  template <typename T> bool validPtr( edm::Ptr<T>& ptr );


  // Electron "images"
  typedef math::XYZVector Vector;
  typedef math::XYZPoint Point;
  BaseParticlePropagator extrapolate_track( const Vector& mom, const Point& pos, int charge,
					    int& reach_ECAL, GlobalPoint& pos_ECAL,
					    int& reach_HCAL, GlobalPoint& pos_HCAL,
					    int& reach_EXIT, GlobalPoint& pos_EXIT );
  //void debug_image( const edm::Event& event, const edm::EventSetup& setup );
  void build_image( const edm::Event& event, const edm::EventSetup& setup );

  inline float adj_eta( float eta, float ref_eta, int charge = 0 ) { 
    return eta-ref_eta; 
    //return eta; //@@ DON'T MODIFY ETA !!!
  }
  inline float adj_phi( float phi, float ref_phi, int charge = 0 ) { 
    return charge == 0 ? reco::deltaPhi(phi,ref_phi) : float(charge)*reco::deltaPhi(phi,ref_phi); 
    //return charge == 0 ? phi : float(charge)*phi; //@@ DON'T SUBTRACT REF_PHI !!!
    //return phi; //@@ DON'T MODIFY PHI !!!
  }
  
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
  float gsfPtThreshold_;
  float gsfEtaThreshold_;
  float tagMuonPtThreshold_;
  float tagMuonEtaThreshold_;
  bool filterNtupleContent_;

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

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > pfToPackedCands_; // MINIAOD
  edm::Handle< edm::Association<pat::PackedCandidateCollection> > pfToPackedCandsH_;
  
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

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPtDepth10_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtDepth10H_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPtDepth11_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtDepth11H_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPtDepth13_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtDepth13H_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPtDepth15_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtDepth15H_;

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

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGammaRetrained_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaRetrainedH_;

  // Conversions

  //@@ const edm::EDGetTokenT< edm::ValueMap<float> > convVtxFitProb_;

  std::vector<ElectronChain> chains_;

  std::vector<reco::TrackPtr> tracks_;
  PdgIds pdgids_;

  //////////////////
  // Obsolete methods
  //////////////////

  // Top-level method that creates ElectronChain objects for low-pT electrons
  void lowPtElectrons( std::set<reco::CandidatePtr>& signal_electrons,
		       std::vector<SigToTrkDR2>& sig2trk,
		       std::vector<SigToTrkDR2>& other_trk,
		       std::vector<SigToGsfDR2>& sig2gsf,
		       std::vector<TrkToGsfDR2>& trk2gsf,
		       std::vector<GsfToGsfDR2>& gsf2pfgsf );

  // Method that creates ElectronChain objects for low-pT electron signal candiates
  void lowPtElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
			      std::vector<SigToTrkDR2>& sig2trk,
			      std::vector<SigToTrkDR2>& other_trk,
			      std::vector<SigToGsfDR2>& sig2gsf,
			      std::vector<GsfToGsfDR2>& gsf2pfgsf,
			      std::vector<GsfToEleDR2>& gsf2ele );

  // Method that creates ElectronChain objects for low-pT electron fake candidates
  void lowPtElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
			     std::vector<TrkToGsfDR2>& trk2gsf,
			     std::vector<GsfToGsfDR2>& gsf2pfgsf,
			     std::vector<GsfToEleDR2>& gsf2ele );
  
  // Top-level method that creates ElectronChain objects for EGamma electrons
  void pfElectrons( std::set<reco::CandidatePtr>& signal_electrons,
		    std::vector<SigToTrkDR2>& sig2trk,
		    std::vector<SigToTrkDR2>& other_trk,
		    std::vector<SigToGsfDR2>& sig2gsf,
		    std::vector<TrkToGsfDR2>& trk2gsf,
		    std::vector<GsfToGsfDR2>& gsf2pfgsf );

  // Method that creates ElectronChain objects for EGamma signal candidates
  void pfElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
			   std::vector<SigToTrkDR2>& sig2trk,
			   std::vector<SigToTrkDR2>& other_trk,
			   std::vector<SigToGsfDR2>& sig2gsf,
			   std::vector<SigToGsfDR2>& sig2pfgsf,
			   std::vector<TrkToGsfDR2>& trk2pfgsf,
			   std::vector<GsfToGsfDR2>& gsf2pfgsf,
			   std::vector<GsfToEleDR2>& pfgsf2ele );
  
  // Method that creates ElectronChain objects for EGamma fake candidates
  void pfElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
			  std::vector<TrkToGsfDR2>& trk2gsf,
			  std::vector<TrkToGsfDR2>& trk2pfgsf,
			  std::vector<GsfToGsfDR2>& gsf2pfgsf,
			  std::vector<GsfToEleDR2>& pfgsf2ele );
  
  // Method that creates ElectronChain objects for EGamma fake candidates
  void pfElectrons_fakes_temp( std::vector<SigToTrkDR2>& other_trk,
			       std::vector<TrkToGsfDR2>& trk2gsf,
			       std::vector<TrkToGsfDR2>& trk2pfgsf,
			       std::vector<GsfToGsfDR2>& gsf2pfgsf,
			       std::vector<TrkToEleDR2>& trk2ele );

  std::map<int,int> pf_pdgids_;

};

#endif // LowPtElectrons_LowPtElectrons_IDNtuplizer
