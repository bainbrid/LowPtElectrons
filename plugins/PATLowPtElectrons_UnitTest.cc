#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class PATLowPtElectrons_UnitTest: public edm::EDAnalyzer {

public:
  
  explicit PATLowPtElectrons_UnitTest( const edm::ParameterSet& );
  
  ~PATLowPtElectrons_UnitTest() {
  
  }

private:

  virtual void analyze( const edm::Event&, const edm::EventSetup& );

  const edm::EDGetTokenT< edm::View<pat::Electron> > patElectrons_;
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaId_;
  bool verbose_;

};

PATLowPtElectrons_UnitTest::PATLowPtElectrons_UnitTest( const edm::ParameterSet& cfg ) :
  patElectrons_{consumes<edm::View<pat::Electron> >(cfg.getParameter<edm::InputTag>("patElectrons"))},
  gsfElectrons_{consumes<edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))},
  mvaId_{consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaId"))},
  verbose_(cfg.getParameter<bool>("verbose"))
  {;}

void PATLowPtElectrons_UnitTest::analyze( const edm::Event& iEvent, 
					  const edm::EventSetup& iSetup )
{

  std::cout << "[PATLowPtElectrons_UnitTest::analyze]" << std::endl;

  edm::Handle<edm::View<pat::Electron> > patElectrons;
  iEvent.getByToken(patElectrons_, patElectrons);
  if ( patElectrons.failedToGet() ) throw cms::Exception("patElectrons.failedToGet()");
  if ( !patElectrons.isValid() ) throw cms::Exception("!patElectrons.isValid()");

  edm::Handle<edm::View<reco::GsfElectron> > gsfElectrons;
  iEvent.getByToken(gsfElectrons_, gsfElectrons);
  if ( gsfElectrons.failedToGet() ) throw cms::Exception("gsfElectrons.failedToGet()");
  if ( !gsfElectrons.isValid() ) throw cms::Exception("!gsfElectrons.isValid()");

  edm::Handle<edm::ValueMap<float> > mvaId;
  iEvent.getByToken(mvaId_, mvaId);
  if ( mvaId.failedToGet() ) throw cms::Exception("mvaId.failedToGet()");
  if ( !mvaId.isValid() ) throw cms::Exception("!mvaId.isValid()");

  std::cout << "patElectrons: " << patElectrons->size()
	    << " gsfElectrons: " << gsfElectrons->size()
	    << std::endl;

  // LOOP THROUGH PAT ELECTRONS
  for ( size_t i = 0; i < patElectrons->size(); ++i ) {

    // GET PAT ELECTRON AND ITS CORE
    edm::Ref<edm::View<pat::Electron> > pat(patElectrons,i);

    // GET ELECTRON CORE, GSF TRACK, SUPERCLUSTER, SEED CALOCLUSTER, (CLOSEST) TRACK
    edm::Ref<std::vector<reco::GsfElectronCore> > core = pat->core();
    edm::Ref<std::vector<reco::GsfTrack> > gsf = pat->gsfTrack();
    edm::Ref<std::vector<reco::SuperCluster> > sc = pat->superCluster();
    edm::Ptr<reco::CaloCluster> seed = pat->seed();
    edm::Ref<std::vector<reco::Track> > trk = pat->closestCtfTrackRef();

    // PRINT OUT
    if ( verbose_ ) 
      std::cout << "index: " << i
	        << std::endl
		<< " PAT:  "
		<< " null: " << pat.isNull()
		<< " id:   " << pat.id()
		<< " key:  " << pat.key()
		<< std::endl
		<< " CORE: "
		<< " null: " << core.isNull()
		<< " id:   " << core.id()
		<< " key:  " << core.key()
		<< std::endl
		<< " GSF:  "
		<< " null: " << gsf.isNull()
		<< " id:   " << gsf.id()
		<< " key:  " << gsf.key()
		<< std::endl
		<< " SC:   "
		<< " null: " << sc.isNull()
		<< " id:   " << sc.id()
		<< " key:  " << sc.key()
		<< std::endl
		<< " SEED: "
		<< " null: " << seed.isNull()
		<< " id:   " << seed.id()
		<< " key:  " << seed.key()
		<< std::endl
		<< " TRK:  "
		<< " null: " << trk.isNull()
		<< " id:   " << trk.id()
		<< " key:  " << trk.key()
		<< std::endl
		<< " ID:   "
	        << " has: " << pat->hasUserFloat("ID")
		<< " val: " << ( pat->hasUserFloat("ID") ? pat->userFloat("ID") : -100. )
                << std::endl;
    
    // ATTEMPT TO MATCH TO GSF ELECTRON VIA CORE
    auto matched_ele = std::find_if( gsfElectrons->begin(), 
				     gsfElectrons->end(), 
				     [core]( const reco::GsfElectron& ele ) { 
				       return ele.core() == core;
				     }
				     );
    if ( matched_ele == gsfElectrons->end() ) {
      throw cms::Exception("NoMatch")
	<< " No match from pat::Electron to reco::GsfElectron via reco::GsfElectronCore!" 
	<< std::endl
	<< " pat.isNonnull():" << pat.isNonnull() 
	<< " pat->core().isNonnull():" << pat->core().isNonnull() 
	<< " pat.id(): " << pat.id() 
	<< " pat.key(): " << pat.key() 
	<< " pat->pt(): " << pat->pt()
	<< " pat->eta(): " << pat->eta()
	<< std::endl;
    }

    // GET ELECTRON CORE, GSF TRACK, SUPERCLUSTER, SEED CALOCLUSTER, (CLOSEST) TRACK
    edm::Ref<edm::View<reco::GsfElectron> > ele(gsfElectrons,std::distance(gsfElectrons->begin(),matched_ele));
    edm::Ref<std::vector<reco::GsfElectronCore> > core_ = ele->core();
    edm::Ref<std::vector<reco::GsfTrack> > gsf_ = ele->gsfTrack();
    edm::Ref<std::vector<reco::SuperCluster> > sc_ = ele->superCluster();
    edm::Ptr<reco::CaloCluster> seed_ = ele->electronCluster(); //@@ pat->seed();
    edm::Ref<std::vector<reco::Track> > trk_ = ele->closestCtfTrackRef();
    float id = mvaId->get( ele.key() );

    // PRINT OUT
    if ( verbose_ ) 
      std::cout	<< " ELE:  "
		<< " null: " << ele.isNull()
		<< " id:   " << ele.id()
		<< " key:  " << ele.key()
		<< std::endl
		<< " CORE: "
		<< " null: " << core_.isNull()
		<< " id:   " << core_.id()
		<< " key:  " << core_.key()
		<< std::endl
		<< " GSF:  "
		<< " null: " << gsf_.isNull()
		<< " id:   " << gsf_.id()
		<< " key:  " << gsf_.key()
		<< std::endl
		<< " SC:   "
		<< " null: " << sc_.isNull()
		<< " id:   " << sc_.id()
		<< " key:  " << sc_.key()
		<< std::endl
		<< " SEED: "
		<< " null: " << seed_.isNull()
		<< " id:   " << seed_.id()
		<< " key:  " << seed_.key()
		<< std::endl
		<< " TRK:  "
		<< " null: " << trk_.isNull()
		<< " id:   " << trk_.id()
		<< " key:  " << trk_.key()
		<< std::endl
		<< " ID:   "
		<< " val: " << id
		<< std::endl;

  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATLowPtElectrons_UnitTest);
