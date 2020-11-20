#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"

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
  const edm::EDGetTokenT< edm::View<reco::GsfTrack> > gsfTracks_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaId_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaUnbiased_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaPtbiased_;
  const edm::InputTag lowPtGsfLinksTag1_;
  edm::EDGetTokenT<edm::Association<pat::PackedCandidateCollection> > lowPtGsfLinks1_;
  const edm::InputTag lowPtGsfLinksTag2_;
  edm::EDGetTokenT<edm::Association<pat::PackedCandidateCollection> > lowPtGsfLinks2_;
  bool verbose_;
  const edm::EDGetTokenT<reco::PFCandidateCollection> pfcands_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> packed_;
  const edm::EDGetTokenT< edm::View<pat::Electron> > patPfElectrons_;
  
};

PATLowPtElectrons_UnitTest::PATLowPtElectrons_UnitTest( const edm::ParameterSet& cfg ) :
  patElectrons_(consumes<edm::View<pat::Electron> >(cfg.getParameter<edm::InputTag>("patElectrons"))),
  gsfElectrons_(consumes<edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
  gsfTracks_(consumes<edm::View<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
  mvaId_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaId"))),
  mvaUnbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaUnbiased"))),
  mvaPtbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaPtbiased"))),
  lowPtGsfLinksTag1_(cfg.getParameter<edm::InputTag>("lowPtGsfLinks1")),
  lowPtGsfLinksTag2_(cfg.getParameter<edm::InputTag>("lowPtGsfLinks2")),
  verbose_(cfg.getParameter<bool>("verbose")),
  pfcands_(consumes<reco::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("PFCandidates"))),
  packed_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
  patPfElectrons_(consumes<edm::View<pat::Electron> >(cfg.getParameter<edm::InputTag>("patPfElectrons")))
{
  lowPtGsfLinks1_ = consumes<edm::Association<pat::PackedCandidateCollection> >(lowPtGsfLinksTag1_);
  lowPtGsfLinks2_ = consumes<edm::Association<pat::PackedCandidateCollection> >(lowPtGsfLinksTag2_);
}

void PATLowPtElectrons_UnitTest::analyze( const edm::Event& event, 
					  const edm::EventSetup& setup )
{
  
   if ( verbose_ ) { std::cout << "[PATLowPtElectrons_UnitTest::analyze]" << std::endl; }

   edm::Handle<edm::View<pat::Electron> > patElectrons;
   event.getByToken(patElectrons_, patElectrons);
   if ( patElectrons.failedToGet() ) edm::LogWarning("patElectrons.failedToGet()");
   if ( !patElectrons.isValid() ) edm::LogWarning("!patElectrons.isValid()");

   edm::Handle<edm::View<reco::GsfElectron> > gsfElectrons;
   event.getByToken(gsfElectrons_, gsfElectrons);
   if ( gsfElectrons.failedToGet() ) edm::LogWarning("gsfElectrons.failedToGet()");
   if ( !gsfElectrons.isValid() ) edm::LogWarning("!gsfElectrons.isValid()");

   edm::Handle<edm::View<reco::GsfTrack> > gsfTracks;
   event.getByToken(gsfTracks_, gsfTracks);
   if ( gsfTracks.failedToGet() ) edm::LogWarning("gsfTracks.failedToGet()");
   if ( !gsfTracks.isValid() ) edm::LogWarning("!gsfTracks.isValid()");

   edm::Handle<edm::ValueMap<float> > mvaId;
   event.getByToken(mvaId_, mvaId);
   if ( mvaId.failedToGet() ) edm::LogWarning("mvaId.failedToGet()");
   if ( !mvaId.isValid() ) edm::LogWarning("!mvaId.isValid()");

   edm::Handle<edm::ValueMap<float> > mvaUnbiased;
   event.getByToken(mvaUnbiased_, mvaUnbiased);
   if ( mvaUnbiased.failedToGet() ) edm::LogWarning("mvaUnbiased.failedToGet()");
   if ( !mvaUnbiased.isValid() ) edm::LogWarning("!mvaUnbiased.isValid()");

   edm::Handle<edm::ValueMap<float> > mvaPtbiased;
   event.getByToken(mvaPtbiased_, mvaPtbiased);
   if ( mvaPtbiased.failedToGet() ) edm::LogWarning("mvaPtbiased.failedToGet()");
   if ( !mvaPtbiased.isValid() ) edm::LogWarning("!mvaPtbiased.isValid()");

   edm::Handle<edm::Association<pat::PackedCandidateCollection> > lowPtGsfLinks1;
   event.getByToken(lowPtGsfLinks1_, lowPtGsfLinks1);
   if ( lowPtGsfLinks1.failedToGet() ) edm::LogWarning("lowPtGsfLinks1.failedToGet()");
   if ( !lowPtGsfLinks1.isValid() ) edm::LogWarning("!lowPtGsfLinks1.isValid()");

   edm::Handle<edm::Association<pat::PackedCandidateCollection> > lowPtGsfLinks2;
   event.getByToken(lowPtGsfLinks2_, lowPtGsfLinks2);
   if ( lowPtGsfLinks2.failedToGet() ) edm::LogWarning("lowPtGsfLinks2.failedToGet()");
   if ( !lowPtGsfLinks2.isValid() ) edm::LogWarning("!lowPtGsfLinks2.isValid()");

   if ( verbose_ ) { 
     std::cout << "patElectrons: " << patElectrons->size()
	       << " gsfElectrons: " << gsfElectrons->size()
	       << " gsfTracks: " << gsfTracks->size()
	       << std::endl;
   }

   // PF COLLECTIONS
   
   edm::Handle<reco::PFCandidateCollection> pfcands;
   event.getByToken(pfcands_, pfcands);
   if ( pfcands.failedToGet() ) edm::LogWarning("pfcands.failedToGet()");
   if ( !pfcands.isValid() ) edm::LogWarning("!pfcands.isValid()");
   
   edm::Handle<pat::PackedCandidateCollection> packed;
   event.getByToken(packed_, packed);
   if ( packed.failedToGet() ) edm::LogWarning("packed.failedToGet()");
   if ( !packed.isValid() ) edm::LogWarning("!packed.isValid()");

   edm::Handle<edm::View<pat::Electron> > patPfElectrons;
   event.getByToken(patPfElectrons_, patPfElectrons);
   if ( patPfElectrons.failedToGet() ) edm::LogWarning("patPfElectrons.failedToGet()");
   if ( !patPfElectrons.isValid() ) edm::LogWarning("!patPfElectrons.isValid()");

   // LOOP THROUGH PAT ELECTRONS
   for ( size_t i = 0; i < patElectrons->size(); ++i ) {

     // GET PAT ELECTRON AND ITS CORE
     edm::Ref<edm::View<pat::Electron> > pat(patElectrons,i);

     // GET ELECTRON CORE, GSF TRACK, SUPERCLUSTER, SEED CALOCLUSTER, (CLOSEST) TRACK
     edm::RefToBase<reco::GsfElectron> base = edm::RefToBase<reco::GsfElectron>(pat);
     const edm::Ptr<reco::Candidate>& orig = pat->originalObjectRef();
     edm::Ref<std::vector<reco::GsfElectronCore> > core = pat->core();
     edm::Ref<std::vector<reco::GsfTrack> > gsf = pat->gsfTrack();
     edm::Ref<std::vector<reco::SuperCluster> > sc = pat->superCluster();
     edm::Ptr<reco::CaloCluster> seed = pat->seed();
     edm::Ref<std::vector<reco::Track> > trk = pat->closestCtfTrackRef();
     const edm::Ptr<pat::PackedCandidate>* link1 = pat->hasUserData("ele2packed") ? 
       pat->userData<edm::Ptr<pat::PackedCandidate> >("ele2packed") : NULL;
     const edm::Ptr<pat::PackedCandidate>* link2 = pat->hasUserData("ele2lost") ? 
       pat->userData<edm::Ptr<pat::PackedCandidate> >("ele2lost") : NULL;
     float id = ( pat->isElectronIDAvailable("ID") ? pat->electronID("ID") : -100. );
     float unbiased = ( pat->isElectronIDAvailable("unbiased") ? pat->electronID("unbiased") : -100. );
     float ptbiased = ( pat->isElectronIDAvailable("ptbiased") ? pat->electronID("ptbiased") : -100. );
     const pat::PFIsolation& iso = pat->miniPFIsolation();

     // Try to match to PF electron
     edm::Ptr<pat::Electron> matched;
     auto matched1 = std::find_if( patPfElectrons->begin(), 
				   patPfElectrons->end(), 
				   [link1]( const pat::Electron& pf ) {
				     for ( auto packed : pf.associatedPackedPFCandidates() ) {
				       if ( *link1 == edm::refToPtr(packed) ) return true;
				     }
				     return false;
				   }
				   );
     auto matched2 = std::find_if( patPfElectrons->begin(), 
				   patPfElectrons->end(), 
				   [link2]( const pat::Electron& pf ) {
				     for ( auto packed : pf.associatedPackedPFCandidates() ) {
				       if ( *link2 == edm::refToPtr(packed) ) return true;
				     }
				     return false;
				   }
				   );
     
     // PRINT OUT
     if ( verbose_ ) {
       
       std::cout << "index: "<<i<<std::endl
		 << " PAT:  "<<" null:"<<pat.isNull() <<" id:"<<pat.id() <<" key:"<<pat.key()<<" pt:"<<pat->pt()<<std::endl
		 << " ELE:  "<<" null:"<<base.isNull()<<" id:"<<base.id()<<" key:"<<base.key()<<" pt:"<<base->pt()<<std::endl
		 << " ORIG: "<<" null:"<<orig.isNull()<<" id:"<<orig.id()<<" key:"<<orig.key()<<" pt:"<<orig->pt()<<std::endl
		 << " CORE: "<<" null:"<<core.isNull()<<" id:"<<core.id()<<" key:"<<core.key()<<std::endl
		 << " gsf1: "<<" null:"<<pat->gsfTrack().isNull()<<" id:"<<pat->gsfTrack().id()<<" key:"<<pat->gsfTrack().key()<<" pt:"<<pat->gsfTrack()->pt()<<std::endl
		 << " gsf2: "<<" null:"<<core->gsfTrack().isNull()<<" id:"<<core->gsfTrack().id()<<" key:"<<core->gsfTrack().key()<<" pt:"<<core->gsfTrack()->pt()<<std::endl
		 << " GSF:  "<<" null:"<<gsf.isNull() <<" id:"<<gsf.id() <<" key:"<<gsf.key()<<" pt:"<<gsf->pt()<<std::endl
		 << " SC:   "<<" null:"<<sc.isNull()  <<" id:"<<sc.id()  <<" key:"<<sc.key()<<" E:"<<sc->energy()<<std::endl
		 << " SEED: "<<" null:"<<seed.isNull()<<" id:"<<seed.id()<<" key:"<<seed.key()<<" E:"<<seed->energy()<<std::endl
		 << " TRK:  "<<" null:"<<trk.isNull() <<" id:"<<trk.id() <<" key:"<<trk.key()<<std::endl
		 << " OBJS: "
		 << " nPat:"<<patElectrons->size()
		 << " nEle:"<<gsfElectrons->size()
		 << " nGsf:"<<gsfTracks->size()
		 << " nLinks1:"<<lowPtGsfLinks1->size()
		 << " nLinks2:"<<lowPtGsfLinks2->size()
		 <<std::endl;
       std::cout << " NAMES:" <<" size:"<<pat->userDataNames().size();
       for ( auto name : pat->userDataNames() ) { std::cout << ", '" << name << "'"; }
       std::cout << std::endl
		 << " LINK1:" 
		 << " null:"<<(*link1).isNull() 
		 << " id:"<<(*link1).id()  
		 << " key:"<<(*link1).key() 
		 << " pt:"<<((*link1).isNull() ?-1.:(*link1)->pt())
		 <<std::endl
		 << " LINK2:" 
		 << " null:"<<(*link2).isNull() 
		 << " id:"<<(*link2).id()  
		 << " key:"<<(*link2).key() 
		 << " pt:"<<((*link2).isNull() ?-1.:(*link2)->pt())
		 <<std::endl
		 << " PF1:  "<<" ok?:"<<(matched1!=patPfElectrons->end())<<" pt:"<<(matched1!=patPfElectrons->end()?matched1->pt():-1.)<<std::endl
		 << " PF2:  "<<" ok?:"<<(matched2!=patPfElectrons->end())<<" pt:"<<(matched2!=patPfElectrons->end()?matched2->pt():-1.)<<std::endl
		 << " ID:   "<<" val:"<<id<<std::endl
		 << " UNB:  "<<" val:"<<unbiased<<std::endl
		 << " PTB:  "<<" val:"<<ptbiased<<std::endl
		 << " ISO:  "
		 <<" ch:"<<iso.chargedHadronIso()
		 <<" nh:"<<iso.neutralHadronIso()
		 <<" ph:"<<iso.photonIso()
		 <<" pu:"<<iso.puChargedHadronIso()
		 << std::endl << std::endl;
     }

     // IF CORE IS NOT EMBEDDED, ATTEMPT TO MATCH PAT TO GSF ELECTRON VIA CORE
     auto matched_ele = std::find_if( gsfElectrons->begin(), 
				      gsfElectrons->end(), 
				      [pat]( const reco::GsfElectron& ele ) {
					return 
					pat.isNull() ||
					pat->core().isNull() ||
					pat->core()->gsfTrack().isNull() ||
					ele.core().isNull() ||
					ele.core()->gsfTrack().isNull() ? 
					false :
					pat->core()->gsfTrack() == ele.core()->gsfTrack();
				      }
				      );
     if ( matched_ele == gsfElectrons->end() ) {
       edm::LogWarning("No PAT-to-GSF match")
	 << " No match from pat::Electron to reco::GsfElectron via reco::GsfTrack!"
	 << std::endl
	 << " pat.isNonnull():" << pat.isNonnull()
	 << " pat->core().isNonnull():" << pat->core().isNonnull()
	 << " pat->core()->gsfTrack().isNonnull():" << pat->core()->gsfTrack().isNonnull()
	 << " pat.id(): " << pat.id()
	 << " pat.key(): " << pat.key()
	 << " pat->pt(): " << pat->pt()
	 << " pat->eta(): " << pat->eta()
	 << std::endl;
       continue;
     }

     // GET ELECTRON CORE, GSF TRACK, SUPERCLUSTER, SEED CALOCLUSTER, (CLOSEST) TRACK
     size_t iele_ = std::distance(gsfElectrons->begin(),matched_ele);
     edm::Ref<edm::View<reco::GsfElectron> > ele_(gsfElectrons,iele_);
     edm::Ref<std::vector<reco::GsfElectronCore> > core_ = ele_->core();
     edm::Ref<std::vector<reco::GsfTrack> > gsf_ = ele_->gsfTrack();
     edm::Ref<std::vector<reco::SuperCluster> > sc_ = ele_->superCluster();
     edm::Ptr<reco::CaloCluster> sd_ = ele_->electronCluster(); //@@ pat->seed();
     edm::Ref<std::vector<reco::Track> > trk_ = ele_->closestCtfTrackRef();
     std::string instance1 = lowPtGsfLinksTag1_.instance();
     std::string instance2 = lowPtGsfLinksTag2_.instance();
     edm::Ptr<pat::PackedCandidate> link1_ = edm::refToPtr( (*lowPtGsfLinks1)[gsf_] );
     edm::Ptr<pat::PackedCandidate> link2_ = edm::refToPtr( (*lowPtGsfLinks2)[gsf_] );
     float id_ = mvaId->get( ele_.key() );
     float unbiased_ = mvaUnbiased->get( ele_->gsfTrack().key() );
     float ptbiased_ = mvaPtbiased->get( ele_->gsfTrack().key() );

     // PRINT OUT
     if ( verbose_ ) 
       std::cout << " ELE:  "<<" null:"<<ele_.isNull() <<" id:"<<ele_.id() <<" key:"<<ele_.key()<<" pt:"<<ele_->pt()
		 << " IELE:"<<iele_<<std::endl
		 << " CORE: "<<" null:"<<core_.isNull()<<" id:"<<core_.id()<<" key:"<<core_.key()<<std::endl
		 << " gsf1: "<<" null:"<<ele_->gsfTrack().isNull()<<" id:"<<ele_->gsfTrack().id()<<" key:"<<ele_->gsfTrack().key()<<" pt:"<<ele_->gsfTrack()->pt()<<std::endl
		 << " gsf2: "<<" null:"<<core_->gsfTrack().isNull()<<" id:"<<core_->gsfTrack().id()<<" key:"<<core_->gsfTrack().key()<<" pt:"<<core_->gsfTrack()->pt()<<std::endl
		 << " GSF:  "<<" null:"<<gsf_.isNull() <<" id:"<<gsf_.id() <<" key:"<<gsf_.key()<<" pt:"<<gsf_->pt()<<std::endl
		 << " SC:   "<<" null:"<<sc_.isNull()  <<" id:"<<sc_.id()  <<" key:"<<sc_.key()<<" E:"<<sc_->energy()<<std::endl
		 << " SEED: "<<" null:"<<sd_.isNull()<<" id:"<<sd_.id()<<" key:"<<sd_.key()<<" E:"<<sd_->energy()<<std::endl
		 << " TRK:  "<<" null:"<<trk_.isNull() <<" id:"<<trk_.id() <<" key:"<<trk_.key()<<std::endl
		 << " LINK1:" 
		 << " null:"<<link1_.isNull() 
		 << " id:"<<link1_.id()  
		 << " key:"<<link1_.key()
		 << " pt:"<<(link1_.isNull()?-1.:link1_->pt())
		 << " name:"<<instance1<<std::endl
		 << " LINK2:" 
		 << " null:"<<link2_.isNull() 
		 << " id:"<<link2_.id()  
		 << " key:"<<link2_.key()
		 << " pt:"<<(link2_.isNull()?-1.:link2_->pt())
		 << " name:"<<instance2<<std::endl
		 << " ID:   "<<" val:  "<<id_<<std::endl
		 << " UNB:  "<<" val:  "<<unbiased_<<std::endl
		 << " PTB:  "<<" val:  "<<ptbiased_<<std::endl
		 << std::endl;
     
     // CHECKS
     if ( core.isNonnull() && core_.isNonnull() && 
	  core->gsfTrack().isNonnull() && core_->gsfTrack().isNonnull() && 
	  std::fabs(core->gsfTrack()->pt() - core_->gsfTrack()->pt()) > 1.e-6 )
       edm::LogWarning("NoMatch") << "reco::GsfElectronCore: " 
				  << " core->gsfTrack()->pt(): " << core->gsfTrack()->pt()
				  << " core_->gsfTrack()->pt(): " << core_->gsfTrack()->pt()
				  << std::endl;
     if ( gsf.isNonnull() && gsf_.isNonnull() && 
	  std::fabs(gsf->pt() - gsf_->pt()) > 1.e-6 )
       edm::LogError("NoMatch") << "reco::GsfTrack" << std::endl;
     if ( sc.isNonnull() && sc_.isNonnull() && 
	  std::fabs(sc->energy() - sc_->energy()) > 1.e-6 )
       edm::LogError("NoMatch") << "reco::SuperCluster" << std::endl;
     if ( seed.isNonnull() && sd_.isNonnull() && 
	  std::fabs(seed->energy() - sd_->energy()) > 1.e-6 )
       edm::LogWarning("NoMatch") << "Seed reco::CaloCluster" << std::endl;
//     if ( trk.isNonnull() && trk_.isNonnull() && 
//	  std::fabs(trk->pt() - trk_->pt()) > 1.e-6 )
//       edm::LogError("NoMatch") << "Closest reco::Track" << std::endl;
     if ( (*link1).isNonnull() && link1_.isNonnull() && 
	  std::fabs((*link1)->pt() - link1_->pt()) > 1.e-6 )
       edm::LogError("NoMatch") << "edm::Ptr for PackedCands" << std::endl;
     if ( (*link2).isNonnull() && link2_.isNonnull() && 
	  std::fabs((*link2)->pt() - link2_->pt()) > 1.e-6 )
       edm::LogError("NoMatch") << "edm::Ptr for LostTracks" << std::endl;
     if ( id != id_ ) 
       edm::LogError("NoMatch") << "ID BDT score" << std::endl;
     if ( unbiased != unbiased_ ) 
       edm::LogError("NoMatch") << "Seed unbiased BDT" << std::endl;
     if ( ptbiased != ptbiased_ ) 
       edm::LogError("NoMatch") << "Seed ptbiased BDT" << std::endl;

  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATLowPtElectrons_UnitTest);
