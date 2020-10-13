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
   const edm::EDGetTokenT< edm::ValueMap<float> > mvaId_;
   const edm::EDGetTokenT< edm::ValueMap<float> > mvaUnbiased_;
   const edm::EDGetTokenT< edm::ValueMap<float> > mvaPtbiased_;
   const edm::InputTag lowPtGsfLinksTag_;
   edm::EDGetTokenT< edm::ValueMap<edm::Ptr<reco::Candidate> > > lowPtGsfLinks_;
   bool verbose_;

 };

 PATLowPtElectrons_UnitTest::PATLowPtElectrons_UnitTest( const edm::ParameterSet& cfg ) :
   patElectrons_(consumes<edm::View<pat::Electron> >(cfg.getParameter<edm::InputTag>("patElectrons"))),
   gsfElectrons_(consumes<edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
   mvaId_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaId"))),
   mvaUnbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaUnbiased"))),
   mvaPtbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaPtbiased"))),
   lowPtGsfLinksTag_(cfg.getParameter<edm::InputTag>("lowPtGsfLinks")),
   verbose_(cfg.getParameter<bool>("verbose"))
 {
   lowPtGsfLinks_ = consumes< edm::ValueMap<edm::Ptr<reco::Candidate> > >(lowPtGsfLinksTag_);
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

   edm::Handle<edm::ValueMap<float> > mvaId;
   event.getByToken(mvaId_, mvaId);
   if ( mvaId.failedToGet() ) edm::LogWarning("mvaId.failedToGet()");
   if ( !mvaId.isValid() ) edm::LogWarning("!mvaId.isValid()");

   edm::Handle< edm::ValueMap<float> > mvaUnbiased;
   event.getByToken(mvaUnbiased_, mvaUnbiased);
   if ( mvaUnbiased.failedToGet() ) edm::LogWarning("mvaUnbiased.failedToGet()");
   if ( !mvaUnbiased.isValid() ) edm::LogWarning("!mvaUnbiased.isValid()");

   edm::Handle< edm::ValueMap<float> > mvaPtbiased;
   event.getByToken(mvaPtbiased_, mvaPtbiased);
   if ( mvaPtbiased.failedToGet() ) edm::LogWarning("mvaPtbiased.failedToGet()");
   if ( !mvaPtbiased.isValid() ) edm::LogWarning("!mvaPtbiased.isValid()");

   edm::Handle< edm::ValueMap<edm::Ptr<reco::Candidate> > > lowPtGsfLinks;
   event.getByToken(lowPtGsfLinks_, lowPtGsfLinks);
   if ( lowPtGsfLinks.failedToGet() ) edm::LogWarning("lowPtGsfLinks.failedToGet()");
   if ( !lowPtGsfLinks.isValid() ) edm::LogWarning("!lowPtGsfLinks.isValid()");

   if ( verbose_ ) { 
     std::cout << "patElectrons: " << patElectrons->size()
	       << " gsfElectrons: " << gsfElectrons->size()
	       << std::endl;
   }

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
     std::string label = lowPtGsfLinksTag_.label();
     edm::Ptr<reco::Candidate> link = pat->hasUserCand(label) ? pat->userCand(label) : edm::Ptr<reco::Candidate>();
     float id = ( pat->hasUserFloat("ID") ? pat->userFloat("ID") : -100. );
     float unbiased = ( pat->isElectronIDAvailable("unbiased") ? pat->electronID("unbiased") : -100. );
     float ptbiased = ( pat->isElectronIDAvailable("ptbiased") ? pat->electronID("ptbiased") : -100. );
     const pat::PFIsolation& iso = pat->miniPFIsolation();
     
     // PRINT OUT
     if ( verbose_ ) 
       std::cout << "index: "<<i<<std::endl
		 << " PAT:  "<<" null: "<<pat.isNull() <<" id: "<<pat.id() <<" key: "<<pat.key()<<" "<<pat->pt()<<std::endl
		 << " CORE: "<<" null: "<<core.isNull()<<" id: "<<core.id()<<" key: "<<core.key()<<std::endl
		 << " GSF:  "<<" null: "<<gsf.isNull() <<" id: "<<gsf.id() <<" key: "<<gsf.key()<<" "<<gsf->pt()<<std::endl
		 << " SC:   "<<" null: "<<sc.isNull()  <<" id: "<<sc.id()  <<" key: "<<sc.key()<<" "<<sc->energy()<<std::endl
		 << " SEED: "<<" null: "<<seed.isNull()<<" id: "<<seed.id()<<" key: "<<seed.key()<<" "<<seed->energy()<<std::endl
		 << " TRK:  "<<" null: "<<trk.isNull() <<" id: "<<trk.id() <<" key: "<<trk.key()<<std::endl
		 << " LINK: " 
		 << " null: "<<link.isNull() 
		 << " id:   "<<link.id()  
		 << " key:  "<<link.key() 
		 << " name: "<<label<<std::endl
		 << " ID:   "<<" val:  "<<id<<std::endl
		 << " UNB:  "<<" val:  "<<unbiased<<std::endl
		 << " PTB:  "<<" val:  "<<ptbiased<<std::endl
		 << " ISO:  "
		 <<" ch: "<<iso.chargedHadronIso()
		 <<" nh: "<<iso.neutralHadronIso()
		 <<" ph: "<<iso.photonIso()
		 <<" pu: "<<iso.puChargedHadronIso()
		 << std::endl;

     // IF CORE IS NOT EMBEDDED, ATTEMPT TO MATCH PAT TO GSF ELECTRON (VIA CORE)
     auto matched_ele = std::find_if( gsfElectrons->begin(), 
				      gsfElectrons->end(), 
				      [core]( const reco::GsfElectron& ele ) { 
					return ele.core() == core;
				      }
				      );
     if ( matched_ele == gsfElectrons->end() ) {
       edm::LogWarning("No PAT-to-GSF match")
	 << " No match from pat::Electron to reco::GsfElectron via reco::GsfElectronCore!"
	 << std::endl
	 << " pat.isNonnull():" << pat.isNonnull()
	 << " pat->core().isNonnull():" << pat->core().isNonnull()
	 << " pat.id(): " << pat.id()
	 << " pat.key(): " << pat.key()
	 << " pat->pt(): " << pat->pt()
	 << " pat->eta(): " << pat->eta()
	 << std::endl;
       continue;
     }

     // GET ELECTRON CORE, GSF TRACK, SUPERCLUSTER, SEED CALOCLUSTER, (CLOSEST) TRACK
     edm::Ref<edm::View<reco::GsfElectron> > ele(gsfElectrons,std::distance(gsfElectrons->begin(),matched_ele));
     edm::Ref<std::vector<reco::GsfElectronCore> > core_ = ele->core();
     edm::Ref<std::vector<reco::GsfTrack> > gsf_ = ele->gsfTrack();
     edm::Ref<std::vector<reco::SuperCluster> > sc_ = ele->superCluster();
     edm::Ptr<reco::CaloCluster> seed_ = ele->electronCluster(); //@@ pat->seed();
     edm::Ref<std::vector<reco::Track> > trk_ = ele->closestCtfTrackRef();
     edm::Ptr<reco::Candidate> link_ = lowPtGsfLinks->get( ele.key() );
     float id_ = mvaId->get( ele.key() );
     float unbiased_ = mvaUnbiased->get( ele->gsfTrack().key() );
     float ptbiased_ = mvaPtbiased->get( ele->gsfTrack().key() );

     // PRINT OUT
     if ( verbose_ ) 
       std::cout << " ELE:  "<<" null: "<<ele.isNull()  <<" id: "<<ele.id()  <<" key: "<<ele.key()<<std::endl
		 << " CORE: "<<" null: "<<core_.isNull()<<" id: "<<core_.id()<<" key: "<<core_.key()<<std::endl
		 << " GSF:  "<<" null: "<<gsf_.isNull() <<" id: "<<gsf_.id() <<" key: "<<gsf_.key()<<std::endl
		 << " SC:   "<<" null: "<<sc_.isNull()  <<" id: "<<sc_.id()  <<" key: "<<sc_.key()<<std::endl
		 << " SEED: "<<" null: "<<seed_.isNull()<<" id: "<<seed_.id()<<" key: "<<seed_.key()<<std::endl
		 << " TRK:  "<<" null: "<<trk_.isNull() <<" id: "<<trk_.id() <<" key: "<<trk_.key()<<std::endl
		 << " LINK: " 
		 << " null: "<<link_.isNull() 
		 << " id:   "<<link_.id()  
		 << " key:  "<<link_.key()<<std::endl
		 << " ID:   "<<" val:  "<<id_<<std::endl
		 << " UNB:  "<<" val:  "<<unbiased_<<std::endl
		 << " PTB:  "<<" val:  "<<ptbiased_<<std::endl
		 << std::endl;
     
     // CHECKS
     if ( !(core == core_ || std::fabs(core->gsfTrack()->pt() == core_->gsfTrack()->pt()) < 1.e-6) )
       edm::LogWarning("NoMatch") << "reco::GsfElectronCore" << std::endl;
     if ( !(gsf == gsf_ || std::fabs(gsf->pt() == gsf_->pt()) < 1.e-6) )
       edm::LogWarning("NoMatch") << "reco::GsfTrack" << std::endl;
     if ( !(sc == sc_ || std::fabs(sc->energy() == sc_->energy()) < 1.e-6) )
       edm::LogWarning("NoMatch") << "reco::SuperCluster" << std::endl;
     if ( !(seed == seed_ || //@@ NOT THE SAME!
	    std::fabs(seed->energy() == seed_->energy()) < 1.e-6) )
       edm::LogWarning("NoMatch") << "Seed reco::CaloCluster" << std::endl;
     if ( !(trk == trk_ || std::fabs(trk->pt() == trk_->pt()) < 1.e-6) )
       edm::LogWarning("NoMatch") << "Closest reco::Track" << std::endl;
     if ( !(link == link_ || std::fabs(link->pt() == link_->pt()) < 1.e-6) )
       edm::LogWarning("NoMatch") << "Linked reco::Track" << std::endl;
     //if ( id != id_ ) edm::LogWarning("NoMatch") << "ID BDT score" << std::endl;
     if ( unbiased != unbiased_ ) edm::LogWarning("NoMatch") << "Seed unbiased BDT" << std::endl;
     if ( ptbiased != ptbiased_ ) edm::LogWarning("NoMatch") << "Seed ptbiased BDT" << std::endl;

  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATLowPtElectrons_UnitTest);
