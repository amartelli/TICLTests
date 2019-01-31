#include "RecoHGCal/TICLTests/test/TICLAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Math/interface/deltaR.h"


using namespace std;
using namespace edm;
using namespace reco;

TICLAnalyzer::TICLAnalyzer(const edm::ParameterSet& iConfig) :
  genToken_(consumes<std::vector<CaloParticle> >(edm::InputTag("mix:MergedCaloTruth"))),
  simclusToken_(consumes<std::vector<SimCluster> >(edm::InputTag("mix:MergedCaloTruth"))),
  tkToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcTracks::RECO2"))),
  mcMIPToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiClusterMIP:MIPMultiClustersFromTracksterByCA"))),
  mcToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiCluster:MultiClustersFromTracksterByCA"))),
  hits_eeToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"))),
  hits_fhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"))),
  hits_bhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits")))
{
  //book some histos here                                                                                                                 
  //histos_["dr"]   = new TH1F("dr",   ";#Delta R;", 50,0,0.5);                                                                           

  histos_["hitsFraction_tk"]   = new TH1F("hitsFraction_tk",   ";fraction of gen-matched hits in hgcTracks;", 101,0,1.1);
  histos_["hitsFraction_mcMIP"]   = new TH1F("hitsFraction_mcMIP",   ";fraction of gen-matched hits in trackstersMIP;", 101,0,1.1);
  histos_["hitsFraction_mc"]   = new TH1F("hitsFraction_mc",   ";fraction of gen-matched hits in tracksters2DCl;", 101,0,1.1);
  histos_["hitsFraction_all"]   = new TH1F("hitsFraction_all",   ";fraction of gen-matched hits in all reco;", 101,0,1.1);

  histos_["EFraction_tk"]   = new TH1F("EFraction_tk",   ";fraction of gen-matched energy in hgcTracks;", 101,0,1.1);
  histos_["EFraction_mcMIP"]   = new TH1F("EFraction_mcMIP",   ";fraction of gen-matched energy in trackstersMIP;", 101,0,1.1);
  histos_["EFraction_mc"]   = new TH1F("EFraction_mc",   ";fraction of gen-matched energy in tracksters2DCl;", 101,0,1.1);
  histos_["EFraction_all"]   = new TH1F("EFraction_all",   ";fraction of gen-matched energy in all reco;", 101,0,1.1);
}



TICLAnalyzer::~TICLAnalyzer() { 
  TFile *outF = TFile::Open("ticl_analysis.root","RECREATE");
  for(std::map<TString,TH1 *>::iterator it=histos_.begin();
      it!=histos_.end();
      it++)
    it->second->Write();
  outF->Close();

}



void TICLAnalyzer::beginRun(const edm::Run& run, 
                            const edm::EventSetup & es) { }


void  TICLAnalyzer::analyze(const Event& iEvent, 
                            const EventSetup& iSetup) {

  edm::Handle<std::vector<CaloParticle> > gpH;
  iEvent.getByToken( genToken_, gpH);

  edm::Handle<std::vector<SimCluster> > scH;
  iEvent.getByToken( simclusToken_, scH);

  edm::Handle<std::vector<reco::CaloCluster> > tkH;
  iEvent.getByToken( tkToken_, tkH);

  edm::Handle<std::vector<reco::HGCalMultiCluster> > mcMIPH,mcH;
  iEvent.getByToken( mcMIPToken_,mcMIPH);
  iEvent.getByToken( mcToken_,mcH);


  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;
  iEvent.getByToken(hits_eeToken_,ee_hits);
  iEvent.getByToken(hits_fhToken_,fh_hits);
  iEvent.getByToken(hits_bhToken_,bh_hits);

  std::map<DetId, const HGCRecHit*> hitmap;
  for(auto const& it: *ee_hits) hitmap[it.detid()] = &it;
  for(auto const& it: *fh_hits) hitmap[it.detid()] = &it;
  for(auto const& it: *bh_hits) hitmap[it.detid()] = &it;
  
  for(auto cp : *gpH) {
    
    float eta=cp.eta();
    float cpEnergy = cp.energy();

    std::vector<std::pair<uint32_t,float> > allHits;
    std::vector<uint32_t> allMatched_tk;
    std::vector<uint32_t> allMatched_mcMIP,allMatched_mc;
    std::set<uint32_t> allMatched_global;

    float allMatched_tkE = 0.;
    float allMatched_mcMIPE = 0.;
    float allMatched_mcE = 0.;
    float allMatched_globalE = 0.;

    //iterate over all the attached sim clusters
    for(CaloParticle::sc_iterator scIt=cp.simCluster_begin();
        scIt!=cp.simCluster_end();
        scIt++) {

      //all hits and energy fractions at sim level
      std::vector<std::pair<uint32_t,float> > hits=scH->at( scIt->key() ).hits_and_fractions();
      allHits.insert(allHits.end(),hits.begin(),hits.end());

      //check which ones are matched by HGC tracking
      std::vector<uint32_t> matched_tk = getTrackedHitsList(eta>0,hits,*tkH);
      allMatched_tk.insert(allMatched_tk.end(),matched_tk.begin(),matched_tk.end());

      //check which ones are matched by the multicluster algorithmb
      std::vector<uint32_t> matched_mcMIP=getClusteredHitsList(eta>0,hits,*mcMIPH);
      allMatched_mcMIP.insert(allMatched_mcMIP.end(),matched_mcMIP.begin(),matched_mcMIP.end());
      std::vector<uint32_t> matched_mc=getClusteredHitsList(eta>0,hits,*mcH);
      allMatched_mc.insert(allMatched_mc.end(),matched_mc.begin(),matched_mc.end());

      for(auto ij : allMatched_tk) {
	//in principle every reco id exists in the collection of recHits                                                                 
	const HGCRecHit *hit = hitmap[ij];
	unsigned int preSize = allMatched_global.size();
	allMatched_global.insert(ij);
	if(allMatched_global.size() != preSize) allMatched_tkE += hit->energy();
      }
      for(auto ij : allMatched_mcMIP){
	const HGCRecHit *hit = hitmap[ij];
	unsigned int preSize = allMatched_global.size();
	allMatched_global.insert(ij);
	if(allMatched_global.size() != preSize) allMatched_mcMIPE += hit->energy();
      }
      for(auto ij : allMatched_mc){
	const HGCRecHit *hit = hitmap[ij];
	unsigned int preSize = allMatched_global.size();
	allMatched_global.insert(ij);
	if(allMatched_global.size() != preSize) allMatched_mcE += hit->energy();
      }
      allMatched_globalE = allMatched_tkE + allMatched_mcMIPE + allMatched_mcE;
    }
    std::cout << cp.pdgId() << " nSimClusters = " << cp.simClusters().size()
	      << " sim " << allHits.size()
	      << " tk " << allMatched_tk.size()
	      << " mcMIP " << allMatched_mcMIP.size() << " mc2D " << allMatched_mc.size() << " global = " << allMatched_global.size()
	      << " cpE = " << cpEnergy << " tkE = " << allMatched_tkE << " mcMIPE = " << allMatched_mcMIPE
	      << " mc2DE = " << allMatched_mcE << " globalE = " << allMatched_globalE << std::endl;

    if(cp.simClusters().size() > 1) continue;
    histos_["hitsFraction_tk"]->Fill(1.*allMatched_tk.size()/allHits.size());
    histos_["hitsFraction_mcMIP"]->Fill(1.*allMatched_mcMIP.size()/allHits.size());
    histos_["hitsFraction_mc"]->Fill(1.*allMatched_mc.size()/allHits.size());
    histos_["hitsFraction_all"]->Fill(1.*allMatched_global.size()/allHits.size());

    histos_["EFraction_tk"]->Fill(1.*allMatched_tkE/cpEnergy);
    histos_["EFraction_mcMIP"]->Fill(1.*allMatched_mcMIPE/cpEnergy);
    histos_["EFraction_mc"]->Fill(1.*allMatched_mcE/cpEnergy);
    histos_["EFraction_all"]->Fill(1.*allMatched_globalE/cpEnergy);
    if(allMatched_globalE/cpEnergy > 1.) std::cout << " problem " << std::endl;

  }
}

//
std::vector<uint32_t> TICLAnalyzer::getClusteredHitsList(bool pos,
                                                         const std::vector<std::pair<uint32_t,float> > &hits,
                                                         const std::vector<reco::HGCalMultiCluster> &mcs) {
  
  std::vector<uint32_t> matchedList;
  
  for(auto mc : mcs) {
    
    //require on the same side
    bool mcPos(mc.eta()>0);
    if(mcPos!=pos) continue;
    
    //loop over all the layer clusters in the multicluster
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){
      const std::vector< std::pair<DetId, float> > &recHits = (*it)->hitsAndFractions();
      std::vector<uint32_t> imatches=getMatched(hits,recHits);
      matchedList.insert(matchedList.end(), imatches.begin(), imatches.end());
    }
    
  }

  return matchedList;
}


//
std::vector<uint32_t> TICLAnalyzer::getTrackedHitsList(bool pos,
                                                       const std::vector<std::pair<uint32_t,float> >  &hits,
                                                       const std::vector<reco::CaloCluster> &ccs) {
  
  std::vector<uint32_t> matchedList;
  
  for(auto cc : ccs) {
    
    //require on the same side
    bool ccPos(cc.eta()>0);
    if(ccPos!=pos) continue;
    
    const std::vector< std::pair<DetId, float> > &recHits =cc.hitsAndFractions();
    std::vector<uint32_t> imatches=getMatched(hits,recHits);
    matchedList.insert(matchedList.end(), imatches.begin(), imatches.end());
    
  }
  
  return matchedList;  
}

//
std::vector<uint32_t> TICLAnalyzer::getMatched(const std::vector<std::pair<uint32_t,float> > &a,
                                               const std::vector<std::pair<DetId,float> > &b){

  std::vector<uint32_t> matchedList;
  for(size_t i=0; i<a.size(); i++) {

    //find first match in second list
    for(size_t j=0; j<b.size(); j++) {
      if(a[i].first!=b[j].first.rawId()) continue;
      matchedList.push_back(a[i].first);
      break;
    }
  }

  return matchedList;
}





DEFINE_FWK_MODULE(TICLAnalyzer);
