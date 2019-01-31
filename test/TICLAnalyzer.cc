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
  eetkToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcTracks:EE"))),
  fhtkToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcTracks:FH"))),
  mcMIPToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiClusterMIP:MIPMultiClustersFromTracksterByCA"))),
  mcToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiCluster:MultiClustersFromTracksterByCA")))
{
  //book some histos here
  //histos_["dr"]   = new TH1F("dr",   ";#Delta R;", 50,0,0.5);
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

  edm::Handle<std::vector<reco::CaloCluster> > eetkH,fhtkH;
  iEvent.getByToken( eetkToken_, eetkH);
  iEvent.getByToken( fhtkToken_, fhtkH);

  edm::Handle<std::vector<reco::HGCalMultiCluster> > mcMIPH,mcH;
  iEvent.getByToken( mcMIPToken_,mcMIPH);
  iEvent.getByToken( mcToken_,mcH);

  
  for(auto cp : *gpH) {
    
    float eta=cp.eta();

    std::vector<std::pair<uint32_t,float> > allHits;
    std::vector<uint32_t> allMatched_eetk, allMatched_fhtk;
    std::vector<uint32_t> allMatched_mcMIP,allMatched_mc;

    //iterate over all the attached sim clusters
    for(CaloParticle::sc_iterator scIt=cp.simCluster_begin();
        scIt!=cp.simCluster_end();
        scIt++) {

      //all hits and energy fractions at sim level
      std::vector<std::pair<uint32_t,float> > hits=scH->at( scIt->key() ).hits_and_fractions();
      allHits.insert(allHits.end(),hits.begin(),hits.end());

      //check which ones are matched by HGC tracking
      std::vector<uint32_t> matched_eetk = getTrackedHitsList(eta>0,hits,*eetkH);
      allMatched_eetk.insert(allMatched_eetk.end(),matched_eetk.begin(),matched_eetk.end());
      std::vector<uint32_t> matched_fhtk = getTrackedHitsList(eta>0,hits,*eetkH);
      allMatched_fhtk.insert(allMatched_fhtk.end(),matched_fhtk.begin(),matched_fhtk.end());

      //check which ones are matched by the multicluster algorithmb
      std::vector<uint32_t> matched_mcMIP=getClusteredHitsList(eta>0,hits,*mcMIPH);
      allMatched_mcMIP.insert(allMatched_mcMIP.end(),matched_mcMIP.begin(),matched_mcMIP.end());
      std::vector<uint32_t> matched_mc=getClusteredHitsList(eta>0,hits,*mcH);
      allMatched_mc.insert(allMatched_mc.end(),matched_mc.begin(),matched_mc.end());
    }


    std::cout << cp.pdgId() << " " << allHits.size() 
              << " " << allMatched_eetk.size() << " " << allMatched_fhtk.size()
              << " " << allMatched_mcMIP.size() << " " << allMatched_mc.size() << std::endl;
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
