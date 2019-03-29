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

  histos_["hitsRecoOKFraction_tk"]   = new TH1F("hitsRecoOKFraction_tk",   ";gen-matched hits over all reco in hgcTracks;", 101,0,1.1);
  histos_["hitsRecoOKFraction_mcMIP"]   = new TH1F("hitsRecoOKFraction_mcMIP",   ";gen-matched hits over all reco in trackstersMIP;", 101,0,1.1);
  histos_["hitsRecoOKFraction_mc"]   = new TH1F("hitsRecoOKFraction_mc",   ";gen-matched hits over all reco in tracksters2DCl;", 101,0,1.1);
  histos_["hitsRecoOKFraction_all"]   = new TH1F("hitsRecoOKFraction_all",   ";gen-matched hits over all reco", 101,0,1.1);


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


  std::cout << " \n \n new evt "<< std::endl;

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

  std::map<uint32_t, const HGCRecHit*> hitmap;
  for(auto const& it: *ee_hits) hitmap[it.detid().rawId()] = &it;
  for(auto const& it: *fh_hits) hitmap[it.detid().rawId()] = &it;
  for(auto const& it: *bh_hits) hitmap[it.detid().rawId()] = &it;
  
  std::cout << " loop over CP " << std::endl;
  for(auto cp : *gpH) {
    //   if(( cp.pdgId() != 22 && cp.simClusters().size() > 1) || (cp.simClusters().size() > 2)) continue;
    
    float eta = cp.eta();
    //float cpEnergy = cp.energy();
    float cpEnergy = 0;

    std::set<std::pair<uint32_t,float> > allHits;
    std::map<uint32_t,float> allHitsM;
    std::set<uint32_t> all_tk;
    std::set<uint32_t> all_mcMIP, all_mc;
    std::set<uint32_t> all_global;

    std::set<uint32_t> allMatched_tk;
    std::set<uint32_t> allMatched_mcMIP,allMatched_mc;
    std::set<uint32_t> allMatched_global;

    float allMatched_tkE = 0.;
    float allMatched_mcMIPE = 0.;
    float allMatched_mcE = 0.;
    float allMatched_globalE = 0.;

    std::cout << " loop over CP simclusters simCluster.size() = " << cp.simClusters().size() << std::endl;
    //iterate over all the attached sim clusters
    for(CaloParticle::sc_iterator scIt=cp.simCluster_begin(); scIt!=cp.simCluster_end(); scIt++){

      //all hits and energy fractions at sim level
      std::vector<std::pair<uint32_t,float> > hits = scH->at( scIt->key() ).hits_and_fractions();
      //allHits.insert(allHits.end(),hits.begin(),hits.end());
      std::cout << " >>> hits.size() = " << hits.size() << std::endl;
      for(auto ij : hits) {
	auto finder = allHitsM.find(ij.first);
	if(finder != allHitsM.end()) std::cout << " previous fraction = " << allHitsM[ij.first] << " new = " << ij.second << std::endl;
        allHitsM[ij.first] = ij.second;
	allHits.insert(ij);

	auto finder2 = hitmap.find(ij.first);
	if(finder2 != hitmap.end()){
	  const HGCRecHit *hit = hitmap[ij.first];
	  cpEnergy += hit->energy() * ij.second;
	}
	//if(ij.second != 1) std::cout << " >>> fraction sim != 1 =  " << ij.second  << std::endl;
      }
    }// calo cluster loop

    std::cout << " now match singles " << std::endl;
    //check which ones are matched by HGC tracking
    std::set<uint32_t> matched_tk = getTrackedHitsList(eta>0, allHits, *tkH, all_tk);
        
    //check which ones are matched by the multicluster algorithmb
    std::set<uint32_t> matched_mcMIP = getClusteredHitsList(eta>0, allHits, *mcMIPH, all_mcMIP);
    
    std::set<uint32_t> matched_mc = getClusteredHitsList(eta>0, allHits, *mcH, all_mc);
    
    std::cout << " post matching conta E " << std::endl;
    for(auto ij : matched_tk) {
      allMatched_tk.insert(ij);
      //in principle every reco id exists in the collection of recHits                                                                 
      const HGCRecHit *hit = hitmap[ij];
      unsigned int preSize = allMatched_global.size();
      allMatched_global.insert(ij);
      if(allMatched_global.size() != preSize){
	float fraction = 1.;
	for(auto kk : allHits){
	  if(kk.first == ij) fraction = kk.second;
	  break;
	}
	if(fraction != 1) std::cout << " >> trks fraction " << fraction << std::endl;
	allMatched_tkE += hit->energy() * fraction;
      }
    }
    for(auto ij : matched_mcMIP){
      allMatched_mcMIP.insert(ij);
      const HGCRecHit *hit = hitmap[ij];
      unsigned int preSize = allMatched_global.size();
      allMatched_global.insert(ij);
      if(allMatched_global.size() != preSize){
	float fraction = 1.;
	for(auto kk : allHits){
	  if(kk.first == ij) fraction = kk.second;
	  break;
	}
	if(fraction != 1) std::cout << " >> MCMIP fraction " << fraction << std::endl;
	allMatched_mcMIPE += hit->energy()*fraction;
      }
    }
    for(auto ij : matched_mc){
      allMatched_mc.insert(ij);
      const HGCRecHit *hit = hitmap[ij];
      unsigned int preSize = allMatched_global.size();
      allMatched_global.insert(ij);
      if(allMatched_global.size() != preSize){
	float fraction = 1.;
	for(auto kk : allHits){
	  if(kk.first == ij) fraction = kk.second;
	  break;
	}
	if(fraction != 1) std::cout << " >> MC2d fraction " << fraction << std::endl;
	allMatched_mcE += hit->energy()*fraction;
      }
    }
    
    allMatched_globalE = allMatched_tkE + allMatched_mcMIPE + allMatched_mcE;
    
    std::cout << cp.pdgId() << " nSimClusters = " << cp.simClusters().size()
	      << " sim " << allHits.size()
	      << " tk " << allMatched_tk.size()
	      << " mcMIP " << allMatched_mcMIP.size() << " mc2D " << allMatched_mc.size() << " global = " << allMatched_global.size()
	      << " cpE = " << cpEnergy << " tkE = " << allMatched_tkE << " mcMIPE = " << allMatched_mcMIPE
	      << " mc2DE = " << allMatched_mcE << " globalE = " << allMatched_globalE << std::endl;
    
    std::cout << " tk matched size = " << allMatched_tk.size() << " all = " << all_tk.size() << std::endl; 
    std::cout << " mc matched size = " << allMatched_mc.size() << " all = " << all_mc.size() << std::endl; 
    std::cout << " mc matchedMIP size = " << allMatched_mcMIP.size() << " all = " << all_mcMIP.size() << std::endl; 

    if(allMatched_mc.size() > all_mc.size()) std::cout << " problem ratio 1 " << std::endl;

    if(allHits.size() != 0){
    histos_["hitsFraction_tk"]->Fill(1.*allMatched_tk.size()/allHits.size());
    histos_["hitsFraction_mcMIP"]->Fill(1.*allMatched_mcMIP.size()/allHits.size());
    histos_["hitsFraction_mc"]->Fill(1.*allMatched_mc.size()/allHits.size());
    histos_["hitsFraction_all"]->Fill(1.*allMatched_global.size()/allHits.size());
    if(1.*allMatched_mcMIP.size()/allHits.size() > 0.9) std::cout << " problem frac" <<  1.*allMatched_mcMIP.size()/allHits.size() << std::endl;
    }

    if(all_tk.size() != 0) histos_["hitsRecoOKFraction_tk"]->Fill(1.*allMatched_tk.size() / all_tk.size());
    if(all_mcMIP.size() != 0) histos_["hitsRecoOKFraction_mcMIP"]->Fill(1.*allMatched_mcMIP.size() / all_mcMIP.size());
    if(all_mc.size() != 0) histos_["hitsRecoOKFraction_mc"]->Fill(1.*allMatched_mc.size() / all_mc.size());

    histos_["EFraction_tk"]->Fill(1.*allMatched_tkE / cpEnergy);
    histos_["EFraction_mcMIP"]->Fill(1.*allMatched_mcMIPE / cpEnergy);
    histos_["EFraction_mc"]->Fill(1.*allMatched_mcE / cpEnergy);
    histos_["EFraction_all"]->Fill(1.*allMatched_globalE / cpEnergy);
    if(allMatched_globalE/cpEnergy > 1.) std::cout << " problem " << std::endl;

  }
}


  std::set<uint32_t> TICLAnalyzer::getClusteredHitsList(bool pos,
							const std::set<std::pair<uint32_t,float> > &hits,
						      const std::vector<reco::HGCalMultiCluster> &mcs, 
						      std::set<uint32_t> &mcHits) {
  
  std::set<uint32_t> matchedList;
  //  std::set<uint32_t> testRh;
  
  for(auto mc : mcs) {
    
    //require on the same side
    bool mcPos(mc.eta()>0);
    if(mcPos!=pos) continue;
    
    //loop over all the layer clusters in the multicluster
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){
      const std::vector< std::pair<DetId, float> > &recHits = (*it)->hitsAndFractions();
      for(auto ij : recHits)
	if(ij.second != 0) mcHits.insert(ij.first.rawId());

      //std::vector<uint32_t> imatches=getMatched(hits,recHits);
      //matchedList.insert(matchedList.end(), imatches.begin(), imatches.end());
      std::set<uint32_t> imatches=getMatched(hits, mcHits);
      for(auto ij : imatches) matchedList.insert(ij);

      // for(auto ij : imatches){
      //  	unsigned int testSize = testRh.size();
      // 	testRh.insert(ij);
      //  	//if(testSize == testRh.size()) std::cout << " duplicated it = " << ij << std::endl;
      //  	if(testSize != testRh.size()) matchedList.push_back(ij);
      // }
    }
    
  }

  return matchedList;
}


//
std::set<uint32_t> TICLAnalyzer::getTrackedHitsList(bool pos,
						    const std::set<std::pair<uint32_t,float> >  &hits,
						    const std::vector<reco::CaloCluster> &ccs, 
						    std::set<uint32_t> &tkHits) {
  
  std::set<uint32_t> matchedList;
  
  for(auto cc : ccs) {
    
    //require on the same side
    bool ccPos(cc.eta()>0);
    if(ccPos!=pos) continue;
    
    const std::vector< std::pair<DetId, float> > &recHits =cc.hitsAndFractions();
    for(auto ij : recHits)
      if(ij.second != 0) tkHits.insert(ij.first.rawId());

    //std::vector<uint32_t> imatches=getMatched(hits,recHits);
    //matchedList.insert(matchedList.end(), imatches.begin(), imatches.end());
    std::set<uint32_t> imatches=getMatched(hits,tkHits);
    for(auto ij : imatches) matchedList.insert(ij);
  }
  
  return matchedList;  
}

//
std::vector<uint32_t> TICLAnalyzer::getMatched(const std::set<std::pair<uint32_t,float> > &a,
                                               const std::vector<std::pair<DetId,float> > &b){

  std::vector<uint32_t> matchedList;
  for(auto ii : a) {

    //find first match in second list
    for(size_t j=0; j<b.size(); j++) {
      if(ii.first != b[j].first.rawId()) continue;
      matchedList.push_back(ii.first);
      break;
    }
  }

  return matchedList;
}


std::set<uint32_t> TICLAnalyzer::getMatched(const std::set<std::pair<uint32_t,float> > &a,
					    const std::set<uint32_t> &b){

  std::set<uint32_t> matchedList;
  for(auto ii : a) {
    //find first match in second list
    for(auto jj : b) {
      if(ii.first != jj) continue;
      matchedList.insert(ii.first);
      break;
    }
  }

  return matchedList;
}



DEFINE_FWK_MODULE(TICLAnalyzer);
