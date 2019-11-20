#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <memory>
#include <cstdio>
#include <utility>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/RecoHGCal/interface/Trackster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TTree.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoParticleFlow/PFClusterProducer/plugins/SimMappers/ComputeClusterTime.h"


using namespace std;
using namespace edm;
using namespace reco;

class CheckResolution : public edm::EDAnalyzer {
public:

  explicit CheckResolution(const edm::ParameterSet&);

  ~CheckResolution();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

private:

  hgcal::RecHitTools         recHitTools;


  std::vector<reco::HGCalMultiCluster> getHighestEnergyMC(std::vector<reco::HGCalMultiCluster>& allmc, 
							  std::vector<reco::HGCalMultiCluster>& allmcMIP, 
							  std::vector<reco::HGCalMultiCluster>& all);

  std::vector<reco::HGCalMultiCluster> cleanTimeMC(const std::vector<reco::HGCalMultiCluster>& allmc);

  std::vector<reco::HGCalMultiCluster> buildTrackerItSameSeed(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs, 
							      std::vector<reco::HGCalMultiCluster>& allmcPhoton);

  void getClusteredHitsList(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs);
  void getTrackedHitsList(bool pos, const std::vector<reco::CaloCluster> &ccs);
  void getTrackedItHitsList(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs);

  void getTimeClusteredHitsList(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs, std::vector<float> &times);
  void getTimeTrackedHitsList(bool pos, const std::vector<reco::CaloCluster> &ccs, std::vector<float> &times);


  std::set<uint32_t> getMatchedHitsList(bool pos,
					const std::set<std::pair<uint32_t,float> >  &hits,
					const std::set<uint32_t>  &mcs,
					std::set<uint32_t> &recoHits);


  std::set<uint32_t> getMatchedClusteredHitsList(bool pos,
						 const std::set<std::pair<uint32_t,float> >  &hits,
						 const std::vector<reco::HGCalMultiCluster> &mcs,
						 std::set<uint32_t> &recoHits);

  std::set<uint32_t> getMatchedTrackedHitsList(bool pos,
					       const std::set<std::pair<uint32_t,float> >  &hits,
					       const std::vector<reco::CaloCluster> &ccs,
					       std::set<uint32_t> &recoHits);

  std::set<uint32_t> getMatched(const std::set<std::pair<uint32_t,float> > &a,
                                const std::set<uint32_t> &b);



  std::set<uint32_t> getMatched(const std::set<uint32_t> &a,
				const std::set<uint32_t> &b);


  std::set<uint32_t> getNOTMatched(const std::set<uint32_t> &a,
				   const std::set<uint32_t> &b);


  edm::EDGetTokenT<std::vector<CaloParticle> > genToken_;

  edm::EDGetTokenT<std::vector<reco::CaloCluster> > tkToken_;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > mcTrkToken_, mcMIPToken_,mcToken_;

  edm::EDGetTokenT<HGCRecHitCollection> hits_eeToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hits_fhToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hits_bhToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> timeMap_;
  edm::Handle<edm::ValueMap<float> > time2DMap;

  std::map<uint32_t, const HGCRecHit*> hitmap;
  
  TH1F* recoBP_genAll_RatioEnergy;
  TH1F* recoAll_genAll_RatioEnergy;
  TH1F* recoBPmatched_genAll_RatioEnergy;
  TH1F* recoAllmatched_genAll_RatioEnergy;

  TH1F* caloEnergyRatio;
  TH1F* dR_bestTr_gen;
  TH1F* deltaT;
  TH1F* tracksterT;

  bool debug;

  bool isPhoton;
  bool isTrkIteration;
  int nEvent;
};



CheckResolution::CheckResolution(const edm::ParameterSet& iConfig) :
  genToken_(consumes<std::vector<CaloParticle> >(edm::InputTag("mix:MergedCaloTruth"))),
  tkToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcTracks::RECO2"))),
  mcTrkToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiClusterTrk:TrkMultiClustersFromTracksterByCA"))),
  mcMIPToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("MultiClustersFromTrackstersMIP:MIPMultiClustersFromTracksterByCA"))),
  mcToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("MultiClustersFromTracksters:MultiClustersFromTracksterByCA"))),
  hits_eeToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"))),
  hits_fhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"))),
  hits_bhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"))),
  timeMap_(consumes<edm::ValueMap<float>>(edm::InputTag("hgcalLayerClusters:timeLayerCluster")))
{
  //book some histos here                                                                                                                 

  edm::Service<TFileService> fs;


  recoBP_genAll_RatioEnergy = fs->make<TH1F>("recoBP_genAll_RatioEnergy", "", 5000, 0., 5.);
  recoAll_genAll_RatioEnergy = fs->make<TH1F>("recoAll_genAll_RatioEnergy", "", 5000, 0., 5.);
  recoBPmatched_genAll_RatioEnergy = fs->make<TH1F>("recoBPmatched_genAll_RatioEnergy", "", 5000, 0., 5.);
  recoAllmatched_genAll_RatioEnergy = fs->make<TH1F>("recoAllmatched_genAll_RatioEnergy", "", 5000, 0., 5.);
  caloEnergyRatio = fs->make<TH1F>("caloEnergyRatio", "", 200, 0., 2.);
  dR_bestTr_gen = fs->make<TH1F>("dR_bestTr_gen", "", 100,0,10.);
  deltaT = fs->make<TH1F>("deltaT", "", 5000,-5.,5.);
  tracksterT = fs->make<TH1F>("tracksterT", "", 5000,-25.,25.);

  debug = false;
  isTrkIteration = false;
  isPhoton = true;
  nEvent = 0;
}



CheckResolution::~CheckResolution() { 
}



void CheckResolution::beginRun(const edm::Run& run, 
                            const edm::EventSetup & es) { }


void  CheckResolution::analyze(const Event& iEvent, 
                            const EventSetup& iSetup) {

  ++nEvent;
  recHitTools.getEventSetup(iSetup);

  if(debug)  std::cout << " \n \n new evt "<< std::endl;

  edm::Handle<std::vector<CaloParticle> > gpH;
  iEvent.getByToken( genToken_, gpH);

  iEvent.getByToken(timeMap_, time2DMap);
  if(debug) std::cout << " time2DMap->size() = " << time2DMap->size() << std::endl;

  edm::Handle<std::vector<reco::HGCalMultiCluster> > mcMIPH, mcH;
  iEvent.getByToken( mcMIPToken_,mcMIPH);
  iEvent.getByToken( mcToken_, mcH);

  //first clean tracksters
  std::vector<reco::HGCalMultiCluster> mcHClean = cleanTimeMC(*mcH);
  std::vector<reco::HGCalMultiCluster> mcMIPHClean = cleanTimeMC(*mcMIPH);


  //now select highest cluster for each side
  std::vector<reco::HGCalMultiCluster> allTracksters;
  std::vector<reco::HGCalMultiCluster> bestTrackster = getHighestEnergyMC(mcHClean, mcMIPHClean, allTracksters);

  if(debug) std::cout << " qui ci sono bestTrackster.size() = " << bestTrackster.size() 
		      << " allTracksters size = " << allTracksters.size() 
		      << std::endl;
  if(debug) std::cout << " mcH->size() = " << mcH->size() << " mcMIPH->size() = " << mcMIPH->size() << std::endl;
  if(debug) std::cout << " mcHClean.size() = " << mcHClean.size() << " mcMIPHClean.size() = " << mcMIPHClean.size() << std::endl;



  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;
  iEvent.getByToken(hits_eeToken_,ee_hits);
  iEvent.getByToken(hits_fhToken_,fh_hits);
  iEvent.getByToken(hits_bhToken_,bh_hits);


  //filling rechit map
  hitmap.clear();
  for(auto const& it: *ee_hits){
    hitmap[it.detid().rawId()] = &it;
  }
  for(auto const& it: *fh_hits){
    hitmap[it.detid().rawId()] = &it;
  }




  if(debug){
    std::cout << " caloParticle size = " << gpH->size() << std::endl;
    std::cout << " loop over CP " << std::endl;
  }


  for(auto cp : *gpH) {
    if(debug) std::cout << " >>> cp.eventId().event() = " << cp.eventId().event() 
			<< " cp.eventId().bunchCrossing() = " << cp.eventId().bunchCrossing() 
			<< " cp.pdgId() = " << cp.pdgId() << " cp.simClusters().size() = " << cp.simClusters().size() << std::endl;
    if(cp.pdgId() != 22 || cp.eventId().event() != 0 || cp.eventId().bunchCrossing() != 0) continue;
    if(cp.pdgId() != 22) continue;
    if(cp.simClusters().size() > 2) continue;
    

    float eta = cp.eta();
    float phi = cp.phi();
    float cpEnergy = 0;

    //all caloP hits
    std::set<uint32_t> all_simHits;
    std::set<std::pair<uint32_t,float> > allHits;

    //all hits in the single collection matched to caloP
    std::set<uint32_t> all_BP;
    std::set<uint32_t> all_mc;

    //all hits in finalP or single collection matched to caloP
    std::set<uint32_t> matched_all;
    std::set<uint32_t> matched_BP;

    //as matched_all but summing among the single collections
    float matchedAll_energy = 0.;
    float matchedBP_energy = 0.;

    //non necessary matched
    float BP_energy = 0.;
    float all_energy = 0.;

    //iterate over all the attached sim clusters
    float genRecoenergy = 0.;

    for(CaloParticle::sc_iterator scIt=cp.simCluster_begin(); scIt!=cp.simCluster_end(); scIt++){

      //all hits and energy fractions at sim level
      std::vector<std::pair<uint32_t,float> > hits = (*scIt)->hits_and_fractions();

      if(debug)  std::cout << " >>> hits.size() = " << hits.size() << std::endl;
      for(auto ij : hits) {

	auto finder2 = hitmap.find(ij.first);
	if(finder2 != hitmap.end()){

	  all_simHits.insert(ij.first);

	  genRecoenergy += ij.second; // * hit->energy();
	} 
      }
    }//sim clusters iterator
  
    caloEnergyRatio->Fill(genRecoenergy / cp.energy());
    
    
    for(auto ij : all_simHits){
      allHits.insert(std::pair<uint32_t,float>(ij, 1));
      auto finder2 = hitmap.find(ij);
      if(finder2 != hitmap.end()){
        const HGCRecHit *hit = hitmap[ij];
	cpEnergy += hit->energy();
      }
    }

    if(debug) std::cout << " cpEnergy = " << cpEnergy << std::endl;

    int recoSide = int(eta > 0);
    // if(recoSide > 0) std::cout << " recoSide = " << recoSide << " eta = " << eta << " bestTrackster[0].eta() = " << bestTrackster[0].eta() 
    // 			       << " bestTrackster[1].eta() = " << bestTrackster[1].eta() << std::endl;
    float dRloc = reco::deltaR(eta, phi, bestTrackster.at(recoSide).eta(), bestTrackster.at(recoSide).phi());

    if(dRloc > 3)std::cout << " eta = " << bestTrackster.at(recoSide).eta() << " cp eta = " << eta 
			<< " phi = " << bestTrackster.at(recoSide).phi() << " cp phi = " << phi << " dRloc = " << dRloc << std::endl;

    if(debug) std::cout << " eta = " << bestTrackster.at(recoSide).eta() << " cp eta = " << eta 
			<< " phi = " << bestTrackster.at(recoSide).phi() << " cp phi = " << phi << " dRloc = " << dRloc << std::endl;
    
    dR_bestTr_gen->Fill(dRloc);


    //all IDs in finalParticles that are matched to allHits = from caloP
    matched_BP = getMatchedClusteredHitsList(eta>0, allHits, bestTrackster, all_BP);
    matched_all = getMatchedClusteredHitsList(eta>0, allHits, allTracksters, all_mc);

    if(debug) std::cout << " now compute energies " << std::endl;

    //all matched to sim
    for(auto ij : matched_all){
      auto finder2 = hitmap.find(ij);
      if(finder2 == hitmap.end()) continue;
      const HGCRecHit *hit = hitmap[ij];
      matchedAll_energy += hit->energy();
    }
    
    //all 
    for(auto ij : all_mc){
      auto finder2 = hitmap.find(ij);
      if(finder2 == hitmap.end()) continue;
      const HGCRecHit *hit = hitmap[ij];
      all_energy += hit->energy();
    }

    //BP matched to sim
    for(auto ij : matched_BP){
      auto finder2 = hitmap.find(ij);
      if(finder2 == hitmap.end()) continue;
      const HGCRecHit *hit = hitmap[ij];
      matchedBP_energy += hit->energy();
    }
    //all 
    for(auto ij : all_BP){
      auto finder2 = hitmap.find(ij);
      if(finder2 == hitmap.end()) continue;
      const HGCRecHit *hit = hitmap[ij];
      BP_energy += hit->energy();
    }

    if(debug)    std::cout << " post matching conta E now all " << std::endl;

    if(debug) std::cout << " BP_energy = " << BP_energy << " matchedBP_energy = " << matchedBP_energy 
			<< " all_energy = " << all_energy << " matchedAll_energy = " << matchedAll_energy
			<< " cpEnergy = " << cpEnergy << std::endl;

    recoBP_genAll_RatioEnergy->Fill(BP_energy/cpEnergy);
    recoAll_genAll_RatioEnergy->Fill(all_energy/cpEnergy);

    recoBPmatched_genAll_RatioEnergy->Fill(matchedBP_energy/cpEnergy);
    recoAllmatched_genAll_RatioEnergy->Fill(matchedAll_energy/cpEnergy);


  }//caloParticles
}

std::vector<reco::HGCalMultiCluster> CheckResolution::cleanTimeMC(const std::vector<reco::HGCalMultiCluster>& allmc){

  if(debug) std::cout << " >>> in cleanTimeMC "  << std::endl;

  std::vector<reco::HGCalMultiCluster> result;

  for(auto mc : allmc) {

    if(debug)    std::cout << " thisMC energy = " << mc.energy() << " size = " << mc.size() << std::endl;

    std::vector<float> localTimeAll;
    std::vector<float> localTime;
    std::vector<float> localEnergy;

    localTimeAll.resize(mc.size());
    unsigned iCount=0;
    //loop over all the layer clusters in the multicluster to compute time
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){

      reco::CaloClusterPtr sClPtr(*it);
      float time2D = (*time2DMap)[sClPtr];
      localTimeAll[iCount] = time2D;
      if(time2D > -50.){
        localTime.push_back(time2D);
        localEnergy.push_back(sClPtr->energy());
      }
      ++iCount;
    }

    if(debug) std::cout << " now compute time localTime.size() = " << localTime.size() << std::endl;
    //use method from Shameena => try for the moment with truncation
    // can try truncation + weighted mean
    float finalT = (localTime.size() >= 3) ? (hgcalsimclustertime::fixSizeHighestDensity(localTime)) : (-99.);
    float finalTW = (localTime.size() >= 3) ? (hgcalsimclustertime::fixSizeHighestDensity(localTime, localEnergy)) : (-99.);

    if(debug) std::cout << " finalT = " << finalT << " finalTW = " << finalTW << std::endl;

    //prepare new multicluster
    reco::HGCalMultiCluster temp;
    double baricenter[3] = {0., 0., 0.};
    double total_weight = 0.;

    if(finalTW == -99) {
      result.push_back(mc);
      continue;
    }
    tracksterT->Fill(finalTW);

    //loop over all the layer clusters in the multicluster to clean wrt time
    for(unsigned ij=0; ij<localTimeAll.size(); ++ij){

      if(localTimeAll.at(ij) > -50.) {
	deltaT->Fill(localTimeAll.at(ij) - finalTW);
	// if(localTimeAll.at(ij) - finalTW > 5.) std::cout << " localTimeAll.at(ij) = " << localTimeAll.at(ij) << " finalTW = " << finalTW
	// 						 << " diff = " << localTimeAll.at(ij) - finalTW << std::endl;
      }
      
      if((std::abs(localTimeAll.at(ij) - finalTW) < 0.210 && localTimeAll.at(ij) > -50.) || (localTimeAll.at(ij) < -50.)){
	reco::HGCalMultiCluster::component_iterator it = mc.begin() + ij;
	reco::CaloClusterPtr sClPtr(*it);
	temp.push_back(sClPtr);
	auto weight = sClPtr->energy();
	total_weight += weight;
	baricenter[0] += sClPtr->x() * weight;
	baricenter[1] += sClPtr->y() * weight;
	baricenter[2] += sClPtr->z() * weight;
      }
    }

    if(total_weight == 0) continue;

    baricenter[0] /= total_weight;
    baricenter[1] /= total_weight;
    baricenter[2] /= total_weight;

    temp.setEnergy(total_weight);
    temp.setPosition(math::XYZPoint(baricenter[0], baricenter[1], baricenter[2]));
    temp.setAlgoId(reco::CaloCluster::hgcal_em);

    if(debug)    std::cout << " thisMC new energy = " << temp.energy() << " size = " << temp.size() << std::endl;

    result.push_back(temp);
  }

  return result;
}


std::vector<reco::HGCalMultiCluster> CheckResolution::getHighestEnergyMC(std::vector<reco::HGCalMultiCluster>& allmc, 
									 std::vector<reco::HGCalMultiCluster>& allmcMIP, 
									 std::vector<reco::HGCalMultiCluster>& all){

  if(debug)    std::cout << " >>> in getHighestEnergyMC initial size mc = " << allmc.size() << " mip = " << allmcMIP.size() << " empty = " << all.size() << std::endl;

  float maxEnergy[2] = {0., 0.};
  float secondEnergy[2] = {0., 0.};
  reco::HGCalMultiCluster maxMC[2];
  reco::HGCalMultiCluster secondMC[2];
  for(auto mc : allmc) {

    int etaSide = int(mc.eta() > 0);
    if(etaSide > 0 && mc.eta() < 0) std::cout << " getHighestEnergyMC etaSide = " << etaSide << " mc.eta() = " << mc.eta() << std::endl;
    if(mc.energy() > maxEnergy[etaSide]){
      secondMC[etaSide] = maxMC[etaSide];
      secondEnergy[etaSide] = maxEnergy[etaSide];
      maxMC[etaSide] = mc;
      maxEnergy[etaSide] = mc.energy();
    }
    else if(mc.energy() > secondEnergy[etaSide]){
      secondMC[etaSide] = mc;
      secondEnergy[etaSide] = mc.energy();
    }

    all.push_back(mc);
  }
  /*
  for(auto mc : allmcMIP) {

    int etaSide = int(mc.eta() > 0);
    if(etaSide > 0 && mc.eta() < 0) std::cout << " getHighestEnergyMC MIP etaSide = " << etaSide << " mc.eta() = " << mc.eta() << std::endl;
    if(mc.energy() > maxEnergy[etaSide]){
      secondMC[etaSide] = maxMC[etaSide];
      secondEnergy[etaSide] = maxEnergy[etaSide];
      maxMC[etaSide] = mc;
      maxEnergy[etaSide] = mc.energy();
    }
    else if(mc.energy() > secondEnergy[etaSide]){
      secondMC[etaSide] = mc;
      secondEnergy[etaSide] = mc.energy();
    }
    all.push_back(mc);
  }
  if(debug) std::cout << " eta < 0 energy = " << maxEnergy[0] << " = " << " maxMC[0].energy() = " << maxMC[0].energy() 
		      << " eta > 0 energy = " << maxEnergy[1] << " = " << " maxMC[1].energy() = " << maxMC[1].energy() 
		      << std::endl;

  */
  std::vector<reco::HGCalMultiCluster> result;
  for(int ij=0; ij<2; ++ij){
    if(maxEnergy[ij] > 0)
      result.push_back(maxMC[ij]);
    // if(secondEnergy[ij] > 0)
    //   result.push_back(secondMC[ij]);
  }
  return (result);
}


std::set<uint32_t> CheckResolution::getMatchedClusteredHitsList(bool pos,
								  const std::set<std::pair<uint32_t,float> > &hits,
								  const std::vector<reco::HGCalMultiCluster> &mcs, 
								  std::set<uint32_t> &mcHits) {
  

  if(debug)  std::cout << " in getMatchedClusteredHitsList mcs.size() = " << mcs.size()  << std::endl;

  for(auto mc : mcs) {
    
    //require on the same side
    bool mcPos(mc.eta()>0);
    if(mcPos!=pos) continue;
    
    //loop over all the layer clusters in the multicluster
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){
      const std::vector< std::pair<DetId, float> > &recHits = (*it)->hitsAndFractions();
      for(auto ij : recHits){
	if(ij.second != 0) 
	  mcHits.insert(ij.first.rawId());
      }
    }
  }

  if(debug)  std::cout << "  mcHits.size() = " <<   mcHits.size() << std::endl;
  std::set<uint32_t> matchedList = getMatched(hits, mcHits);
  //  for(auto ij : imatches) matchedList.insert(ij);

  if(debug)  std::cout << " looping in  getMatchedClusteredHitsList matchedList.size() = " <<  matchedList.size() << std::endl;
  return matchedList;
}



std::set<uint32_t> CheckResolution::getMatchedHitsList(bool pos,
						      const std::set<std::pair<uint32_t,float> > &hits,
						      const std::set<uint32_t> &mcs,
						      std::set<uint32_t> &mcHits) {
  

  if(debug)  std::cout << " in getMatchedHitsList" << std::endl;

  std::set<uint32_t> matchedList;
  
  std::set<uint32_t> imatches=getMatched(hits, mcs);
  for(auto ij : imatches) matchedList.insert(ij);
  

  if(debug)  std::cout << " >>> looping in getMatchedHitsList matchedList.size() = " <<  matchedList.size() << std::endl;
  return matchedList;
}



//
std::set<uint32_t> CheckResolution::getMatchedTrackedHitsList(bool pos,
							     const std::set<std::pair<uint32_t,float> >  &hits,
							     const std::vector<reco::CaloCluster> &ccs, 
							     std::set<uint32_t> &tkHits) {
  
  if(debug)  std::cout << " in getMatchedTrackedHitsList" << std::endl;

  std::set<uint32_t> matchedList;

  for(auto cc : ccs) {
    
    //require on the same side
    bool ccPos(cc.eta()>0);
    if(ccPos!=pos) continue;
    
    const std::vector< std::pair<DetId, float> > &recHits =cc.hitsAndFractions();
    for(auto ij : recHits){
      if(ij.second != 0) tkHits.insert(ij.first.rawId());
    }
    //std::vector<uint32_t> imatches=getMatched(hits,recHits);
    //matchedList.insert(matchedList.end(), imatches.begin(), imatches.end());
    std::set<uint32_t> imatches=getMatched(hits,tkHits);
    for(auto ij : imatches) matchedList.insert(ij);
  }
  
  if(debug)  std::cout << " >>> looping  getMatchedTrackedHitsList matchedList.size() = " <<  matchedList.size() << std::endl;
  return matchedList;  
}




std::set<uint32_t> CheckResolution::getMatched(const std::set<std::pair<uint32_t,float> > &a,
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


std::set<uint32_t> CheckResolution::getMatched(const std::set<uint32_t> &a,
					      const std::set<uint32_t> &b){

  std::set<uint32_t> matchedList;
  for(auto ii : a) {
    //find first match in second list
    for(auto jj : b) {
      if(ii != jj) continue;
      matchedList.insert(ii);
      break;
    }
  }

  return matchedList;
}


//in 1st not in 2nd
std::set<uint32_t> CheckResolution::getNOTMatched(const std::set<uint32_t> &a,
						 const std::set<uint32_t> &b){

  std::set<uint32_t> matchedList;
  for(auto ii : a) {
    bool isNew = true;
    //find first match in second list
    for(auto jj : b) {
      if(ii == jj) {
	isNew = false;
	break;
      }
    }
    if(isNew) matchedList.insert(ii);
  }

  return matchedList;
}



DEFINE_FWK_MODULE(CheckResolution);
