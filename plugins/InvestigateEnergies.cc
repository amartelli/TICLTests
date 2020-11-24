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
#include <numeric>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TTree.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/Common/interface/ValueMap.h"



using namespace std;
using namespace edm;
using namespace reco;

class InvestigateEnergies : public edm::EDAnalyzer {
public:

  explicit InvestigateEnergies(const edm::ParameterSet&);

  ~InvestigateEnergies();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

  virtual void beginJob() override;
  virtual void endJob() override;

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
  void get2DclAssociated(std::map<uint32_t, std::vector<uint32_t> >& mappa, std::map<uint32_t, std::vector<uint32_t> >& mappa_ni, 
			 const std::vector<CaloParticle>& cP, const std::vector<reco::CaloCluster>& lC);


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


  std::vector<size_t> decrease_sorted_indices(const std::vector<float>& v);

  std::pair<float, float> fixSizeHighestDensity(std::vector<float>& time, std::vector<float> weight,
                                                unsigned int minNhits = 3, float deltaT=0.210, float timeWidthBy=0.5);




  edm::EDGetTokenT<std::vector<CaloParticle> > genToken_;
  edm::EDGetTokenT<float> tGenToken_;

  //  edm::EDGetTokenT<std::vector<reco::CaloCluster> > tkToken_;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > mcTrkToken_, mcMIPToken_,mcToken_;

  edm::EDGetTokenT<std::vector<reco::CaloCluster> > caloClToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hits_eeToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hits_fhToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hits_bhToken_;

  edm::EDGetTokenT<edm::ValueMap<pair<float,float> > > timeMap_;
  edm::Handle<edm::ValueMap<pair<float, float> > > time2DMap;

  std::map<uint32_t, const HGCRecHit*> hitmap;
 
  std::map<uint32_t, std::vector<uint32_t> > caloP_2dCl_map;
  std::map<uint32_t, std::vector<uint32_t> > caloP_2dCl_map_ni;

  TH2F* layerClenergy_vsCaloPpt;
  TH2F* layerClenergy_vsCaloPenergy;
  TH2F* layerClenergySum_vsCaloPenergy;
  TH2F* layerClenergy_vslayerClenergySum;

  TH2F* fr_layerClenergy_vsCaloPenergy;
  TH2F* fr_layerClenergySum_vsCaloPenergy;
  TH2F* fr_layerClenergy_vslayerClenergySum;

  ///// nn interacting
  TH2F* layerClenergy_vsCaloPenergy_ni;
  TH2F* layerClenergySum_vsCaloPenergy_ni;
  TH2F* layerClenergy_vslayerClenergySum_ni;


  TH2F* dR_vs_layer;
  TH2F* dEta_vs_layer;
  TH2F* dPhi_vs_layer;

  
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



InvestigateEnergies::InvestigateEnergies(const edm::ParameterSet& iConfig) :
  genToken_(consumes<std::vector<CaloParticle> >(edm::InputTag("mix:MergedCaloTruth"))),
  tGenToken_(consumes<float>(edm::InputTag("genParticles:t0"))),
  //tkToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcTracks::RECO2"))),
  mcTrkToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("ticlMultiClustersFromTrackstersTrk::RECO"))),
  mcMIPToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("ticlMultiClustersFromTrackstersMIP::RECO"))),
  mcToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("ticlMultiClustersFromTrackstersMerge::RECO"))),
  caloClToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcalLayerClusters::RECO"))),
  hits_eeToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"))),
  hits_fhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"))),
  hits_bhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"))),
  timeMap_(consumes<edm::ValueMap<pair<float,float> > > (edm::InputTag("hgcalLayerClusters:timeLayerCluster")))
{
  //book some histos here                                                                                                                 

  edm::Service<TFileService> fs;

  layerClenergy_vsCaloPpt = fs->make<TH2F>("layerClenergy_vsCaloPpt", "", 5000, 0., 20., 500, 0., 10.);
  layerClenergy_vsCaloPenergy = fs->make<TH2F>("layerClenergy_vsCaloPenergy", "", 5000, 0., 20., 500, 0., 10.);
  layerClenergySum_vsCaloPenergy = fs->make<TH2F>("layerClenergySum_vsCaloPenergy", "", 5000, 0., 20., 500, 0., 10.);
  layerClenergy_vslayerClenergySum = fs->make<TH2F>("layerClenergy_vslayerClenergySum", "", 5000, 0., 20., 500, 0., 10.);

  fr_layerClenergy_vsCaloPenergy = fs->make<TH2F>("fr_layerClenergy_vsCaloPenergy", "", 5000, 0., 20., 50, 0., 1.);
  fr_layerClenergySum_vsCaloPenergy = fs->make<TH2F>("fr_layerClenergySum_vsCaloPenergy", "", 5000, 0., 20., 50, 0., 1.);
  fr_layerClenergy_vslayerClenergySum = fs->make<TH2F>("fr_layerClenergy_vslayerClenergySum", "", 5000, 0., 20., 50, 0., 1.);

  ///non interactiong
  layerClenergy_vsCaloPenergy_ni = fs->make<TH2F>("layerClenergy_vsCaloPenergy_ni", "", 5000, 0., 20., 500, 0., 10.);
  layerClenergySum_vsCaloPenergy_ni = fs->make<TH2F>("layerClenergySum_vsCaloPenergy_ni", "", 5000, 0., 20., 500, 0., 10.);
  layerClenergy_vslayerClenergySum_ni = fs->make<TH2F>("layerClenergy_vslayerClenergySum_ni", "", 5000, 0., 20., 500, 0., 10.);


  dR_vs_layer = fs->make<TH2F>("dR_vs_layer", "", 100, -50., 50., 500, -10., 10.);
  dEta_vs_layer = fs->make<TH2F>("dEta_vs_layer", "", 100, -50., 50., 500, -10., 10.);
  dPhi_vs_layer = fs->make<TH2F>("dPhi_vs_layer", "", 100, -50., 50., 500, -10., 10.);

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



InvestigateEnergies::~InvestigateEnergies() { 
}



void InvestigateEnergies::beginRun(const edm::Run& run, 
			       const edm::EventSetup & es) { }


void  InvestigateEnergies::analyze(const Event& iEvent, 
			       const EventSetup& iSetup) {

  ++nEvent;
  edm::ESHandle<CaloGeometry> geom;
  iSetup.get<CaloGeometryRecord>().get(geom);
  recHitTools.setGeometry(*geom);
  
  //  recHitTools.getEventSetup(iSetup);
  
  if(debug)  std::cout << " \n \n new evt "<< std::endl;
  
  edm::Handle<std::vector<CaloParticle> > gpH;
  iEvent.getByToken( genToken_, gpH);

  iEvent.getByToken(timeMap_, time2DMap);
  if(debug) std::cout << " time2DMap->size() = " << time2DMap->size() << std::endl;

  edm::Handle<std::vector<reco::HGCalMultiCluster> > mcMIPH, mcH;
  iEvent.getByToken( mcMIPToken_,mcMIPH);
  iEvent.getByToken( mcToken_, mcH);

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

  edm::Handle<std::vector<reco::CaloCluster>> cluster_h;
  iEvent.getByToken(caloClToken_, cluster_h);
  const auto &layerClusters = *cluster_h;

  caloP_2dCl_map.clear();
  caloP_2dCl_map_ni.clear();
  //map caloParticle-2dCl with more than 90% energy
  get2DclAssociated(caloP_2dCl_map, caloP_2dCl_map_ni, *gpH, layerClusters);
  //  std::cout << " tornato size = " << caloP_2dCl_map.size()  << std::endl;

  //  std::cout << " >>> now filling histo " << std::endl;

  /*
  for(auto ij:caloP_2dCl_map){
    float CPenergy = (*gpH)[ij.first].energy();
    if(CPenergy <= 0) std::cout << " >> CPenergy = " << CPenergy << " evento = " << nEvent << std::endl;
    float sum2Denergy = 0.;
    for(auto ijH:ij.second){
      float layerCenergy = layerClusters[ijH].energy();
      if(layerCenergy < 0) std::cout << " >> layerCenergy = " << layerCenergy << std::endl;
      layerClenergy_vsCaloPenergy->Fill(CPenergy, layerCenergy);    
      sum2Denergy += layerCenergy;
    }
    if(sum2Denergy != 0.) {
      layerClenergySum_vsCaloPenergy->Fill(CPenergy, sum2Denergy);
      for(auto ijH:ij.second){
	float layerCenergy = layerClusters[ijH].energy();
	layerClenergy_vslayerClenergySum->Fill(sum2Denergy, layerCenergy);    
      }
    }
  }

  std::cout << " nuovo evento " << std::endl;
  */
  /*
  int firstBin = 0;
  //  std::cout << " layerClenergy_vsCaloPenergy->GetNbinsX() = " << layerClenergy_vsCaloPenergy->GetNbinsX() << std::endl;
  for(int iB = 1; iB < layerClenergy_vsCaloPenergy->GetNbinsX()+1; ++iB){
    if(firstBin != 0 && iB > firstBin  + 5) break;
    float xVal = (iB-1) * 0.004 + 0.002;
    //    std::cout << " 2D binX = " << xVal << std::endl;
    auto dummy = (TH1D*)layerClenergy_vsCaloPenergy->ProjectionY("dummy", iB, iB+1);
    int allEvt = dummy->GetEntries();
    float localEvt = 0.;

    auto dummySumCP = (TH1D*)layerClenergySum_vsCaloPenergy->ProjectionY("dummySumCP", iB, iB+1);
    int allEvtSumCP = dummySumCP->GetEntries();
    float localEvtSumCP = 0.;

    auto dummySum = (TH1D*)layerClenergy_vslayerClenergySum->ProjectionY("dummySum", iB, iB+1);
    int allEvtSum = dummySum->GetEntries();
    float localEvtSum = 0.;

    if(allEvt == 0 || allEvtSumCP == 0 || allEvtSum == 0) continue;
    if( firstBin==0) firstBin = iB;
    for(int ij = 1; ij<dummy->GetNbinsX()+1; ++ij){
      float yVal = (ij-1) * 0.02 + 0.01;
      //      std::cout << " 2D binY = " << yVal << " dummy->GetNbinsX() = " << dummy->GetNbinsX() << " ij = " << ij << " " << (ij-1 * 0.02) + 0.01 << std::endl;
      localEvt += dummy->GetBinContent(ij);
      fr_layerClenergy_vsCaloPenergy->Fill(xVal, yVal/xVal, (allEvt - localEvt) / allEvt);
      std::cout << " x = " << xVal << " y = " << yVal/xVal << " val = " << (allEvt - localEvt) / allEvt << std::endl;

      localEvtSumCP += dummySumCP->GetBinContent(ij);
      fr_layerClenergySum_vsCaloPenergy->Fill(xVal, yVal/xVal, (allEvtSumCP - localEvtSumCP) / allEvtSumCP);

      localEvtSum += dummySum->GetBinContent(ij);
      fr_layerClenergy_vslayerClenergySum->Fill(xVal, yVal/xVal, (allEvtSum - localEvtSum) / allEvtSum);
    }
    delete dummy;
    delete dummySum;
    delete dummySumCP;

  }
  */

  /*
  // non interacting
  for(auto ij:caloP_2dCl_map_ni){
    float CPenergy = (*gpH)[ij.first].energy();
    if(CPenergy < 0) std::cout << " >> CPenergy = " << CPenergy << " evento = " << nEvent << std::endl;
    float sum2Denergy = 0.;
    for(auto ijH:ij.second){
      float layerCenergy = layerClusters[ijH].energy();
      if(layerCenergy < 0) std::cout << " >> layerCenergy = " << layerCenergy << std::endl;
      layerClenergy_vsCaloPenergy_ni->Fill(CPenergy, layerCenergy);    
      sum2Denergy += layerCenergy;
    }
    if(sum2Denergy != 0.) {
      layerClenergySum_vsCaloPenergy_ni->Fill(CPenergy, sum2Denergy);
      for(auto ijH:ij.second){
	float layerCenergy = layerClusters[ijH].energy();
	layerClenergy_vslayerClenergySum_ni->Fill(sum2Denergy, layerCenergy);    
      }
    }
  }

  */

  return;

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

std::vector<reco::HGCalMultiCluster> InvestigateEnergies::cleanTimeMC(const std::vector<reco::HGCalMultiCluster>& allmc){

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
      float time2D = ((*time2DMap)[sClPtr]).first;
      float timeE = ((*time2DMap)[sClPtr]).second;
      localTimeAll[iCount] = time2D;
      if(time2D > -50.){
        localTime.push_back(time2D);
        //localEnergy.push_back(sClPtr->energy());
	localEnergy.push_back(1. / pow(timeE, 2));
      }
      ++iCount;
    }

    if(debug) std::cout << " now compute time localTime.size() = " << localTime.size() << std::endl;
    //use method from Shameena => try for the moment with truncation
    // can try truncation + weighted mean
    float finalT = (localTime.size() >= 3) ? (fixSizeHighestDensity(localTime, localEnergy)).first : (-99);
    float finalTW = (localTime.size() >= 3) ? (fixSizeHighestDensity(localTime, localEnergy)).first : (-99);

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


std::vector<reco::HGCalMultiCluster> InvestigateEnergies::getHighestEnergyMC(std::vector<reco::HGCalMultiCluster>& allmc, 
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


std::set<uint32_t> InvestigateEnergies::getMatchedClusteredHitsList(bool pos,
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



std::set<uint32_t> InvestigateEnergies::getMatchedHitsList(bool pos,
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
std::set<uint32_t> InvestigateEnergies::getMatchedTrackedHitsList(bool pos,
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




std::set<uint32_t> InvestigateEnergies::getMatched(const std::set<std::pair<uint32_t,float> > &a,
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


std::set<uint32_t> InvestigateEnergies::getMatched(const std::set<uint32_t> &a,
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
std::set<uint32_t> InvestigateEnergies::getNOTMatched(const std::set<uint32_t> &a,
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

std::vector<size_t> InvestigateEnergies::decrease_sorted_indices(const std::vector<float>& v) {
  // initialize original index locations                                                       
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indices based on comparing values in v (decreasing order)                            
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  return idx;
};

std::pair<float, float> InvestigateEnergies::fixSizeHighestDensity(std::vector<float>& time, std::vector<float> weight,
                                                               unsigned int minNhits, float deltaT, float timeWidthBy) {
  if (time.size() < minNhits)
    return std::pair<float, float>(-99., -1.);

  if (weight.empty())
    weight.resize(time.size(), 1.);

  std::vector<float> t(time.size(), 0.);
  std::vector<float> w(time.size(), 0.);
  std::vector<size_t> sortedIndex = decrease_sorted_indices(time);
  for (std::size_t i = 0; i < sortedIndex.size(); ++i) {
    t[i] = time[sortedIndex[i]];
    w[i] = weight[sortedIndex[i]];
  }

  int max_elements = 0;
  int start_el = 0;
  int end_el = 0;
  float timeW = 0.f;
  float tolerance = 0.05f;

  for (auto start = t.begin(); start != t.end(); ++start) {
    const auto startRef = *start;
    int c = count_if(start, t.end(), [&](float el) { return el - startRef <= deltaT + tolerance; });
    if (c > max_elements) {
      max_elements = c;
      auto last_el = find_if_not(start, t.end(), [&](float el) { return el - startRef <= deltaT + tolerance; });
      auto valTostartDiff = *(--last_el) - startRef;
      if (std::abs(deltaT - valTostartDiff) < tolerance) {
        tolerance = std::abs(deltaT - valTostartDiff);
      }
      start_el = distance(t.begin(), start);
      end_el = distance(t.begin(), last_el);
      timeW = valTostartDiff;
    }
  }

  // further adjust time width around the chosen one based on the hits density
  // proved to improve the resolution: get as many hits as possible provided they are close in time
  float HalfTimeDiff = timeW * timeWidthBy;
  float sum = 0.;
  float num = 0;
  int totSize = t.size();

  for (int ij = 0; ij <= start_el; ++ij) {
    if (t[ij] > (t[start_el] - HalfTimeDiff)) {
      for (int kl = ij; kl < totSize; ++kl) {
        if (t[kl] < (t[end_el] + HalfTimeDiff)) {
          sum += t[kl] * w[kl];
          num += w[kl];
        } else
          break;
      }
      break;
    }
  }

  if (num == 0) {
    return std::pair<float, float>(-99., -1.);
  }
  return std::pair<float, float>(sum / num, 1. / sqrt(num));
}


void InvestigateEnergies::get2DclAssociated(std::map<uint32_t, std::vector<uint32_t> >& mappa, 
					    std::map<uint32_t, std::vector<uint32_t> >& mappa_ni, 
					    const std::vector<CaloParticle>& cP, 
					    const std::vector<reco::CaloCluster>& lC){


  std::vector<int> used;
  used.resize(lC.size());

  if(debug) std::cout << " now double loop for match " << std::endl; 
  for(auto itCP : cP){
    float ptCP = itCP.pt();
    float etaCP = itCP.eta();
    if(std::abs(itCP.eta()) < 1.4 || itCP.pt() < 2.) continue;

    float energyCP = itCP.energy();
    if(debug) std::cout << " itCP.energy = " << itCP.energy() << std::endl;
    if(itCP.energy() < 0) std::cout << " NEGGGG  itCP.energy = " << itCP.energy() << std::endl;
    float phiCP = itCP.phi();
    int nIntCP = int(itCP.simClusters().size());
    std::vector<uint32_t> cpHits;
    cpHits.clear();
    for(CaloParticle::sc_iterator scIt=itCP.simCluster_begin(); scIt!=itCP.simCluster_end(); scIt++){
      //all hits and energy fractions at sim level 
      std::vector<std::pair<uint32_t,float> > hits = (*scIt)->hits_and_fractions();
      if(debug) std::cout << " >>> hits.size() = " << hits.size() << std::endl;
      for(auto ij : hits) {
        auto finder2 = hitmap.find(ij.first);
        if(finder2 != hitmap.end()){
          if(debug) std::cout << " trovato nel gen " << std::endl;
          //hitsInCaloP[counter].push_back(ij.first);
	  cpHits.push_back(ij.first);
        }
      }
    }//loop over CP clusters

    float sum2Denergy = 0.;
    std::vector<float> clEnergy;
    clEnergy.clear();
    unsigned int counter2D = 0;
    for(const reco::CaloCluster &sCl : lC){
      if(used[counter2D] == 1) continue;
      float energy2D = sCl.energy();
      if(energy2D == 0) continue;
      if(debug) std::cout << " energy2D = " << energy2D << std::endl;
      float eta2D = sCl.eta();
      float phi2D = sCl.phi();
      if(debug) std::cout << " 2D energy = " << energy2D << " eta = " << eta2D << " phi = " << phi2D<< std::endl;
      if(eta2D * etaCP < 0.) continue;

      const HGCalDetId hitid = sCl.hitsAndFractions()[0].first;
      float z2D = recHitTools.getPosition(hitid).z();
      int layer2D = recHitTools.getLayerWithOffset(hitid);
      //if(layer2D == 1) continue;
      float dEta = (eta2D - etaCP);
      float dPhi = reco::deltaPhi(phi2D, phiCP);
      dEta_vs_layer->Fill(z2D > 0 ? layer2D : -1.*layer2D, dEta);
      dPhi_vs_layer->Fill(z2D > 0 ? layer2D : -1.*layer2D, dPhi);
      float deltaR = sqrt(dEta*dEta + dPhi*dPhi);
      if(debug) std::cout << " cp energy = " << energyCP << " eta = " << etaCP << " phi = " << phiCP << " deltaR = " << deltaR << std::endl;
      layer2D = z2D > 0 ? layer2D : -1.*layer2D;
      dR_vs_layer->Fill(layer2D, deltaR);
      //std::cout << " layer2D = " << layer2D << std::endl;

      float energy2D_local = 0.;
      for(auto ij2D : sCl.hitsAndFractions()){
	for(auto ijHitsCP : cpHits){
	  if(ij2D.first == ijHitsCP){
	    if(debug) std::cout << " hit in 2Dcl matchata a CP " << std::endl;
	    const HGCRecHit *hit = hitmap[ij2D.first];
	    energy2D_local += hit->energy();
	    break;
	  }
	}
      }
      if(energy2D_local > energy2D * 0.8){
	if(debug) std::cout << " inserito" << std::endl;
	used[counter2D] = 1;
	layerClenergy_vsCaloPpt->Fill(ptCP, energy2D);
	layerClenergy_vsCaloPenergy->Fill(energyCP, energy2D);
	sum2Denergy += energy2D;
        clEnergy.push_back(energy2D);

	if(nIntCP == 1){
          //mappa_ni[ijCP.first].push_back(ijCl.first);
          layerClenergy_vsCaloPenergy_ni->Fill(energyCP, energy2D);
        }
      }
      ++counter2D;
    } //2Dcl loop
    
    if(sum2Denergy != 0.){
      layerClenergySum_vsCaloPenergy->Fill(energyCP, sum2Denergy);
      if(nIntCP == 1) layerClenergySum_vsCaloPenergy_ni->Fill(energyCP, sum2Denergy);
      
      for(auto iE : clEnergy){
	layerClenergy_vslayerClenergySum->Fill(sum2Denergy, iE);
        if(nIntCP == 1) layerClenergy_vslayerClenergySum_ni->Fill(sum2Denergy, iE);
      }
    }
  }// loop CP



  /*
  for(auto ijCP : hitsInCaloP){
    float energyCP = cP[ijCP.first].energy();
    float etaCP = cP[ijCP.first].eta();
    float phiCP = cP[ijCP.first].phi();
    int nIntCP = int(cP[ijCP.first].simClusters().size());

    float sum2Denergy = 0.;  
    std::vector<float> clEnergy;
    clEnergy.clear();

    for(auto ijCl : hitsIn2Dcl){
      float energy2D = lC[ijCl.first].energy();
      float eta2D = lC[ijCl.first].eta();
      float phi2D = lC[ijCl.first].phi();
      const HGCalDetId hitid = lC[ijCl.first].hitsAndFractions()[0].first;
      float z2D = recHitTools.getPosition(hitid).z();
      int layer2D = recHitTools.getLayerWithOffset(hitid);

      if(debug) std::cout << " 2D energy = " << energy2D << " eta = " << eta2D << " phi = " << phi2D<< std::endl;
      if(eta2D * etaCP < 0.) continue;

      float dEta = (eta2D - etaCP);
      float dPhi = reco::deltaPhi(phi2D, phiCP);
      dEta_vs_layer->Fill(z2D > 0 ? layer2D : -1.*layer2D, dEta);
      dPhi_vs_layer->Fill(z2D > 0 ? layer2D : -1.*layer2D, dPhi);
      float deltaR = sqrt(dEta*dEta + dPhi*dPhi);
      if(debug) std::cout << " cp energy = " << energyCP << " eta = " << etaCP << " phi = " << phiCP << " deltaR = " << deltaR << std::endl; 
      dR_vs_layer->Fill(z2D > 0 ? layer2D : -1.*layer2D, deltaR);

      float energy2D_local_v2 = 0.;
      for(auto ijHits2D : ijCl.second){
      	for(auto ijHitsCP : ijCP.second){
	  if(ijHits2D.first == ijHitsCP){
	    energy2D_local_v2 += ijHits2D.second;
	    if(debug) std::cout << " matchato caloP 2Dcl => ratio local/all = " << energy2D_local_v2 / energy2D << std::endl;
	    break;
	  }
	}
      }
      if(energy2D_local_v2 > energy2D * 0.8){
	mappa[ijCP.first].push_back(ijCl.first);
	layerClenergy_vsCaloPenergy->Fill(energyCP, energy2D);
	sum2Denergy += energy2D;
	clEnergy.push_back(energy2D);

	if(nIntCP == 1){
	  mappa_ni[ijCP.first].push_back(ijCl.first);
	  layerClenergy_vsCaloPenergy_ni->Fill(energyCP, energy2D);
	}
	if(debug) std::cout << " inserito" << std::endl;
      }
      //      else if(energy2D_local_v2 != 0) std::cout << " typic E ratio = " << energy2D_local_v2/energy2D << std::endl;
    }//hits in 2Dcl
    if(sum2Denergy != 0.){
      layerClenergySum_vsCaloPenergy->Fill(energyCP, sum2Denergy);
      if(nIntCP == 1) layerClenergySum_vsCaloPenergy_ni->Fill(energyCP, sum2Denergy);

      for(auto iE : clEnergy){
	layerClenergy_vslayerClenergySum->Fill(sum2Denergy, iE);
	if(nIntCP == 1) layerClenergy_vslayerClenergySum_ni->Fill(sum2Denergy, iE);
      }
    }
  }//hits in CP



  //caloparticle level map
  std::map<uint32_t, std::vector<uint32_t>> hitsInCaloP;
  unsigned int counter = 0;
  for(auto itCP : cP){
    //    std::cout << " CP energy = " << itCP.energy() << " pt = " << itCP.pt() << std::endl;
    if(std::abs(itCP.eta()) < 1.4 || itCP.pt() < 2.) continue;
    if(debug) std::cout << " itCP.energy = " << itCP.energy() << std::endl;
    if(itCP.energy() < 0) std::cout << " NEGGGG  itCP.energy = " << itCP.energy() << std::endl;
    for(CaloParticle::sc_iterator scIt=itCP.simCluster_begin(); scIt!=itCP.simCluster_end(); scIt++){
      
      //all hits and energy fractions at sim level
      std::vector<std::pair<uint32_t,float> > hits = (*scIt)->hits_and_fractions();
      if(debug) std::cout << " >>> hits.size() = " << hits.size() << std::endl;
      for(auto ij : hits) {
	
	auto finder2 = hitmap.find(ij.first);
	if(finder2 != hitmap.end()){
	  
	  if(debug) std::cout << " trovato nel gen " << std::endl;
	  hitsInCaloP[counter].push_back(ij.first);
	}
      }
    }
    ++counter;
  }

  //layer cluster level
  std::map<uint32_t, std::vector<std::pair<uint32_t, float>>> hitsIn2Dcl;
  unsigned int counter2D = 0;
  for(const reco::CaloCluster &sCl : lC){
    float energy2D = sCl.energy();  
    if(debug) std::cout << " energy2D = " << energy2D << std::endl;
    float energy2D_local = 0.;
    for(auto ij2D : sCl.hitsAndFractions()){
      auto finder2D = hitmap.find(ij2D.first);
      if(finder2D != hitmap.end()){
	const HGCRecHit *hit = hitmap[ij2D.first];
	energy2D_local += hit->energy();
	if(debug) std::cout << " trovata reco hit energy = " << hit->energy() << std::endl;
	if(hit->energy() < 0) std::cout << " NEGGGG  hit->energy() = " << hit->energy() << std::endl;
	if(debug) std::cout << " trovato nel 2D energy = " << hit->energy() << " energy2D_local = " << energy2D_local << " ratio to 2D = " << energy2D_local/energy2D << std::endl;
	hitsIn2Dcl[counter2D].push_back(std::pair<uint32_t, float>(ij2D.first, hit->energy()));
	if(debug) std::cout << " counter2D = " << counter2D << " hitsIn2Dcl.size() = " << hitsIn2Dcl.size() << " hitsIn2Dcl[counter2D].size = " << hitsIn2Dcl[counter2D].size() << std::endl;
	//delete hit;
	if(debug) std::cout << " fine giro 1 " << std::endl;
      }
      if(debug) std::cout << " fine giro 2 " << std::endl;
    }
    if(debug) std::cout << " fine giro 3 " << std::endl;
    ++counter2D;
  }
  
  return;
  if(debug) std::cout << " now double loop for match " << std::endl; 
  for(auto ijCP : hitsInCaloP){
    float energyCP = cP[ijCP.first].energy();
    float etaCP = cP[ijCP.first].eta();
    float phiCP = cP[ijCP.first].phi();
    int nIntCP = int(cP[ijCP.first].simClusters().size());

    float sum2Denergy = 0.;  
    std::vector<float> clEnergy;
    clEnergy.clear();

    for(auto ijCl : hitsIn2Dcl){
      float energy2D = lC[ijCl.first].energy();
      float eta2D = lC[ijCl.first].eta();
      float phi2D = lC[ijCl.first].phi();
      const HGCalDetId hitid = lC[ijCl.first].hitsAndFractions()[0].first;
      float z2D = recHitTools.getPosition(hitid).z();
      int layer2D = recHitTools.getLayerWithOffset(hitid);

      if(debug) std::cout << " 2D energy = " << energy2D << " eta = " << eta2D << " phi = " << phi2D<< std::endl;
      if(eta2D * etaCP < 0.) continue;

      float dEta = (eta2D - etaCP);
      float dPhi = reco::deltaPhi(phi2D, phiCP);
      dEta_vs_layer->Fill(z2D > 0 ? layer2D : -1.*layer2D, dEta);
      dPhi_vs_layer->Fill(z2D > 0 ? layer2D : -1.*layer2D, dPhi);
      float deltaR = sqrt(dEta*dEta + dPhi*dPhi);
      if(debug) std::cout << " cp energy = " << energyCP << " eta = " << etaCP << " phi = " << phiCP << " deltaR = " << deltaR << std::endl; 
      dR_vs_layer->Fill(z2D > 0 ? layer2D : -1.*layer2D, deltaR);

      float energy2D_local_v2 = 0.;
      for(auto ijHits2D : ijCl.second){
      	for(auto ijHitsCP : ijCP.second){
	  if(ijHits2D.first == ijHitsCP){
	    energy2D_local_v2 += ijHits2D.second;
	    if(debug) std::cout << " matchato caloP 2Dcl => ratio local/all = " << energy2D_local_v2 / energy2D << std::endl;
	    break;
	  }
	}
      }
      if(energy2D_local_v2 > energy2D * 0.8){
	mappa[ijCP.first].push_back(ijCl.first);
	layerClenergy_vsCaloPenergy->Fill(energyCP, energy2D);
	sum2Denergy += energy2D;
	clEnergy.push_back(energy2D);

	if(nIntCP == 1){
	  mappa_ni[ijCP.first].push_back(ijCl.first);
	  layerClenergy_vsCaloPenergy_ni->Fill(energyCP, energy2D);
	}
	if(debug) std::cout << " inserito" << std::endl;
      }
      //      else if(energy2D_local_v2 != 0) std::cout << " typic E ratio = " << energy2D_local_v2/energy2D << std::endl;
    }//hits in 2Dcl
    if(sum2Denergy != 0.){
      layerClenergySum_vsCaloPenergy->Fill(energyCP, sum2Denergy);
      if(nIntCP == 1) layerClenergySum_vsCaloPenergy_ni->Fill(energyCP, sum2Denergy);

      for(auto iE : clEnergy){
	layerClenergy_vslayerClenergySum->Fill(sum2Denergy, iE);
	if(nIntCP == 1) layerClenergy_vslayerClenergySum_ni->Fill(sum2Denergy, iE);
      }
    }
  }//hits in CP
  */

  if(debug) std::cout << " mappa.size() = " << mappa.size() << std::endl;
  return;
}





void
InvestigateEnergies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
InvestigateEnergies::endJob()
{
  int firstBin = 0;
  //  std::cout << " layerClenergy_vsCaloPenergy->GetNbinsX() = " << layerClenergy_vsCaloPenergy->GetNbinsX() << std::endl;                                            
  for(int iB = 1; iB < layerClenergy_vsCaloPenergy->GetNbinsX()+1; ++iB){
    //    if(firstBin != 0 && iB > firstBin  + 5) break;
    float xVal = (iB-1) * 0.004 + 0.002;
    //    std::cout << " 2D binX = " << xVal << std::endl;                                                                                                             
    auto dummy = (TH1D*)layerClenergy_vsCaloPenergy->ProjectionY("dummy", iB, iB+1);
    int allEvt = dummy->GetEntries();
    float localEvt = 0.;
    auto dummySumCP = (TH1D*)layerClenergySum_vsCaloPenergy->ProjectionY("dummySumCP", iB, iB+1);
    int allEvtSumCP = dummySumCP->GetEntries();
    float localEvtSumCP = 0.;
    auto dummySum = (TH1D*)layerClenergy_vslayerClenergySum->ProjectionY("dummySum", iB, iB+1);
    int allEvtSum = dummySum->GetEntries();
    float localEvtSum = 0.;

    if(allEvt == 0 || allEvtSumCP == 0 || allEvtSum == 0) continue;

    if( firstBin==0) firstBin = iB;
    for(int ij = 1; ij<dummy->GetNbinsX()+1; ++ij){
      float yVal = (ij-1) * 0.02 + 0.01;
      //      std::cout << " 2D binY = " << yVal << " dummy->GetNbinsX() = " << dummy->GetNbinsX() << " ij = " << ij << " " << (ij-1 * 0.02) + 0.01 << std::endl;      
      int xBin = fr_layerClenergy_vsCaloPenergy->GetXaxis()->FindBin(xVal);
      int yBin = fr_layerClenergy_vsCaloPenergy->GetYaxis()->FindBin(yVal/xVal);

      localEvt += dummy->GetBinContent(ij);
      fr_layerClenergy_vsCaloPenergy->SetBinContent(xBin, yBin, (allEvt - localEvt) / allEvt);
      //      std::cout << " x = " << xVal << " y = " << yVal/xVal << " val = " << (allEvt - localEvt) / allEvt << std::endl;
	    
      localEvtSumCP += dummySumCP->GetBinContent(ij);
      fr_layerClenergySum_vsCaloPenergy->SetBinContent(xBin, yBin, (allEvtSumCP - localEvtSumCP) / allEvtSumCP);

      localEvtSum += dummySum->GetBinContent(ij);
      fr_layerClenergy_vslayerClenergySum->SetBinContent(xBin, yBin, (allEvtSum - localEvtSum) / allEvtSum);
    }
    delete dummy;
    delete dummySum;
    delete dummySumCP;
  }

}




DEFINE_FWK_MODULE(InvestigateEnergies);

