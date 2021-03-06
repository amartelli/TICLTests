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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "RecoHGCal/TICL/interface/Trackster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "TSystem.h"
#include "TH1F.h"
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

class BuildParticles : public edm::EDAnalyzer {
public:

  explicit BuildParticles(const edm::ParameterSet&);

  ~BuildParticles();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

  //  usesResource("TFileService");

private:

  hgcal::RecHitTools         recHitTools;

  void getClusteredHitsList(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs);
  void getTrackedHitsList(bool pos, const std::vector<reco::CaloCluster> &ccs);

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
  
  std::vector<uint32_t> getMatched(const std::set<std::pair<uint32_t,float> > &a,
                                   const std::vector<std::pair<DetId,float> > &b);
  
  std::set<uint32_t> getMatched(const std::set<std::pair<uint32_t,float> > &a,
                                const std::set<uint32_t> &b);



  std::set<uint32_t> getMatched(const std::set<uint32_t> &a,
				const std::set<uint32_t> &b);

  edm::EDGetTokenT<std::vector<CaloParticle> > genToken_;
  //edm::EDGetTokenT<std::vector<SimCluster> > simclusToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster> > tkToken_;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > mcMIPToken_,mcToken_;

  edm::EDGetTokenT<HGCRecHitCollection> hits_eeToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hits_fhToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hits_bhToken_;


  edm::EDGetTokenT<reco::CaloClusterCollection> layerToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> timeMap_;

  edm::Handle<reco::CaloClusterCollection> layerClusters;  
  edm::Handle<edm::ValueMap<float> > time2DMap;
  std::map<uint32_t, const HGCRecHit*> hitmap;

  std::map<TString,TH1 *> histos_;
  TH1F* caloEnergyRatio;
  TH1F* caloEnergyRatioPos;
  TH1F* caloEnergyRatioNeg;

  std::vector<std::set<uint32_t>> finalParticles;
  std::vector<std::vector<float>> finalParticles_etaPhi;
  std::vector<float> finalParticles_time;
  std::vector<float> tkParticles_time;
  std::vector<float> mcParticles_time;
  std::vector<float> mipParticles_time;
  bool debug;
};



BuildParticles::BuildParticles(const edm::ParameterSet& iConfig) :
  genToken_(consumes<std::vector<CaloParticle> >(edm::InputTag("mix:MergedCaloTruth"))),
  //  simclusToken_(consumes<std::vector<SimCluster> >(edm::InputTag("mix:MergedCaloTruth"))),
  tkToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcTracks::RECO2"))),
  mcMIPToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiClusterMIP:MIPMultiClustersFromTracksterByCA"))),
  mcToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiCluster:MultiClustersFromTracksterByCA"))),
  hits_eeToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"))),
  hits_fhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"))),
  hits_bhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"))),
  layerToken_(consumes<reco::CaloClusterCollection>(edm::InputTag("hgcalLayerClusters::RECO"))),
  timeMap_(consumes<edm::ValueMap<float>>(edm::InputTag("hgcalLayerClusters:timeLayerCluster")))
{
  //book some histos here                                                                                                                 
  //histos_["dr"]   = new TH1F("dr",   ";#Delta R;", 50,0,0.5);                                                                           

  edm::Service<TFileService> fs;

  caloEnergyRatio = fs->make<TH1F>("caloEnergyRatio", "", 200, 0., 2.);
  caloEnergyRatioPos = fs->make<TH1F>("caloEnergyRatioPos", "", 200, 0., 2.);
  caloEnergyRatioNeg = fs->make<TH1F>("caloEnergyRatioNeg", "", 200, 0., 2.);

  histos_["hitsFraction_tk"]   = fs->make<TH1F>("hitsFraction_tk",   ";fraction of gen-matched hits in hgcTracks;", 101,0,1.1);
  histos_["hitsFraction_mcMIP"]   = fs->make<TH1F>("hitsFraction_mcMIP",   ";fraction of gen-matched hits in trackstersMIP;", 101,0,1.1);
  histos_["hitsFraction_mc"]   = fs->make<TH1F>("hitsFraction_mc",   ";fraction of gen-matched hits in tracksters2DCl;", 101,0,1.1);
  histos_["hitsFraction_all"]   = fs->make<TH1F>("hitsFraction_all",   ";fraction of gen-matched hits in all reco;", 101,0,1.1);
  histos_["hitsFraction_allV2"]   = fs->make<TH1F>("hitsFraction_allV2",   ";fraction of gen-matched hits in allV2 reco;", 101,0,1.1);
  histos_["hitsFraction_allV2_cpAsreco"]   = fs->make<TH1F>("hitsFraction_allV2_cpAsreco",   ";fraction of gen-matched hits in allV2 reco;", 101,0,1.1);

  histos_["hitsRecoOKFraction_tk"]   = fs->make<TH1F>("hitsRecoOKFraction_tk",   ";gen-matched hits over all reco in hgcTracks;", 101,0,1.1);
  histos_["hitsRecoOKFraction_mcMIP"]   = fs->make<TH1F>("hitsRecoOKFraction_mcMIP",   ";gen-matched hits over all reco in trackstersMIP;", 101,0,1.1);
  histos_["hitsRecoOKFraction_mc"]   = fs->make<TH1F>("hitsRecoOKFraction_mc",   ";gen-matched hits over all reco in tracksters2DCl;", 101,0,1.1);
  histos_["hitsRecoOKFraction_all"]   = fs->make<TH1F>("hitsRecoOKFraction_all",   ";gen-matched hits over all reco", 101,0,1.1);
  histos_["hitsRecoOKFraction_allV2"]   = fs->make<TH1F>("hitsRecoOKFraction_allV2",   ";gen-matched hits over all reco V2", 101,0,1.1);

  histos_["EFraction_tk"]   = fs->make<TH1F>("EFraction_tk",   ";fraction of gen-matched energy in hgcTracks;", 101,0,1.1);
  histos_["EFraction_mcMIP"]   = fs->make<TH1F>("EFraction_mcMIP",   ";fraction of gen-matched energy in trackstersMIP;", 101,0,1.1);
  histos_["EFraction_mc"]   = fs->make<TH1F>("EFraction_mc",   ";fraction of gen-matched energy in tracksters2DCl;", 101,0,1.1);
  histos_["EFraction_all"]   = fs->make<TH1F>("EFraction_all",   ";fraction of gen-matched energy in all reco;", 101,0,1.1);
  histos_["EFraction_allV2"]   = fs->make<TH1F>("EFraction_allV2",   ";fraction of gen-matched energy in allV2 reco;", 101,0,1.1);
  histos_["EFraction_allV2_cpAsreco"]   = fs->make<TH1F>("EFraction_allV2_cpAsreco",   ";fraction of gen-matched energy in allV2 reco;", 101,0,1.1);
  debug = false;
}



BuildParticles::~BuildParticles() { 
  // TFile *outF = TFile::Open("ticl_analysis_BuildP.root","RECREATE");
  // for(std::map<TString,TH1 *>::iterator it=histos_.begin();
  //     it!=histos_.end();
  //     it++)
  //   it->second->Write();
  // outF->Close();



}



void BuildParticles::beginRun(const edm::Run& run, 
                            const edm::EventSetup & es) { }


void  BuildParticles::analyze(const Event& iEvent, 
                            const EventSetup& iSetup) {



  recHitTools.getEventSetup(iSetup);
  finalParticles_etaPhi.clear();
  finalParticles_time.clear();
  tkParticles_time.clear();
  mcParticles_time.clear();
  mipParticles_time.clear();
  finalParticles.clear();

  if(debug)  std::cout << " \n \n new evt "<< std::endl;

  edm::Handle<std::vector<CaloParticle> > gpH;
  iEvent.getByToken( genToken_, gpH);

  // edm::Handle<std::vector<SimCluster> > scH;
  // iEvent.getByToken( simclusToken_, scH);

  edm::Handle<std::vector<reco::CaloCluster> > tkH;
  iEvent.getByToken( tkToken_, tkH);

  edm::Handle<std::vector<reco::HGCalMultiCluster> > mcMIPH, mcH;
  iEvent.getByToken( mcMIPToken_,mcMIPH);
  iEvent.getByToken( mcToken_,mcH);


  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;
  iEvent.getByToken(hits_eeToken_,ee_hits);
  iEvent.getByToken(hits_fhToken_,fh_hits);
  iEvent.getByToken(hits_bhToken_,bh_hits);


  iEvent.getByToken(layerToken_, layerClusters);
  iEvent.getByToken(timeMap_, time2DMap);

  std::map<uint32_t, const HGCRecHit*> hitmapPos;
  std::map<uint32_t, const HGCRecHit*> hitmapNeg;
  hitmap.clear();
  hitmapPos.clear();
  hitmapNeg.clear();
  for(auto const& it: *ee_hits){
    hitmap[it.detid().rawId()] = &it;
    if (recHitTools.getEta(it.detid().rawId()) > 0.) hitmapPos[it.detid().rawId()] = &it;
    else hitmapNeg[it.detid().rawId()] = &it;
  }
  for(auto const& it: *fh_hits){
    hitmap[it.detid().rawId()] = &it;
    if (recHitTools.getEta(it.detid().rawId()) > 0.) hitmapPos[it.detid().rawId()] = &it;
    else hitmapNeg[it.detid().rawId()] = &it;
  }
  for(auto const& it: *bh_hits){
    hitmap[it.detid().rawId()] = &it;
    if (recHitTools.getEta(it.detid().rawId()) > 0.) hitmapPos[it.detid().rawId()] = &it;
    else hitmapNeg[it.detid().rawId()] = &it;
  }

  //loop over reco tracks
  for(unsigned int ij=0; ij< tkH->size(); ++ij){
    const reco::CaloCluster cc = tkH->at(ij);
    const std::vector< std::pair<DetId, float> > &recHits = cc.hitsAndFractions();
    if(debug)std::cout << " recHits.size = " << recHits.size() << std::endl;
  }
   
  finalParticles.resize(tkH->size());
  finalParticles_etaPhi.resize(tkH->size());
  finalParticles_time.resize(tkH->size());
  tkParticles_time.resize(tkH->size());
  mcParticles_time.resize(mcH->size());
  mipParticles_time.resize(mcMIPH->size());
  if(debug)  std::cout << " pure tracks size = " << finalParticles.size() << std::endl;

  getTimeTrackedHitsList(-1, *tkH, tkParticles_time);
  getTimeClusteredHitsList(-1, *mcMIPH, mipParticles_time);
  getTimeClusteredHitsList(-1, *mcH, mcParticles_time);

  if(debug)  std::cout << " pre " << tkH->size() << " " << mcMIPH->size() << " " << mcH->size() << std::endl;
  getTrackedHitsList(-1, *tkH);
  getClusteredHitsList(-1, *mcMIPH);
  getClusteredHitsList(-1, *mcH);
  if(debug)  std::cout << " post " << tkH->size() << " " << mcMIPH->size() << " " << mcH->size() << std::endl;

  if(debug){
    std::cout << " >> finalParticles size = " << finalParticles.size() 
	      << " caloParticle size = " << gpH->size() << std::endl;
    std::cout << " loop over CP " << std::endl;
  }

  //  std::cout << " cp size = " << gpH->size() << std::endl;
  for(auto cp : *gpH) {
    
    float eta = cp.eta();
    float phi = cp.phi();
    float cpEnergy = 0;

    std::set<std::pair<uint32_t,float> > allHits;
    std::map<uint32_t,float> allHitsM;

    std::set<uint32_t> all_tk;
    std::set<uint32_t> all_mcMIP, all_mc;
    std::set<uint32_t> all_global;

    std::set<uint32_t> allMatched;
    std::set<uint32_t> allMatched_tk;
    std::set<uint32_t> allMatched_mcMIP,allMatched_mc;
    std::set<uint32_t> allMatched_global;

    float allMatched_E = 0.;
    float allMatched_tkE = 0.;
    float allMatched_mcMIPE = 0.;
    float allMatched_mcE = 0.;
    float allMatched_globalE = 0.;

    //iterate over all the attached sim clusters
    float genRecoenergy = 0.;
    float genRecoenergyPos = 0.;
    float genRecoenergyNeg = 0.;
    if(cp.simClusters().size() > 1) continue;
    float cpEnergyNofraction = 0;
    for(CaloParticle::sc_iterator scIt=cp.simCluster_begin(); scIt!=cp.simCluster_end(); scIt++){
      /*
      std::cout << " size = " << cp.simClusters().size() << "energy Sc = " << (*scIt)->energy() << std::endl;
      std::vector<std::pair<uint32_t,float> > hitsA = (*scIt)->hits_and_fractions();
      std::cout << " hits size = " << hitsA.size() << std::endl;
      for(auto ij : hitsA) {
	auto finder2 = hitmap.find(ij.first);
	if(finder2 != hitmap.end()){
        const HGCRecHit *hit = hitmap[ij.first];
        cpEnergy += hit->energy() * ij.second;
	cpEnergyNofraction += hit->energy();
	}
      }

      std::cout << " >>>> cpEnergy = " << cpEnergy << " cp.energy() = " << cp.energy() << " cpEnergyNofraction = " << cpEnergyNofraction
		<< " size = " << cp.simClusters().size() 
		<< std::endl;
      */

      //all hits and energy fractions at sim level
      std::vector<std::pair<uint32_t,float> > hits = (*scIt)->hits_and_fractions();
      //allHits.insert(allHits.end(),hits.begin(),hits.end());
      if(debug)  std::cout << " >>> hits.size() = " << hits.size() << std::endl;
      for(auto ij : hits) {
	//auto finder = allHitsM.find(ij.first);
	//if(finder != allHitsM.end()) std::cout << " previous fraction = " << allHitsM[ij.first] << " new = " << ij.second << std::endl;

	auto finder2 = hitmap.find(ij.first);
	if(finder2 != hitmap.end()){
	  const HGCRecHit *hit = hitmap[ij.first];
	  genRecoenergy += ij.second * hit->energy();
	  if(cp.eta() > 0.) {
	    auto finder2P = hitmapPos.find(ij.first);
	    if(finder2P != hitmapPos.end()){
	      const HGCRecHit *hitPos = hitmapPos[ij.first];
	      genRecoenergyPos += ij.second * hitPos->energy();
	    }
	  }
	  else{
	    auto finder2N = hitmapNeg.find(ij.first);
	    if(finder2N != hitmapNeg.end()){
	      const HGCRecHit *hitNeg = hitmapNeg[ij.first];
	      genRecoenergyNeg += ij.second * hitNeg->energy();
	    }
	  }
	} 
	allHitsM[ij.first] += ij.second;
	allHits.insert(std::pair<uint32_t,float>(ij.first, allHitsM[ij.first])); //save the cumulative fraction
	//if(!(allHitsM[ij.first] <= 1.)) std::cout << " total fraction > 1 " << allHitsM[ij.first] << std::endl;
      }
    }//sim clusters iterator

    caloEnergyRatio->Fill(genRecoenergy / cp.energy());
    if(cp.eta() > 0.) caloEnergyRatioPos->Fill(genRecoenergyPos / cp.energy());
    else caloEnergyRatioNeg->Fill(genRecoenergyNeg / cp.energy());

    // std::cout << " >>>> end simCP loop cpEnergy = " << cpEnergy << " cp.energy() = " << cp.energy() << " cpEnergyNofraction = " << cpEnergyNofraction
    // 	      << " size = " << cp.simClusters().size() 
    // 	      << " eta = " << cp.eta() << " phi = " << cp.phi() << std::endl;

 
    for(auto ij : allHits){
      auto finder2 = hitmap.find(ij.first);
      if(finder2 != hitmap.end()){
	const HGCRecHit *hit = hitmap[ij.first];
	//if(!(ij.second <= 1.)) std::cout << " fraction > 1 while computing energy " << ij.second << " hit->energy() = " << hit->energy() << std::endl;
	cpEnergy += hit->energy() ; // * ij.second;
	//	if(cpEnergy >= cp.energy()) std::cout << " cpEnergy becomes = " << cpEnergy << " cp.energy() =  " << cp.energy() << std::endl;
	// eta += recHitTools.getEta(ij.first) * hit->energy();
	// eta += recHitTools.getPhi(ij.first) * hit->energy();
      }
    } // vector of simHits with cumulative fraction for this given caloParticle
  
    // eta = eta/cpEnergy;
    // phi = phi/cpEnergy;

    //    if(cpEnergy > cp.energy()) cpEnergy = cp.energy();

    int bestReco = -1;
    float dR = 99;
    for(unsigned int ij=0; ij<finalParticles_etaPhi.size(); ++ij){
      float dRloc = reco::deltaR(eta, phi, finalParticles_etaPhi.at(ij)[0], finalParticles_etaPhi.at(ij)[1]);
      if(dRloc < dR){
	dR = dRloc;
	bestReco = ij;
      }

      if(debug){
	std::cout << " dRloc = " << dRloc << " ij = " << ij
		  << " eta = " << finalParticles_etaPhi.at(ij)[0]
		  << " phi = " << finalParticles_etaPhi.at(ij)[1] << std::endl;
      }
    }

    if(debug){
      std::cout << " caloparticle eta = " << eta << " phi = " << phi << " pdgId = " << cp.pdgId() << " energy = " << cpEnergy
		<< " best matched = " << bestReco << " dR = " << dR ;
      if(bestReco != -1) std::cout << " eta = " << finalParticles_etaPhi.at(bestReco)[0];
      if(bestReco != -1) std::cout << " phi = " << finalParticles_etaPhi.at(bestReco)[1];
      std::cout << " " << std::endl;
    }

    std::set<uint32_t> matched_all;
    if(bestReco != -1)
      matched_all = getMatchedHitsList(eta>0, allHits, finalParticles[bestReco], all_tk);
    std::set<uint32_t> matched_tk = getMatchedTrackedHitsList(eta>0, allHits, *tkH, all_tk);
    std::set<uint32_t> matched_mcMIP = getMatchedClusteredHitsList(eta>0, allHits, *mcMIPH, all_mcMIP);
    std::set<uint32_t> matched_mc = getMatchedClusteredHitsList(eta>0, allHits, *mcH, all_mc);

   if( bestReco == 0 && getMatched(allHits, finalParticles[bestReco]).size() < matched_tk.size()) {
     std::cout << " big problem " << std::endl; 
     std::cout << " matched_tk.size() = " << matched_tk.size() << std::endl; 
     std::cout << " getMatched(allHits, finalParticles[bestReco]).size() = " << getMatched(allHits, finalParticles[bestReco]).size() << std::endl; 
   }


    
    if(debug)    std::cout << " post matching conta E " << std::endl;
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
	if(fraction != 1 && debug) std::cout << " >> trks fraction " << fraction << std::endl;
	allMatched_tkE += hit->energy() * fraction;
	//allMatched_tkE += hit->energy();
      }
    }
    if(debug)    std::cout << " post matching conta E now MIP " << std::endl;
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
	if(fraction != 1 && debug) std::cout << " >> MCMIP fraction " << fraction << std::endl;
	allMatched_mcMIPE += hit->energy()*fraction;
	//allMatched_mcMIPE += hit->energy();
      }
    }
    if(debug)    std::cout << " post matching conta E now all " << std::endl;
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
	//if(fraction != 1 && debug) std::cout << " >> MC2d fraction " << fraction << std::endl;
	allMatched_mcE += hit->energy()*fraction;
	//allMatched_mcE += hit->energy();
      }
    }
    
    allMatched_globalE = allMatched_tkE + allMatched_mcMIPE + allMatched_mcE;

    if(debug)    std::cout << " post matching conta E now all V2 " << std::endl;
    for(auto ij : matched_all){
      allMatched.insert(ij);
      const HGCRecHit *hit = hitmap[ij];
      float fraction = 1.;
      for(auto kk : allHits){
	if(kk.first == ij) fraction = kk.second;
	break;
      }
      if(fraction != 1 && debug) std::cout << " >> MC2d fraction " << fraction << std::endl;
      allMatched_E += hit->energy()*fraction;
      //allMatched_E += hit->energy();
    }

   
    if(!(cpEnergy <= cp.energy())) std::cout << " AHHHHHHHHHHHHHHHHHH cpEnergy = " << cpEnergy << " cp.energy() = " << cp.energy() << std::endl;

    if(debug){
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
    }

    if(allHits.size() != 0){
    histos_["hitsFraction_tk"]->Fill(1.*allMatched_tk.size()/allHits.size());
    histos_["hitsFraction_mcMIP"]->Fill(1.*allMatched_mcMIP.size()/allHits.size());
    histos_["hitsFraction_mc"]->Fill(1.*allMatched_mc.size()/allHits.size());
    histos_["hitsFraction_all"]->Fill(1.*allMatched_global.size()/allHits.size());
    histos_["hitsFraction_allV2"]->Fill(1.*allMatched.size()/allHits.size());
    if(finalParticles.size() == gpH->size()) histos_["hitsFraction_allV2_cpAsreco"]->Fill(1.*allMatched.size()/allHits.size());
    if(debug && 1.*allMatched_mcMIP.size()/allHits.size() > 0.9) std::cout << " problem frac" <<  1.*allMatched_mcMIP.size()/allHits.size() << std::endl;
    }

    if(all_tk.size() != 0) histos_["hitsRecoOKFraction_tk"]->Fill(1.*allMatched_tk.size() / all_tk.size());
    if(all_mcMIP.size() != 0) histos_["hitsRecoOKFraction_mcMIP"]->Fill(1.*allMatched_mcMIP.size() / all_mcMIP.size());
    if(all_mc.size() != 0) histos_["hitsRecoOKFraction_mc"]->Fill(1.*allMatched_mc.size() / all_mc.size());
    if(matched_mc.size() != 0) histos_["hitsRecoOKFraction_allV2"]->Fill(1.*allMatched.size() / finalParticles[bestReco].size());
    
    histos_["EFraction_tk"]->Fill(1.*allMatched_tkE / cpEnergy);
    histos_["EFraction_mcMIP"]->Fill(1.*allMatched_mcMIPE / cpEnergy);
    histos_["EFraction_mc"]->Fill(1.*allMatched_mcE / cpEnergy);
    histos_["EFraction_all"]->Fill(1.*allMatched_globalE / cpEnergy);
    histos_["EFraction_allV2"]->Fill(1.*allMatched_E / cpEnergy);
    if(finalParticles.size() == gpH->size())     histos_["EFraction_allV2_cpAsreco"]->Fill(1.*allMatched_E / cpEnergy);
    if(debug && allMatched_globalE/cpEnergy > 1.) std::cout << " problem " << std::endl;

  }//caloParticles
}


//// add some get truckster time
void BuildParticles::getTimeClusteredHitsList(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs, std::vector<float>& times){

  std::cout << " >>> found n multiclusters = " << mcs.size()  << std::endl;
  for(auto mc : mcs) {

    std::vector<float> localTime;
    //loop over all the layer clusters in the multicluster
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){

      reco::CaloClusterPtr sClPtr(*it);
      float time2D = (*time2DMap)[sClPtr];
      if(time2D > -1) localTime.push_back(time2D);
      std::cout << " >>> void getTimeHitList mc time = " << time2D  << std::endl;
    }

    //use method from Shameena => try for the moment with truncation
    // can try truncation + weighted mean 
    float finalT = (localTime.size() >= 3) ? (hgcalsimclustertime::fixSizeHighestDensity(localTime)) : (-1.);
    times.push_back(finalT);
    std::cout << " >>> void getTimeHitList final time = " << finalT  << std::endl;
  }
  return;
}



void BuildParticles::getTimeTrackedHitsList(bool pos, const std::vector<reco::CaloCluster> &ccs, std::vector<float>& times){

  //vector => 1 position 1 track
  for(auto cc : ccs) {
    std::vector<float> localTime;

    //1 track = collection hits = caloCluster
    //compute time of caloCluster from hits

    const std::vector< std::pair<DetId, float> > &recHits =cc.hitsAndFractions();
    for(auto ij : recHits){

      const HGCRecHit *hit = hitmap[ij.first];
      float timeHit = hit->time();
      std::cout << " >>> void getTimeHitList trk hit time = " << timeHit  << std::endl;
      if(timeHit > -1) localTime.push_back(timeHit);
    }
    //se recHit singola da capire

    float finalT = (localTime.size() >= 3) ? (hgcalsimclustertime::fixSizeHighestDensity(localTime)) : (-1.);
    times.push_back(finalT);
    std::cout << " >>> void getTimeHitList final time = " << finalT  << std::endl;  
  }
  return;
}




////////////////////
void BuildParticles::getClusteredHitsList(bool pos,
					  const std::vector<reco::HGCalMultiCluster> &mcs){
  

  if(debug)  std::cout << " in getClusteredHitsList" << std::endl;
  for(auto mc : mcs) {
    
    if(debug)    std::cout << " >>> add multi cluster eta = " << mc.eta() << " phi = " << mc.phi() << std::endl;

    int assignToFinalParticle = -1;
    int counter = 0;
    std::set<uint32_t> testRh;

    //loop over all the layer clusters in the multicluster
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){
      const std::vector< std::pair<DetId, float> > &recHits = (*it)->hitsAndFractions();
      
      for(auto ij : recHits){
	if(ij.second != 0) {
	  testRh.insert(ij.first.rawId());
	}
      }
    }

    for(auto sw : finalParticles){
      std::set<uint32_t> imatches=getMatched(testRh, sw);
      if(imatches.size() > 0) {
	assignToFinalParticle = counter;
	break;
      }
      ++counter;
    }

    if(debug) std::cout << " matching found assignToFinalParticle = " << assignToFinalParticle << " presize = " << finalParticles[assignToFinalParticle].size() << std::endl;

    if(assignToFinalParticle != -1){
      for(auto ij : testRh) finalParticles[assignToFinalParticle].insert(ij);
    }
    else {
      finalParticles.push_back(testRh);
      std::vector<float> dummy;
      dummy.push_back(mc.eta());
      dummy.push_back(mc.phi());
      finalParticles_etaPhi.push_back(dummy);
    }
    
    if(debug && assignToFinalParticle != -1) std::cout << " loop into getClusteredHitsList size = " << finalParticles[assignToFinalParticle].size() << std::endl;    
    if(debug) std::cout << " loop into getClusteredHitsList size = " << finalParticles[finalParticles.size()-1].size() << std::endl;    
  }
  return;
}



void BuildParticles::getTrackedHitsList(bool pos,
					const std::vector<reco::CaloCluster> &ccs){
						      
  if(debug)  std::cout << " in getTrackedHitsList" << std::endl;
  int ipos = 0;
  for(auto cc : ccs) {

    if(debug)    std::cout << " >>> new trk cluster eta = " << cc.eta() << " phi = " << cc.phi() << std::endl;

    const std::vector< std::pair<DetId, float> > &recHits =cc.hitsAndFractions();
    for(auto ij : recHits){
      if(ij.second != 0) 
	finalParticles[ipos].insert(ij.first.rawId());
    }

    finalParticles_etaPhi[ipos].push_back(cc.eta());
    finalParticles_etaPhi[ipos].push_back(cc.phi());

    if(debug){
    std::cout << " >>> new trk cluster eta = " << cc.eta() << " phi = " << cc.phi() << " and "
	      << " finalParticles_etaPhi[ipos][0] = " << finalParticles_etaPhi[ipos][0] << " " << finalParticles_etaPhi[ipos][1] << std::endl;
    }

    ++ipos;

    if(debug) std::cout << " loop into getTrackedHitsList size = " << finalParticles[ipos-1].size() << std::endl;    
  }
  
  return;
}





std::set<uint32_t> BuildParticles::getMatchedClusteredHitsList(bool pos,
							       const std::set<std::pair<uint32_t,float> > &hits,
							       const std::vector<reco::HGCalMultiCluster> &mcs, 
							       std::set<uint32_t> &mcHits) {
  

  if(debug)  std::cout << " in getMatchedClusteredHitsList" << std::endl;

  std::set<uint32_t> matchedList;
  
  for(auto mc : mcs) {
    
    //require on the same side
    bool mcPos(mc.eta()>0);
    if(mcPos!=pos) continue;
    
    //loop over all the layer clusters in the multicluster
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){
      const std::vector< std::pair<DetId, float> > &recHits = (*it)->hitsAndFractions();
      std::set<uint32_t> testRh;
      for(auto ij : recHits){
	if(ij.second != 0) 
	  mcHits.insert(ij.first.rawId());

      //std::vector<uint32_t> imatches=getMatched(hits,recHits);
      //matchedList.insert(matchedList.end(), imatches.begin(), imatches.end());
      std::set<uint32_t> imatches=getMatched(hits, mcHits);
      for(auto ij : imatches) matchedList.insert(ij);

      }
    }
  }

  if(debug)  std::cout << " looping in  getMatchedClusteredHitsList matchedList.size() = " <<  matchedList.size() << std::endl;
  return matchedList;
}



std::set<uint32_t> BuildParticles::getMatchedHitsList(bool pos,
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
std::set<uint32_t> BuildParticles::getMatchedTrackedHitsList(bool pos,
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

//
std::vector<uint32_t> BuildParticles::getMatched(const std::set<std::pair<uint32_t,float> > &a,
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


std::set<uint32_t> BuildParticles::getMatched(const std::set<std::pair<uint32_t,float> > &a,
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


std::set<uint32_t> BuildParticles::getMatched(const std::set<uint32_t> &a,
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



DEFINE_FWK_MODULE(BuildParticles);
