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

class BuildParticles : public edm::EDAnalyzer {
public:

  explicit BuildParticles(const edm::ParameterSet&);

  ~BuildParticles();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

private:

  hgcal::RecHitTools         recHitTools;

  std::vector<reco::HGCalMultiCluster> buildTrackerItSameSeed(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs);

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


  edm::EDGetTokenT<reco::CaloClusterCollection> layerToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> timeMap_;

  edm::Handle<reco::CaloClusterCollection> layerClusters;  
  edm::Handle<edm::ValueMap<float> > time2DMap;
  std::map<uint32_t, const HGCRecHit*> hitmap;
  std::map<uint32_t, const HGCRecHit*> hitmap_below37;

  std::map<TString,TH1 *> histos_;
  TH1F* caloEnergyRatio;
  TH1F* caloEnergyRatioPos;
  TH1F* caloEnergyRatioNeg;
  TH1F* deltaTtrack;

  TH2F* recoLayer2D_RvsL;
  TH2F* recoLayer2D_RvsZ;
  TH2F* recoLayer2D_YvsX_pos;
  TH2F* recoLayer2D_YvsX_neg;

  TH2F* recoRecHit_EvsL_neg;
  TH2F* recoLayer2D_EvsL_neg;
  TH2F* recoSimHit_EvsL_neg;
  TH2F* recoRecHit_EvsL_pos;
  TH2F* recoLayer2D_EvsL_pos;
  TH2F* recoSimHit_EvsL_pos;

  TH1F* recoNOTAll_genAll_RatioEnergy_below37;
  TH2F* nonRecoHits_Rvslayer_below37;
  TH1F* recoAll_genAll_RatioEnergy_below37;
  TH2F* recoHits_Rvslayer_below37;

  TH2F* recoHits_YvsX;
  TH2F* recoHits_YvsX_below37;
  TH2F* nonRecoHits_YvsX;
  TH2F* nonRecoHits_YvsX_below37;

  TH1F* recoNOTAll_genAll_RatioEnergy;
  TH2F* nonRecoHits_Rvslayer;
  TH2F* recoHits_Rvslayer;
  TH2F* nonRecoHits_Rvslayer_pos;
  TH2F* recoHits_Rvslayer_pos;
  TH2F* nonRecoHits_Rvslayer_neg;
  TH2F* recoHits_Rvslayer_neg;

  TH1F* recoAll_genAll_RatioEnergy;
  TH2F* ratioEnergy_vsEta;
  TH1F* genEnergy_forRatioEnergy;
  TH1F* recoEnergy_forRatioEnergy;
  TH1F* genEnergy_forRatioEnergyBelow0p6;
  TH1F* recoEnergy_forRatioEnergyBelow0p6;
  TH2F* genEnergy_vsEta;
  TH2F* recoEnergy_vsEta;
  TH2F* genVsrecoEnergy;

  std::vector<std::set<uint32_t>> finalParticles;
  std::vector<std::vector<float>> finalParticles_etaPhi;
  std::vector<float> finalParticles_time;
  std::vector<float> tkParticles_time; //can be both from hgcal-Tracks and TrkIteration
  std::vector<float> mcParticles_time;
  std::vector<float> mipParticles_time;
  bool debug;

  bool isTrkIteration;
  int nEvent;
};



BuildParticles::BuildParticles(const edm::ParameterSet& iConfig) :
  genToken_(consumes<std::vector<CaloParticle> >(edm::InputTag("mix:MergedCaloTruth"))),
  tkToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcTracks::RECO2"))),
  mcTrkToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiClusterTrk:TrkMultiClustersFromTracksterByCA"))),
  mcMIPToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiClusterMIP:MIPMultiClustersFromTracksterByCA"))),
  mcToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("TrackstersToMultiCluster:MultiClustersFromTracksterByCA"))),
  hits_eeToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"))),
  hits_fhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"))),
  hits_bhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"))),
  layerToken_(consumes<reco::CaloClusterCollection>(edm::InputTag("hgcalLayerClusters::RECO"))),
  timeMap_(consumes<edm::ValueMap<float>>(edm::InputTag("hgcalLayerClusters:timeLayerCluster")))
{
  //book some histos here                                                                                                                 

  edm::Service<TFileService> fs;

  recoRecHit_EvsL_neg = fs->make<TH2F>("recoRecHit_EvsL_neg", "", 60, 0., 60., 200, 0, 100.);
  recoLayer2D_EvsL_neg = fs->make<TH2F>("recoLayer2D_EvsL_neg", "", 60, 0., 60., 200, 0, 100.);
  recoSimHit_EvsL_neg = fs->make<TH2F>("recoSimHit_EvsL_neg", "", 60, 0., 60., 200, 0, 100.);
  recoRecHit_EvsL_pos = fs->make<TH2F>("recoRecHit_EvsL_pos", "", 60, 0., 60., 200, 0, 100.);
  recoLayer2D_EvsL_pos = fs->make<TH2F>("recoLayer2D_EvsL_pos", "", 60, 0., 60., 200, 0, 100.);
  recoSimHit_EvsL_pos = fs->make<TH2F>("recoSimHit_EvsL_pos", "", 60, 0., 60., 200, 0, 100.);

  recoLayer2D_RvsL = fs->make<TH2F>("recoLayer2D_RvsL", "", 60, 0., 60., 300, 0, 300.);
  recoLayer2D_RvsZ = fs->make<TH2F>("recoLayer2D_RvsZ", "", 1200, -600., 600., 300, 0, 300.);
  recoLayer2D_YvsX_pos = fs->make<TH2F>("recoLayer2D_YvsX_pos", "", 600, -300., 300., 600, -300, 300.);
  recoLayer2D_YvsX_neg = fs->make<TH2F>("recoLayer2D_YvsX_neg", "", 600, -300., 300., 600, -300, 300.);


  recoNOTAll_genAll_RatioEnergy_below37 = fs->make<TH1F>("NOTrecoAll_genAll_RatioEnergy_below37", "", 101,0,1.1);
  nonRecoHits_Rvslayer_below37 = fs->make<TH2F>("nonRecoHits_Rvslayer_below37", "", 60, 0., 60., 300, 0, 300.);
  nonRecoHits_YvsX_below37 = fs->make<TH2F>("nonRecoHits_YvsX_below37", "", 600, -300., 300., 600, -300, 300.);
  recoAll_genAll_RatioEnergy_below37 = fs->make<TH1F>("recoAll_genAll_RatioEnergy_below37", "", 101,0,1.1);
  recoHits_Rvslayer_below37 = fs->make<TH2F>("recoHits_Rvslayer_below37", "", 60, 0., 60., 300, 0, 300.);
  recoHits_YvsX_below37 = fs->make<TH2F>("recoHits_YvsX_below37", "", 600, -300., 300., 600, -300, 300.);


  recoNOTAll_genAll_RatioEnergy = fs->make<TH1F>("NOTrecoAll_genAll_RatioEnergy", "", 101,0,1.1);
  nonRecoHits_Rvslayer = fs->make<TH2F>("nonRecoHits_Rvslayer", "", 60, 0., 60., 300, 0, 300.);
  nonRecoHits_YvsX = fs->make<TH2F>("nonRecoHits_YvsX", "", 600, -300., 300., 600, -300, 300.);
  recoHits_Rvslayer = fs->make<TH2F>("recoHits_Rvslayer", "", 60, 0., 60., 300, 0, 300.);
  recoHits_YvsX = fs->make<TH2F>("recoHits_YvsX", "", 600, -300., 300., 600, -300, 300.);

  nonRecoHits_Rvslayer_pos = fs->make<TH2F>("nonRecoHits_Rvslayer_pos", "", 60, 0., 60., 300, 0, 300.);
  recoHits_Rvslayer_pos = fs->make<TH2F>("recoHits_Rvslayer_pos", "", 60, 0., 60., 300, 0, 300.);
  nonRecoHits_Rvslayer_neg = fs->make<TH2F>("nonRecoHits_Rvslayer_neg", "", 60, 0., 60., 300, 0, 300.);
  recoHits_Rvslayer_neg = fs->make<TH2F>("recoHits_Rvslayer_neg", "", 60, 0., 60., 300, 0, 300.);

  recoAll_genAll_RatioEnergy = fs->make<TH1F>("recoAll_genAll_RatioEnergy", "", 101,0,1.1);
  ratioEnergy_vsEta = fs->make<TH2F>("ratioEnergy_vsEta", "", 6, -3., 3., 101, 0., 1.1);
  genEnergy_forRatioEnergy = fs->make<TH1F>("genEnergy_forRatioEnergy", "", 200, 0., 200);
  recoEnergy_forRatioEnergy = fs->make<TH1F>("recoEnergy_forRatioEnergy", "", 200, 0., 200);

  genEnergy_forRatioEnergyBelow0p6 = fs->make<TH1F>("genEnergy_forRatioEnergyBelow0p6", "", 200, 0., 200);
  recoEnergy_forRatioEnergyBelow0p6 = fs->make<TH1F>("recoEnergy_forRatioEnergyBelow0p6", "", 200, 0., 200);

  genEnergy_vsEta = fs->make<TH2F>("genEnergy_vsEta", "", 6, -3., 3., 200, 0., 200);
  recoEnergy_vsEta = fs->make<TH2F>("recoEnergy_vsEta", "", 6, -3., 3., 200, 0., 200);
  genVsrecoEnergy = fs->make<TH2F>("genVsrecoEnergy", "", 200, 0., 200., 200, 0., 200.);

  caloEnergyRatio = fs->make<TH1F>("caloEnergyRatio", "", 200, 0., 2.);
  caloEnergyRatioPos = fs->make<TH1F>("caloEnergyRatioPos", "", 200, 0., 2.);
  caloEnergyRatioNeg = fs->make<TH1F>("caloEnergyRatioNeg", "", 200, 0., 2.);

  deltaTtrack = fs->make<TH1F>("deltaTtrack", "", 1000, 0., 25);

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
  isTrkIteration = false;
  nEvent = 0;
}



BuildParticles::~BuildParticles() { 
}



void BuildParticles::beginRun(const edm::Run& run, 
                            const edm::EventSetup & es) { }


void  BuildParticles::analyze(const Event& iEvent, 
                            const EventSetup& iSetup) {

  ++nEvent;
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

  edm::Handle<std::vector<reco::CaloCluster> > tkH;
  if(!isTrkIteration)  iEvent.getByToken( tkToken_, tkH);

  auto mcTrkH = std::make_unique<std::vector<reco::HGCalMultiCluster>>();
  edm::Handle<std::vector<reco::HGCalMultiCluster> > mcTrkH_co, mcMIPH, mcH;
  if(isTrkIteration) iEvent.getByToken(mcTrkToken_, mcTrkH_co);
  iEvent.getByToken( mcMIPToken_,mcMIPH);
  iEvent.getByToken( mcToken_,mcH);

  if(isTrkIteration)  *mcTrkH = buildTrackerItSameSeed(-1, *mcTrkH_co);
  if(debug && isTrkIteration){
    std::cout << " mcTrkH->size() = " << mcTrkH->size() << std::endl;
    for(auto mc : *mcTrkH) {
      std::cout << " eta = " << mc.eta() << std::endl;
    }
  }


  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;
  iEvent.getByToken(hits_eeToken_,ee_hits);
  iEvent.getByToken(hits_fhToken_,fh_hits);
  iEvent.getByToken(hits_bhToken_,bh_hits);


  iEvent.getByToken(layerToken_, layerClusters);
  iEvent.getByToken(timeMap_, time2DMap);

  std::vector<float> E2Dcl_pos(60, 0);
  std::vector<float> E2Dcl_neg(60, 0);
  for(const reco::CaloCluster &sCl : *layerClusters){
    float clX = sCl.x();
    float clY = sCl.y();
    float clZ = sCl.z();

    const HGCalDetId clSeedId = sCl.hitsAndFractions()[0].first;
    int rhSeedL = recHitTools.getLayerWithOffset(clSeedId);

    float Radius = sqrt(clX*clX + clY*clY);
    recoLayer2D_RvsL->Fill(rhSeedL, Radius);
    recoLayer2D_RvsZ->Fill(clZ, Radius);

    if(clZ > 0) recoLayer2D_YvsX_pos->Fill(clY, clX);
    else recoLayer2D_YvsX_neg->Fill(clY, clX);

    if(clZ > 0) E2Dcl_pos[rhSeedL] += sCl.energy();
    else E2Dcl_neg[rhSeedL] += sCl.energy();
  }


  std::map<uint32_t, const HGCRecHit*> hitmapPos;
  std::map<uint32_t, const HGCRecHit*> hitmapNeg;
  hitmap_below37.clear();
  hitmap.clear();
  hitmapPos.clear();
  hitmapNeg.clear();

  std::vector<float> ERH_pos(60, 0);
  std::vector<float> ERH_neg(60, 0);
  for(auto const& it: *ee_hits){
    hitmap[it.detid().rawId()] = &it;

    int rhL = recHitTools.getLayerWithOffset(it.detid().rawId());
    if(rhL < 37) hitmap_below37[it.detid().rawId()] = &it;

    float rhZ = recHitTools.getPosition(it.detid().rawId()).z();
    if(rhZ > 0) ERH_pos[rhL] += it.energy();
    else ERH_neg[rhL] += it.energy();

    if (recHitTools.getEta(it.detid().rawId()) > 0.) hitmapPos[it.detid().rawId()] = &it;
    else hitmapNeg[it.detid().rawId()] = &it;
  }
  for(auto const& it: *fh_hits){
    hitmap[it.detid().rawId()] = &it;

    int rhL = recHitTools.getLayerWithOffset(it.detid().rawId());
    if(rhL < 37) hitmap_below37[it.detid().rawId()] = &it;

    float rhZ = recHitTools.getPosition(it.detid().rawId()).z();
    if(rhZ > 0) ERH_pos[rhL] += it.energy();
    else ERH_neg[rhL] += it.energy();

    if (recHitTools.getEta(it.detid().rawId()) > 0.) hitmapPos[it.detid().rawId()] = &it;
    else hitmapNeg[it.detid().rawId()] = &it;
  }
  //RA excluding BH
  /*
  for(auto const& it: *bh_hits){
    hitmap[it.detid().rawId()] = &it;
    if (recHitTools.getEta(it.detid().rawId()) > 0.) hitmapPos[it.detid().rawId()] = &it;
    else hitmapNeg[it.detid().rawId()] = &it;
  }
  */
  //loop over reco tracks
  /*
  for(unsigned int ij=0; ij< tkH->size(); ++ij){
    const reco::CaloCluster cc = tkH->at(ij);
    const std::vector< std::pair<DetId, float> > &recHits = cc.hitsAndFractions();
    if(debug)std::cout << " recHits.size = " << recHits.size() << std::endl;
  }
  */
  //  return ;

  if(!isTrkIteration){ 
    finalParticles.resize(tkH->size());
    finalParticles_etaPhi.resize(tkH->size());
    tkParticles_time.resize(tkH->size());
    finalParticles_time.resize(tkH->size());
  }
  else{
    finalParticles.resize(mcTrkH->size());
    finalParticles_etaPhi.resize(mcTrkH->size());
    tkParticles_time.resize(mcTrkH->size());
    finalParticles_time.resize(mcTrkH->size());
  }
  mcParticles_time.resize(mcH->size());
  mipParticles_time.resize(mcMIPH->size());
  if(debug)  std::cout << " pure tracks size = " << finalParticles.size() << std::endl;

  if(!isTrkIteration)
    getTimeTrackedHitsList(-1, *tkH, tkParticles_time);
  else getTimeClusteredHitsList(-1, *mcTrkH, tkParticles_time);
  getTimeClusteredHitsList(-1, *mcMIPH, mipParticles_time);
  getTimeClusteredHitsList(-1, *mcH, mcParticles_time);

  if(debug && !isTrkIteration)  std::cout << " pre " << tkH->size()  << " " << mcMIPH->size() << " " << mcH->size() << std::endl;
  if(debug && isTrkIteration)  std::cout << " pre " << mcTrkH->size() << " " << mcMIPH->size() << " " << mcH->size() << std::endl;
  if(!isTrkIteration)
    getTrackedHitsList(-1, *tkH);
  else getTrackedItHitsList(-1, *mcTrkH);
  getClusteredHitsList(1, *mcMIPH);
  getClusteredHitsList(0, *mcH);
  if(debug && !isTrkIteration)  std::cout << " post " << tkH->size()  << " " << mcMIPH->size() << " " << mcH->size() << std::endl;
  if(debug && isTrkIteration)  std::cout << " post " << mcTrkH->size() << " " << mcMIPH->size() << " " << mcH->size() << std::endl;

  //  if(debug)  std::cout << " post " << tkH->size() << " " << mcTrkH->size() << " " << mcMIPH->size() << " " << mcH->size() << std::endl;

  if(debug){
    std::cout << " >> finalParticles size = " << finalParticles.size() 
	      << " caloParticle size = " << gpH->size() << std::endl;
    std::cout << " loop over CP " << std::endl;
  }


  for(auto cp : *gpH) {
    //    if(cp.simClusters().size() > 1) continue;
    
    std::vector<float> ESC_pos(60, 0);
    std::vector<float> ESC_neg(60, 0);

    float eta = cp.eta();
    float phi = cp.phi();
    float cpEnergy = 0;
    float cpEnergy_below37 = 0;

    //all caloP hits
    std::set<uint32_t> all_simHits_below37;
    std::set<uint32_t> all_simHits;
    std::set<std::pair<uint32_t,float> > allHits;
    std::map<uint32_t,float> allHitsM;

    //all hits in the single collection matched to caloP
    std::set<uint32_t> all_recHits_below37;
    std::set<uint32_t> all_recHits;
    std::set<uint32_t> all_tk;
    std::set<uint32_t> all_tkBP;
    std::set<uint32_t> all_mcMIP, all_mc;

    //all hits in finalP or single collection matched to caloP
    std::set<uint32_t> matched_all;
    std::set<uint32_t> matched_tk;
    std::set<uint32_t> matched_mcMIP; 
    std::set<uint32_t> matched_mc;

    //as matched_all but summing among the single collections
    std::set<uint32_t> allMatched_global;

    //as previous but counting the corresponding energy
    float allMatched_E = 0.;
    float allMatched_tkE = 0.;
    float allMatched_mcMIPE = 0.;
    float allMatched_mcE = 0.;
    float allMatched_globalE = 0.;

    //iterate over all the attached sim clusters
    float genRecoenergy = 0.;
    float genRecoenergyPos = 0.;
    float genRecoenergyNeg = 0.;   
    float cpEnergyNofraction = 0;


    for(CaloParticle::sc_iterator scIt=cp.simCluster_begin(); scIt!=cp.simCluster_end(); scIt++){

      //all hits and energy fractions at sim level
      std::vector<std::pair<uint32_t,float> > hits = (*scIt)->hits_and_fractions();

      if(debug)  std::cout << " >>> hits.size() = " << hits.size() << std::endl;
      for(auto ij : hits) {

	auto finder37 = hitmap_below37.find(ij.first);
        if(finder37 != hitmap_below37.end()){
	  all_simHits_below37.insert(ij.first);
	}

	auto finder2 = hitmap.find(ij.first);
	if(finder2 != hitmap.end()){

	  int rhL = recHitTools.getLayerWithOffset(ij.first);
	  float rhZ = recHitTools.getPosition(ij.first).z();
	  if(rhZ > 0) ESC_pos[rhL] += hitmap[ij.first]->energy();
	  else ESC_neg[rhL] += hitmap[ij.first]->energy();

	  all_simHits.insert(ij.first);

	  //	  const HGCRecHit *hit = hitmap[ij.first];
	  genRecoenergy += ij.second; // * hit->energy();
	  /*
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
	  */
	  allHitsM[ij.first] += ij.second;
	} 
      }
    }//sim clusters iterator
  

    for(unsigned int ij=1; ij<E2Dcl_pos.size(); ++ij)
      if(E2Dcl_pos[ij] != 0)
	recoLayer2D_EvsL_pos->Fill(ij, E2Dcl_pos[ij]);
    for(unsigned int ij=1; ij<ERH_pos.size(); ++ij)
            if(ERH_pos[ij] != 0) recoRecHit_EvsL_pos->Fill(ij, ERH_pos[ij]);
    for(unsigned int ij=1; ij<ESC_pos.size(); ++ij)
      if(ESC_pos[ij] != 0)
	recoSimHit_EvsL_pos->Fill(ij, ESC_pos[ij]);
    
    for(unsigned int ij=1; ij<E2Dcl_neg.size(); ++ij)
      if(E2Dcl_neg[ij] != 0)      
	recoLayer2D_EvsL_neg->Fill(ij, E2Dcl_neg[ij]);
    for(unsigned int ij=1; ij<ERH_neg.size(); ++ij)
      if(ERH_neg[ij] != 0)
	recoRecHit_EvsL_neg->Fill(ij, ERH_neg[ij]);
    for(unsigned int ij=1; ij<ESC_neg.size(); ++ij)
      if(ESC_neg[ij] != 0) recoSimHit_EvsL_neg->Fill(ij, ESC_neg[ij]);
    
    
    //    return;
    
    caloEnergyRatio->Fill(genRecoenergy / cp.energy());
    
    
    for(auto ij : all_simHits_below37){
      auto finder2 = hitmap_below37.find(ij);
      if(finder2 != hitmap_below37.end()){
        const HGCRecHit *hit = hitmap_below37[ij];
	cpEnergy_below37 += hit->energy();
      }
    }



    for(auto ij : all_simHits){
      allHits.insert(std::pair<uint32_t,float>(ij, 1));
      auto finder2 = hitmap.find(ij);
      if(finder2 != hitmap.end()){
        const HGCRecHit *hit = hitmap[ij];
	cpEnergy += hit->energy();
      }
    }

    /*
    for(auto ij : allHitsM){
      auto finder2 = hitmap.find(ij.first);
      if(finder2 != hitmap.end()){
	const HGCRecHit *hit = hitmap[ij.first];

	cpEnergy += hit->energy() ; // * ij.second;

	// eta += recHitTools.getEta(ij.first) * hit->energy();
	// eta += recHitTools.getPhi(ij.first) * hit->energy();
	allHits.insert(std::pair<uint32_t,float>(ij.first, ij.second)); 
      }
    } // vector of simHits with cumulative fraction for this given caloParticle
    // eta = eta/cpEnergy;
    // phi = phi/cpEnergy;
    */

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

    //all IDs in finalParticles that are matched to allHits = from caloP
    if(bestReco != -1)
      matched_all = getMatchedHitsList(eta>0, allHits, finalParticles[bestReco], all_tkBP);

    if(!isTrkIteration) matched_tk = getMatchedTrackedHitsList(eta>0, allHits, *tkH, all_tk);
    else matched_tk = getMatchedClusteredHitsList(eta>0, allHits, *mcTrkH, all_tk);
    matched_mcMIP = getMatchedClusteredHitsList(eta>0, allHits, *mcMIPH, all_mcMIP);
    matched_mc = getMatchedClusteredHitsList(eta>0, allHits, *mcH, all_mc);


    /*
    //just take the inclusive of all recHits    
    for(auto ij:all_tk){
      all_recHits.insert(ij);
      if(recHitTools.getLayerWithOffset(ij) < 37) all_recHits_below37.insert(ij);
    }
    for(auto ij:all_mcMIP){
      all_recHits.insert(ij);
      if(recHitTools.getLayerWithOffset(ij) < 37) all_recHits_below37.insert(ij);
    }
    for(auto ij:all_mc){
      all_recHits.insert(ij);
      if(recHitTools.getLayerWithOffset(ij) < 37) all_recHits_below37.insert(ij);
    }


    std::set<uint32_t> all_recHits_matchedToSim = getMatched(all_recHits, all_simHits);
    std::set<uint32_t> all_recHits_NOTmatchedToSim = getNOTMatched(all_simHits, all_recHits);
    // std::set<uint32_t>::iterator it_all;
    // it_all = std::set_difference(all_simHits.begin(), all_simHits.end(), all_recHits.begin(), all_recHits.end(), all_recHits_NOTmatchedToSim.begin());
    // all_recHits_NOTmatchedToSim.resize(it_all - all_recHits_NOTmatchedToSim.begin());


    std::set<uint32_t> all_recHits_matchedToSim_below37 = getMatched(all_recHits_below37, all_simHits_below37);
    std::set<uint32_t> all_recHits_NOTmatchedToSim_below37 = getNOTMatched(all_simHits_below37, all_recHits_below37);

    /////////
    float NOTrecoAllEnergy_below37 = 0.;
    for(auto ij : all_recHits_NOTmatchedToSim_below37){
      //const HGCalDetId hitid = all_recHits_NOTmatchedToSim.at(ij);
      //      const HGCalDetId hitid = ij;
      auto hitid = ij;
      auto finder2 = hitmap_below37.find(hitid);
      if(finder2 != hitmap_below37.end()){
	
        const HGCRecHit *hit = hitmap_below37[hitid];
        NOTrecoAllEnergy_below37 += hit->energy();
	int rhL = recHitTools.getLayerWithOffset(hitid);
	float rhX = recHitTools.getPosition(hitid).x();
	float rhY = recHitTools.getPosition(hitid).y();
	float radius = sqrt(rhX*rhX + rhY*rhY);
	nonRecoHits_Rvslayer_below37->Fill(rhL, radius);
	nonRecoHits_YvsX_below37->Fill(rhY, rhX);
      }
    }
    if(NOTrecoAllEnergy_below37 < 0) std::cout << " nEvent = " << nEvent << " NOTrecoAllEnergy_below37 = " << NOTrecoAllEnergy_below37 
					       << " cpEnergy_below37 = " << cpEnergy_below37 << std::endl;

    recoNOTAll_genAll_RatioEnergy_below37->Fill(NOTrecoAllEnergy_below37/cpEnergy_below37);

    float recoAllEnergy_below37 = 0.;
    for(auto ij : all_recHits_matchedToSim_below37){
      auto finder2 = hitmap_below37.find(ij);
      if(finder2 != hitmap_below37.end()){
        const HGCRecHit *hit = hitmap_below37[ij];
        recoAllEnergy_below37 += hit->energy();
	int rhL = recHitTools.getLayerWithOffset(ij);
	float rhX = recHitTools.getPosition(ij).x();
	float rhY = recHitTools.getPosition(ij).y();
	float radius = sqrt(rhX*rhX + rhY*rhY);
	recoHits_Rvslayer_below37->Fill(rhL, radius);
	recoHits_YvsX_below37->Fill(rhY, rhX);
      }
    }
    recoAll_genAll_RatioEnergy_below37->Fill(recoAllEnergy_below37/cpEnergy_below37);

    if(recoAllEnergy_below37 < 0) std::cout << " nEvent = " << nEvent << " recoAllEnergy_below37 = " << recoAllEnergy_below37
					    << " cpEnergy_below37 = " << cpEnergy_below37 << std::endl;

    /////////

    float NOTrecoAllEnergy = 0.;
    for(auto ij : all_recHits_NOTmatchedToSim){
      //const HGCalDetId hitid = all_recHits_NOTmatchedToSim.at(ij);
      //      const HGCalDetId hitid = ij;
      auto hitid = ij;
      auto finder2 = hitmap.find(hitid);
      if(finder2 != hitmap.end()){
        const HGCRecHit *hit = hitmap[hitid];
        NOTrecoAllEnergy += hit->energy();
	int rhL = recHitTools.getLayerWithOffset(hitid);
	float rhX = recHitTools.getPosition(hitid).x();
	float rhY = recHitTools.getPosition(hitid).y();
	float radius = sqrt(rhX*rhX + rhY*rhY);
	nonRecoHits_Rvslayer->Fill(rhL, radius);
	nonRecoHits_YvsX->Fill(rhY, rhX);
	if(recHitTools.getPosition(hitid).z() > 0) nonRecoHits_Rvslayer_pos->Fill(rhL, radius);
	else nonRecoHits_Rvslayer_neg->Fill(rhL, radius);
      }
    }
    recoNOTAll_genAll_RatioEnergy->Fill(NOTrecoAllEnergy/cpEnergy);
*/



    float recoAllEnergy = 0.;
    //for(auto ij : all_recHits_matchedToSim){
    for(auto ij : matched_all){
      auto finder2 = hitmap.find(ij);
      if(finder2 != hitmap.end()){
        const HGCRecHit *hit = hitmap[ij];
        recoAllEnergy += hit->energy();
	int rhL = recHitTools.getLayerWithOffset(ij);
	float rhX = recHitTools.getPosition(ij).x();
	float rhY = recHitTools.getPosition(ij).y();
	float radius = sqrt(rhX*rhX + rhY*rhY);
	recoHits_Rvslayer->Fill(rhL, radius);
	recoHits_YvsX->Fill(rhY, rhX);
	if(recHitTools.getPosition(ij).z() > 0) recoHits_Rvslayer_pos->Fill(rhL, radius);
	else recoHits_Rvslayer_neg->Fill(rhL, radius);

      }
    }

    recoAll_genAll_RatioEnergy->Fill(recoAllEnergy/cpEnergy);
    ratioEnergy_vsEta->Fill(cp.eta(), recoAllEnergy/cpEnergy);
    genEnergy_forRatioEnergy->Fill(cpEnergy);
    recoEnergy_forRatioEnergy->Fill(recoAllEnergy);
    if(recoAllEnergy/cpEnergy < 0.6){
      genEnergy_forRatioEnergyBelow0p6->Fill(cpEnergy);
      recoEnergy_forRatioEnergyBelow0p6->Fill(recoAllEnergy);
    }
    genEnergy_vsEta->Fill(cp.eta(), cpEnergy);
    recoEnergy_vsEta->Fill(cp.eta(), recoAllEnergy);
    genVsrecoEnergy->Fill(recoAllEnergy, cpEnergy);


   if( bestReco == 0 && getMatched(allHits, finalParticles[bestReco]).size() < matched_tk.size()) {
     std::cout << " big problem " << std::endl; 
     std::cout << " matched_tk.size() = " << matched_tk.size() << std::endl; 
     std::cout << " getMatched(allHits, finalParticles[bestReco]).size() = " << getMatched(allHits, finalParticles[bestReco]).size() << std::endl; 
   }


    
    if(debug)    std::cout << " post matching conta E " << std::endl;

    for(auto ij : matched_tk) {
      
      //in principle every reco id exists in the collection of recHits                                                                 
      const HGCRecHit *hit = hitmap[ij];
      unsigned int preSize = allMatched_global.size();
      allMatched_global.insert(ij);
      if(allMatched_global.size() != preSize){
	/*
	  float fraction = 1.;
	  for(auto kk : allHits){
	  if(kk.first == ij) fraction = kk.second;
	  break;
	  }
	  if(fraction != 1 && debug) std::cout << " >> trks fraction " << fraction << std::endl;
	  allMatched_tkE += hit->energy() * fraction;
	*/
	allMatched_tkE += hit->energy();
      }
    }
    if(debug)    std::cout << " post matching conta E now MIP " << std::endl;
    for(auto ij : matched_mcMIP){

      const HGCRecHit *hit = hitmap[ij];
      unsigned int preSize = allMatched_global.size();
      allMatched_global.insert(ij);
      if(allMatched_global.size() != preSize){
	/*
	 float fraction = 1.;
	 for(auto kk : allHits){
	   if(kk.first == ij) fraction = kk.second;
	   break;
	 }
	 if(fraction != 1 && debug) std::cout << " >> MCMIP fraction " << fraction << std::endl;
	 allMatched_mcMIPE += hit->energy()*fraction;
	*/
	allMatched_mcMIPE += hit->energy();
      }
    }
    if(debug)    std::cout << " post matching conta E now all " << std::endl;
    for(auto ij : matched_mc){

      const HGCRecHit *hit = hitmap[ij];
      unsigned int preSize = allMatched_global.size();
      allMatched_global.insert(ij);
      if(allMatched_global.size() != preSize){
	/*
	  float fraction = 1.;
	  for(auto kk : allHits){
	  if(kk.first == ij) fraction = kk.second;
	  break;
	  }
	  if(fraction != 1 && debug) std::cout << " >> MC2d fraction " << fraction << std::endl;
	  allMatched_mcE += hit->energy()*fraction;
	*/
	allMatched_mcE += hit->energy();
      }
    }
    
    allMatched_globalE = allMatched_tkE + allMatched_mcMIPE + allMatched_mcE;


    if(debug)    std::cout << " post matching conta E now all V2 " << std::endl;
    for(auto ij : matched_all){

      const HGCRecHit *hit = hitmap[ij];
      /*
	float fraction = 1.;
	for(auto kk : allHits){
       	if(kk.first == ij) fraction = kk.second;
       	break;
	}
	if(fraction != 1 && debug) std::cout << " >> MC2d fraction " << fraction << std::endl;
	allMatched_E += hit->energy()*fraction;
      */
      allMatched_E += hit->energy();
    }

   

    if(debug){
      std::cout << cp.pdgId() << " nSimClusters = " << cp.simClusters().size() 
		<< " \n nHits matched to caloP for each collection = " 
		<< " sim " << allHits.size()
		<< " tk " << matched_tk.size()
		<< " mcMIP " << matched_mcMIP.size() << " mc2D " << matched_mc.size() 
		<< " global from singleC = " << allMatched_global.size() << " from finalP = " << matched_all.size()
		<< " \n corresponding energy = " 
		<< " cpE = " << cpEnergy << " tkE = " << allMatched_tkE 
		<< " mcMIPE = " << allMatched_mcMIPE << " mc2DE = " << allMatched_mcE 
		<< " globalE from singleC = " << allMatched_globalE << " from finalP = " << allMatched_E << std::endl;
    
      std::cout << " tk matched size = " << matched_tk.size() << " all = " << all_tk.size() << std::endl; 
      std::cout << " mc matched size = " << matched_mc.size() << " all = " << all_mc.size() << std::endl; 
      std::cout << " mc matchedMIP size = " << matched_mcMIP.size() << " all = " << all_mcMIP.size() << std::endl; 

      if(matched_mc.size() > all_mc.size()) std::cout << " problem ratio 1 " << std::endl;
    }

    if(allHits.size() != 0){
      histos_["hitsFraction_tk"]->Fill(1. * matched_tk.size()/allHits.size());
      histos_["hitsFraction_mcMIP"]->Fill(1. * matched_mcMIP.size()/allHits.size());
      histos_["hitsFraction_mc"]->Fill(1. * matched_mc.size()/allHits.size());
      histos_["hitsFraction_all"]->Fill(1. * allMatched_global.size()/allHits.size());
      histos_["hitsFraction_allV2"]->Fill(1. * matched_all.size()/allHits.size());
      if(finalParticles.size() == gpH->size()) histos_["hitsFraction_allV2_cpAsreco"]->Fill(1. * matched_all.size()/allHits.size());
      if(debug && 1. * matched_mcMIP.size()/allHits.size() > 0.9) 
	std::cout << " problem frac" <<  1. * matched_mcMIP.size()/allHits.size() << std::endl;
    }

    if(all_tk.size() != 0) histos_["hitsRecoOKFraction_tk"]->Fill(1. * matched_tk.size() / all_tk.size());
    if(all_mcMIP.size() != 0) histos_["hitsRecoOKFraction_mcMIP"]->Fill(1. * matched_mcMIP.size() / all_mcMIP.size());
    if(all_mc.size() != 0) histos_["hitsRecoOKFraction_mc"]->Fill(1. * matched_mc.size() / all_mc.size());
    if(matched_mc.size() != 0) histos_["hitsRecoOKFraction_allV2"]->Fill(1. * matched_all.size() / finalParticles[bestReco].size());
    
    histos_["EFraction_tk"]->Fill(1. * allMatched_tkE / cpEnergy);
    histos_["EFraction_mcMIP"]->Fill(1. * allMatched_mcMIPE / cpEnergy);
    histos_["EFraction_mc"]->Fill(1. * allMatched_mcE / cpEnergy);
    histos_["EFraction_all"]->Fill(1. * allMatched_globalE / cpEnergy);
    histos_["EFraction_allV2"]->Fill(1. * allMatched_E / cpEnergy);
    if(finalParticles.size() == gpH->size()) histos_["EFraction_allV2_cpAsreco"]->Fill(1. * allMatched_E / cpEnergy);
    if(debug && allMatched_globalE/cpEnergy > 1.) std::cout << " problem " << std::endl;

  }//caloParticles
}

std::vector<reco::HGCalMultiCluster> BuildParticles::buildTrackerItSameSeed(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs){

  if(debug)    std::cout << " >>> build trackerIteration merging objects with same seed initial size = " << mcs.size() << std::endl;

  std::vector<reco::HGCalMultiCluster> tempMC;
  reco::HGCalMultiCluster temp[2];
  //  edm::PtrVector<reco::BasicCluster> clusterPtrs[2];

  double baricenter[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
  double total_weight[2] = {0., 0.};

  //  int etaStartSide = int(mcs[0].eta() > 0); // 1 eta > 0 ; 0 eta < 0
  for(auto mc : mcs) {
    
    int etaSide = int(mc.eta() > 0);
    
    if(debug)    std::cout << " thisMC eta = " << mc.eta() << " phi = " << mc.phi() << " energy = " << mc.energy() << std::endl;
    //loop over all the layer clusters in the multicluster
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){
	
      reco::CaloClusterPtr sClPtr(*it);

      //      clusterPtrs[etaSide].push_back(sClPtr);
      temp[etaSide].push_back(sClPtr);
      auto weight = sClPtr->energy();
      total_weight[etaSide] += weight;
      baricenter[etaSide][0] += sClPtr->x() * weight;
      baricenter[etaSide][1] += sClPtr->y() * weight;
      baricenter[etaSide][2] += sClPtr->z() * weight;
    }
  }
  
  for(int ij=0; ij<2; ++ij){
    if(debug){
      std::cout << " ij = " << ij << " baricenter[ij][0] = " << baricenter[ij][0] 
		<< " baricenter[ij][1] = " << baricenter[ij][1]
		<< " baricenter[ij][2] = " << baricenter[ij][2] << " total_weight[ij] = " << total_weight[ij] << std::endl;
    }

    baricenter[ij][0] /= total_weight[ij];
    baricenter[ij][1] /= total_weight[ij];
    baricenter[ij][2] /= total_weight[ij];
    if(debug){
      std::cout << " ij = " << ij << " baricenter[ij][0] = " << baricenter[ij][0] 
		<< " baricenter[ij][1] = " << baricenter[ij][1]
		<< " baricenter[ij][2] = " << baricenter[ij][2] << " total_weight[ij] = " << total_weight[ij] << std::endl;
    }
    temp[ij].setEnergy(total_weight[ij]);
    temp[ij].setPosition(math::XYZPoint(baricenter[ij][0], baricenter[ij][1], baricenter[ij][2]));
    temp[ij].setAlgoId(reco::CaloCluster::hgcal_em);

    if(debug)    std::cout << " eta built = " << temp[ij].eta() << " phi = " << temp[ij].phi() << " energy = " << temp[ij].energy() << std::endl;

    tempMC.push_back(temp[ij]);
  }
  
  if(debug) std::cout << " >>> final size size = " << tempMC.size() << std::endl;

  return *(&(tempMC));
}


//// add some get truckster time
void BuildParticles::getTimeClusteredHitsList(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs, std::vector<float>& times){

  //  std::cout << " >>> found n multiclusters = " << mcs.size()  << std::endl;
  int counter = 0;
  for(auto mc : mcs) {

    std::vector<float> localTime;
    std::vector<std::pair<float,float>> localTimeW;
    //loop over all the layer clusters in the multicluster
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){

      reco::CaloClusterPtr sClPtr(*it);
      float time2D = (*time2DMap)[sClPtr];
      if(time2D >= 0){
	localTime.push_back(time2D);
	localTimeW.push_back(std::pair<float,float> (time2D, sClPtr->energy()));
      }
    }

    //use method from Shameena => try for the moment with truncation
    // can try truncation + weighted mean 
    float finalT = (localTime.size() >= 3) ? (hgcalsimclustertime::fixSizeHighestDensity(localTime)) : (-1.);
    float finalTW = (localTime.size() >= 3) ? (hgcalsimclustertime::fixSizeHighestDensity(localTimeW)) : (-1.);
    times[counter] = (finalT == -99) ? -1 : finalT;
    //    std::cout << " >>> void getTimeClusteredHitsList final time = " << finalT << " weighted = " << finalTW << std::endl;
    ++counter;
  }
  return;
}



void BuildParticles::getTimeTrackedHitsList(bool pos, const std::vector<reco::CaloCluster> &ccs, std::vector<float>& times){

  //vector => 1 position 1 track
  int counter = 0;
  for(auto cc : ccs) {
    std::vector<float> localTime;
    std::vector<std::pair<float,float>> localTimeW;
    //1 track = collection hits = caloCluster
    //compute time of caloCluster from hits

    const std::vector< std::pair<DetId, float> > &recHits =cc.hitsAndFractions();
    for(auto ij : recHits){

      const HGCRecHit *hit = hitmap[ij.first];
      float timeHit = hit->time();

      if(timeHit >= 0){
	localTime.push_back(timeHit);
	localTimeW.push_back(std::pair<float,float> (timeHit, hit->energy()));
      }
    }
    //se recHit singola da capire

    float finalT = (localTime.size() >= 3) ? (hgcalsimclustertime::fixSizeHighestDensity(localTime)) : (-1.);
    float finalTW = (localTime.size() >= 3) ? (hgcalsimclustertime::fixSizeHighestDensity(localTimeW)) : (-1.);
    times[counter] = (finalT == -99 || finalT == -1) ? -1 : (finalT -5.); // offset is already taken off for clusters
    ++counter;

  }
  return;
}




////////////////////
void BuildParticles::getClusteredHitsList(bool pos,
					  const std::vector<reco::HGCalMultiCluster> &mcs){
  

  if(debug)  std::cout << " in getClusteredHitsList" << std::endl;
  int mcCounter = 0;
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

    float prevDeltaT = 0.03;
    float clTime = pos ? mipParticles_time[mcCounter] : mcParticles_time[mcCounter];
    for(auto sw : finalParticles){
      float trkTime = finalParticles_time[counter];
      float deltaT = (clTime != -1 && trkTime != -1) ? std::abs(trkTime - clTime) : -1;
      if(deltaT != -1) deltaTtrack->Fill(deltaT);

      std::set<uint32_t> imatches=getMatched(testRh, sw);
      if(imatches.size() > 0) {
	assignToFinalParticle = counter;
	//	std::cout << " matchato " << std::endl;
	break;
      }

      
      if(deltaT < prevDeltaT){
	prevDeltaT = deltaT;
	assignToFinalParticle = counter;
	//break;
      } 
       
      ++counter;
    }

    if(debug) std::cout << " matching found assignToFinalParticle = " << assignToFinalParticle << " presize = " << finalParticles[assignToFinalParticle].size() << std::endl;

    if(assignToFinalParticle != -1){
      for(auto ij : testRh) finalParticles[assignToFinalParticle].insert(ij);
    }
    else{
      finalParticles.push_back(testRh);
      std::vector<float> dummy;
      dummy.push_back(mc.eta());
      dummy.push_back(mc.phi());
      finalParticles_etaPhi.push_back(dummy);

      finalParticles_time.push_back(clTime);
    }
    
    if(debug && assignToFinalParticle != -1) std::cout << " loop into getClusteredHitsList size = " << finalParticles[assignToFinalParticle].size() << std::endl;    
    if(debug) std::cout << " loop into getClusteredHitsList size = " << finalParticles[finalParticles.size()-1].size() << std::endl;    

    ++mcCounter;
  }
  return;
}


void BuildParticles::getTrackedItHitsList(bool pos, const std::vector<reco::HGCalMultiCluster> &mcs){

  if(debug)  std::cout << " in getTrackedItHitsList mcs size = " << mcs.size() << std::endl;
  int mcCounter = 0;
  for(auto mc : mcs) {
    if(debug)    std::cout << " >>> add multi cluster eta = " << mc.eta() << " phi = " << mc.phi() << std::endl;

    //loop over all the layer clusters in the multicluster                                                        
    for(reco::HGCalMultiCluster::component_iterator it = mc.begin(); it!=mc.end(); it++){
      const std::vector< std::pair<DetId, float> > &recHits = (*it)->hitsAndFractions();
      for(auto ij : recHits){
	if(ij.second != 0) {
	  finalParticles[mcCounter].insert(ij.first.rawId());
	}
      }
    }

    finalParticles_etaPhi[mcCounter].push_back(mc.eta());
    finalParticles_etaPhi[mcCounter].push_back(mc.phi());
    finalParticles_time[mcCounter] = tkParticles_time[mcCounter];

    if(debug){
      std::cout << " >>> new trk cluster eta = " << mc.eta() << " phi = " << mc.phi() << " and "
		<< " finalParticles_etaPhi[mcCounter][0] = " << finalParticles_etaPhi[mcCounter][0] << " " << finalParticles_etaPhi[mcCounter][1] << std::endl;
    }

    ++mcCounter;

    if(debug) std::cout << " loop into getTrackedItHitsList size = " << finalParticles[mcCounter-1].size() << std::endl;
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

    finalParticles_time[ipos] = tkParticles_time[ipos];

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


//in 1st not in 2nd
std::set<uint32_t> BuildParticles::getNOTMatched(const std::set<uint32_t> &a,
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



DEFINE_FWK_MODULE(BuildParticles);
