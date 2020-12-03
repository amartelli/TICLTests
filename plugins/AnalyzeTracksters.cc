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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "DataFormats/Math/interface/GeantUnits.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "TSystem.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TTree.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/Common/interface/ValueMap.h"

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>> XYZPointF;


constexpr double m_pi = 0.13957018;
constexpr double m_pi_inv2 = 1.0 / m_pi / m_pi;
constexpr double m_k = 0.493677;
constexpr double m_k_inv2 = 1.0 / m_k / m_k;
constexpr double m_p = 0.9382720813;
constexpr double m_p_inv2 = 1.0 / m_p / m_p;
constexpr double c_cm_ns = geant_units::operators::convertMmToCm(CLHEP::c_light);  // [mm/ns] -> [cm/ns]
constexpr double c_inv = 1.0 / c_cm_ns;


using namespace std;
using namespace edm;
using namespace reco;
using namespace ticl;

class AnalyzeTracksters : public edm::EDAnalyzer {
public:

  explicit AnalyzeTracksters(const edm::ParameterSet&);

  ~AnalyzeTracksters();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

private:
  hgcal::RecHitTools         recHitTools;

  float GetStraightTrackLength(const GlobalPoint &startingPoint,
				 const GlobalPoint &endPoint);

  float GetPropagatedTrackLength(const GlobalPoint &startingPoint,
				 const GlobalPoint &endPoint,
				 const reco::Track &track);

  void layerIntersection(std::array<double,3> &to,
			 const std::array<double,3> &from,
			 const std::array<double,3> &vtx);

  void layerIntersection(std::array<double,3> &to,
			 const GlobalPoint &from,
			 const GlobalPoint &vtx);

  std::vector<size_t> decrease_sorted_indices(const std::vector<float>& v);

  std::pair<float, float> fixSizeHighestDensity(std::vector<float>& time, std::vector<float> weight,
                                                unsigned int minNhits = 3, float deltaT=0.210, float timeWidthBy=0.5);


  const edm::EDGetTokenT<std::vector<CaloParticle> > genToken_;
  const edm::EDGetTokenT<XYZPointF> genVtxPositionToken_;
  const edm::EDGetTokenT<float> genVtxTimeToken_;
  const edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vtxsToken_;

  const edm::EDGetTokenT<std::vector<Trackster> > tracksters_mergeToken_;

  const edm::EDGetTokenT<HGCRecHitCollection> hits_eeToken_;
  const edm::EDGetTokenT<HGCRecHitCollection> hits_fhToken_;
  const edm::EDGetTokenT<HGCRecHitCollection> hits_bhToken_;

  const edm::EDGetTokenT<std::vector<reco::CaloCluster> > layerClustersToken_;
  const edm::EDGetTokenT<edm::ValueMap<pair<float,float> > > timeMapToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
  edm::ESHandle<MagneticField> bfield_;

  std::map<uint32_t, const HGCRecHit*> hitmap;
  

  TH1F* h_diff_hgc;
  TH1F* h_error_hgc;
  TProfile* h_diff_vsE_hgc;
  TProfile* h_error_vsE_hgc;
  TH1F* h_pull_hgc;
  TH2F* h_diff_vs_dR_hgc;
  TH2F* h_pull_vs_dR_hgc;



  TH1F* h_diff_hgc_Ecut;
  TH1F* h_pull_hgc_Ecut;
  TH1F* h_diff_hgc_dRcut;
  TH1F* h_pull_hgc_dRcut;

  TH1F* trackster_Energy;

  //stLC
  TH1F* h_diff_stLC;
  TH1F* h_pull_stLC;
  TH2F* h_diff_vs_dR_stLC;
  TH2F* h_pull_vs_dR_stLC;

  //seedLC
  TH1F* h_diff_seedLC;
  TH1F* h_pull_seedLC;
  TH2F* h_diff_vs_dR_seedLC;
  TH2F* h_pull_vs_dR_seedLC;


  //time seed 
  TH1F* h_diff_seedLCT;
  TH1F* h_pull_seedLCT;
  TH2F* h_pull_vs_dR_seedLCT;

  TH2F* h_diff_vs_dR_seedLCT;
  TProfile* tp_diff_vs_dR_seedLCT;
  TH2F* h_diff_vs_E_seedLCT;
  TProfile* tp_diff_vs_E_seedLCT;
  TH2F* h_diff_vs_EdRcut_seedLCT;
  TProfile* tp_diff_vs_EdRcut_seedLCT;

  TH2F* h_diff_seed_lc_vsR;
  TProfile* tp_diff_seed_lc_vsR;
  TH2F* h_diff_seed_lc_vsE;
  TProfile* tp_diff_seed_lc_vsE;
  TH2F* h_diff_seed_lc_vsERcut;
  TProfile* tp_diff_seed_lc_vsERcut;

  TH2F* h2_energy_vs_dR;
  TH2F* h2_time_vs_energy;
  TProfile* tp_time_vs_energy;
  TH1F* hTime_2Dcl;

  bool debug;

  bool isPhoton;
  bool isTrkIteration;
  int nEvent;
};



AnalyzeTracksters::AnalyzeTracksters(const edm::ParameterSet& iConfig) :
  genToken_(consumes<std::vector<CaloParticle> >(edm::InputTag("mix:MergedCaloTruth"))),
  genVtxPositionToken_(consumes<XYZPointF>(edm::InputTag("genParticles:xyz0"))),
  genVtxTimeToken_(consumes<float>(edm::InputTag("genParticles:t0"))),
  tracksToken_(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks"))),
  vtxsToken_(consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices4DWithBS"))),
  tracksters_mergeToken_(consumes<std::vector<Trackster>>(edm::InputTag("ticlTrackstersMerge"))),
  //  tracksters_mergeToken_(consumes<std::vector<Trackster>>(edm::InputTag("ticlTrackstersEM"))),
  hits_eeToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"))),
  hits_fhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"))),
  hits_bhToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"))),
  layerClustersToken_(consumes<std::vector<reco::CaloCluster> >(edm::InputTag("hgcalLayerClusters::RECO"))),
  timeMapToken_(consumes<edm::ValueMap<pair<float,float> > > (edm::InputTag("hgcalLayerClusters:timeLayerCluster")))
{

  //  auto sumes = consumesCollector();
  //  bfield_token_ = sumes.esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>();
  bfield_token_ = esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>();

  //book some histos here
  edm::Service<TFileService> fs;

  trackster_Energy = fs->make<TH1F>("trackster_Energy", "", 2000, 0., 500.);

  h_diff_hgc = fs->make<TH1F>("h_diff_hgc", "", 200, -0.1, 0.1);
  h_error_hgc = fs->make<TH1F>("h_error_hgc", "", 100, 0., 0.1);
  h_diff_vsE_hgc = fs->make<TProfile>("h_diff_vsE_hgc", "", 500, 0., 500.);
  h_error_vsE_hgc = fs->make<TProfile>("h_error_vsE_hgc", "", 500, 0., 500.);
  h_pull_hgc = fs->make<TH1F>("h_pull_hgc", "", 5000, -5., 5.);
  h_diff_vs_dR_hgc = fs->make<TH2F>("h_diff_vs_dR_hgc", "", 100, 0., 10., 200, -0.1, 0.1);
  h_pull_vs_dR_hgc = fs->make<TH2F>("h_pull_vs_dR_hgc", "", 100, 0., 10., 5000, -5., 5.);
  //
  h_diff_hgc_Ecut = fs->make<TH1F>("h_diff_hgc_Ecut", "", 200, -0.1, 0.1);
  h_pull_hgc_Ecut = fs->make<TH1F>("h_pull_hgc_Ecut", "", 5000, -5., 5.);
  //
  h_diff_hgc_dRcut = fs->make<TH1F>("h_diff_hgc_dRcut", "", 200, -0.1, 0.1);
  h_pull_hgc_dRcut = fs->make<TH1F>("h_pull_hgc_dRcut", "", 5000, -5., 5.);


  h_diff_stLC = fs->make<TH1F>("h_diff_stLC", "", 200, -0.1, 0.1);
  h_pull_stLC = fs->make<TH1F>("h_pull_stLC", "", 5000, -5., 5.);
  h_diff_vs_dR_stLC = fs->make<TH2F>("h_diff_vs_dR_stLC", "", 100, 0., 10., 200, -0.1, 0.1);
  h_pull_vs_dR_stLC = fs->make<TH2F>("h_pull_vs_dR_stLC", "", 100, 0., 10., 5000, -5., 5.);

  h_diff_seedLC = fs->make<TH1F>("h_diff_seedLC", "", 200, -0.1, 0.1);
  h_pull_seedLC = fs->make<TH1F>("h_pull_seedLC", "", 5000, -5., 5.);
  h_diff_vs_dR_seedLC = fs->make<TH2F>("h_diff_vs_dR_seedLC", "", 100, 0., 10., 200, -0.1, 0.1);
  h_pull_vs_dR_seedLC = fs->make<TH2F>("h_pull_vs_dR_seedLC", "", 100, 0., 10., 5000, -5., 5.);

  //time seed
  h_diff_seedLCT = fs->make<TH1F>("h_diff_seedLCT", "", 200, -0.1, 0.1);
  h_pull_seedLCT = fs->make<TH1F>("h_pull_seedLCT", "", 5000, -5., 5.);
  h_pull_vs_dR_seedLCT = fs->make<TH2F>("h_pull_vs_dR_seedLCT", "", 100, 0., 10., 5000, -5., 5.);

  h_diff_vs_dR_seedLCT = fs->make<TH2F>("h_diff_vs_dR_seedLCT", "", 100, 0., 10., 200, -0.1, 0.1);
  tp_diff_vs_dR_seedLCT = fs->make<TProfile>("tp_diff_vs_dR_seedLCT", "", 100, 0., 10.);
  h_diff_vs_E_seedLCT = fs->make<TH2F>("h_diff_vs_E_seedLCT", "", 1000, 0., 100., 200, -0.1, 0.1);
  tp_diff_vs_E_seedLCT = fs->make<TProfile>("tp_diff_vs_E_seedLCT", "", 1000, 0., 100.);
  h_diff_vs_EdRcut_seedLCT = fs->make<TH2F>("h_diff_vs_EdRcut_seedLCT", "", 1000, 0., 100., 200, -0.1, 0.1);
  tp_diff_vs_EdRcut_seedLCT = fs->make<TProfile>("tp_diff_vs_EdRcut_seedLCT", "", 1000, 0., 100.);

  h_diff_seed_lc_vsR = fs->make<TH2F>("h_diff_seed_lc_vsR", "", 100, 0., 10., 200, -0.1, 0.1);
  tp_diff_seed_lc_vsR = fs->make<TProfile>("tp_diff_seed_lc_vsR", "", 1000, 0., 10.);
  h_diff_seed_lc_vsE = fs->make<TH2F>("h_diff_seed_lc_vsE", "", 1000, 0., 5., 200, -0.1, 0.1);
  //h_diff_seed_lc_vsE = fs->make<TH2F>("h_diff_seed_lc_vsE", "", 1000, 0., 5., 200, -100., 100.);
  tp_diff_seed_lc_vsE = fs->make<TProfile>("tp_diff_seed_lc_vsE", "", 1000, 0., 5.);
  h_diff_seed_lc_vsERcut = fs->make<TH2F>("h_diff_seed_lc_vsERcut", "", 1000, 0., 5., 200, -0.1, 0.1);
  tp_diff_seed_lc_vsERcut = fs->make<TProfile>("tp_diff_seed_lc_vsERcut", "", 1000, 0., 5.);

  h2_energy_vs_dR = fs->make<TH2F>("h2_energy_vs_dR", "", 100, 0., 10., 2000, 0., 20.);
  h2_time_vs_energy = fs->make<TH2F>("h2_time_vs_energy", "", 2000, 0., 20., 200, -0.1, 0.1);
  tp_time_vs_energy = fs->make<TProfile>("tp_time_vs_energy", "", 2000, 0., 20.);

  hTime_2Dcl = fs->make<TH1F>("hTime_2Dcl", "", 60000, -5., 25.);

  debug = false;
  isTrkIteration = false;
  isPhoton = true;
  nEvent = 0;
}



AnalyzeTracksters::~AnalyzeTracksters() { 
}



void AnalyzeTracksters::beginRun(const edm::Run& run, 
			       const edm::EventSetup & es) { 

  bfield_ = es.getHandle(bfield_token_);
}


void  AnalyzeTracksters::analyze(const Event& iEvent, 
			       const EventSetup& iSetup) {

  ++nEvent;
  edm::ESHandle<CaloGeometry> geom;
  iSetup.get<CaloGeometryRecord>().get(geom);
  recHitTools.setGeometry(*geom);
  
  if(debug)  std::cout << " \n \n new evt "<< std::endl;
  
  edm::Handle<std::vector<CaloParticle> > gpH;
  iEvent.getByToken( genToken_, gpH);

  edm::Handle<std::vector<reco::Track>> track_h;
  iEvent.getByToken(tracksToken_, track_h);
  const auto &tracks = *track_h;

  //reco vtx                                                                                                                                                        
  edm::Handle<std::vector<reco::Vertex>> vtx_h;
  iEvent.getByToken(vtxsToken_, vtx_h);
  const auto &vtxs = *vtx_h;

  //gen vtx pos
  const auto& genVtxPosition_h = iEvent.getHandle(genVtxPositionToken_);
  const auto &genVtxPos = *genVtxPosition_h;

  //gen vtx t                                                                                                                                     
  const auto& genVtxTime_h = iEvent.getHandle(genVtxTimeToken_);
  const auto &genVtxTime = *genVtxTime_h;


  edm::Handle<std::vector<Trackster>> trackstersMergeH;
  iEvent.getByToken(tracksters_mergeToken_, trackstersMergeH);
  const auto &trackstersM = *trackstersMergeH;

  edm::Handle<std::vector<reco::CaloCluster>> clusterHandle;
  iEvent.getByToken(layerClustersToken_, clusterHandle);
  const auto &layerClusters = *clusterHandle;

  edm::Handle<edm::ValueMap<std::pair<float, float>>> clustersTime_h;
  iEvent.getByToken(timeMapToken_, clustersTime_h);
  const auto& layerClustersTimes = *clustersTime_h;
  if(debug) std::cout << " time2DMap->size() = " << layerClustersTimes.size() << std::endl;

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


  //loop over tracksters (multiclusters for the moment)
  for(auto iT : trackstersM){

    size_t N = iT.vertices().size();
    if(N == 0) continue;
    
    if(debug) std::cout << " time = " << iT.time() << " energy = " << iT.raw_energy() << " size = " << N << std::endl;
    if(iT.timeError() == -1) continue;
    
    for(auto cp : *gpH) {
      if(debug) std::cout << " >>> cp.eventId().event() = " << cp.eventId().event() 
			  << " cp.eventId().bunchCrossing() = " << cp.eventId().bunchCrossing() 
			  << " cp.pdgId() = " << cp.pdgId() << " cp.simClusters().size() = " << cp.simClusters().size() << std::endl;

      bool EventOk = false;
      if(cp.pdgId() == 22 && cp.simClusters().size() <= 2) EventOk = true;
      if(abs(cp.pdgId()) == 211 && cp.simClusters().size() <= 1) EventOk = true;
      if(abs(cp.pdgId()) == 211 && cp.simClusters().size() <= 1) EventOk = true;
      if(!EventOk) continue;
    
      float etaGen = cp.eta();
      float phiGen = cp.phi();
      float xGen = cp.momentum().x();
      float yGen = cp.momentum().y();
      float zGen = cp.momentum().z();

      if(iT.barycenter().eta() * etaGen < 0.) continue;

      size_t seedLC = -1;
      float maxE = 0.;
      float minTLC = 99.;
      int minL = 100;
      size_t firstLC = -1;
      size_t earliestLC = -1;

      for (size_t i = 0; i < N; ++i) {
        const reco::CaloCluster &cluster = layerClusters[iT.vertices(i)];
        int j = recHitTools.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
	float lcT = layerClustersTimes.get(iT.vertices(i)).first;
	float lcTE = layerClustersTimes.get(iT.vertices(i)).second;
	if(lcTE == -1.) continue;
        if(j < minL){
          minL = j;
          firstLC = i;
        }
        float energyLC = cluster.energy();
        if(energyLC > maxE){
          maxE = energyLC;
          seedLC = i;
	}
	if(lcT < minTLC){
	  minTLC = lcT;
	  earliestLC = i;
	}
      }
      const reco::CaloCluster &firstcl = layerClusters[iT.vertices(firstLC)];
      const reco::CaloCluster &seedcl = layerClusters[iT.vertices(seedLC)];
      const reco::CaloCluster &earlycl = layerClusters[iT.vertices(earliestLC)];

      GlobalPoint origin(0., 0., 0.);
      GlobalPoint genVtx(genVtxPos.x(), genVtxPos.y(), genVtxPos.z());
      GlobalPoint hgcPoint(iT.barycenter().x(), iT.barycenter().y(), iT.barycenter().z());
      GlobalPoint firstLCPoint(firstcl.x(), firstcl.y(), firstcl.z());
      GlobalPoint seedLCPoint(seedcl.x(), seedcl.y(), seedcl.z());
      GlobalPoint earlyLCPoint(earlycl.x(), earlycl.y(), earlycl.z());

      float trackL0_hgc = GetStraightTrackLength(origin, hgcPoint);
      float trackL_hgc = GetStraightTrackLength(genVtx, hgcPoint);
      float dT_hgc = (trackL0_hgc - trackL_hgc) * c_inv;

      float trackL0_stLC = GetStraightTrackLength(origin, firstLCPoint);
      float trackL_stLC = GetStraightTrackLength(genVtx, firstLCPoint);
      float dT_stLC = (trackL0_stLC - trackL_stLC) * c_inv;

      float trackL0_seedLC = GetStraightTrackLength(origin, seedLCPoint);
      float trackL_seedLC = GetStraightTrackLength(genVtx, seedLCPoint);
      float dT_seedLC = (trackL0_seedLC - trackL_seedLC) * c_inv;


      //if charged
      if(iT.seedID() != edm::ProductID()) {
	auto trackIdx = iT.seedIndex();
	auto const &track = tracks[trackIdx];
	//trackP = track.p();
	trackL_hgc = GetPropagatedTrackLength(genVtx, hgcPoint, track);
	trackL_stLC = GetPropagatedTrackLength(genVtx, firstLCPoint, track);
	trackL_seedLC = GetPropagatedTrackLength(genVtx, seedLCPoint, track);

	float gammasq_pi = 1. + track.p2() * m_pi_inv2;
	float beta_pi = std::sqrt(1. - 1. / gammasq_pi);
	float dt_pi = (trackL0_hgc - trackL_hgc / beta_pi) * c_inv;
	dT_hgc = dt_pi;

	dT_stLC = (trackL0_stLC - trackL_stLC / beta_pi) * c_inv;
	dT_seedLC = (trackL0_seedLC - trackL_seedLC / beta_pi) * c_inv;


      }

      float correctedT_hgc = iT.time() + dT_hgc;
      float tDiff_recoGen_hgc = correctedT_hgc - genVtxTime;
      float tPull_recoGen_hgc = tDiff_recoGen_hgc / iT.timeError();

      //st_LC
      float correctedT_stLC = iT.time() + dT_stLC;
      float tDiff_recoGen_stLC = correctedT_stLC - genVtxTime;
      float tPull_recoGen_stLC = tDiff_recoGen_stLC / iT.timeError();
      //seed_LC
      float correctedT_seedLC = iT.time() + dT_seedLC;
      float tDiff_recoGen_seedLC = correctedT_seedLC - genVtxTime;
      float tPull_recoGen_seedLC = tDiff_recoGen_seedLC / iT.timeError();

      std::array<double,3> to{ {0., 0., hgcPoint.z()} };
      layerIntersection(to, GlobalPoint(xGen, yGen, zGen), genVtx);   
      float dR_hgc = sqrt(pow(hgcPoint.x() - to[0], 2) + pow(hgcPoint.y() - to[1], 2));

      std::array<double,3> to_stLC{ {0., 0., firstLCPoint.z()} };
      layerIntersection(to_stLC, GlobalPoint(xGen, yGen, zGen), genVtx);   
      float dR_stLC = sqrt(pow(firstLCPoint.x() - to_stLC[0], 2) + pow(firstLCPoint.y() - to_stLC[1], 2));
      //
      std::array<double,3> to_seedLC{ {0., 0., seedLCPoint.z()} };
      layerIntersection(to_seedLC, GlobalPoint(xGen, yGen, zGen), genVtx);   
      float dR_seedLC = sqrt(pow(seedLCPoint.x() - to_seedLC[0], 2) + pow(seedLCPoint.y() - to_seedLC[1], 2));

      //hgc
      h_error_hgc->Fill(iT.timeError());
      h_diff_hgc->Fill(tDiff_recoGen_hgc);
      h_error_vsE_hgc->Fill(iT.raw_pt(), iT.timeError());
      h_diff_vsE_hgc->Fill(iT.raw_pt(), tDiff_recoGen_hgc);
      h_pull_hgc->Fill(tPull_recoGen_hgc);
      h_diff_vs_dR_hgc->Fill(dR_hgc, tDiff_recoGen_hgc);
      h_pull_vs_dR_hgc->Fill(dR_hgc, tPull_recoGen_hgc);
      //
      if(dR_stLC < 2){
	h_diff_hgc_dRcut->Fill(tDiff_recoGen_hgc);
	h_pull_hgc_dRcut->Fill(tPull_recoGen_hgc);
      }
      if(iT.raw_pt() > 10.){
	h_diff_hgc_Ecut->Fill(tDiff_recoGen_hgc);
	h_pull_hgc_Ecut->Fill(tPull_recoGen_hgc);
      }


      //st_LC
      h_diff_stLC->Fill(tDiff_recoGen_stLC);
      h_pull_stLC->Fill(tPull_recoGen_stLC);
      h_diff_vs_dR_stLC->Fill(dR_stLC, tDiff_recoGen_stLC);
      h_pull_vs_dR_stLC->Fill(dR_stLC, tPull_recoGen_stLC);

      //seed_LC
      h_diff_seedLC->Fill(tDiff_recoGen_seedLC);
      h_pull_seedLC->Fill(tPull_recoGen_seedLC);
      h_diff_vs_dR_seedLC->Fill(dR_seedLC, tDiff_recoGen_seedLC);
      h_pull_vs_dR_seedLC->Fill(dR_seedLC, tPull_recoGen_seedLC);

      float tSeed = layerClustersTimes.get(iT.vertices(seedLC)).first;
      float tESeed = layerClustersTimes.get(iT.vertices(seedLC)).second;
      float tDiff_recoGen_stLCT = tSeed + dT_seedLC - genVtxTime;
      h_diff_seedLCT->Fill(tDiff_recoGen_stLCT);
      h_pull_seedLCT->Fill(tDiff_recoGen_stLCT/tESeed);
      h_pull_vs_dR_seedLCT->Fill(dR_seedLC, tDiff_recoGen_stLCT/tESeed);
      h_diff_vs_dR_seedLCT->Fill(dR_seedLC, tSeed + dT_seedLC - genVtxTime);
      tp_diff_vs_dR_seedLCT->Fill(dR_seedLC, tSeed + dT_seedLC - genVtxTime);
      h_diff_vs_E_seedLCT->Fill(seedcl.energy(), tSeed + dT_seedLC - genVtxTime);
      tp_diff_vs_E_seedLCT->Fill(seedcl.energy(), tSeed + dT_seedLC - genVtxTime);
      if(dR_seedLC < 2.){
	h_diff_vs_EdRcut_seedLCT->Fill(seedcl.energy(), (tSeed + dT_seedLC - genVtxTime));
	tp_diff_vs_EdRcut_seedLCT->Fill(seedcl.energy(), (tSeed + dT_seedLC - genVtxTime));
      }
      //time diff wrt seed layer cluster
      hTime_2Dcl->Fill(tSeed);
      for (size_t i = 0; i < N; ++i) {
	if(seedLC == i) continue;
        const reco::CaloCluster &cluster = layerClusters[iT.vertices(i)];
        float lcT = layerClustersTimes.get(iT.vertices(i)).first;
	float lcE = layerClustersTimes.get(iT.vertices(i)).second;
	if(lcE == -1 || lcT == -99.) continue;

	std::array<double,3> to_lc{ {0., 0., seedLCPoint.z()} };
	layerIntersection(to_lc, GlobalPoint(xGen, yGen, zGen), genVtx);
	float dR_lc = sqrt(pow(seedLCPoint.x() - to_lc[0], 2) + pow(seedLCPoint.y() - to_lc[1], 2));


	float tdiff = lcT - tSeed;
	if(tdiff > 10.) std::cout << " lcT = " << lcT << " tSeed = " << tSeed << std::endl;
	h_diff_seed_lc_vsR->Fill(dR_lc, tdiff);
	tp_diff_seed_lc_vsR->Fill(dR_lc, tdiff);
	h_diff_seed_lc_vsE->Fill(cluster.energy(), tdiff);
	tp_diff_seed_lc_vsE->Fill(cluster.energy(), tdiff);

	float combE = sqrt(tESeed*tESeed + lcE*lcE);
	// h_pull_seed_lc_vsR->Fill(dR_lc, tdiff/combE);
	// tp_pull_seed_lc_vsR->Fill(dR_lc, tdiff/combE);
	if(dR_lc < 2.){
	  h_diff_seed_lc_vsERcut->Fill(cluster.energy(), tdiff);
	  tp_diff_seed_lc_vsERcut->Fill(cluster.energy(), tdiff);
	  h2_time_vs_energy->Fill(cluster.energy(), tdiff);
	  tp_time_vs_energy->Fill(cluster.energy(), tdiff);

	}
	h2_energy_vs_dR->Fill(dR_lc, cluster.energy());
	hTime_2Dcl->Fill(layerClustersTimes.get(iT.vertices(i)).first);
      }
      

      trackster_Energy->Fill(iT.raw_pt());
    }//caloParticles

  }//loop over tracksters
}

float AnalyzeTracksters::GetStraightTrackLength(const GlobalPoint &startingPoint,
					      const GlobalPoint &endPoint){

  return sqrt(pow(startingPoint.x() - endPoint.x(), 2) + pow(startingPoint.y() - endPoint.y(), 2) + pow(startingPoint.z() - endPoint.z(), 2));
}


float AnalyzeTracksters::GetPropagatedTrackLength(const GlobalPoint &startingPoint,
						const GlobalPoint &endPoint,
						const reco::Track &track){

  GlobalVector startingMomentum(track.px(), track.py(), track.pz());
  auto magField = bfield_.product();
  FreeTrajectoryState trajectory(startingPoint, startingMomentum, track.charge(), magField);
  SteppingHelixPropagator propagator(magField);
  float propagatedTrackL_ = propagator.propagateWithPath(trajectory, endPoint).second;
  return propagatedTrackL_;
}

void AnalyzeTracksters::layerIntersection(std::array<double,3> &to,
					const GlobalPoint &from,
					const GlobalPoint &vtx){

  std::array<double,3> fromV{ {from.x(), from.y(), from.z()} };
  std::array<double,3> vtxV{ {vtx.x(), vtx.y(), vtx.z()} };
  return layerIntersection(to, fromV, vtxV);
}


void AnalyzeTracksters::layerIntersection(std::array<double,3> &to,
					const std::array<double,3> &from,
					const std::array<double,3> &vtx){
  to[0]=(from[0]-vtx[0]) / (from[2] - vtx[2]) * (to[2] - from[2]) + from[0];
  to[1]=(from[1]-vtx[1]) / (from[2] - vtx[2]) * (to[2] - from[2]) + from[1];
  return;
}


std::vector<size_t> AnalyzeTracksters::decrease_sorted_indices(const std::vector<float>& v) {
  // initialize original index locations                                                       
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indices based on comparing values in v (decreasing order)                            
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  return idx;
};

std::pair<float, float> AnalyzeTracksters::fixSizeHighestDensity(std::vector<float>& time, std::vector<float> weight,
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



DEFINE_FWK_MODULE(AnalyzeTracksters);

