#ifndef _ticlanalyzer_h_
#define _ticlanalyzer_h_

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

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

#include "TSystem.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TTree.h"


class TICLAnalyzer : public edm::EDAnalyzer {
 public:

  explicit TICLAnalyzer(const edm::ParameterSet&);

  ~TICLAnalyzer();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

 private:
  
  std::vector<uint32_t> getClusteredHitsList(bool pos,
                                             const std::vector<std::pair<uint32_t,float> >  &hits,
                                             const std::vector<reco::HGCalMultiCluster> &mcs);
  std::vector<uint32_t> getTrackedHitsList(bool pos,
                                           const std::vector<std::pair<uint32_t,float> >  &hits,
                                           const std::vector<reco::CaloCluster> &ccs);
  std::vector<uint32_t> getMatched(const std::vector<std::pair<uint32_t,float> > &a,
                                   const std::vector<std::pair<DetId,float> > &b);
  

  edm::EDGetTokenT<std::vector<CaloParticle> > genToken_;
  edm::EDGetTokenT<std::vector<SimCluster> > simclusToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster> > eetkToken_,fhtkToken_;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > mcMIPToken_,mcToken_;

  std::map<TString,TH1 *> histos_;

};

#endif
