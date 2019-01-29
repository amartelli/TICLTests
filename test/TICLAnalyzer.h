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
  
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > mcToken_;

  std::map<TString,TH1 *> histos_;

};

#endif
