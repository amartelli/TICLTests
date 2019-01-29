#include "RecoHGCal/TICLTests/test/TICLAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

using namespace std;
using namespace edm;
using namespace reco;

TICLAnalyzer::TICLAnalyzer(const edm::ParameterSet& iConfig) :
  mcToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("hgcalMultiClusters")))
{
  histos_["nmc"] = new TH1F("nmc", ";Multicluster multiplicity;",                 10,   0 , 10 );
  histos_["nlc"] = new TH1F("nlc", ";Layer clusters/multicluster multiplicity;",  100,  0 , 100 );
  histos_["en"]  = new TH1F("en",  ";Energy [GeV];",  100,  0 , 100 );  
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
  
  edm::Handle<std::vector<reco::HGCalMultiCluster> > mcH;
  iEvent.getByToken( mcToken_,mcH);
  histos_["nmc"]->Fill(mcH->size());
  for(auto m : *mcH) {
    histos_["nlc"]->Fill(m.size());
    histos_["en"]->Fill(m.energy());
  }

}

DEFINE_FWK_MODULE(TICLAnalyzer);
