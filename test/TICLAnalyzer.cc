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
  genToken_(consumes<std::vector<reco::GenParticle> >(edm::InputTag("genParticles"))),
  mcToken_(consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("hgcalMultiClusters")))
{
  histos_["dr"]   = new TH1F("dr",   ";#Delta R;", 50,0,0.5);
  histos_["nlc"]  = new TH1F("nlc",  ";Layer clusters/multicluster multiplicity;",  100,  0 , 100 );
  histos_["resp"] = new TH1F("resp", ";E(multicluster)/E(gen);",  50,  0 , 2 );  
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

  edm::Handle<std::vector<reco::GenParticle> > gpH;
  iEvent.getByToken( genToken_, gpH);

  edm::Handle<std::vector<reco::HGCalMultiCluster> > mcH;
  iEvent.getByToken( mcToken_,mcH);

  for(size_t i=0; i<gpH->size(); i++) {
    const reco::GenParticle &g=gpH->at(i);

    //closest multicluster
    int bestMatchIdx(-1);
    float minDR2(9999.);
    for(size_t j=0; j<mcH->size(); j++) {
      const reco::HGCalMultiCluster &m=mcH->at(j);

      float dR2=deltaR2(g,m);
      if(dR2>minDR2) continue;
      minDR2=dR2;
      bestMatchIdx=j;
    }
    if(bestMatchIdx<0) continue;
    
    //control dists
    const reco::HGCalMultiCluster &m=mcH->at(bestMatchIdx);
    histos_["nlc"]->Fill(m.size());
    histos_["dr"]->Fill(sqrt(minDR2));
    histos_["resp"]->Fill(m.energy()/g.energy());
  }
}

DEFINE_FWK_MODULE(TICLAnalyzer);
