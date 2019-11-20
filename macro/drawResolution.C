#include "TLegend.h"
#include "TLatex.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include "TPaveStats.h"


/*
float getXmax(TH1F* histo, float& YMax){

  float yVal = 0.;
  int xBin = 1;
  for(int iB=1; iB<histo->GetNbinsX(); ++iB){
    if(histo->GetBinContent(iB) > yVal){
      xBin = iB;
      yVal = histo->GetBinContent(iB);
      YMax = yVal;
      if(yVal > 0 && histo->GetBinContent(iB) < yVal) break;
    }
  }

  std::cout << histo->GetName() << " " << histo->GetBinCenter(xBin) << std::endl; 
  return histo->GetBinCenter(xBin);
}

*/

void drawResolution(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  //  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  //  gROOT->LoadMacro("~/public/myStyle.C");
  //  gStyle->SetOptStat(1);
  //gStyle->SetOptFit(1);

  std::cout << " inizio ci sono " << std::endl; 


  bool doAllTheFits = true;
  
  //  int iColors[16] = {kRed, kOrange+4, kOrange-3, kOrange-2, kBlue, kBlue-9, kAzure-9, kAzure+10, kCyan, kGreen+1, kCyan-2, kYellow+2}; //kGray+1};
  int iColors[4] = {kBlue, kRed+1, kGreen+2, kMagenta};
  int iStyle[3] = {20, 22}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  

  std::vector<std::string> ptName;
  ptName.push_back("10");
  ptName.push_back("60");
  ptName.push_back("150");

  std::vector<std::string> etaName;
  etaName.push_back("LowEta");
  etaName.push_back("HighEta");

  std::vector<float> ptVal;
  ptVal.push_back(10.);
  ptVal.push_back(60.);
  ptVal.push_back(150.);

  std::vector<float> etaVal;
  etaVal.push_back(1.7);
  etaVal.push_back(2.7);


  TGraphErrors* scaleVsPt_recoBP[2];
  TGraphErrors* resolutionVsPt_recoBP[2];
  TGraphErrors* scaleVsPt_recoAll[2];
  TGraphErrors* resolutionVsPt_recoAll[2];
  TGraphErrors* scaleVsPt_recoBPmatched[2];
  TGraphErrors* resolutionVsPt_recoBPmatched[2];
  TGraphErrors* scaleVsPt_recoAllmatched[2];
  TGraphErrors* resolutionVsPt_recoAllmatched[2];

  TFile* inF[3][2];
  for(int ij=0; ij<3; ++ij){
    for(int kl=0; kl<2; ++kl){
      //      inF[ij][kl] = TFile::Open(("../test/TICL0PU/ticl_22Pt"+ptName.at(ij)+"_"+etaName.at(kl)+"_timeV2.root").c_str());
      inF[ij][kl] = TFile::Open(("../test/TICL0PU_timeClean_2best/ticl_22Pt"+ptName.at(ij)+"_"+etaName.at(kl)+"_timeV2_timeClean.root").c_str());

      if(ij == 0){
	scaleVsPt_recoBP[kl] = new TGraphErrors();
	scaleVsPt_recoBP[kl]->SetName(Form("scaleVsPt_recoBP_%d", kl));
	scaleVsPt_recoBP[kl]->SetMarkerStyle(iStyle[kl]);
	scaleVsPt_recoBP[kl]->SetMarkerColor(iColors[0]);

	resolutionVsPt_recoBP[kl] = new TGraphErrors();
	resolutionVsPt_recoBP[kl]->SetName(Form("resolutionVsPt_recoBP_%d", kl));
	resolutionVsPt_recoBP[kl]->SetMarkerStyle(iStyle[kl]);
	resolutionVsPt_recoBP[kl]->SetMarkerColor(iColors[0]);
	//
	scaleVsPt_recoAll[kl] = new TGraphErrors();
	scaleVsPt_recoAll[kl]->SetName(Form("scaleVsPt_recoAll_%d", kl));
	scaleVsPt_recoAll[kl]->SetMarkerStyle(iStyle[kl]);
	scaleVsPt_recoAll[kl]->SetMarkerColor(iColors[1]);

	resolutionVsPt_recoAll[kl] = new TGraphErrors();
	resolutionVsPt_recoAll[kl]->SetName(Form("resolutionVsPt_recoAll_%d", kl));
	resolutionVsPt_recoAll[kl]->SetMarkerStyle(iStyle[kl]);
	resolutionVsPt_recoAll[kl]->SetMarkerColor(iColors[1]);
	//
	scaleVsPt_recoBPmatched[kl] = new TGraphErrors();
	scaleVsPt_recoBPmatched[kl]->SetName(Form("scaleVsPt_recoBPmatched_%d", kl));
	scaleVsPt_recoBPmatched[kl]->SetMarkerStyle(iStyle[kl]);
	scaleVsPt_recoBPmatched[kl]->SetMarkerColor(iColors[2]);

	resolutionVsPt_recoBPmatched[kl] = new TGraphErrors();
	resolutionVsPt_recoBPmatched[kl]->SetName(Form("resolutionVsPt_recoBPmatched_%d", kl));
	resolutionVsPt_recoBPmatched[kl]->SetMarkerStyle(iStyle[kl]);
	resolutionVsPt_recoBPmatched[kl]->SetMarkerColor(iColors[2]);
	//
	scaleVsPt_recoAllmatched[kl] = new TGraphErrors();
	scaleVsPt_recoAllmatched[kl]->SetName(Form("scaleVsPt_recoAllmatched_%d", kl));
	scaleVsPt_recoAllmatched[kl]->SetMarkerStyle(iStyle[kl]);
	scaleVsPt_recoAllmatched[kl]->SetMarkerColor(iColors[3]);

	resolutionVsPt_recoAllmatched[kl] = new TGraphErrors();
	resolutionVsPt_recoAllmatched[kl]->SetName(Form("resolutionVsPt_recoAllmatched_%d", kl));
	resolutionVsPt_recoAllmatched[kl]->SetMarkerStyle(iStyle[kl]);
	resolutionVsPt_recoAllmatched[kl]->SetMarkerColor(iColors[3]);
      }
    }
  }
  std::cout << " >>> fatto = presi " << std::endl;

  std::string folder = "plots";

  TH1F* recoBP_genAll_RatioEnergy[3][2];
  TH1F* recoAll_genAll_RatioEnergy[3][2];
  TH1F* recoBPmatched_genAll_RatioEnergy[3][2];
  TH1F* recoAllmatched_genAll_RatioEnergy[3][2];

  TF1* hfit = new TF1("hfit", "gaus", 0., 1.1);

  for(int ij=0; ij<3; ++ij){
    for(int kl=0; kl<2; ++kl){

      recoBP_genAll_RatioEnergy[ij][kl] = (TH1F*)inF[ij][kl]->Get("ticlAnalyzer/recoBP_genAll_RatioEnergy");
      recoAll_genAll_RatioEnergy[ij][kl] = (TH1F*)inF[ij][kl]->Get("ticlAnalyzer/recoAll_genAll_RatioEnergy");
      recoBPmatched_genAll_RatioEnergy[ij][kl] = (TH1F*)inF[ij][kl]->Get("ticlAnalyzer/recoBPmatched_genAll_RatioEnergy");
      recoAllmatched_genAll_RatioEnergy[ij][kl] = (TH1F*)inF[ij][kl]->Get("ticlAnalyzer/recoAllmatched_genAll_RatioEnergy");

      recoBP_genAll_RatioEnergy[ij][kl]->SetLineColor(iColors[0]);
      recoAll_genAll_RatioEnergy[ij][kl]->SetLineColor(iColors[1]);
      recoBPmatched_genAll_RatioEnergy[ij][kl]->SetLineColor(iColors[2]); 
      recoAllmatched_genAll_RatioEnergy[ij][kl]->SetLineColor(iColors[3]);

      recoBP_genAll_RatioEnergy[ij][kl]->SetLineWidth(2);
      recoAll_genAll_RatioEnergy[ij][kl]->SetLineWidth(2);
      recoBPmatched_genAll_RatioEnergy[ij][kl]->SetLineWidth(2); 
      recoAllmatched_genAll_RatioEnergy[ij][kl]->SetLineWidth(2);

      //
      hfit->SetParameter(1, recoBP_genAll_RatioEnergy[ij][kl]->GetMean());
      hfit->SetParameter(2, recoBP_genAll_RatioEnergy[ij][kl]->GetRMS());
      hfit->SetRange(recoBP_genAll_RatioEnergy[ij][kl]->GetMean() - recoBP_genAll_RatioEnergy[ij][kl]->GetRMS(), 
		     recoBP_genAll_RatioEnergy[ij][kl]->GetMean() + recoBP_genAll_RatioEnergy[ij][kl]->GetRMS());

      recoBP_genAll_RatioEnergy[ij][kl]->Scale(1./recoBP_genAll_RatioEnergy[ij][kl]->Integral());
      recoBP_genAll_RatioEnergy[ij][kl]->GetXaxis()->SetRangeUser(0.5, 1.2);
      recoBP_genAll_RatioEnergy[ij][kl]->Fit("hfit", "Q");
      scaleVsPt_recoBP[kl]->SetPoint(scaleVsPt_recoBP[kl]->GetN(), ptVal.at(ij), hfit->GetParameter(1));
      scaleVsPt_recoBP[kl]->SetPointError(scaleVsPt_recoBP[kl]->GetN()-1, 0, hfit->GetParError(1));
      resolutionVsPt_recoBP[kl]->SetPoint(resolutionVsPt_recoBP[kl]->GetN(), ptVal.at(ij), hfit->GetParameter(2)/hfit->GetParameter(1));
      resolutionVsPt_recoBP[kl]->SetPointError(resolutionVsPt_recoBP[kl]->GetN()-1, 0, hfit->GetParError(2));
      std::cout << " recoBP ij = " << ij << " kl = " << kl << " reso = " << hfit->GetParameter(2)/hfit->GetParameter(1) << std::endl;
      //
      hfit->SetParameter(1, recoAll_genAll_RatioEnergy[ij][kl]->GetMean());
      hfit->SetParameter(2, recoAll_genAll_RatioEnergy[ij][kl]->GetRMS());
      hfit->SetRange(recoAll_genAll_RatioEnergy[ij][kl]->GetMean() - recoAll_genAll_RatioEnergy[ij][kl]->GetRMS(), 
		     recoAll_genAll_RatioEnergy[ij][kl]->GetMean() + recoAll_genAll_RatioEnergy[ij][kl]->GetRMS());

      recoAll_genAll_RatioEnergy[ij][kl]->Scale(1./recoAll_genAll_RatioEnergy[ij][kl]->Integral());
      recoAll_genAll_RatioEnergy[ij][kl]->GetXaxis()->SetRangeUser(0.5, 1.2);
      recoAll_genAll_RatioEnergy[ij][kl]->Fit("hfit", "Q");
      scaleVsPt_recoAll[kl]->SetPoint(scaleVsPt_recoAll[kl]->GetN(), ptVal.at(ij), hfit->GetParameter(1));
      scaleVsPt_recoAll[kl]->SetPointError(scaleVsPt_recoAll[kl]->GetN()-1, 0, hfit->GetParError(1));
      resolutionVsPt_recoAll[kl]->SetPoint(resolutionVsPt_recoAll[kl]->GetN(), ptVal.at(ij), hfit->GetParameter(2)/hfit->GetParameter(1));
      resolutionVsPt_recoAll[kl]->SetPointError(resolutionVsPt_recoAll[kl]->GetN()-1, 0, hfit->GetParError(2));
      std::cout << " recoAll ij = " << ij << " kl = " << kl << " reso = " << hfit->GetParameter(2)/hfit->GetParameter(1) 
		<< " hfit->GetParameter(2) = " << hfit->GetParameter(2) << " hfit->GetParameter(1) = " << hfit->GetParameter(1)<< std::endl;
      //
      hfit->SetParameter(1, recoBPmatched_genAll_RatioEnergy[ij][kl]->GetMean());
      hfit->SetParameter(2, recoBPmatched_genAll_RatioEnergy[ij][kl]->GetRMS());
      hfit->SetRange(recoBPmatched_genAll_RatioEnergy[ij][kl]->GetMean() - recoBPmatched_genAll_RatioEnergy[ij][kl]->GetRMS(), 
		     recoBPmatched_genAll_RatioEnergy[ij][kl]->GetMean() + recoBPmatched_genAll_RatioEnergy[ij][kl]->GetRMS());

      recoBPmatched_genAll_RatioEnergy[ij][kl]->Scale(1./recoBPmatched_genAll_RatioEnergy[ij][kl]->Integral());
      recoBPmatched_genAll_RatioEnergy[ij][kl]->GetXaxis()->SetRangeUser(0.5, 1.2);
      recoBPmatched_genAll_RatioEnergy[ij][kl]->Fit("hfit", "Q");
      scaleVsPt_recoBPmatched[kl]->SetPoint(scaleVsPt_recoBPmatched[kl]->GetN(), ptVal.at(ij), hfit->GetParameter(1));
      scaleVsPt_recoBPmatched[kl]->SetPointError(scaleVsPt_recoBPmatched[kl]->GetN()-1, 0, hfit->GetParError(1));
      resolutionVsPt_recoBPmatched[kl]->SetPoint(resolutionVsPt_recoBPmatched[kl]->GetN(), ptVal.at(ij), hfit->GetParameter(2)/hfit->GetParameter(1));
      resolutionVsPt_recoBPmatched[kl]->SetPointError(resolutionVsPt_recoBPmatched[kl]->GetN()-1, 0, hfit->GetParError(2));
      std::cout << " recoBPmatched ij = " << ij << " kl = " << kl << " reso = " << hfit->GetParameter(2)/hfit->GetParameter(1) << std::endl;
      //
      hfit->SetParameter(1, recoAllmatched_genAll_RatioEnergy[ij][kl]->GetMean());
      hfit->SetParameter(2, recoAllmatched_genAll_RatioEnergy[ij][kl]->GetRMS());
      hfit->SetRange(recoAllmatched_genAll_RatioEnergy[ij][kl]->GetMean() - recoAllmatched_genAll_RatioEnergy[ij][kl]->GetRMS(), 
		     recoAllmatched_genAll_RatioEnergy[ij][kl]->GetMean() + recoAllmatched_genAll_RatioEnergy[ij][kl]->GetRMS());

      recoAllmatched_genAll_RatioEnergy[ij][kl]->Scale(1./recoAllmatched_genAll_RatioEnergy[ij][kl]->Integral());
      recoAllmatched_genAll_RatioEnergy[ij][kl]->GetXaxis()->SetRangeUser(0.5, 1.2);
      recoAllmatched_genAll_RatioEnergy[ij][kl]->Fit("hfit", "Q");
      scaleVsPt_recoAllmatched[kl]->SetPoint(scaleVsPt_recoAllmatched[kl]->GetN(), ptVal.at(ij), hfit->GetParameter(1));
      scaleVsPt_recoAllmatched[kl]->SetPointError(scaleVsPt_recoAllmatched[kl]->GetN()-1, 0, hfit->GetParError(1));
      resolutionVsPt_recoAllmatched[kl]->SetPoint(resolutionVsPt_recoAllmatched[kl]->GetN(), ptVal.at(ij), hfit->GetParameter(2)/hfit->GetParameter(1));
      resolutionVsPt_recoAllmatched[kl]->SetPointError(resolutionVsPt_recoAllmatched[kl]->GetN()-1, 0, hfit->GetParError(2));
      std::cout << " recoAllmatched ij = " << ij << " kl = " << kl << " reso = " << hfit->GetParameter(2)/hfit->GetParameter(1) << std::endl;
    }
  }

  std::cout << " ci sono ora stampo " << std::endl;


  TLegend *legPt = new TLegend(0.3,0.4,0.45,0.6,NULL,"brNDC");
  legPt->SetTextFont(42);
  legPt->SetTextSize(0.04);
  legPt->SetFillColor(kWhite);
  legPt->SetLineColor(kWhite);
  legPt->SetShadowColor(kWhite);

  legPt->AddEntry(recoBP_genAll_RatioEnergy[0][0], "best", "l");
  legPt->AddEntry(recoAll_genAll_RatioEnergy[0][0], "all", "l");
  legPt->AddEntry(recoBPmatched_genAll_RatioEnergy[0][0], "matched in best", "l");
  legPt->AddEntry(recoAllmatched_genAll_RatioEnergy[0][0], "matched in all", "l");

  TLegend *legEta = new TLegend(0.65,0.45,0.8,0.55,NULL,"brNDC");
  legEta->SetTextFont(42);
  legEta->SetTextSize(0.04);
  legEta->SetFillColor(kWhite);
  legEta->SetLineColor(kWhite);
  legEta->SetShadowColor(kWhite);
  legEta->AddEntry(scaleVsPt_recoAll[0], etaName.at(0).c_str(), "p");
  legEta->AddEntry(scaleVsPt_recoAll[1], etaName.at(1).c_str(), "p");



  std::cout << " legends ok  " << std::endl;

  //gStyle->SetOptStat("emr");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas* ch_Ratio[3][2];
  for(int ij=0; ij<3; ++ij){
    for(int kl=0; kl<2; ++kl){
      ch_Ratio[ij][kl] = new TCanvas();
      ch_Ratio[ij][kl]->cd();

      recoBP_genAll_RatioEnergy[ij][kl]->GetYaxis()->SetRangeUser(0., 0.4);
      recoBP_genAll_RatioEnergy[ij][kl]->GetXaxis()->SetTitle("energy fraction");
      recoBP_genAll_RatioEnergy[ij][kl]->Draw("h");
      recoAll_genAll_RatioEnergy[ij][kl]->Draw("h, same");
      recoBPmatched_genAll_RatioEnergy[ij][kl]->Draw("h, same");
      recoAllmatched_genAll_RatioEnergy[ij][kl]->Draw("h, same");

      legPt->SetHeader((ptName.at(ij)+" "+etaName.at(kl)).c_str());
      legPt->Draw("same");
      
      ch_Ratio[ij][kl]->Print((folder+"/Ratio_pt"+ptName.at(ij)+"_eta"+etaName.at(kl)+".png").c_str(), "png");
    }
  }


  TCanvas* ch_Scale = new TCanvas();
  scaleVsPt_recoBP[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
  scaleVsPt_recoBP[0]->GetYaxis()->SetTitle("energy fraction");
  scaleVsPt_recoBP[0]->GetYaxis()->SetRangeUser(0.5, 1.1);
  scaleVsPt_recoBP[0]->Draw("ap");
  for(int ij=0; ij<2; ++ij){
    scaleVsPt_recoBP[ij]->Draw("p, same");
    scaleVsPt_recoAll[ij]->Draw("p, same");
    scaleVsPt_recoBPmatched[ij]->Draw("p, same");
    scaleVsPt_recoAllmatched[ij]->Draw("p, same");
  }
  legPt->SetHeader("");
  legPt->Draw("same");
  legEta->Draw("same");
  ch_Scale->Print((folder+"/Scale.png").c_str(), "png");


  TCanvas* ch_Resolution = new TCanvas();
  resolutionVsPt_recoBP[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
  resolutionVsPt_recoBP[0]->GetYaxis()->SetTitle("#sigma(E)/E");
  resolutionVsPt_recoBP[0]->GetYaxis()->SetRangeUser(0., 0.1);
  resolutionVsPt_recoBP[0]->Draw("ap");
  for(int ij=0; ij<2; ++ij){
    resolutionVsPt_recoBP[ij]->Draw("p, same");
    resolutionVsPt_recoAll[ij]->Draw("p, same");
    resolutionVsPt_recoBPmatched[ij]->Draw("p, same");
    resolutionVsPt_recoAllmatched[ij]->Draw("p, same");
  }
  legPt->SetHeader("");
  legPt->Draw("same");
  legEta->Draw("same");
  ch_Resolution->Print((folder+"/Resolution.png").c_str(), "png");


  return;

}
