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

void drawHistos(){
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
  int iColors[3] = {kBlue, kRed+1, kGreen+2};
  int iStyle[3] = {1, 2, 5}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  
  std::vector<std::string> fileName;
  fileName.push_back("noTime_all");
  fileName.push_back("noTime_default_all");
  //  fileName.push_back("withTime_all");
  fileName.push_back("withTime_loose_noBH");
  //fileName.push_back("withTime_default_all");
  fileName.push_back("withTime_default_noBH");
  fileName.push_back("noTime_default_trkIt");
  //  fileName.push_back("withTime_default_trkIt");
  fileName.push_back("withTime_trkIt_noBH");


  std::vector<std::string> legNameHead;
  legNameHead.push_back("tracking (KF) [loose CA]");
  legNameHead.push_back("tracking (KF)");
  legNameHead.push_back("tracking (KF) [loose CA + time]");
  legNameHead.push_back("tracking (KF) [time]");
  legNameHead.push_back("tracking (CA)");
  legNameHead.push_back("tracking (CA) [time]");


  TFile* inF[6];
  for(int ij=0; ij<6; ++ij){
    inF[ij] = TFile::Open(("../test/ticl_analysis_BuildP_211_"+fileName.at(ij)+".root").c_str());
  }
  std::cout << " >>> fatto = presi " << std::endl;


  std::vector<std::string> legName;
  legName.push_back("all");
  legName.push_back("best reco (dR)");
  legName.push_back("best reco (dR && N)");

  std::string folder = "plots";


  TH1F* hEFraction_all[6];
  TH1F* hEFraction_allV2[6];
  TH1F* hEFraction_allV2_cpAsreco[6];


  for(int ij=0; ij<6; ++ij){
    hEFraction_all[ij] = (TH1F*)inF[ij]->Get("ticlAnalyzer/EFraction_all");
    //hEFraction_all[ij]->SetName(Form("EFraction_all_%d", ij));
    hEFraction_all[ij]->SetLineColor(iColors[0]);
    hEFraction_all[ij]->SetLineStyle(iStyle[0]);
    hEFraction_all[ij]->SetLineWidth(2);

    hEFraction_allV2[ij] = (TH1F*)inF[ij]->Get("ticlAnalyzer/EFraction_allV2");
    // hEFraction_allV2[ij]->SetName(Form("EFraction_allV2_%d", ij));
    hEFraction_allV2[ij]->SetLineColor(iColors[1]);
    hEFraction_allV2[ij]->SetLineStyle(iStyle[1]);
    hEFraction_allV2[ij]->SetLineWidth(2);

    hEFraction_allV2_cpAsreco[ij] = (TH1F*)inF[ij]->Get("ticlAnalyzer/EFraction_allV2_cpAsreco");
    //hEFraction_allV2_cpAsreco[ij]->SetName(Form("EFraction_allV2_cpAsreco_%d", ij));
    hEFraction_allV2_cpAsreco[ij]->SetLineColor(iColors[2]);
    hEFraction_allV2_cpAsreco[ij]->SetLineStyle(iStyle[2]);
    hEFraction_allV2_cpAsreco[ij]->SetLineWidth(2);
  }


  std::cout << " ci sono ora stampo " << std::endl;


  TLegend *legTGM = new TLegend(0.32,0.50,0.5,0.7,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.04);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  //  for(int iF=0; iF<3; ++iF){
  legTGM->AddEntry(hEFraction_all[0], legName.at(0).c_str(), "l");
  legTGM->AddEntry(hEFraction_allV2[0], legName.at(1).c_str(), "l");
  legTGM->AddEntry(hEFraction_allV2_cpAsreco[0], legName.at(2).c_str(), "l");


  std::cout << " legends ok  " << std::endl;

  gStyle->SetOptStat("emr");

  TCanvas* ch_Histo[6];
  for(int ij=0; ij<6; ++ij){
    ch_Histo[ij] = new TCanvas();
    ch_Histo[ij]->cd();
    //gPad->SetLogy();

    hEFraction_all[ij]->GetXaxis()->SetTitle("energy fraction (reco/gen)");
    hEFraction_all[ij]->GetYaxis()->SetRangeUser(0., 800.);
    //    hTime_Eta_dRadius[ij]->GetYaxis()->SetRangeUser(0.1, 1.e+5);
    hEFraction_all[ij]->Draw("h");
    hEFraction_allV2[ij]->Draw("h, sames");
    hEFraction_allV2_cpAsreco[ij]->Draw("h,  sames");
    gPad->Update();
    ch_Histo[ij]->Update();

    
    TPaveStats* p1 = (TPaveStats*)hEFraction_all[ij]->FindObject("stats");
    p1->SetStatFormat(".2f");
    p1->SetX1NDC(0.22);
    p1->SetX2NDC(0.4);
    p1->SetY1NDC(0.78);
    p1->SetY2NDC(0.93);
    p1->SetTextColor(iColors[0]);
    p1->SetLineColor(iColors[0]);
    p1->SetLineStyle(iStyle[0]);
    p1->SetLineWidth(2);
    p1->SetTextSize(0.03);
    
    
    TPaveStats* p2 = (TPaveStats*)hEFraction_allV2[ij]->GetListOfFunctions()->FindObject("stats");
    p2->SetStatFormat(".2f");
    p2->SetX1NDC(0.4);
    p2->SetX2NDC(0.58);
    p2->SetY1NDC(0.78);
    p2->SetY2NDC(0.93);
    p2->SetTextColor(iColors[1]);
    p2->SetLineColor(iColors[1]);
    p2->SetLineStyle(iStyle[1]);
    p2->SetLineWidth(2);
    p2->SetTextSize(0.03);

    TPaveStats* p3 = (TPaveStats*)hEFraction_allV2_cpAsreco[ij]->FindObject("stats");
    p3->SetStatFormat(".2f");
    p3->SetX1NDC(0.58);
    p3->SetX2NDC(0.76);
    p3->SetY1NDC(0.78);
    p3->SetY2NDC(0.93);
    p3->SetTextColor(iColors[2]);
    p3->SetLineColor(iColors[2]);
    p3->SetLineStyle(iStyle[2]);
    p3->SetLineWidth(2);
    p3->SetTextSize(0.03);


    legTGM->SetHeader(legNameHead.at(ij).c_str(), "ticlAnalyzer");    
    legTGM->Draw("same");

    ch_Histo[ij]->Print((folder+"h_Histo"+fileName.at(ij)+".png").c_str(), "png");
    ch_Histo[ij]->Print((folder+"h_Histo"+fileName.at(ij)+".pdf").c_str(), "pdf");
    // ch_Histo[ij]->Print((folder+"h_Histo"+fileName.at(ij)+".png").c_str(), "png");
  }


  return;

}
