#include <iostream>
#include "TString.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"

TString particle="Bplus";
const int nbinsRFB = 2;
Double_t xbinsRFB[nbinsRFB]={0.5,1.465};
Double_t exlRFB[nbinsRFB]={0.5,0.465};
Double_t yRFBSystTotHighRel[nbinsRFB]= {0.162,0.196};
Double_t yRFBSystTotLowRel[nbinsRFB]= {0.162,0.196};
/*
const int nbins = 5;
Double_t xbins[nbins]={-2.4,-1.470,-0.535,0.465,1.465};
Double_t exl[nbins]={0.465,0.465,0.5,0.5,0.465};
*/
const int nbins = 5;
Double_t xbins[nbins]={-2.4,-1.465,-0.5,0.5,1.465};
Double_t exl[nbins]={0.465,0.465,0.5,0.5,0.465};

Double_t yPercSigmapPbSystTotHigh[nbins]={0.146,0.143,0.143,0.141,0.141};
Double_t yPercSigmapPbSystTotLow[nbins]={0.146,0.143,0.143,0.141,0.141};

Double_t commonErrorP = TMath::Sqrt(0.035*0.035+0.031*0.031);
Double_t commonErrorN = TMath::Sqrt(0.035*0.035+0.031*0.031);
Double_t tagandprobcorrection[nbins]={1.054,1.044,1.037,1.034,1.036};

void NuclearModificationY_v2(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
/*
  gStyle->SetPadRightMargin(0.010);//###0.020
  gStyle->SetPadLeftMargin(0.165);
  gStyle->SetPadTopMargin(0.075);
  gStyle->SetPadBottomMargin(0.145);
*/

    
  //TFile*filePPReference=new TFile(Form("../../../fonll/output%sY.root",particle.Data()));  
  TFile*filePPReference=new TFile("../../../fonll/Code/Compdsigmady_ScaleCor_v3.root");  


  TGraphAsymmErrors*gaeBplusReference=(TGraphAsymmErrors*)filePPReference->Get(Form("gaeSigmaDecayv2%s",particle.Data()));
  gaeBplusReference->SetName(Form("gae%sReference",particle.Data()));
  
  TH1D* hbin_ATL7 = (TH1D*)filePPReference->Get("hbin_ATL7");
  TH1D* hbin_CMS7 = (TH1D*)filePPReference->Get("hbin_CMS7");

  TFile*filepPb=new TFile(Form("../Results%s_y/Sigma%s.root",particle.Data(),particle.Data()));
  TH1F*hSigmapPbStat=(TH1F*)filepPb->Get("hPtSigma");  
  
  double yvalue,xvalue,yerrorhigh,yerrorlow;
  double scalingfactor=1.0e-6;
    
  for (int i=0;i<nbins;i++){
    hSigmapPbStat->SetBinContent(i+1,scalingfactor*hSigmapPbStat->GetBinContent(i+1));
    hSigmapPbStat->SetBinError(i+1,scalingfactor*hSigmapPbStat->GetBinError(i+1));
    
    yvalue=-1.;
    xvalue=-1.;
    yerrorhigh=-1.;
    yerrorlow=-1.;
  } 
  
  for (int i=0;i<nbins;i++){
    hSigmapPbStat->SetBinContent(i+1,(1./tagandprobcorrection[i])*(hSigmapPbStat->GetBinContent(i+1)));
    hSigmapPbStat->SetBinError(i+1,(1./tagandprobcorrection[i])*(hSigmapPbStat->GetBinError(i+1)));
  } 

  
  Double_t yRefPP[nbins];                        //value y reference
  Double_t xRefPP[nbins];                        //value x reference
  Double_t yPPsystFONLLhigh[nbins];              //y err syst FONLL high
  Double_t yPPsystFONLLlow[nbins];               //y err syst FONLL low
  Double_t yPercPPsystFONLLhigh[nbins];          //y percentuale err syst FONLL high
  Double_t yPercPPsystFONLLlow[nbins];           //y percentuale err syst FONLL low
    
  Double_t ySigmapPb[nbins];                     //value y pPb 
  Double_t xSigmapPb[nbins];                     //value x pPb
  Double_t ySigmapPbStat[nbins];                 //y err stat pPb
  Double_t yPercSigmapPbStat[nbins];             //y err stat pPb
  
  Double_t yFONLL[nbins];                        //1
  Double_t yRpA[nbins];                          //value y RpA 
  Double_t yRpAStat[nbins];                      //y err stat RpA 
  Double_t yRpAsystFONLLhigh[nbins];             //y err syst FONLL RpA high
  Double_t yRpAsystFONLLlow[nbins];              //y err syst FONLL RpA lzow
  Double_t yPercRpAsystFONLLhigh[nbins];         //y percentuale err syst FONLL RpA high
  Double_t yPercRpAsystFONLLlow[nbins];          //y percentuale err syst FONLL RpA low
    
  Double_t ySigmapPbSystTotHigh[nbins];              //y percentuale err syst pPb TOT
  Double_t ySigmapPbSystTotLow[nbins];              //y percentuale err syst pPb TOT

  Double_t yPercRpPbSystTotHigh[nbins];          //y percentuale err syst RpPb TOT
  Double_t yPercRpPbSystTotLow[nbins];          //y percentuale err syst RpPb TOT
  
  Double_t yRpPbSystTotHigh[nbins];              //y percentuale err syst RpPb TOT
  Double_t yRpPbSystTotLow[nbins];              //y percentuale err syst RpPb TOT

  Double_t yRFB[nbinsRFB];
  Double_t yRFBStat[nbinsRFB];
  Double_t yRFBSystTotHigh[nbinsRFB];
  Double_t yRFBSystTotLow[nbinsRFB];

//CMS
  Double_t yRefPP_CMS7[nbins];                        //value y reference
  //Double_t xRefPP_CMS7[nbins];                        //value x reference
  Double_t yPPsystFONLLhigh_CMS7[nbins];              //y err syst FONLL high
  Double_t yPPsystFONLLlow_CMS7[nbins];               //y err syst FONLL low
  Double_t yPercPPsystFONLLhigh_CMS7[nbins];          //y percentuale err syst FONLL high
  Double_t yPercPPsystFONLLlow_CMS7[nbins];           //y percentuale err syst FONLL low
    
  Double_t yFONLL_CMS7[nbins];                        //1
  Double_t yRpA_CMS7[nbins];                          //value y RpA 
  Double_t yRpAStat_CMS7[nbins];                      //y err stat RpA 
  Double_t yRpAsystFONLLhigh_CMS7[nbins];             //y err syst FONLL RpA high
  Double_t yRpAsystFONLLlow_CMS7[nbins];              //y err syst FONLL RpA lzow
  Double_t yPercRpAsystFONLLhigh_CMS7[nbins];         //y percentuale err syst FONLL RpA high
  Double_t yPercRpAsystFONLLlow_CMS7[nbins];          //y percentuale err syst FONLL RpA low
    
  Double_t yRpPbSystTotHigh_CMS7[nbins];              //y percentuale err syst RpPb TOT
  Double_t yRpPbSystTotLow_CMS7[nbins];              //y percentuale err syst RpPb TOT
//ATLAS
  Double_t yRefPP_ATL7[nbins];                        //value y reference
  //Double_t xRefPP_ATL7[nbins];                        //value x reference
  Double_t yPPsystFONLLhigh_ATL7[nbins];              //y err syst FONLL high
  Double_t yPPsystFONLLlow_ATL7[nbins];               //y err syst FONLL low
  Double_t yPercPPsystFONLLhigh_ATL7[nbins];          //y percentuale err syst FONLL high
  Double_t yPercPPsystFONLLlow_ATL7[nbins];           //y percentuale err syst FONLL low
    
  Double_t yFONLL_ATL7[nbins];                        //1
  Double_t yRpA_ATL7[nbins];                          //value y RpA 
  Double_t yRpAStat_ATL7[nbins];                      //y err stat RpA 
  Double_t yRpAsystFONLLhigh_ATL7[nbins];             //y err syst FONLL RpA high
  Double_t yRpAsystFONLLlow_ATL7[nbins];              //y err syst FONLL RpA lzow
  Double_t yPercRpAsystFONLLhigh_ATL7[nbins];         //y percentuale err syst FONLL RpA high
  Double_t yPercRpAsystFONLLlow_ATL7[nbins];          //y percentuale err syst FONLL RpA low
    
  Double_t yRpPbSystTotHigh_ATL7[nbins];              //y percentuale err syst RpPb TOT
  Double_t yRpPbSystTotLow_ATL7[nbins];              //y percentuale err syst RpPb TOT


  double x,y;
// FONLL
  for (Int_t i=0;i<nbins;i++) {
    gaeBplusReference->GetPoint(i,xRefPP[i],yRefPP[i]);
    yPPsystFONLLhigh[i]=gaeBplusReference->GetEYhigh()[i];
    yPPsystFONLLlow[i]=gaeBplusReference->GetEYlow()[i];
    yPercPPsystFONLLhigh[i]=yPPsystFONLLhigh[i]/yRefPP[i];
    yPercPPsystFONLLlow[i]=yPPsystFONLLlow[i]/yRefPP[i];
  }
  
  for(Int_t i=0;i<nbins;i++) {
    ySigmapPb[i]=hSigmapPbStat->GetBinContent(i+1);
    ySigmapPbStat[i]=hSigmapPbStat->GetBinError(i+1);
    yPercSigmapPbStat[i]=ySigmapPbStat[i]/ySigmapPb[i];
    ySigmapPbSystTotHigh[i]=yPercSigmapPbSystTotHigh[i]*ySigmapPb[i];
    ySigmapPbSystTotLow[i]=yPercSigmapPbSystTotLow[i]*ySigmapPb[i];
  }
  
  for(Int_t i=0;i<nbins;i++) {
    yRpA[i]=ySigmapPb[i]/yRefPP[i];
    yRpAStat[i]=ySigmapPbStat[i]/yRefPP[i];
    yFONLL[i]=yRpA[i];
    yPercRpAsystFONLLhigh[i]=(yPercPPsystFONLLlow[i]/(1-yPercPPsystFONLLlow[i]));
    yPercRpAsystFONLLlow[i]=(yPercPPsystFONLLhigh[i]/(1+yPercPPsystFONLLhigh[i]));
    yRpAsystFONLLhigh[i]=yPercRpAsystFONLLhigh[i]*yRpA[i];
    yRpAsystFONLLlow[i]=yPercRpAsystFONLLlow[i]*yRpA[i];
    yRpPbSystTotHigh[i]=yPercSigmapPbSystTotHigh[i]*yRpA[i];
    yRpPbSystTotLow[i]=yPercSigmapPbSystTotLow[i]*yRpA[i];
  }
// CMS
  for (Int_t i=0;i<nbins;i++) {
    yRefPP_CMS7[i]=hbin_CMS7->GetBinContent(i+1);
    yPPsystFONLLhigh_CMS7[i]=hbin_CMS7->GetBinError(i+1);
    yPPsystFONLLlow_CMS7[i]=hbin_CMS7->GetBinError(i+1);
    yPercPPsystFONLLhigh_CMS7[i]=yPPsystFONLLhigh_CMS7[i]/yRefPP_CMS7[i];
    yPercPPsystFONLLlow_CMS7[i]=yPPsystFONLLlow_CMS7[i]/yRefPP_CMS7[i];
  }
  
  for(Int_t i=0;i<nbins;i++) {
    yRpA_CMS7[i]=ySigmapPb[i]/yRefPP_CMS7[i];
    yRpAStat_CMS7[i]=ySigmapPbStat[i]/yRefPP_CMS7[i];
    yFONLL_CMS7[i]=yRpA_CMS7[i];
    yPercRpAsystFONLLhigh_CMS7[i]=(yPercPPsystFONLLlow_CMS7[i]/(1-yPercPPsystFONLLlow_CMS7[i]));
    yPercRpAsystFONLLlow_CMS7[i]=(yPercPPsystFONLLhigh_CMS7[i]/(1+yPercPPsystFONLLhigh_CMS7[i]));
    yRpAsystFONLLhigh_CMS7[i]=yPercRpAsystFONLLhigh_CMS7[i]*yRpA_CMS7[i];
    yRpAsystFONLLlow_CMS7[i]=yPercRpAsystFONLLlow_CMS7[i]*yRpA_CMS7[i];
    yRpPbSystTotHigh_CMS7[i]=yPercSigmapPbSystTotHigh[i]*yRpA_CMS7[i];
    yRpPbSystTotLow_CMS7[i]=yPercSigmapPbSystTotLow[i]*yRpA_CMS7[i];
  }
// ATLAS
  for (Int_t i=0;i<nbins;i++) {
    yRefPP_ATL7[i]=hbin_ATL7->GetBinContent(i+1);
    yPPsystFONLLhigh_ATL7[i]=hbin_ATL7->GetBinError(i+1);
    yPPsystFONLLlow_ATL7[i]=hbin_ATL7->GetBinError(i+1);
    yPercPPsystFONLLhigh_ATL7[i]=yPPsystFONLLhigh_ATL7[i]/yRefPP_ATL7[i];
    yPercPPsystFONLLlow_ATL7[i]=yPPsystFONLLlow_ATL7[i]/yRefPP_ATL7[i];
  }
  
  for(Int_t i=0;i<nbins;i++) {
    yRpA_ATL7[i]=ySigmapPb[i]/yRefPP_ATL7[i];
    yRpAStat_ATL7[i]=ySigmapPbStat[i]/yRefPP_ATL7[i];
    yFONLL_ATL7[i]=yRpA_ATL7[i];
    yPercRpAsystFONLLhigh_ATL7[i]=(yPercPPsystFONLLlow_ATL7[i]/(1-yPercPPsystFONLLlow_ATL7[i]));
    yPercRpAsystFONLLlow_ATL7[i]=(yPercPPsystFONLLhigh_ATL7[i]/(1+yPercPPsystFONLLhigh_ATL7[i]));
    yRpAsystFONLLhigh_ATL7[i]=yPercRpAsystFONLLhigh_ATL7[i]*yRpA_ATL7[i];
    yRpAsystFONLLlow_ATL7[i]=yPercRpAsystFONLLlow_ATL7[i]*yRpA_ATL7[i];
    yRpPbSystTotHigh_ATL7[i]=yPercSigmapPbSystTotHigh[i]*yRpA_ATL7[i];
    yRpPbSystTotLow_ATL7[i]=yPercSigmapPbSystTotLow[i]*yRpA_ATL7[i];
  }



  //RFB stuff
  //central value
  yRFB[0] =  ySigmapPb[3]/ySigmapPb[2];
  yRFB[1] =  ySigmapPb[4]/ySigmapPb[1];
  //Stat.
  yRFBStat[0] =  yRFB[0]*sqrt(pow(ySigmapPbStat[2]/ySigmapPb[2],2)+pow(ySigmapPbStat[3]/ySigmapPb[3],2));
  yRFBStat[1] =  yRFB[1]*sqrt(pow(ySigmapPbStat[4]/ySigmapPb[4],2)+pow(ySigmapPbStat[1]/ySigmapPb[1],2));
  //Syst.

  yRFBSystTotHigh[0] = yRFB[0]*yRFBSystTotHighRel[0];
  yRFBSystTotLow[0] = yRFB[0]*yRFBSystTotLowRel[0];
  yRFBSystTotHigh[1] = yRFB[1]*yRFBSystTotHighRel[1];
  yRFBSystTotLow[1] = yRFB[1]*yRFBSystTotLowRel[1];
    
  TGraphAsymmErrors *gSigmasyst = new TGraphAsymmErrors(nbins,xbins,ySigmapPb,exl,exl,ySigmapPbSystTotLow,ySigmapPbSystTotHigh);
  gSigmasyst->SetTitle("Sigma syst uncertainty from pPb");
  gSigmasyst->SetMarkerColor(1);
  gSigmasyst->SetLineColor(1);
  gSigmasyst->SetLineWidth(1);   
  gSigmasyst->SetMarkerStyle(21);
  gSigmasyst->SetMarkerColor(1);

  TGraphAsymmErrors *gSigmastat = new TGraphAsymmErrors(nbins,xbins,ySigmapPb,exl,exl,ySigmapPbStat,ySigmapPbStat);
  gSigmastat->SetTitle("Sigma stat uncertainty from pPb");
  gSigmastat->SetMarkerColor(1);
  gSigmastat->SetLineColor(1);
  gSigmastat->SetLineWidth(1);   
  gSigmastat->SetMarkerStyle(21);
  gSigmastat->SetMarkerColor(1);
  
  gSigmastat->SetFillColor(0);
  gSigmastat->SetFillStyle(0);


  TCanvas *canvasSigma=new TCanvas("canvasSigma","canvasSigma",600,500);   
  canvasSigma->cd();
  canvasSigma->Range(-1.989924,-0.2917772,25.49622,2.212202);
  canvasSigma->SetFillColor(0);
  canvasSigma->SetBorderMode(0);
  canvasSigma->SetBorderSize(2);
  canvasSigma->SetLeftMargin(0.200);
  canvasSigma->SetRightMargin(0.025);
  canvasSigma->SetTopMargin(0.080);
  canvasSigma->SetBottomMargin(0.150);
  canvasSigma->SetFrameBorderMode(0);
  
  TH2F* hempty=new TH2F("hempty","",10,-3.5,2.5,10,0,650.);  
  hempty->GetXaxis()->SetTitle("y_{CM}");
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetYaxis()->SetTitle("d#sigma / dy_{CM}(#mub GeV^{-1}c)");
  hempty->GetXaxis()->SetTitleOffset(0.90);
  hempty->GetYaxis()->SetTitleOffset(1.34);
  hempty->GetXaxis()->SetTitleSize(0.070);//###0.045
  hempty->GetYaxis()->SetTitleSize(0.070);//###0.045
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.060);//###0.037
  hempty->GetYaxis()->SetLabelSize(0.060);//###0.037  
  hempty->SetMaximum(2);
  hempty->SetMinimum(0.);
  hempty->Draw();
    
  
  gaeBplusReference->SetMarkerColor(1);
  gaeBplusReference->SetMarkerStyle(25);  
  gaeBplusReference->SetFillColor(5);
  gaeBplusReference->SetFillStyle(1001);
  gaeBplusReference->SetLineColor(5);
  gaeBplusReference->SetLineWidth(1);
  gaeBplusReference->Draw("2psame");

  
  gSigmastat->SetMarkerColor(1);
  gSigmastat->SetLineColor(1);
  gSigmastat->SetLineWidth(2);   
  gSigmastat->SetMarkerStyle(21);
  gSigmastat->SetMarkerColor(1);
  gSigmastat->SetFillColor(0);
  gSigmastat->Draw("epsame");
  
  
  //###TLegend *legendSigma=new TLegend(0.233871,0.6701903,0.5322581,0.7526427,"");
  TLegend *legendSigma=new TLegend(0.50,0.60,0.80,0.68,"");
  legendSigma->SetBorderSize(0);
  legendSigma->SetLineColor(0);
  legendSigma->SetFillColor(0);
  legendSigma->SetFillStyle(1001);
  legendSigma->SetTextFont(42);
  legendSigma->SetTextSize(0.055);//###0.045

  TBox *c = new TBox(3,1-commonErrorN,7,1+commonErrorP);
  c->SetLineColor(5);
  c->SetFillColor(5);
  c->Draw();
  
  TLegendEntry *ent_SigmapPb=legendSigma->AddEntry(gSigmastat,"pPb","pf");
  ent_SigmapPb->SetTextFont(42);
  ent_SigmapPb->SetLineColor(1);
  ent_SigmapPb->SetMarkerColor(1);

  TLegendEntry *ent_Sigmapp=legendSigma->AddEntry(gaeBplusReference,"FONLL pp ref.","fp");
  ent_Sigmapp->SetTextFont(42);
  ent_Sigmapp->SetLineColor(5);
  ent_Sigmapp->SetMarkerColor(1);
  ent_Sigmapp->SetMarkerStyle(25);

  legendSigma->Draw("psame");
  gSigmasyst->SetFillColor(0);
  gSigmasyst->SetFillStyle(0);
  gSigmasyst->Draw("2same");
  
  TBox *d = new TBox(3,1-commonErrorN,7,1+commonErrorP);
  d->SetLineColor(1);
  d->SetFillColor(0);
  d->Draw();

  //###TLatex * tlatex1=new TLatex(0.1612903,0.8625793,"CMS                 pPb #sqrt{s_{NN}}= 5.02 TeV");
  TLatex * tlatex1=new TLatex(0.23,0.84,"CMS");
  tlatex1->SetNDC();
  tlatex1->SetTextColor(1);
  tlatex1->SetTextFont(62);//###42
  tlatex1->SetTextSize(0.07);//###0.045
  tlatex1->Draw();
  

  TString mypar="B^{+}";
  
  //###TLatex * tlatexlumi=new TLatex(0.671371,0.7801268,"L = 34.8 nb^{-1}");
  TLatex * tlatexlumi=new TLatex(0.40,0.94,"34.6 nb^{-1} (pPb 5.02 TeV)");
  tlatexlumi->SetNDC();
  tlatexlumi->SetTextColor(1);
  tlatexlumi->SetTextFont(42);
  tlatexlumi->SetTextSize(0.06);//###0.045
  tlatexlumi->Draw();

  //###double xpos=0.8528226;
  //###double ypos=0.6849894;
  double xpos=0.85;
  double ypos=0.75;
  
  TLatex * tlatex3=new TLatex(xpos,ypos,mypar.Data());
   tlatex3->SetNDC();
  tlatex3->SetTextColor(1);
  tlatex3->SetTextFont(42);
  tlatex3->SetTextSize(0.07);//###0.06
  tlatex3->Draw();
  
  TLatex * tlatex4=new TLatex(0.51,0.84,"10 < p_{T}^{B} < 60 GeV/c");
  tlatex4->SetNDC();
  tlatex4->SetTextColor(1);
  tlatex4->SetTextFont(42);
  tlatex4->SetTextSize(0.06);
  tlatex4->Draw();
 
  canvasSigma->SaveAs(Form("../Results%s_y/canvasSigma%s.pdf",particle.Data(),particle.Data()));  
  
  TGraphAsymmErrors *gRpAstat = new TGraphAsymmErrors(nbins,xbins,yRpA,exl,exl,yRpAStat,yRpAStat);
  gRpAstat->SetTitle("RpA stat uncertainty from pPb");
  gRpAstat->SetMarkerStyle(21);
  gRpAstat->SetMarkerColor(1);
  gRpAstat->SetLineColor(1);
  gRpAstat->SetLineWidth(2);  
  gRpAstat->SetFillColor(0);
  
  TGraphAsymmErrors *gRpAsyst = new TGraphAsymmErrors(nbins,xbins,yRpA,exl,exl,yRpPbSystTotLow,yRpPbSystTotHigh);
  gRpAsyst->SetTitle("RpA syst uncertainty from pPb");
  gRpAsyst->SetFillColor(0);
  gRpAsyst->SetMarkerSize(0);
  gRpAsyst->SetLineColor(1);
  gRpAsyst->SetLineWidth(2);
  gRpAsyst->SetFillStyle(0);

   
  TGraphAsymmErrors *gRpAsystFONLL = new TGraphAsymmErrors(nbins,xbins,yFONLL,exl,exl,yRpAsystFONLLlow,yRpAsystFONLLhigh);
  gRpAsystFONLL->SetTitle("RpA syst uncertainty from FONLL reference");
  gRpAsystFONLL->SetTitle("RpA syst uncertainty from FONLL reference");
  gRpAsystFONLL->SetFillColor(5);
  gRpAsystFONLL->SetLineColor(5);//kAzure-3);
  gRpAsystFONLL->SetMarkerColor(4);//kAzure-3);

//CMS
  TGraphAsymmErrors *gRpAstat_CMS7 = new TGraphAsymmErrors(nbins,xbins,yRpA_CMS7,exl,exl,yRpAStat_CMS7,yRpAStat_CMS7);
  gRpAstat_CMS7->SetTitle("RpA stat uncertainty from pPb");
  gRpAstat_CMS7->SetMarkerStyle(21);
  gRpAstat_CMS7->SetMarkerColor(1);
  gRpAstat_CMS7->SetLineColor(1);
  gRpAstat_CMS7->SetLineWidth(2);  
  gRpAstat_CMS7->SetFillColor(0);
  
  TGraphAsymmErrors *gRpAsyst_CMS7 = new TGraphAsymmErrors(nbins,xbins,yRpA_CMS7,exl,exl,yRpPbSystTotLow_CMS7,yRpPbSystTotHigh_CMS7);
  gRpAsyst_CMS7->SetTitle("RpA syst uncertainty from pPb");
  gRpAsyst_CMS7->SetFillColor(0);
  gRpAsyst_CMS7->SetMarkerSize(0);
  gRpAsyst_CMS7->SetLineColor(1);
  gRpAsyst_CMS7->SetLineWidth(2);
  gRpAsyst_CMS7->SetFillStyle(0);

   
  TGraphAsymmErrors *gRpAsystFONLL_CMS7 = new TGraphAsymmErrors(nbins,xbins,yFONLL_CMS7,exl,exl,yRpAsystFONLLlow_CMS7,yRpAsystFONLLhigh_CMS7);
  gRpAsystFONLL_CMS7->SetTitle("RpA syst uncertainty from FONLL reference");
  gRpAsystFONLL_CMS7->SetTitle("RpA syst uncertainty from FONLL reference");
  gRpAsystFONLL_CMS7->SetFillColor(5);
  gRpAsystFONLL_CMS7->SetLineColor(5);//kAzure-3);
  gRpAsystFONLL_CMS7->SetMarkerColor(4);//kAzure-3);

//ATLAS
  TGraphAsymmErrors *gRpAstat_ATL7 = new TGraphAsymmErrors(nbins,xbins,yRpA_ATL7,exl,exl,yRpAStat_ATL7,yRpAStat_ATL7);
  gRpAstat_ATL7->SetTitle("RpA stat uncertainty from pPb");
  gRpAstat_ATL7->SetMarkerStyle(21);
  gRpAstat_ATL7->SetMarkerColor(1);
  gRpAstat_ATL7->SetLineColor(1);
  gRpAstat_ATL7->SetLineWidth(2);  
  gRpAstat_ATL7->SetFillColor(0);
  
  TGraphAsymmErrors *gRpAsyst_ATL7 = new TGraphAsymmErrors(nbins,xbins,yRpA_ATL7,exl,exl,yRpPbSystTotLow_ATL7,yRpPbSystTotHigh_ATL7);
  gRpAsyst_ATL7->SetTitle("RpA syst uncertainty from pPb");
  gRpAsyst_ATL7->SetFillColor(0);
  gRpAsyst_ATL7->SetMarkerSize(0);
  gRpAsyst_ATL7->SetLineColor(1);
  gRpAsyst_ATL7->SetLineWidth(2);
  gRpAsyst_ATL7->SetFillStyle(0);

   
  TGraphAsymmErrors *gRpAsystFONLL_ATL7 = new TGraphAsymmErrors(nbins,xbins,yFONLL_ATL7,exl,exl,yRpAsystFONLLlow_ATL7,yRpAsystFONLLhigh_ATL7);
  gRpAsystFONLL_ATL7->SetTitle("RpA syst uncertainty from FONLL reference");
  gRpAsystFONLL_ATL7->SetTitle("RpA syst uncertainty from FONLL reference");
  gRpAsystFONLL_ATL7->SetFillColor(5);
  gRpAsystFONLL_ATL7->SetLineColor(5);//kAzure-3);
  gRpAsystFONLL_ATL7->SetMarkerColor(4);//kAzure-3);




  TCanvas *canvasRpA=new TCanvas("canvasRpA","canvasRpA",600,500);   
  
  canvasRpA->Range(-1.989924,-0.2917772,25.49622,2.212202);
  canvasRpA->SetFillColor(0);
  canvasRpA->SetBorderMode(0);
  canvasRpA->SetBorderSize(2);
  canvasRpA->SetLeftMargin(0.200);
  canvasRpA->SetRightMargin(0.025);
  canvasRpA->SetTopMargin(0.080);
  canvasRpA->SetBottomMargin(0.150);
  canvasRpA->SetFrameBorderMode(0);
  canvasRpA->SetFrameBorderMode(0);
  
  TLegend *legendRpA=new TLegend(0.23,0.56,0.57,0.75,"");
  legendRpA->SetBorderSize(0);
  legendRpA->SetLineColor(0);
  legendRpA->SetFillColor(0);
  legendRpA->SetFillStyle(1001);
  legendRpA->SetTextFont(42);
  legendRpA->SetTextSize(0.055);//###0.04

  //TH2F* hempty=new TH2F("hempty","",4,-3.5,2.7,10.,0.,3.5);  
//  hempty=new TH2F("hempty","",4,-3.5,2.7,10.,0.,3.5);  
  hempty=new TH2F("hempty","",4,-3.5,2.7,10.,0.,5.0);  



  hempty->GetXaxis()->SetTitle("y_{CM}");
  if(particle=="Bplus") hempty->GetYaxis()->SetTitle("R^{FONLL}_{pA}  (B^{+})");
  if(particle=="Bzero") hempty->GetYaxis()->SetTitle("R^{FONLL}_{pA}  (B^{0})");
  if(particle=="Bs") hempty->GetYaxis()->SetTitle("R^{FONLL}_{pA}  (B_{s})");
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(0.90);
  hempty->GetYaxis()->SetTitleOffset(1.34);
  hempty->GetXaxis()->SetTitleSize(0.070);//###0.045
  hempty->GetYaxis()->SetTitleSize(0.070);//###0.045
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.060);//###0.040
  hempty->GetYaxis()->SetLabelSize(0.060);//###0.040  
  hempty->SetMaximum(3);
  hempty->SetMinimum(0.);
  hempty->Draw();
  
  TLine *l = new TLine(0,1,70,1);
  l->SetLineStyle(2);
  
  legendRpA->Draw();
  gRpAstat->SetMarkerStyle(21);
  gRpAstat->SetLineColor(1);
  gRpAstat->SetMarkerColor(1);

  gRpAsystFONLL->Draw("2same");
  gRpAsyst->Draw("2esame");
  gRpAstat->Draw("psame");
  

  TBox *b = new TBox(-3.3,1-commonErrorN,-2.9,1+commonErrorP);
  b->SetLineColor(1);
  b->SetFillColor(kGray);
  b->Draw();
  
  
  TLegendEntry *ent_RpAstat=legendRpA->AddEntry(gRpAstat,"R^{FONLL}_{pA}","pf");
  ent_RpAstat->SetTextFont(42);
  ent_RpAstat->SetLineColor(2);
  ent_RpAstat->SetMarkerColor(2);
  
  TLegendEntry *ent_RpAsystData=legendRpA->AddEntry(b,"Syst. Lumi+BR","f");
  ent_RpAsystData->SetTextFont(42);
  ent_RpAsystData->SetLineColor(2);
  ent_RpAsystData->SetMarkerColor(2);
  
  TLegendEntry *ent_RpAsystFONLL=legendRpA->AddEntry(gRpAsystFONLL,"Syst. err. from FONLL pp ref.","f");
  ent_RpAsystFONLL->SetTextFont(42);
  ent_RpAsystFONLL->SetLineColor(5);
  ent_RpAsystFONLL->SetLineStyle(1);
  ent_RpAsystFONLL->SetMarkerColor(5);
  
  tlatex1->Draw();
  tlatexlumi->Draw();
  tlatex3->Draw();
  tlatex4->Draw();
  canvasRpA->SaveAs(Form("../Results%s_y/canvasRpA%s_FONLL.pdf",particle.Data(),particle.Data()));  
  canvasRpA->Clear();
   hempty->Draw();
  gRpAsystFONLL_CMS7->Draw("2same");
  gRpAsyst_CMS7->Draw("2esame");
  gRpAstat_CMS7->Draw("psame");
  canvasRpA->SaveAs(Form("../Results%s_y/canvasRpA%s_CMS7.pdf",particle.Data(),particle.Data()));  
   canvasRpA->Clear();
  hempty->Draw();
 gRpAsystFONLL_ATL7->Draw("2same");
  gRpAsyst_ATL7->Draw("2esame");
  gRpAstat_ATL7->Draw("psame");
  canvasRpA->SaveAs(Form("../Results%s_y/canvasRpA%s_ATL7.pdf",particle.Data(),particle.Data()));  


/*
  TCanvas *canvasRFB=new TCanvas("canvasRFB","canvasRFB",500,500);   
  canvasRFB->cd();
  canvasRFB->Range(-1.989924,-0.2917772,25.49622,2.212202);
  canvasRFB->SetFillColor(0);
  canvasRFB->SetBorderMode(0);
  canvasRFB->SetBorderSize(2);
  canvasRFB->SetLeftMargin(0.1451613);
  canvasRFB->SetRightMargin(0.05443548);
  canvasRFB->SetTopMargin(0.08474576);
  canvasRFB->SetBottomMargin(0.1165254);
  canvasRFB->SetFrameBorderMode(0);
  canvasRFB->SetFrameBorderMode(0);
 
  TLegend *legendRFB=new TLegend(0.2560484,0.6490486,0.5483871,0.7716702,"");
  legendRFB->SetBorderSize(0);
  legendRFB->SetLineColor(0);
  legendRFB->SetFillColor(0);
  legendRFB->SetFillStyle(1001);
  legendRFB->SetTextFont(42);
  legendRFB->SetTextSize(0.045);

  TH2F* hempty=new TH2F("hempty","",4,-0.1,2.,10.,0.,2.5);  
  hempty->GetXaxis()->SetTitle("|y_{CM}|");
  hempty->GetYaxis()->SetTitle("R_{FB}");
  hempty->GetXaxis()->SetTitleOffset(1.1);
  hempty->GetYaxis()->SetTitleOffset(1.3);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.04);  
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->SetMaximum(3);
  hempty->SetMinimum(0.);
  hempty->Draw();
  legendRFB->Draw();

  TGraphAsymmErrors *gRFBsyst = new TGraphAsymmErrors(nbinsRFB,xbinsRFB,yRFB,exlRFB,exlRFB,yRFBSystTotLow,yRFBSystTotHigh);
  gRFBsyst->SetTitle("gRFBsyst");
   gRFBsyst->SetName("gRFBsyst");
  gRFBsyst->SetMarkerColor(4);
  gRFBsyst->SetLineColor(4);
  gRFBsyst->SetLineWidth(2);   
  gRFBsyst->SetMarkerStyle(0);
  gRFBsyst->SetMarkerColor(1);


  


  TGraphAsymmErrors *gRFBstat = new TGraphAsymmErrors(nbinsRFB,xbinsRFB,yRFB,exlRFB,exlRFB,yRFBStat, yRFBStat);
  gRFBstat->SetTitle("gRFBstat");
  gRFBstat->SetName("gRFBstat");
  gRFBstat->SetMarkerColor(1);
  gRFBstat->SetLineColor(1);
  gRFBstat->SetLineWidth(2);   
  gRFBstat->SetMarkerStyle(21);
  gRFBstat->SetMarkerColor(1);
  gRFBstat->SetFillColor(0);
  gRFBstat->SetFillStyle(0);

  gRFBsyst->SetFillColor(0);
  gRFBsyst->SetFillStyle(0);
  
  gRFBsyst->Draw("2same");
  gRFBstat->Draw("psame");



  TLegendEntry *ent_RFB=legendRFB->AddEntry(gRFBstat,"R_{FB} stat.","pl");
  ent_RFB->SetTextFont(42);
  ent_RFB->SetLineColor(2);
  ent_RFB->SetMarkerColor(2);

  TLegendEntry *ent_RFB=legendRFB->AddEntry(gRFBsyst,"R_{FB} syst.","l");
  ent_RFB->SetTextFont(42);
  ent_RFB->SetLineColor(4);
  // ent_RFB->SetMarkerColor(2);

  TLatex * tlatex1=new TLatex(0.1612903,0.8625793,"CMS                 pPb #sqrt{s_{NN}}= 5.02 TeV");
  tlatex1->SetNDC();
  tlatex1->SetTextColor(1);
  tlatex1->SetTextFont(42);
  tlatex1->SetTextSize(0.045);
  tlatex1->Draw();

    
  TLatex * tlatexlumi=new TLatex(0.671371,0.7801268,"L = 34.8 nb^{-1}");
  tlatexlumi->SetNDC();
  tlatexlumi->SetTextColor(1);
  tlatexlumi->SetTextFont(42);
  tlatexlumi->SetTextSize(0.045);
  tlatexlumi->Draw();

  
  //double xpos1=0.6528226;
  // double ypos1=0.6849894;
  TLatex * tlatex3=new TLatex(xpos,ypos,mypar.Data());
   tlatex3->SetNDC();
  tlatex3->SetTextColor(1);
  tlatex3->SetTextFont(42);
  tlatex3->SetTextSize(0.06);
  tlatex3->Draw();
  canvasRFB->SaveAs(Form("../Results%s_y/CanvasRDFB.pdf",particle.Data(),particle.Data()));

  */


//  l->Draw();  
//  canvasRpA->SaveAs(Form("../Results%sY/canvasRpA%s.pdf",particle.Data(),particle.Data()));  
//  canvasRpA->SaveAs(Form("../Results%s_y/canvasRpA%s.pdf",particle.Data(),particle.Data()));  
  
//  TFile *fout=new TFile(Form("../Results%sY/fileRpA%s.root",particle.Data(),particle.Data()),"recreate");  
  TFile *fout=new TFile(Form("../Results%s_y/fileRpA%s.root",particle.Data(),particle.Data()),"recreate");  
  fout->cd();
  gSigmasyst->SetName("gSigmasyst");
  gSigmasyst->Write();
  gaeBplusReference->Write();
  hSigmapPbStat->Write();






}
