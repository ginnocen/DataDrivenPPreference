#include "TH1F.h"
#include <cmath>
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>

#define BIN_NUM 460 //pPb_pt:220,pPb_y:40,pp_pt:222,pp_y:45
#define REBIN 5     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
#define REBINp 6    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5
#define HMIN 5      //pPb_pt:5,pPb_y:-2,pp_pt:9,pp_y:0
#define HMAX 120     //pPb_pt:55,pPb_y:2,pp_pt:120,pp_y:2.25

int Bplusdsigmadpt_pp_5p02TeV_PPbBin()
{
  TString infile,outfile;
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  infile="../FONLLInputs/fo_Curve_pp_pt_rap24_5p02TeV.dat";
  outfile="Rootf/outputBplus_pp_pt_rap24_5p02TeV_PPbBin.root";

  ifstream getdata(infile.Data());

  if(!getdata.is_open())
    {
      cout<<"Opening the file fails"<<endl;
    }

  float central[BIN_NUM];
  float min_all[BIN_NUM],max_all[BIN_NUM],min_sc[BIN_NUM],max_sc[BIN_NUM],min_mass[BIN_NUM],max_mass[BIN_NUM];
  float min_pdf[BIN_NUM],max_pdf[BIN_NUM];
  float fr_05_05[BIN_NUM],fr_20_20[BIN_NUM],fr_20_10[BIN_NUM],fr_10_20[BIN_NUM],fr_10_05[BIN_NUM],fr_05_10[BIN_NUM];

  int i;
  float tem;
  for(i=0;i<BIN_NUM;i++)
    {
      getdata>>tem;
      getdata>>central[i];
      getdata>>min_all[i];
      getdata>>max_all[i];
      getdata>>min_sc[i];
      getdata>>max_sc[i];
      getdata>>min_mass[i];
      getdata>>max_mass[i];
	  getdata>>min_pdf[i];
	  getdata>>max_pdf[i];
	  getdata>>fr_05_05[i];
	  getdata>>fr_20_20[i];
	  getdata>>fr_20_10[i];
	  getdata>>fr_10_20[i];
	  getdata>>fr_10_05[i];
	  getdata>>fr_05_10[i];

    }
  
  TH1F* hpt = new TH1F("hpt","",BIN_NUM,HMIN,HMAX);
  TH1F* hminall = new TH1F("hminall","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxall = new TH1F("hmaxall","",BIN_NUM,HMIN,HMAX);
  TH1F* hminsc = new TH1F("hminsc","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxsc = new TH1F("hmaxsc","",BIN_NUM,HMIN,HMAX);
  TH1F* hminmass = new TH1F("hminmass","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxmass = new TH1F("hmaxmass","",BIN_NUM,HMIN,HMAX);
  TH1F* hminpdf = new TH1F("hminpdf","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxpdf = new TH1F("hmaxpdf","",BIN_NUM,HMIN,HMAX);
  TH1F* hfr_05_05 = new TH1F("hfr_05_05","",BIN_NUM,HMIN,HMAX);
  TH1F* hfr_20_20 = new TH1F("hfr_20_20","",BIN_NUM,HMIN,HMAX);
  TH1F* hfr_20_10 = new TH1F("hfr_20_10","",BIN_NUM,HMIN,HMAX);
  TH1F* hfr_10_20 = new TH1F("hfr_10_20","",BIN_NUM,HMIN,HMAX);
  TH1F* hfr_10_05 = new TH1F("hfr_10_05","",BIN_NUM,HMIN,HMAX);
  TH1F* hfr_05_10 = new TH1F("hfr_05_10","",BIN_NUM,HMIN,HMAX);


  for(i=0;i<BIN_NUM;i++)
    {
      hpt->SetBinContent(i+1,central[i]);
      hminall->SetBinContent(i+1,min_all[i]);
      hmaxall->SetBinContent(i+1,max_all[i]);
      hminsc->SetBinContent(i+1,min_sc[i]);
      hmaxsc->SetBinContent(i+1,max_sc[i]);
      hminmass->SetBinContent(i+1,min_mass[i]);
      hmaxmass->SetBinContent(i+1,max_mass[i]);
      hminpdf->SetBinContent(i+1,min_pdf[i]);
      hmaxpdf->SetBinContent(i+1,max_pdf[i]);
	  hminpdf->SetBinContent(i+1,min_pdf[i]);
	  hmaxpdf->SetBinContent(i+1,max_pdf[i]);
	  hfr_05_05->SetBinContent(i+1,fr_05_05[i]);
	  hfr_20_20->SetBinContent(i+1,fr_20_20[i]);
	  hfr_20_10->SetBinContent(i+1,fr_20_10[i]);
	  hfr_10_20->SetBinContent(i+1,fr_10_20[i]);
	  hfr_05_10->SetBinContent(i+1,fr_05_10[i]);
	  hfr_10_05->SetBinContent(i+1,fr_10_05[i]);

    }
  //Rebin Edge
  double rebin[REBINp] = {10,15,20,25,30,60};

  TH1F* hpt_rebin = (TH1F*)hpt->Rebin(REBIN,"hpt_rebin",rebin);
  TH1F* hminall_rebin = (TH1F*)hminsc->Rebin(REBIN,"hminall_rebin",rebin);
  TH1F* hmaxall_rebin = (TH1F*)hmaxsc->Rebin(REBIN,"hmaxall_rebin",rebin);
  TH1F* hminsc_rebin = (TH1F*)hminsc->Rebin(REBIN,"hminsc_rebin",rebin);
  TH1F* hmaxsc_rebin = (TH1F*)hmaxsc->Rebin(REBIN,"hmaxsc_rebin",rebin);
  TH1F* hminmass_rebin = (TH1F*)hminmass->Rebin(REBIN,"hminmass_rebin",rebin);
  TH1F* hmaxmass_rebin = (TH1F*)hmaxmass->Rebin(REBIN,"hmaxmass_rebin",rebin);
  TH1F* hminpdf_rebin = (TH1F*)hminpdf->Rebin(REBIN,"hminpdf_rebin",rebin);
  TH1F* hmaxpdf_rebin = (TH1F*)hmaxpdf->Rebin(REBIN,"hmaxpdf_rebin",rebin);

  //bin middle
  double apt[REBIN] = {12.5,17.5,22.5,27.5,45};//pPb_pt
  //bin half width
  double aptl[REBIN] = {2.5,2.5,2.5,2.5,15};//pPb_pt
  double asigma[REBIN],aminall[REBIN],amaxall[REBIN],aminsc[REBIN],amaxsc[REBIN],aminmass[REBIN],amaxmass[REBIN],aminpdf[REBIN],amaxpdf[REBIN],aerrorl[REBIN],aerrorh[REBIN];

  //number of every rebined bin
  double bin_num[REBIN] = {20,20,20,20,120};//pPb_pt
  int j;
  double norm=1.;
  
  for(j=0;j<REBIN;j++)
    {

      tem = hpt_rebin->GetBinContent(j+1);
      asigma[j] = tem*norm/bin_num[j];

      tem = hminall_rebin->GetBinContent(j+1);
      aminall[j] = tem*norm/bin_num[j];

      tem = hmaxsc_rebin->GetBinContent(j+1);
      amaxall[j] = tem*norm/bin_num[j];

      tem = hminsc_rebin->GetBinContent(j+1);
      aminsc[j] = tem*norm/bin_num[j];

      tem = hmaxsc_rebin->GetBinContent(j+1);
      amaxsc[j] = tem*norm/bin_num[j];

      tem = hminmass_rebin->GetBinContent(j+1);
      aminmass[j] = tem*norm/bin_num[j];

      tem = hmaxmass_rebin->GetBinContent(j+1);
      amaxmass[j] = tem*norm/bin_num[j];

      tem = hminpdf_rebin->GetBinContent(j+1);
      aminpdf[j] = tem*norm/bin_num[j];

      tem = hmaxpdf_rebin->GetBinContent(j+1);
      amaxpdf[j] = tem*norm/bin_num[j];

      aerrorl[j] = asigma[j]-aminall[j];//all,sc,mass,pdf
      aerrorh[j] = amaxall[j]-asigma[j];//all,sc,mass,pdf
    }

  cout<<"------- pp_7------"<<endl;
  cout<<endl;
 
  TGraphAsymmErrors* gae = new TGraphAsymmErrors(REBIN, apt, asigma, aptl, aptl, aerrorl, aerrorh);
  gae->SetTitle(";p_{T}(GeV/c);d#sigma (B admix) /dp_{T}(pb c/GeV)");
  gae->SetFillColor(2);
  gae->SetFillStyle(3001);

  TCanvas* cr = new TCanvas("cr","cr",600,500);
  cr->SetLogy();
  TH2F* hempty=new TH2F("hempty","",10,5,60.,10.,10,100000000);  
  hempty->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hempty->GetYaxis()->SetTitle("d#sigma(B admix)/dp_{T}(pb/GeV)");
  hempty->GetXaxis()->SetTitleOffset(1.);
  hempty->GetYaxis()->SetTitleOffset(.9);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.04);  
  hempty->Draw();
  hminall->SetLineColor(2);
  hmaxall->SetLineColor(2);
  hpt->SetLineColor(2);
  hminall->Draw("same");
  hmaxall->Draw("same");
  hpt->Draw("same");
  gae->SetLineWidth(3);
  gae->Draw("psame");
  
  TLatex * tlatex=new TLatex(0.18,0.85,"pp collisions at 5.02 from FONLL, |y_{LAB}|<2.4");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.04);
  tlatex->Draw();
  TLatex * tlatex=new TLatex(0.18,0.80,"Total syst uncertainties shown");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.04);
  tlatex->Draw();
  cr->SaveAs("Plots/cBmesonPredFONLLBplusBinning_pp_pt_rap24_5p02TeV_PPbBin.pdf");
  TGraphAsymmErrors* gaeSigmaDecay=(TGraphAsymmErrors*)gae->Clone();
  gaeSigmaDecay->SetName("gaeSigmaDecay");
  //double BRchain=6.09604e-5;
  double BRchain=1.;
  double Fraction=0.401;
  
  double norm=1.;
  double BRFraction=BRchain*Fraction;
  
  for (int i=0;i<gaeSigmaDecay->GetN();i++){
    gaeSigmaDecay->GetY()[i] *= BRFraction*norm;
    gaeSigmaDecay->SetPointEYhigh(i,gaeSigmaDecay->GetErrorYhigh(i)*BRFraction*norm);
    gaeSigmaDecay->SetPointEYlow(i,gaeSigmaDecay->GetErrorYlow(i)*BRFraction*norm); 
  }
  
  gaeSigmaDecay->SetFillColor(2);
  gaeSigmaDecay->SetFillStyle(3001); 
  gaeSigmaDecay->SetTitle(";p_{T}(GeV/c);d#sigma/dp_{T} (B^{+}) #times A (GeV^{-1}c)");
   
  TH2F* hempty=new TH2F("hempty","",10,0,70.,10.,100.,50000000);  
  hempty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hempty->GetXaxis()->SetTitleOffset(1.);
  hempty->GetYaxis()->SetTitleOffset(.9);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.035);  
  hempty->GetYaxis()->SetTitle("d#sigma/dp_{T}(B^{+}) #times A (pb c/GeV)");

  TCanvas*canvas=new TCanvas("canvas","canvas",600,500);
  canvas->SetLogy();
  hempty->Draw();
  gaeSigmaDecay->Draw("psame");
  
  //TLatex * tlatex=new TLatex(0.2,0.85,"B^{+}=40.1%, |y|<1.93, BR unc not shown");
  TLatex * tlatex=new TLatex(0.2,0.85,"B^{+},|y_{LAB}|<2.4, BR unc not shown");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.05);
  tlatex->Draw();
  
  gae->SetName("gaeBplus");
  gaeSigmaDecay->SetName("gaeSigmaDecayBplus");
  canvas->SaveAs("Plots/canvasBplus_pp_pt_rap24_5p02TeV_PPbBin.pdf");

  TFile*foutput=new TFile(outfile.Data(),"recreate");
  foutput->cd();
  gae->Write();
  gaeSigmaDecay->Write();
  hpt->Write();
  hfr_05_05->Write();
  hfr_20_20->Write();
  hfr_20_10->Write();
  hfr_10_20->Write();
  hfr_10_05->Write();
  hfr_05_10->Write();
  hminall->Write();
  hmaxall->Write();
  hminsc->Write();
  hmaxsc->Write();
  hminmass->Write();
  hmaxmass->Write();
  hminpdf->Write();
  hmaxpdf->Write();
  
	TCanvas* c1 = new TCanvas("c1","",500,500);
	c1->SetLogy(1);

	TH1F* hbase = new TH1F("hbase","",13,5.0,65.0);
	hbase->GetXaxis()->SetTitle("B^{+} p_{T}");
	hbase->GetYaxis()->SetRangeUser(100.0,1000000000.0);
	hbase->GetXaxis()->CenterTitle();
	hbase->GetYaxis()->CenterTitle();
	hbase->GetXaxis()->SetTitleOffset(1.0);
	hbase->GetYaxis()->SetTitleOffset(1.4);
	hbase->SetLineColor(0);
	hbase->GetXaxis()->SetTitleSize(0.045);
	hbase->GetYaxis()->SetTitleSize(0.045);
	hbase->GetXaxis()->SetTitleFont(42);
	hbase->GetYaxis()->SetTitleFont(42);
	hbase->GetXaxis()->SetLabelFont(42);
	hbase->GetYaxis()->SetLabelFont(42);
	hbase->GetXaxis()->SetLabelSize(0.04);
	hbase->GetYaxis()->SetLabelSize(0.04);
	hbase->Draw("");

	hpt->SetLineColor(kRed);
	hpt->SetLineWidth(3);
	hfr_05_05->SetLineColor(kOrange+7);
	hfr_20_20->SetLineColor(kYellow-6);
	hfr_20_10->SetLineColor(kTeal+3);
	hfr_10_20->SetLineColor(kAzure+1);
	hfr_05_10->SetLineColor(kBlue+4);
	hfr_10_05->SetLineColor(kViolet+4);
	hpt->Draw("same");
	hfr_05_05->Draw("same");
	hfr_20_20->Draw("same");
	hfr_20_10->Draw("same");
	hfr_10_20->Draw("same");
	hfr_05_10->Draw("same");
	hfr_10_05->Draw("same");

	TLegend* leg = new TLegend(0.68,0.51,0.86,0.90);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.040);
	leg->AddEntry(hbase,"#mu_{F}/#mu_{0} , #mu_{R}/#mu_{0}","l");
	leg->AddEntry(hpt,"  1.0   ,   1.0","l");
	leg->AddEntry(hfr_05_05,"  0.5   ,   0.5","l");
	leg->AddEntry(hfr_20_20,"  2.0   ,   2.0","l");
	leg->AddEntry(hfr_20_10,"  2.0   ,   1.0","l");
	leg->AddEntry(hfr_10_20,"  1.0   ,   2.0","l");
	leg->AddEntry(hfr_10_05,"  1.0   ,   0.5","l");
	leg->AddEntry(hfr_05_10,"  0.5   ,   1.0","l");
	leg->Draw("");
	c1->SaveAs("../ResultsBplus/CompFONLLdsigmadpt_5p02TeV.pdf");


}
