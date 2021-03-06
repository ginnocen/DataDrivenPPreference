#include "TH1F.h"
#include <cmath>
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>

#define BIN_NUM 460 //pPb_pt:220,pPb_y:40,pp_pt:222,pp_y:45
#define REBIN 8     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
#define REBINp 9    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5
#define HMIN 5      //pPb_pt:5,pPb_y:-2,pp_pt:9,pp_y:0
#define HMAX 120     //pPb_pt:55,pPb_y:2,pp_pt:120,pp_y:2.25

int Bplusdsigmadpt_pp_7TeV_AtlasBin()
{
  TString infile,outfile;
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  infile="../FONLLInputs/fo_pp_pt_rap225_7TeV.dat";
  outfile="Rootf/outputBplus_pp_pt_rap225_7TeV.root";
  
  ifstream getdata(infile.Data());

  if(!getdata.is_open())
    {
      cout<<"Opening the file fails"<<endl;
    }

  float central[BIN_NUM];
  float min_all[BIN_NUM],max_all[BIN_NUM],min_sc[BIN_NUM],max_sc[BIN_NUM],min_mass[BIN_NUM],max_mass[BIN_NUM],min_pdf[BIN_NUM],max_pdf[BIN_NUM];
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
  TH1F* hratio = new TH1F("hratio","",BIN_NUM,HMIN,HMAX);


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
    }

  double data[REBIN] = {103.4,36.03,15.33,6.056,1.814,0.3477,0.06244,0.006099};
  double staterror[REBIN] = {4.,0.8,0.3,0.1,0.03,0.008,0.003,0.0006};
  double syserror[REBIN] = {8.,2.3,1.0,0.4,0.12,0.028,0.005,.0007};

  //Rebin Edge
  double rebin[REBINp] = {9.,13.,16.,20.,25.,35.,50.,70.,120.};

  TH1F* hpt_rebin = (TH1F*)hpt->Rebin(REBIN,"hpt_rebin",rebin);
  TH1F* hminall_rebin = (TH1F*)hminsc->Rebin(REBIN,"hminall_rebin",rebin);
  TH1F* hmaxall_rebin = (TH1F*)hmaxsc->Rebin(REBIN,"hmaxall_rebin",rebin);
  TH1F* hminsc_rebin = (TH1F*)hminsc->Rebin(REBIN,"hminsc_rebin",rebin);
  TH1F* hmaxsc_rebin = (TH1F*)hmaxsc->Rebin(REBIN,"hmaxsc_rebin",rebin);
  TH1F* hminmass_rebin = (TH1F*)hminmass->Rebin(REBIN,"hminmass_rebin",rebin);
  TH1F* hmaxmass_rebin = (TH1F*)hmaxmass->Rebin(REBIN,"hmaxmass_rebin",rebin);
  TH1F* hminpdf_rebin = (TH1F*)hminpdf->Rebin(REBIN,"hminpdf_rebin",rebin);
  TH1F* hmaxpdf_rebin = (TH1F*)hmaxpdf->Rebin(REBIN,"hmaxpdf_rebin",rebin);

  TH1F* hratio_rebin = (TH1F*)hratio->Rebin(REBIN,"hratio_rebin",rebin);


  //bin middle
  double apt[REBIN] = {11.,14.5,18.,22.5,30.,42.5,60.,95.};//pPb_pt
  //bin half width
  double aptl[REBIN] = {2.,1.5,2.,2.5,5.,7.5,10.,25.};//pPb_pt
  double asigma[REBIN],aminall[REBIN],amaxall[REBIN],aminsc[REBIN],amaxsc[REBIN],aminmass[REBIN],amaxmass[REBIN],aminpdf[REBIN],amaxpdf[REBIN],aerrorl[REBIN],aerrorh[REBIN];

  //number of every rebined bin
  double bin_num[REBIN] = {16.,12,16.,20.,40.,60.,80.,200.};//pPb_pt
  
  int j;
  double norm=1.;
  double inte=1.;

  for(j=0;j<REBIN;j++)
    {

      tem = hpt_rebin->GetBinContent(j+1);
      asigma[j] = tem*norm*inte/bin_num[j];

      tem = hminall_rebin->GetBinContent(j+1);
      aminall[j] = tem*norm*inte/bin_num[j];

      tem = hmaxsc_rebin->GetBinContent(j+1);
      amaxall[j] = tem*norm*inte/bin_num[j];

      tem = hminsc_rebin->GetBinContent(j+1);
      aminsc[j] = tem*norm*inte/bin_num[j];

      tem = hmaxsc_rebin->GetBinContent(j+1);
      amaxsc[j] = tem*norm*inte/bin_num[j];

      tem = hminmass_rebin->GetBinContent(j+1);
      aminmass[j] = tem*norm*inte/bin_num[j];

      tem = hmaxmass_rebin->GetBinContent(j+1);
      amaxmass[j] = tem*norm*inte/bin_num[j];

      tem = hminpdf_rebin->GetBinContent(j+1);
      aminpdf[j] = tem*norm*inte/bin_num[j];

      tem = hmaxpdf_rebin->GetBinContent(j+1);
      amaxpdf[j] = tem*norm*inte/bin_num[j];

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
  
  TLatex * tlatex=new TLatex(0.18,0.85,"pp collisions at 7TeV from FONLL, |y_{LAB}|<2.25");
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
  cr->SaveAs("Plots/cBmesonPredFONLLBplusBinning_pp_pt_rap225_7TeV_AtlasBin.pdf");
  TGraphAsymmErrors* gaeSigmaDecay=(TGraphAsymmErrors*)gae->Clone();
  gaeSigmaDecay->SetName("gaeSigmaDecay");
  //double BRchain=6.09604e-5;
  double BRchain=1.;
  double Fraction=0.401;
  double BR=0.001016*0.0593;

  double norm=1.;
  double BRFraction=BRchain*Fraction;

  double syserry[REBIN],syserrey[REBIN];
  
  for (int i=0;i<gaeSigmaDecay->GetN();i++){
    gaeSigmaDecay->GetY()[i] *= BRFraction*norm/1000000;
    gaeSigmaDecay->SetPointEYhigh(i,gaeSigmaDecay->GetErrorYhigh(i)*BRFraction*norm/1000000);
    gaeSigmaDecay->SetPointEYlow(i,gaeSigmaDecay->GetErrorYlow(i)*BRFraction*norm/1000000);
    hratio_rebin->SetBinContent(i+1,data[i]/(gaeSigmaDecay->GetY()[i]*1000000*BR));
    hratio_rebin->SetBinError(i+1,staterror[i]/(gaeSigmaDecay->GetY()[i]*1000000*BR));
    syserry[i] = data[i]/(gaeSigmaDecay->GetY()[i]*1000000*BR);
    syserrey[i] = syserror[i]/(gaeSigmaDecay->GetY()[i]*1000000*BR);
  }

  gaeSigmaDecay->SetFillColor(2);
  gaeSigmaDecay->SetFillStyle(3001); 
  gaeSigmaDecay->SetTitle(";p_{T}(GeV/c);d#sigma/dp_{T} (B^{+}) #times A (GeV^{-1}c)");
   
  TH2F* hempty=new TH2F("hempty","",10,0,120.,10.,1.e-6,10);  
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
  hempty->GetYaxis()->SetTitle("d#sigma/dp_{T}(B^{+}) #times A (#mub c/GeV)");

  TCanvas*canvas=new TCanvas("canvas","canvas",600,500);
  canvas->SetLogy();
  hempty->Draw();
  gaeSigmaDecay->Draw("psame");
  
  //TLatex * tlatex=new TLatex(0.2,0.85,"B^{+}=40.1%, |y|<1.93, BR unc not shown");
  TLatex * tlatex=new TLatex(0.2,0.85,"B^{+},|y_{LAB}|<2.25, BR unc not shown");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.05);
  tlatex->Draw();
  
  gae->SetName("gaeBplus");
  gaeSigmaDecay->SetName("gaeSigmaDecayBplus");
  canvas->SaveAs("Plots/canvasBplus_pp_pt_rap225_7TeV_AtlasBin.pdf");

  TCanvas*cratio=new TCanvas("cratio","cratio",500,500);
  hratio_rebin->SetMaximum(2.);
  hratio_rebin->SetMinimum(0.5);
  hratio_rebin->SetXTitle("p_{T}(GeV/c)");
  hratio_rebin->SetYTitle("#sigma / #sigma(FONLL)");
  hratio_rebin->SetTitleOffset(1.2,"Y");
  hratio_rebin->SetMarkerStyle(8);
  hratio_rebin->SetStats(0);
  hratio_rebin->Draw("lep");

  TGraphErrors* gsyserror = new TGraphErrors(REBIN,apt,syserry,aptl,syserrey);  
  gsyserror->SetMarkerColor(1);
  gsyserror->SetLineColor(1);
  gsyserror->SetLineWidth(1);   
  gsyserror->SetMarkerStyle(21);
  gsyserror->SetMarkerColor(1);
  gsyserror->SetFillColor(0);
  gsyserror->SetFillStyle(0);
  gsyserror->Draw("2same");

  TLegend *leg = new TLegend(0.5,0.75,0.9,0.9);
  leg->AddEntry((TObject*)0,"ATLAS pp 7TeV","");
  leg->AddEntry((TObject*)0,"|y_{LAB}|<2.25","");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("same");

  TLegend *leg2 = new TLegend(0.55,0.60,0.9,0.75);
  leg2->AddEntry(gsyserror,"Sys uncert","f");
  leg2->AddEntry(hratio_rebin,"Stat uncert","ple");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw("same");

  TLine* lin0=new TLine(9,1,120,1);
  lin0->SetLineStyle(2);
  lin0->SetLineColor(1);
  lin0->SetLineWidth(3);
  lin0->Draw("same");

  cratio->SaveAs("Plots/cratio_pp_pt_rap225_7TeV_AtlasBin.pdf");

  TFile*foutput=new TFile(outfile.Data(),"recreate");
  foutput->cd();
  gae->Write();
  gaeSigmaDecay->Write();
  
}
