//Made by Gian Michele, modified by Hyunchul Kim

#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"

#include <cmath>
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TF1.h"

#include <fstream>
#include <iostream>

//#define BIN_NUM 220 //pPb_pt:220,pPb_y:40,pp_pt:222,pp_y:45
//#define REBIN 55     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
//#define REBINp 56    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5
//#define HMIN 5      //pPb_pt:5,pPb_y:-2,pp_pt:9,pp_y:0
//#define HMAX 60     //pPb_pt:55,pPb_y:2,pp_pt:120,pp_y:2.25

//#define BIN_NUM 480 //pPb_pt:220,pPb_y:40,pp_pt:222,pp_y:45
#define BIN_NUM 600 //pPb_pt:220,pPb_y:40,pp_pt:222,pp_y:45

//#define REBIN 48     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
//#define REBINp 49    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5

//#define REBIN 60     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
//#define REBINp 61    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5

//#define HMINp -2.865      //pPb_pt:5,pPb_y:-2,pp_pt:9,pp_y:0
//#define HMAXp 1.935     //pPb_pt:55,pPb_y:2,pp_pt:120,pp_y:2.25

#define HMINp -3.0      //pPb_pt:5,pPb_y:-2,pp_pt:9,pp_y:0
#define HMAXp 3.0     //pPb_pt:55,pPb_y:2,pp_pt:120,pp_y:2.25

#define REBIN 5     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
#define REBINp 6    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5


void BplusdsigmadyUnbinned_all_v5(int option, bool isBinned)
{

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  std::string ename;
  ifstream getdata;
  if (option==5) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_10_60_all_5TeV.dat");ename="5TeV_rap30_pt1060";}
  if (option==7) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_10_60_all_7TeV.dat");ename="7TeV_rap30_pt1060";}
  if (option==71) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_5_120_all_7TeV.dat");ename="7TeV_rap30_pt5120";}
  if (option==72) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_9_120_all_7TeV.dat");ename="7TeV_rap30_pt9120";}

  //TFile*foutput=new TFile(Form("../outputBplusyUnbinned_%s.root",ename.c_str()),"recreate");
  TFile*foutput=new TFile(Form("../outputBplusyUnbinned_%s_Binned.root",ename.c_str()),"recreate");

  if(!getdata.is_open())
  {
    cout<<"Opening the file fails"<<endl;
  }

  float central[BIN_NUM],y[BIN_NUM];
  float min_all[BIN_NUM],max_all[BIN_NUM],min_sc[BIN_NUM],max_sc[BIN_NUM],min_mass[BIN_NUM],max_mass[BIN_NUM],min_pdf[BIN_NUM],max_pdf[BIN_NUM];
  float fr_05_05[BIN_NUM],fr_20_20[BIN_NUM],fr_20_10[BIN_NUM],fr_10_20[BIN_NUM],fr_10_05[BIN_NUM],fr_05_10[BIN_NUM];
  float tem;
  int i;
  for(i=0;i<BIN_NUM;i++)
  {
    getdata>>tem;
    // all the value is 
    if (option==5) y[i]=tem-0.465; else y[i]=tem;
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
  double HMIN,HMAX;

  if (option==5) {HMIN=HMINp-0.465;HMAX=HMAXp-0.465;}
  else {HMIN=HMINp;HMAX=HMAXp;}

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
    hfr_05_05->SetBinContent(i+1,fr_05_05[i]);
    hfr_20_20->SetBinContent(i+1,fr_20_20[i]);
    hfr_20_10->SetBinContent(i+1,fr_20_10[i]);
    hfr_10_20->SetBinContent(i+1,fr_10_20[i]);
    hfr_05_10->SetBinContent(i+1,fr_05_10[i]);
    hfr_10_05->SetBinContent(i+1,fr_10_05[i]);
  }

  TCanvas* c1 = new TCanvas("c1","",500,500);
  //c1->SetLogy(1);
  hpt->SetLineColor(kRed);
  hpt->GetYaxis()->SetRangeUser(0,7000000);
  hfr_05_05->SetLineColor(kOrange+7);
  hfr_20_20->SetLineColor(kYellow-6);
  hfr_20_10->SetLineColor(kTeal+3);
  hfr_10_20->SetLineColor(kAzure+1);
  hfr_05_10->SetLineColor(kBlue+4);
  hfr_10_05->SetLineColor(kViolet+4);
  hpt->Draw("");
  hfr_05_05->Draw("same");
  hfr_20_20->Draw("same");
  hfr_20_10->Draw("same");
  hfr_10_20->Draw("same");
  hfr_05_10->Draw("same");
  hfr_10_05->Draw("same");
  c1->SaveAs(Form("CompFONLL_Unbinned_%s.pdf",ename.c_str()));

  // REBIN here

  double rebiny[REBIN+1] = {-2.4,-1.470,-0.535,0.465,1.465,2.4};//pPb_lab

  Double_t rebin[REBINp];
  Float_t apt[REBIN];
  Float_t aptl[REBIN];

  Float_t asigma[REBIN],aminall[REBIN],amaxall[REBIN],aminsc[REBIN],amaxsc[REBIN],aminmass[REBIN],amaxmass[REBIN],aminpdf[REBIN],amaxpdf[REBIN],aerrorl[REBIN],aerrorh[REBIN];
  Float_t bin_num[REBIN];

  if (isBinned) {
    for (i=0;i<REBIN;i++) {
      if (option==5) {
	apt[i]=(rebiny[i+1]+rebiny[i])/2;
      }
      else {
	apt[i]=((rebiny[i+1]-0.465)+(rebiny[i]-0.465))/2;
      }
      aptl[i]=(rebiny[i+1]-rebiny[i])/2;
      rebin[i]=rebiny[i];
    }
  }
  else {
    for (i=0;i<48;i++) {
      if (option==5) {apt[i]=-2.4+i*0.1;}
      else {apt[i]=-2.865+i*0.1;}
      aptl[i]=0.005;
      bin_num[i]=10;
    }
  }
  if (option==5) rebin[REBIN] = 2.4;
  else rebin[REBIN] = 1.935;

  TH1F* hpt_rebin = (TH1F*)hpt->Rebin(REBIN,"hpt_rebin",rebin);
  TH1F* hminall_rebin = (TH1F*)hminsc->Rebin(REBIN,"hminall_rebin",rebin);
  TH1F* hmaxall_rebin = (TH1F*)hmaxsc->Rebin(REBIN,"hmaxall_rebin",rebin);
  TH1F* hminsc_rebin = (TH1F*)hminsc->Rebin(REBIN,"hminsc_rebin",rebin);
  TH1F* hmaxsc_rebin = (TH1F*)hmaxsc->Rebin(REBIN,"hmaxsc_rebin",rebin);
  TH1F* hminmass_rebin = (TH1F*)hminmass->Rebin(REBIN,"hminmass_rebin",rebin);
  TH1F* hmaxmass_rebin = (TH1F*)hmaxmass->Rebin(REBIN,"hmaxmass_rebin",rebin);
  TH1F* hminpdf_rebin = (TH1F*)hminpdf->Rebin(REBIN,"hminpdf_rebin",rebin);
  TH1F* hmaxpdf_rebin = (TH1F*)hmaxpdf->Rebin(REBIN,"hmaxpdf_rebin",rebin);

  int j;
  double norm=1.;

  for(j=0;j<REBIN;j++)
  {
    double binwidth=(rebiny[j+1]-rebiny[j])/0.01;
    tem = hpt_rebin->GetBinContent(j+1);
    asigma[j] = tem*norm/binwidth;

    tem = hminall_rebin->GetBinContent(j+1);
    aminall[j] = tem*norm/binwidth;

    tem = hmaxall_rebin->GetBinContent(j+1);
    amaxall[j] = tem*norm/binwidth;

    tem = hminsc_rebin->GetBinContent(j+1);
    aminsc[j] = tem*norm/binwidth;

    tem = hmaxsc_rebin->GetBinContent(j+1);
    amaxsc[j] = tem*norm/binwidth;

    tem = hminmass_rebin->GetBinContent(j+1);
    aminmass[j] = tem*norm/binwidth;

    tem = hmaxmass_rebin->GetBinContent(j+1);
    amaxmass[j] = tem*norm/binwidth;

    tem = hminpdf_rebin->GetBinContent(j+1);
    aminpdf[j] = tem*norm/binwidth;

    tem = hmaxpdf_rebin->GetBinContent(j+1);
    amaxpdf[j] = tem*norm/binwidth;

    aerrorl[j] = asigma[j]-aminall[j];//all,sc,mass,pdf
    aerrorh[j] = amaxall[j]-asigma[j];//all,sc,mass,pdf
  }

  cout<<"------- pPb_5.02------"<<endl;
  //cout<<"------- pp 7------"<<endl;

  cout<<endl;

  TGraphAsymmErrors* gae = new TGraphAsymmErrors(REBIN, apt, asigma, aptl, aptl, aerrorl, aerrorh);

  gae->SetTitle(";y_{CM};d#sigma (B admix) / dy (pb)");
  gae->SetFillColor(2);
  gae->SetFillStyle(3001);

  TCanvas* cr = new TCanvas("cr","cr",600,500);
  //cr->SetLogy();
  //TH2F* hempty=new TH2F("hempty","",10,5,60.,10.,10,100000000);  

  TH2F* hempty;
  /*
     if (option==7 || option==71 || option==72) hempty=new TH2F("hempty","",50,-2.5,2.5,300.,0,30000000);
  //hempty->SetAxisRange(0.,22000000.,"Y");
  else hempty=new TH2F("hempty","",50,-2.965,2.035,70.,0,7000000);
   */
  hempty=new TH2F("hempty","",50,-2.965,2.035,70.,0,7000000);

  hempty->GetXaxis()->SetTitle("y_{CM}");
  hempty->GetYaxis()->SetTitle("d#sigma(B admix)/dy (pb)");
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
  gae->SetLineWidth(1);

  TF1* fitft = new TF1("fitft","[0]+[1]*pow((x-0.465),2)",-3.0,3.0);

  std::cout << "##################################################################" << std::endl;
  gae->Fit(fitft,"","",-2.865,1.935);
  std::cout << "##################################################################" << std::endl;

  gae->Draw("psame");

  
  TLatex * tlatex;
  std::string tlatexremt;
  int tlatexptmin, tlatexptmax;

  if (option==5) {tlatexremt="5TeV";tlatexptmin=10;tlatexptmax=60;}
  if (option==7) {tlatexremt="7TeV";tlatexptmin=10;tlatexptmax=60;}
  if (option==71) {tlatexremt="7TeV";tlatexptmin=5;tlatexptmax=120;}
  if (option==72) {tlatexremt="7TeV";tlatexptmin=9;tlatexptmax=120;}

  tlatex=new TLatex(0.18,0.85,Form("pp collisions at %s from FONLL, -2.865<y_{CM}<1.935, %i<p_{T}<%i",tlatexrem.c_str(),tlatexptmin,tlatexptmax));

  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.04);
  tlatex->Draw();
  tlatex=new TLatex(0.18,0.80,"Total syst uncertainties shown");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.04);
  tlatex->Draw();
  //cr->SaveAs("Plots/cBmesonPredFONLLBplusBinning_7TeV.eps");
  cr->SaveAs("Plots/cBmesonPredFONLLBplusyBinning.eps");
  cr->SaveAs(Form("Plots/cBmesonPredFONLLBplusyBinning_%s.pdf",ename.c_str()));

  TGraphAsymmErrors* gaeSigmaDecay=(TGraphAsymmErrors*)gae->Clone();
  gaeSigmaDecay->SetName("gaeSigmaDecay");

  TGraphAsymmErrors* gaeSigmaDecayv2=(TGraphAsymmErrors*)gae->Clone();
  gaeSigmaDecayv2->SetName("gaeSigmaDecayv2");

  double BRchain=6.09604e-5;
  double Fraction=0.401;

  //double norm=208.;
  norm=208.;
  double BRFraction=BRchain*Fraction;

  for (i=0;i<gaeSigmaDecay->GetN();i++){
    gaeSigmaDecay->GetY()[i] *= BRFraction*norm;
    //gaeSigmaDecay->SetPoint(i,gaeSigmaDecay->GetX()[i],gaeSigmaDecay->GetY()[i]*BRFraction*norm);
    //gaeSigmaDecay->SetPoint(i,gaeSigmaDecay->GetX()[i],gaeSigmaDecay->GetY()[i]*norm*Fraction);


    gaeSigmaDecay->SetPointEYhigh(i,gaeSigmaDecay->GetErrorYhigh(i)*BRFraction*norm);
    gaeSigmaDecay->SetPointEYlow(i,gaeSigmaDecay->GetErrorYlow(i)*BRFraction*norm); 

    //gaeSigmaDecayv2->GetY()[i] *= Fraction*norm*(1.0e-6);
    gaeSigmaDecayv2->SetPoint(i,gaeSigmaDecayv2->GetX()[i],gaeSigmaDecayv2->GetY()[i]*Fraction*norm*(1.0e-6));
    gaeSigmaDecayv2->SetPointEYhigh(i,gaeSigmaDecayv2->GetErrorYhigh(i)*Fraction*norm*(1.0e-6));
    gaeSigmaDecayv2->SetPointEYlow(i,gaeSigmaDecayv2->GetErrorYlow(i)*Fraction*norm*(1.0e-6)); 


  }

  TH1F* hpt_rebinv2=(TH1F*)hpt_rebin->Clone();
  TH1F* hminall_rebinv2=(TH1F*)hminall_rebin->Clone();
  TH1F* hmaxall_rebinv2=(TH1F*)hmaxall_rebin->Clone();
  hpt_rebinv2->SetName("hpt_rebinv2");
  hminall_rebinv2->SetName("hminall_rebinv2");
  hmaxall_rebinv2->SetName("hmaxall_rebinv2");



  for (i=0;i<5;i++) {
    hpt_rebinv2->SetBinContent(i+1,hpt_rebin->GetBinContent(i+1)*Fraction*norm*(1.0e-6));
    hminall_rebinv2->SetBinContent(i+1,hminall_rebin->GetBinContent(i+1)*Fraction*norm*(1.0e-6));
    hmaxall_rebinv2->SetBinContent(i+1,hmaxall_rebin->GetBinContent(i+1)*Fraction*norm*(1.0e-6));


  }

  std::cout << "GetY: " << gaeSigmaDecay->GetY()[0] << std::endl; 
  //###TF1* fitft2 = new TF1("fitft2","[0]+[1]*x*x",-3.0,3.0);
  TF1* fitft2 = new TF1("fitft2","[0]+[1]*pow((x-0.465),2)",-3.0,3.0);


  std::cout << "####### Fit with FONLL 7 TeV, pT 10,60, |y|<3.0 #########" << std::endl;
  /*    
	if (option==7 || option==71 || option==72) gaeSigmaDecay->Fit(fitft2,"","",-2.4,2.4);
	else gaeSigmaDecay->Fit(fitft2,"","",-2.865,1.935);
   */
  /*
     if (option==7 || option==71 || option==72) gaeSigmaDecayv2->Fit(fitft2,"","",-3.0,3.0);
     else gaeSigmaDecayv2->Fit(fitft2,"","",-2.865,1.935);
   */


  gaeSigmaDecay->SetFillColor(2);
  gaeSigmaDecay->SetFillStyle(3001); 
  if (option==5 || option==7 || option==71 || option==72) gaeSigmaDecay->SetTitle(";y_{LAB};d#sigma(B^{+} full chain)/dy #times A");
  else gaeSigmaDecay->SetTitle(";y_{CM};d#sigma(B^{+} full chain)/dy #times A");

  //#TH2F* hempty=new TH2F("hempty","",10,0,70.,10.,1.,500000);  
  //hempty=new TH2F("hempty","",10,0,70.,10.,1.,500000);  

  //TH2F* hempty;
  /*
     if (option==7 || option==71 || option==72) hempty=new TH2F("hempty","",50,-2.5,2.5,150.,0,150000);
  //hempty->SetAxisRange(0.,22000000.,"Y");
  else hempty=new TH2F("hempty","",50,-2.965,2.035,15.,0,15000);
   */
  if (option==5 || option==7 || option==71 || option==72) hempty=new TH2F("hempty","",60,-3.0,3.0,150.,0,150000);
  //hempty->SetAxisRange(0.,22000000.,"Y");
  else hempty=new TH2F("hempty","",50,-2.965,2.035,15.,0,15000);


  if (option==5 || option==7 || option==71 || option==72) hempty->GetXaxis()->SetTitle("y_{LAB}");
  else hempty->GetXaxis()->SetTitle("y_{CM}");
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
  //hempty->GetYaxis()->SetTitle("d#sigma (B^{+} full chain)/dp_{T} #times A (pb c/GeV)");
  hempty->GetYaxis()->SetTitle("d#sigma/dy(B^{+}) #times A (pb)");

  TCanvas*canvas=new TCanvas("canvas","canvas",600,500);
  //canvas->SetLogy();
  hempty->Draw();
  gaeSigmaDecay->Draw("psame");

  //#TLatex * tlatex=new TLatex(0.2,0.85,"B^{+}=40.1%, |y_{LAB}|<2.4, BR unc not shown");
  tlatex=new TLatex(0.2,0.85,"B^{+}=40.1%, 10<p_{T}<60GeV/c, BR unc not shown");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.04);
  tlatex->Draw();

  gae->SetName("gaeBplus");
  gaeSigmaDecay->SetName("gaeSigmaDecayBplus");
  gaeSigmaDecayv2->SetName("gaeSigmaDecayv2Bplus");


  //canvas->SaveAs("Plots/canvasBplusFinePt_7TeV.eps");
  //TFile*foutput=new TFile("../outputBplusFinePt_7TeV.root","recreate");
  canvas->SaveAs("Plots/canvasBplusyFinePt.eps");
  canvas->SaveAs(Form("Plots/canvasBplusyFinePt_%s.pdf",ename.c_str()));

  hmaxall_rebinv2->GetYaxis()->SetRangeUser(0,500);
  hmaxall_rebinv2->Draw();
  hpt_rebinv2->Draw("same");
  hminall_rebinv2->Draw("same");
  canvas->SaveAs(Form("Plots/canvasBplusyFinePt_%s_Binned.pdf",ename.c_str()));



  //TFile*foutput=new TFile("../outputBplusyFinePt.root","recreate");

  foutput->cd();
  gae->Write();
  gaeSigmaDecay->Write();
  gaeSigmaDecayv2->Write();

  /*
     hpt_rebin->Write();
     hminall_rebin->Write();
     hmaxall_rebin->Write();
     hminsc_rebin->Write();
     hmaxsc_rebin->Write();
     hminmass_rebin->Write();
     hmaxmass_rebin->Write();
     hminpdf_rebin->Write();
     hmaxpdf_rebin->Write();
   */
  hpt->Write();
  hminall->Write();
  hmaxall->Write();
  hminsc->Write();
  hmaxsc->Write();
  hminmass->Write();
  hmaxmass->Write();
  hminpdf->Write();
  hmaxpdf->Write();
  hfr_05_05->Write();
  hfr_20_20->Write();
  hfr_20_10->Write();
  hfr_10_20->Write();
  hfr_10_05->Write();
  hfr_05_10->Write();
  hpt_rebin->Write();
  hminall_rebin->Write();
  hmaxall_rebin->Write();


  hpt_rebinv2->Write();
  hminall_rebinv2->Write();
  hmaxall_rebinv2->Write();

  //foutput->Write();

}
