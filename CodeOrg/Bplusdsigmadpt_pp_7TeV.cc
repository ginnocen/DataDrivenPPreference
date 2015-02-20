#include "TH1F.h"
#include <cmath>
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>

#define BIN_NUM 460 //pPb_pt:220,pPb_y:40,pp_pt:222,pp_y:45
#define REBINa 8     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
#define REBINpa 9    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5
#define HMIN 5      //pPb_pt:5,pPb_y:-2,pp_pt:9,pp_y:0
#define HMAX 120     //pPb_pt:55,pPb_y:2,pp_pt:120,pp_y:2.25
#define REBIN 5     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
#define REBINp 6    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5

int Bplusdsigmadpt_pp_7TeV()
{
  TString infilea;
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  infilea="../FONLLInputs/fo_pp_pt_rap225_7TeV.dat";
  
  ifstream getdataa(infilea.Data());

  if(!getdataa.is_open())
    {
      cout<<"Opening the file fails"<<endl;
    }

  float centrala[BIN_NUM];
  float min_alla[BIN_NUM],max_alla[BIN_NUM],min_sca[BIN_NUM],max_sca[BIN_NUM],min_massa[BIN_NUM],max_massa[BIN_NUM],min_pdfa[BIN_NUM],max_pdfa[BIN_NUM];
  int i;
  float tem;
  for(i=0;i<BIN_NUM;i++)
    {
      getdataa>>tem;
      getdataa>>centrala[i];
      getdataa>>min_alla[i];
      getdataa>>max_alla[i];
      getdataa>>min_sca[i];
      getdataa>>max_sca[i];
      getdataa>>min_massa[i];
      getdataa>>max_massa[i];
      getdataa>>min_pdfa[i];
      getdataa>>max_pdfa[i];
    }
  
  TH1F* hpta = new TH1F("hpta","",BIN_NUM,HMIN,HMAX);
  TH1F* hminalla = new TH1F("hminalla","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxalla = new TH1F("hmaxalla","",BIN_NUM,HMIN,HMAX);
  TH1F* hminsca = new TH1F("hminsca","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxsca = new TH1F("hmaxsca","",BIN_NUM,HMIN,HMAX);
  TH1F* hminmassa = new TH1F("hminmassa","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxmassa = new TH1F("hmaxmassa","",BIN_NUM,HMIN,HMAX);
  TH1F* hminpdfa = new TH1F("hminpdfa","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxpdfa = new TH1F("hmaxpdfa","",BIN_NUM,HMIN,HMAX);
  TH1F* hratioa = new TH1F("hratioa","",BIN_NUM,HMIN,HMAX);


  for(i=0;i<BIN_NUM;i++)
    {
      hpta->SetBinContent(i+1,centrala[i]);
      hminalla->SetBinContent(i+1,min_alla[i]);
      hmaxalla->SetBinContent(i+1,max_alla[i]);
      hminsca->SetBinContent(i+1,min_sca[i]);
      hmaxsca->SetBinContent(i+1,max_sca[i]);
      hminmassa->SetBinContent(i+1,min_massa[i]);
      hmaxmassa->SetBinContent(i+1,max_massa[i]);
      hminpdfa->SetBinContent(i+1,min_pdfa[i]);
      hmaxpdfa->SetBinContent(i+1,max_pdfa[i]);
    }

  double dataa[REBINa] = {103.4,36.03,15.33,6.056,1.814,0.3477,0.06244,0.006099};
  double staterrora[REBINa] = {4.,0.8,0.3,0.1,0.03,0.008,0.003,0.0006};
  double syserrora[REBINa] = {8.,2.3,1.0,0.4,0.12,0.028,0.005,.0007};

  //Rebin Edge
  double rebina[REBINpa] = {9.,13.,16.,20.,25.,35.,50.,70.,120.};

  TH1F* hpt_rebina = (TH1F*)hpta->Rebin(REBINa,"hpt_rebina",rebina);
  TH1F* hminall_rebina = (TH1F*)hminsca->Rebin(REBINa,"hminall_rebina",rebina);
  TH1F* hmaxall_rebina = (TH1F*)hmaxsca->Rebin(REBINa,"hmaxall_rebina",rebina);
  TH1F* hminsc_rebina = (TH1F*)hminsca->Rebin(REBINa,"hminsc_rebina",rebina);
  TH1F* hmaxsc_rebina = (TH1F*)hmaxsca->Rebin(REBINa,"hmaxsc_rebina",rebina);
  TH1F* hminmass_rebina = (TH1F*)hminmassa->Rebin(REBINa,"hminmass_rebina",rebina);
  TH1F* hmaxmass_rebina = (TH1F*)hmaxmassa->Rebin(REBINa,"hmaxmass_rebina",rebina);
  TH1F* hminpdf_rebina = (TH1F*)hminpdfa->Rebin(REBINa,"hminpdf_rebina",rebina);
  TH1F* hmaxpdf_rebina = (TH1F*)hmaxpdfa->Rebin(REBINa,"hmaxpdf_rebina",rebina);

  TH1F* hratio_rebina = (TH1F*)hratioa->Rebin(REBINa,"hratio_rebina",rebina);


  //bin middle
  double apta[REBINa] = {11.,14.5,18.,22.5,30.,42.5,60.,95.};//pPb_pt
  //bin half width
  double aptla[REBINa] = {2.,1.5,2.,2.5,5.,7.5,10.,25.};//pPb_pt
  double asigmaa[REBINa],aminalla[REBINa],amaxalla[REBINa],aminsca[REBINa],amaxsca[REBINa],aminmassa[REBINa],amaxmassa[REBINa],aminpdfa[REBINa],amaxpdfa[REBINa],aerrorla[REBINa],aerrorha[REBINa];

  //number of every rebined bin
  double bin_numa[REBINa] = {16.,12,16.,20.,40.,60.,80.,200.};//pPb_pt
  
  int j;
  double norm=1.;
  double inte=1.;

  for(j=0;j<REBINa;j++)
    {

      tem = hpt_rebina->GetBinContent(j+1);
      asigmaa[j] = tem*norm*inte/bin_numa[j];

      tem = hminall_rebina->GetBinContent(j+1);
      aminalla[j] = tem*norm*inte/bin_numa[j];

      tem = hmaxsc_rebina->GetBinContent(j+1);
      amaxalla[j] = tem*norm*inte/bin_numa[j];

      tem = hminsc_rebina->GetBinContent(j+1);
      aminsca[j] = tem*norm*inte/bin_numa[j];

      tem = hmaxsc_rebina->GetBinContent(j+1);
      amaxsca[j] = tem*norm*inte/bin_numa[j];

      tem = hminmass_rebina->GetBinContent(j+1);
      aminmassa[j] = tem*norm*inte/bin_numa[j];

      tem = hmaxmass_rebina->GetBinContent(j+1);
      amaxmassa[j] = tem*norm*inte/bin_numa[j];

      tem = hminpdf_rebina->GetBinContent(j+1);
      aminpdfa[j] = tem*norm*inte/bin_numa[j];

      tem = hmaxpdf_rebina->GetBinContent(j+1);
      amaxpdfa[j] = tem*norm*inte/bin_numa[j];

      aerrorla[j] = asigmaa[j]-aminalla[j];//all,sc,mass,pdf
      aerrorha[j] = amaxalla[j]-asigmaa[j];//all,sc,mass,pdf
    }

  cout<<"------- pp_7------"<<endl;
  cout<<endl;
 
  TGraphAsymmErrors* gaea = new TGraphAsymmErrors(REBINa, apta, asigmaa, aptla, aptla, aerrorla, aerrorha);
  gaea->SetTitle(";p_{T}(GeV/c);d#sigma (B admix) /dp_{T}(pb c/GeV)");
  gaea->SetFillColor(2);
  gaea->SetFillStyle(3001);

  TGraphAsymmErrors* gaeSigmaDecaya=(TGraphAsymmErrors*)gaea->Clone();
  gaeSigmaDecaya->SetName("gaeSigmaDecaya");
  double BRchain=1.;
  double Fraction=0.401;
  double BR=0.001016*0.0593;
  double norm=1.;
  double BRFraction=BRchain*Fraction;

  double syserrya[REBINa],syserreya[REBINa];
  
  for (int i=0;i<gaeSigmaDecaya->GetN();i++)
    {
      gaeSigmaDecaya->GetY()[i] *= BRFraction*norm/1000000;
      gaeSigmaDecaya->SetPointEYhigh(i,gaeSigmaDecaya->GetErrorYhigh(i)*BRFraction*norm/1000000);
      gaeSigmaDecaya->SetPointEYlow(i,gaeSigmaDecaya->GetErrorYlow(i)*BRFraction*norm/1000000);
      hratio_rebina->SetBinContent(i+1,dataa[i]/(gaeSigmaDecaya->GetY()[i]*1000000*BR));
      hratio_rebina->SetBinError(i+1,staterrora[i]/(gaeSigmaDecaya->GetY()[i]*1000000*BR));
      syserrya[i] = dataa[i]/(gaeSigmaDecaya->GetY()[i]*1000000*BR);
      syserreya[i] = syserrora[i]/(gaeSigmaDecaya->GetY()[i]*1000000*BR);
    }

  gaeSigmaDecaya->SetFillColor(2);
  gaeSigmaDecaya->SetFillStyle(3001); 
  gaeSigmaDecaya->SetTitle(";p_{T}(GeV/c);d#sigma/dp_{T} (B^{+}) #times A (GeV^{-1}c)");

  TCanvas* cratioa = new TCanvas("cratioa","cratioa",500,500);

  hratio_rebina->SetMaximum(2.);
  hratio_rebina->SetMinimum(0.5);
  hratio_rebina->SetXTitle("p_{T}(GeV/c)");
  hratio_rebina->SetYTitle("#sigma / #sigma(FONLL)");
  hratio_rebina->SetTitleOffset(1.2,"Y");
  hratio_rebina->SetMarkerStyle(22);
  hratio_rebina->SetLineWidth(2);
  hratio_rebina->SetMarkerColor(kRed+2);
  hratio_rebina->SetLineColor(kRed+2);
  hratio_rebina->SetStats(0);
  hratio_rebina->Draw("lep");
  
  TGraphErrors* gsyserrora = new TGraphErrors(REBINa,apta,syserrya,aptla,syserreya);  
  gsyserrora->SetMarkerColor(1);
  gsyserrora->SetLineColor(kRed+2);
  gsyserrora->SetLineWidth(2);
  gsyserrora->SetMarkerStyle(22);
  gsyserrora->SetMarkerColor(kRed+2);
  gsyserrora->SetFillColor(0);
  gsyserrora->SetFillStyle(0);
  gsyserrora->Draw("2same");

  TLine* lin0a=new TLine(9,1,120,1);
  lin0a->SetLineStyle(2);
  lin0a->SetLineColor(1);
  lin0a->SetLineWidth(3);
  lin0a->Draw("same");

  
  /////////////////////////////////////////////////////////////////////////////////////////////////

  TString infile, outfile;
  infile="../FONLLInputs/fo_pp_pt_rap24_7TeV.dat";
  outfile="Rootf/outputBplus_pp_pt_rap24_7TeV_CmsBin.root";
  
  ifstream getdata(infile.Data());

  if(!getdata.is_open())
    {
      cout<<"Opening the file fails"<<endl;
    }

  float central[BIN_NUM];
  float min_all[BIN_NUM],max_all[BIN_NUM],min_sc[BIN_NUM],max_sc[BIN_NUM],min_mass[BIN_NUM],max_mass[BIN_NUM],min_pdf[BIN_NUM],max_pdf[BIN_NUM];
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

  double data[REBIN] = {4.07,1.47,0.412,0.181,0.042}; //#mub
  double staterror[REBIN] = {0.47,0.13,0.041,0.015,0.007};
  double syserror[REBIN] = {0.31,0.09,0.026,0.012,0.004};

  //Rebin Edge
  double rebin[REBINp] = {5.,10.,13.,17.,24.,30.};

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
  double apt[REBIN] = {7.5,11.5,15,20.5,27};
  //bin half width
  double aptl[REBIN] = {2.5,1.5,2.,3.5,3.};
  double asigma[REBIN],aminall[REBIN],amaxall[REBIN],aminsc[REBIN],amaxsc[REBIN],aminmass[REBIN],amaxmass[REBIN],aminpdf[REBIN],amaxpdf[REBIN],aerrorl[REBIN],aerrorh[REBIN];

  //number of every rebined bin
  double bin_num[REBIN] = {20.,12.,16.,28.,24.};
  
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

  TGraphAsymmErrors* gaeSigmaDecay=(TGraphAsymmErrors*)gae->Clone();
  gaeSigmaDecay->SetName("gaeSigmaDecay");
  double syserry[REBIN],syserrey[REBIN];
  BR=1;
  for (int i=0;i<gaeSigmaDecay->GetN();i++)
    {
      gaeSigmaDecay->GetY()[i] *= BRFraction*norm/1000000;
      gaeSigmaDecay->SetPointEYhigh(i,gaeSigmaDecay->GetErrorYhigh(i)*BRFraction*norm/1000000);
      gaeSigmaDecay->SetPointEYlow(i,gaeSigmaDecay->GetErrorYlow(i)*BRFraction*norm/1000000);
      hratio_rebin->SetBinContent(i+1,data[i]/(gaeSigmaDecay->GetY()[i]*BR));
      hratio_rebin->SetBinError(i+1,staterror[i]/(gaeSigmaDecay->GetY()[i]*BR));
      syserry[i] = data[i]/(gaeSigmaDecay->GetY()[i]*BR);
      syserrey[i] = syserror[i]/(gaeSigmaDecay->GetY()[i]*BR); 
    }

  gaeSigmaDecay->SetFillColor(2);
  gaeSigmaDecay->SetFillStyle(3001); 
  gaeSigmaDecay->SetTitle(";p_{T}(GeV/c);d#sigma/dp_{T} (B^{+}) #times A (GeV^{-1}c)");
   
  hratio_rebin->SetMaximum(2.);
  hratio_rebin->SetMinimum(0.5);
  hratio_rebin->SetXTitle("p_{T}(GeV/c)");
  hratio_rebin->SetYTitle("#sigma / #sigma(FONLL)");
  hratio_rebin->SetTitleOffset(1.2,"Y");
  hratio_rebin->SetMarkerStyle(8);
  hratio_rebin->SetMarkerColor(kAzure-6);
  hratio_rebin->SetLineColor(kAzure-6);
  hratio_rebin->SetLineWidth(2);
  hratio_rebin->SetStats(0);
  hratio_rebin->Draw("lep same");

  TGraphErrors* gsyserror = new TGraphErrors(REBIN,apt,syserry,aptl,syserrey);  
  gsyserror->SetMarkerColor(1);
  gsyserror->SetLineColor(kAzure-6);
  gsyserror->SetLineWidth(2);   
  gsyserror->SetMarkerStyle(8);
  gsyserror->SetMarkerColor(kAzure-6);
  gsyserror->SetFillColor(0);
  gsyserror->SetFillStyle(0);
  gsyserror->Draw("2same");

  TLegend *lega = new TLegend(0.6,0.82,0.9,0.9);
  lega->AddEntry((TObject*)0,"pp 7 TeV","");
  lega->SetBorderSize(0);
  lega->SetFillStyle(0);
  lega->Draw("same");

  TLegend *leg2a = new TLegend(0.63,0.65,0.98,0.82);
  TLegendEntry *entatlas=leg2a->AddEntry(gsyserrora,"ATLAS","pf");
  TLegendEntry *entcms=leg2a->AddEntry(gsyserror,"CMS","pf");

  leg2a->SetBorderSize(0);
  leg2a->SetFillStyle(0);
  leg2a->Draw("same");

  cratioa->SaveAs("Plots/cratio_pp_pt_rap24_7TeV.pdf");
}
