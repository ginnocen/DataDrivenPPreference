#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include "TGraphAsymmErrors.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>

#define APb 208
#define NUM_ORI 5

void GetBinned_w2760GeVpp(int codeNum, int codeDen)
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	gStyle->SetPadTopMargin(0.095);
	gStyle->SetPadLeftMargin(0.130);
	gStyle->SetPadRightMargin(0.045);

	//TFile* fin_FONLL=new TFile("../ResultsBplus_pp/dSigmadpt_FONLL.root");
	//TFile* fin_FONLL=new TFile("../ResultsBplus/dSigmadpt_FONLL_pp2760GeV.root");
	TFile* fin_FONLL=new TFile("../ResultsBplus_2760GeVpp/dSigmadpt_Bplus_2760GeVpp.root");

  TGraphAsymmErrors* gaeDatastat = (TGraphAsymmErrors*)fin_FONLL->Get("gSigmastat");
	TGraphAsymmErrors* gaeDatasyst = (TGraphAsymmErrors*)fin_FONLL->Get("gSigmasyst");

	std::string codeName[3]={"5TeV","7TeV","2760GeV"};

	TFile* fout = new TFile(Form("../ResultsBplus/Estimatedpp%s_with%s.root",codeName[codeNum].c_str(),codeName[codeDen].c_str()),"RECREATE");

	double ptbin[6]={10.0,15.0,20.0,25.0,30.0,60.0};

	double scalefac[6];
	double scalefac_plerr[6];
	double scalefac_mierr[6];

	TFile* fin_Val;
	TFile* fin_Norm;
	TH1D* hfrval;
	TH1D* hfrmaxerr;
	TH1D* hfrminerr;

	if (codeNum!=codeDen) {
		fin_Val=new TFile(Form("../ResultsBplus/CompFONLL_Bplus_Binned_Val_%svs%s.root",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
		fin_Norm=new TFile(Form("../ResultsBplus/CompFONLL_Bplus_Binned_Norm_%svs%s.root",codeName[codeNum].c_str(),codeName[codeDen].c_str()));

		hfrval=(TH1D*)fin_Val->Get("hfrval");
		hfrmaxerr=(TH1D*)fin_Norm->Get("hfrmaxerr");
		hfrminerr=(TH1D*)fin_Norm->Get("hfrminerr");
	}

	for (int i=0;i<6;i++) {
		if (codeNum!=codeDen) {
			scalefac[i]=hfrval->GetBinContent(i+1);
			scalefac_plerr[i]=hfrmaxerr->GetBinContent(i+1);
			scalefac_mierr[i]=hfrminerr->GetBinContent(i+1);
		}
		else {
			scalefac[i]=1.0;
			scalefac_plerr[i]=0.0;
			scalefac_mierr[i]=0.0;
		}
	}

	double apt[6],aptl[6],asigma[6],aerrorh[6],aerrorl[6];

  double asigmaNoSca[6];
  double aerrorhStat[6],aerrorlStat[6];
  double aerrorhSyst[6],aerrorlSyst[6];

	//double selerr = 0.192;//sqrt(0.167^2+0.095^2)
  //double selerr;
  double selerrl, selerrh;

	// recalculated with our binning
	std::cout << "recalculated dsigma/dp_{T} for " << codeName[codeDen].c_str() << ", not times A=208" << std::endl;
	for (int i=0;i<5;i++) {
    apt[i] = (ptbin[i]+ptbin[i+1])/2;//centered point
		aptl[i] = (ptbin[i+1]-ptbin[i])/2;//half bin width

		asigma[i] = gaeDatastat->GetY()[i]*scalefac[i+1];

    asigmaNoSca[i] = gaeDatastat->GetY()[i];
    aerrorhStat[i] = gaeDatastat->GetEYhigh()[i];
    aerrorlStat[i] = gaeDatastat->GetEYlow()[i];
    aerrorhSyst[i] = gaeDatasyst->GetEYhigh()[i];
    aerrorlSyst[i] = gaeDatasyst->GetEYlow()[i];

		selerrl = sqrt(pow(gaeDatastat->GetEYlow()[i]/gaeDatastat->GetY()[i],2)+pow(gaeDatasyst->GetEYlow()[i]/gaeDatasyst->GetY()[i],2));//recalculated errorPerc from pp data, not times A
		selerrh = sqrt(pow(gaeDatastat->GetEYhigh()[i]/gaeDatastat->GetY()[i],2)+pow(gaeDatasyst->GetEYhigh()[i]/gaeDatasyst->GetY()[i],2));//recalculated errorPerc from pp data, not times A
		aerrorl[i] = sqrt(pow(scalefac_mierr[i+1],2)+pow(selerrl,2))*(asigma[i]);
		aerrorh[i] = sqrt(pow(scalefac_plerr[i+1],2)+pow(selerrh,2))*(asigma[i]);
	}
/*
	TGraphAsymmErrors* gaeEstimatedpp5TeV = new TGraphAsymmErrors(6,apt,asigma,aptl,aptl,aerrorl,aerrorh);
	gaeEstimatedpp5TeV->SetName("gaeEstimatedpp5TeV");
  TGraphAsymmErrors* gaeStatNoSca2760GeV = new TGraphAsymmErrors(6,apt,asigmaNoSca,aptl,aptl,aerrorlStat,aerrorhStat);
  gaeStatNoSca2760GeV->SetName("gaeStatNoSca2760GeV");
  TGraphAsymmErrors* gaeSystNoSca2760GeV = new TGraphAsymmErrors(6,apt,asigmaNoSca,aptl,aptl,aerrorlSyst,aerrorhSyst);
  gaeSystNoSca2760GeV->SetName("gaeSystNoSca2760GeV");
*/
	TGraphAsymmErrors* gaeEstimatedpp5TeV = new TGraphAsymmErrors(5,apt,asigma,aptl,aptl,aerrorl,aerrorh);
	gaeEstimatedpp5TeV->SetName("gaeEstimatedpp5TeV");
  TGraphAsymmErrors* gaeStatNoSca2760GeV = new TGraphAsymmErrors(5,apt,asigmaNoSca,aptl,aptl,aerrorlStat,aerrorhStat);
  gaeStatNoSca2760GeV->SetName("gaeStatNoSca2760GeV");
  TGraphAsymmErrors* gaeSystNoSca2760GeV = new TGraphAsymmErrors(5,apt,asigmaNoSca,aptl,aptl,aerrorlSyst,aerrorhSyst);
  gaeSystNoSca2760GeV->SetName("gaeSystNoSca2760GeV");

	TCanvas *c1 = new TCanvas("c1","",500,500);
	c1->SetLogy(1);

	gaeEstimatedpp5TeV->SetFillColor(kPink-4);
	gaeEstimatedpp5TeV->SetFillStyle(3001);
	gaeEstimatedpp5TeV->SetMarkerStyle(25);
	gaeEstimatedpp5TeV->SetMarkerSize(1.5);
	gaeEstimatedpp5TeV->SetMarkerColor(kViolet+1);
	gaeEstimatedpp5TeV->SetTitle(";p_{T}(GeV/c);d#sigma(B^{+} full chain)/dp_{T} (#mub GeV^{-1}c)");

	TH2F* hempty=new TH2F("hempty","",10,0,70.,10.,0.0001,100000.0);
	hempty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hempty->GetXaxis()->SetTitleOffset(1.0);
	hempty->GetYaxis()->SetTitleOffset(1.3);
	hempty->GetXaxis()->SetTitleSize(0.045);
	hempty->GetYaxis()->SetTitleSize(0.045);
	hempty->GetXaxis()->SetTitleFont(42);
	hempty->GetYaxis()->SetTitleFont(42);
	hempty->GetXaxis()->SetLabelFont(42);
	hempty->GetYaxis()->SetLabelFont(42);
	hempty->GetXaxis()->SetLabelSize(0.04);
	hempty->GetYaxis()->SetLabelSize(0.04);
	hempty->GetYaxis()->SetTitle("d#sigma/dp_{T}(B^{+}) (#mub c/GeV)");
	//hempty->GetYaxis()->SetRangeUser(1.0,1000000.0);
	hempty->GetYaxis()->SetLimits(0.0001,100.0);
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();

	hempty->Draw("");
	gaeEstimatedpp5TeV->Draw("samee2");
	gaeEstimatedpp5TeV->Draw("samep");

	//TLegend* leg = new TLegend(0.15,0.70,0.45,0.90,"");
	TLegend* leg = new TLegend(0.15,0.80,0.45,0.90,"");

	leg->SetBorderSize(0);
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	leg->SetTextFont(42);
	leg->SetTextSize(0.040);

	leg->AddEntry(gaeEstimatedpp5TeV,Form("Estimated CMS %s pp data",codeName[codeNum].c_str()),"lpf");
	leg->Draw("");

	c1->SaveAs(Form("../ResultsBplus/gaeEstimatedpp%s_from%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));

	fout->cd(); 
	gaeEstimatedpp5TeV->Write();
  gaeStatNoSca2760GeV->Write();
  gaeSystNoSca2760GeV->Write();
	fout->Write();
}


