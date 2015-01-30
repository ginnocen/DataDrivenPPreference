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

void GetBinned (int codeNum, int codeDen)
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	gStyle->SetPadTopMargin(0.095);
	gStyle->SetPadLeftMargin(0.130);
	gStyle->SetPadRightMargin(0.045);

	std::string codeName[3]={"5TeV","7TeV","2760GeV"};

	TFile* fout = new TFile(Form("../ResultsBplus/Estimatedpp%s_with%s.root",codeName[codeNum].c_str(),codeName[codeDen].c_str()),"RECREATE");

	TF1* fitft = new TF1 ("fitft","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);

	TFile* fin_fitpar = new TFile(Form("../ResultsBplus/ScaleCor_reweighted%sdata.root",codeName[codeDen].c_str()));
	TH1D* hparafit = (TH1D*)fin_fitpar->Get("hparafit");
	fitft->SetParameter(0,hparafit->GetBinContent(1));
	fitft->SetParameter(1,hparafit->GetBinContent(2));
	fitft->SetParameter(2,hparafit->GetBinContent(3));

	double ptbin[7]={5.0,10.0,15.0,20.0,25.0,30.0,60.0};

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
	double selerr;

	// Estimated error before consider A
	double eststa_befA[6];
	double estsys_befA[6];

	double apt_wgt[6],aptl_wgt[6],apth_wgt[6];

	TFile* fin_ErrPer = new TFile(Form("../ResultsBplus/CalcErrPerc_%s.root",codeName[codeDen].c_str()));
	TH1D* hstaerrPerc = (TH1D*)fin_ErrPer->Get("hstaerrPrec");
	TH1D* hsyserrPerc = (TH1D*)fin_ErrPer->Get("hsyserrPrec");
	TH1D* hWgtcen = (TH1D*)fin_ErrPer->Get("hWgtcen");

	for (int i=0;i<6;i++) {
		eststa_befA[i]=hstaerrPerc->GetBinContent(i+1);
		estsys_befA[i]=hsyserrPerc->GetBinContent(i+1);
	}

	// recalculated with our binning
	std::cout << "recalculated dsigma/dp_{T} for " << codeName[codeDen].c_str() << ", not times A=208" << std::endl;
	for (int i=0;i<6;i++) {
		apt[i] = (ptbin[i]+ptbin[i+1])/2;//centered point
		aptl[i] = (ptbin[i+1]-ptbin[i])/2;//half bin width
		apt_wgt[i] = hWgtcen->GetBinContent(i+1);
		aptl_wgt[i] = apt_wgt[i]-ptbin[i];//left bin width
		apth_wgt[i] = ptbin[i+1]-apt_wgt[i];//right bin width

		asigma[i] = (fitft->Integral(ptbin[i],ptbin[i+1]))/(ptbin[i+1]-ptbin[i]);//recalculated y, not times A
		asigmaNoSca[i] = (fitft->Integral(ptbin[i],ptbin[i+1]))/(ptbin[i+1]-ptbin[i]);//recalculated y, not times A
		aerrorhStat[i] = asigmaNoSca[i]*eststa_befA[i];
		aerrorlStat[i] = asigmaNoSca[i]*eststa_befA[i];
		aerrorhSyst[i] = asigmaNoSca[i]*estsys_befA[i];
		aerrorlSyst[i] = asigmaNoSca[i]*estsys_befA[i];

		std :: cout << "# " << i << " : " << asigma[i] << std::endl;
		selerr = sqrt(pow(eststa_befA[i],2)+pow(estsys_befA[i],2));//recalculated errorPerc from pp data, not times A
		std :: cout << std::string(10,'-') << " error from data - stat. : " << eststa_befA[i] << " , syst. : " << estsys_befA[i] << "total : " << selerr << std::endl;

		asigma[i] = (fitft->Integral(ptbin[i],ptbin[i+1]))/(ptbin[i+1]-ptbin[i])*scalefac[i];//scaled for 2.76TeV
		aerrorh[i] = sqrt(pow(scalefac_plerr[i],2)+pow(selerr,2))*(asigma[i]);
		aerrorl[i] = sqrt(pow(scalefac_mierr[i],2)+pow(selerr,2))*(asigma[i]);


	}
	TGraphAsymmErrors* gaeEstimatedpp5TeV = new TGraphAsymmErrors(6,apt,asigma,aptl,aptl,aerrorl,aerrorh);
	gaeEstimatedpp5TeV->SetName("gaeEstimatedpp5TeV");

	TGraphAsymmErrors* gaeStatNoSca7TeV = new TGraphAsymmErrors(6,apt,asigmaNoSca,aptl,aptl,aerrorlStat,aerrorhStat);
	gaeStatNoSca7TeV->SetName("gaeStatNoSca7TeV");
	TGraphAsymmErrors* gaeSystNoSca7TeV = new TGraphAsymmErrors(6,apt,asigmaNoSca,aptl,aptl,aerrorlSyst,aerrorhSyst);
	gaeSystNoSca7TeV->SetName("gaeSystNoSca7TeV");




	TGraphAsymmErrors* gaeEstimatedpp5TeV_wgt = new TGraphAsymmErrors(6,apt_wgt,asigma,aptl_wgt,apth_wgt,aerrorl,aerrorh);
	gaeEstimatedpp5TeV_wgt->SetName("gaeEstimatedpp5TeV_wgt");

	// CMS 7TeV, pp->(B+)X, unit : micro barn
	double ptbin_ORI[NUM_ORI+1]={5.,10.,13.,17.,24.,30.};
	double sigmapt_ORI[NUM_ORI]={4.07,1.47,0.412,0.181,0.042};
	double sigmapt_sta_ORI[NUM_ORI]={0.47,0.13,0.041,0.015,0.007};
	double asigma_Recal_pp7[NUM_ORI];

	TH1D* hsigmapt_ORI = new TH1D("hsigmapt_ORI","",NUM_ORI,ptbin_ORI);
	TH1D* hsigmapt_Recal = new TH1D("hsigmapt_Recal","",NUM_ORI,ptbin_ORI);
	TH1D* hsigmapt_Comp = new TH1D("hsigmapt_Comp","",NUM_ORI,ptbin_ORI);

	hsigmapt_ORI->Sumw2();
	hsigmapt_Recal->Sumw2();
	hsigmapt_Comp->Sumw2();

	std::cout << "recalculated dsigma/dp_{T} for " << codeName[codeDen].c_str() << ", not times A=208 with CMS pp 7 TeV binning" << std::endl;
	for (int i=0;i<5;i++) {
		asigma_Recal_pp7[i] = (fitft->Integral(ptbin_ORI[i],ptbin_ORI[i+1]))/(ptbin_ORI[i+1]-ptbin_ORI[i]);//recalculated y, not times A
		double rat=asigma_Recal_pp7[i]/sigmapt_ORI[i];
		hsigmapt_ORI->SetBinContent(i+1,sigmapt_ORI[i]);
		hsigmapt_ORI->SetBinError(i+1,sigmapt_sta_ORI[i]);
		hsigmapt_Recal->SetBinContent(i+1,asigma_Recal_pp7[i]);
		hsigmapt_Recal->SetBinError(i+1,sigmapt_sta_ORI[i]);
		std :: cout << "# " << i << " Origin : " << sigmapt_ORI[i] << " --- Recalculated : " << asigma_Recal_pp7[i] << " Recalculated/Origin : " << rat << std::endl;
	}
	hsigmapt_Comp->Divide(hsigmapt_Recal,hsigmapt_ORI,1,1,"B");

	TCanvas *c1 = new TCanvas("c1","",500,500);
	c1->SetLogy(1);

	gaeEstimatedpp5TeV->SetFillColor(kPink-4);
	gaeEstimatedpp5TeV->SetFillStyle(3001);
	gaeEstimatedpp5TeV->SetMarkerStyle(25);
	gaeEstimatedpp5TeV->SetMarkerSize(1.5);
	gaeEstimatedpp5TeV->SetMarkerColor(kViolet+1);
	gaeEstimatedpp5TeV->SetTitle(";p_{T}(GeV/c);d#sigma(B^{+} full chain)/dp_{T} (#mub GeV^{-1}c)");

	gaeEstimatedpp5TeV_wgt->SetFillColor(kPink-4);
	gaeEstimatedpp5TeV_wgt->SetFillStyle(3001);
	//gaeEstimatedpp5TeV_wgt->SetFillColor(0);
	//gaeEstimatedpp5TeV_wgt->SetFillStyle(0);
	gaeEstimatedpp5TeV_wgt->SetMarkerStyle(33);
	gaeEstimatedpp5TeV_wgt->SetMarkerSize(2);
	gaeEstimatedpp5TeV_wgt->SetMarkerColor(kViolet+1);
	gaeEstimatedpp5TeV_wgt->SetLineColor(kViolet+1);
	gaeEstimatedpp5TeV_wgt->SetLineWidth(1);
	gaeEstimatedpp5TeV_wgt->SetTitle(";p_{T}(GeV/c);d#sigma(B^{+} full chain)/dp_{T} (#mub GeV^{-1}c)");

	TH2F* hempty=new TH2F("hempty","",10,0,70.,10.,0.000001,1000000.0);
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
	hempty->GetYaxis()->SetLimits(0.00001,50.0);
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();

	hempty->Draw("");
	gaeEstimatedpp5TeV->Draw("samee2");
	gaeEstimatedpp5TeV_wgt->Draw("samee2");
	gaeEstimatedpp5TeV->Draw("samep");
	gaeEstimatedpp5TeV_wgt->Draw("samep");

	TLegend* leg = new TLegend(0.15,0.70,0.45,0.90,"");

	leg->SetBorderSize(0);
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	leg->SetTextFont(42);
	leg->SetTextSize(0.040);

	leg->AddEntry(gaeEstimatedpp5TeV,Form("Estimated CMS %s pp data",codeName[codeNum].c_str()),"lpf");
	leg->AddEntry(gaeEstimatedpp5TeV_wgt,Form("#splitline{Estimated CMS %s pp data}{with weighted center}",codeName[codeNum].c_str()),"lpf");
	leg->Draw("");

	c1->SaveAs(Form("../ResultsBplus/gaeEstimatedpp%s_from%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));

	c1->Clear();
	c1->SetLogy(1);

	hsigmapt_ORI->SetMarkerStyle(24);
	hsigmapt_ORI->SetMarkerSize(1.5);
	hsigmapt_ORI->SetMarkerColor(kRed);
	hsigmapt_ORI->SetLineColor(kRed);
	hsigmapt_ORI->SetLineWidth(2);

	hsigmapt_Recal->SetMarkerStyle(20);
	hsigmapt_Recal->SetMarkerSize(1.5);
	hsigmapt_Recal->SetMarkerColor(kRed-9);
	hsigmapt_Recal->SetLineColor(kRed-9);
	hsigmapt_Recal->SetLineWidth(2);

	hempty=new TH2F("hempty","",35,0,35.,10.,0.001,1000.0);
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
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();
	hempty->Draw("");
	hsigmapt_ORI->Draw("samepe1");
	hsigmapt_Recal->Draw("samepe1");

	leg->DeleteEntry();
	leg->DeleteEntry();
	leg->AddEntry(hsigmapt_ORI,Form("Published CMS %s pp data",codeName[codeDen].c_str()),"lp");
	leg->AddEntry(hsigmapt_Recal,Form("Recalculated CMS %s pp data with FONLL",codeName[codeDen].c_str()),"lp");
	leg->Draw("");

	c1->SaveAs(Form("../ResultsBplus/CompRecal%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
	c1->Clear();
	c1->SetLogy(0);
	hsigmapt_Comp->SetMarkerStyle(34);
	hsigmapt_Comp->SetMarkerSize(1.5);
	hsigmapt_Comp->SetMarkerColor(kRed-9);
	hsigmapt_Comp->SetLineColor(kRed-9);
	hsigmapt_Comp->SetLineWidth(2);

	hempty=new TH2F("hempty","",35,0,35.,10.,0.8,1.3);
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
	hempty->GetYaxis()->SetTitle("Ratio of d#sigma/dp_{T}(B^{+}), Recalculated./Estimated.");
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();

	hempty->Draw("");
	hsigmapt_Comp->Draw("samepe1");
	c1->SaveAs(Form("../ResultsBplus/CompRecalRatio%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));

	fout->cd(); 
	gaeEstimatedpp5TeV->Write();
	gaeEstimatedpp5TeV_wgt->Write();
	hsigmapt_ORI->Write();
	hsigmapt_Recal->Write();
	hsigmapt_Comp->Write();
	gaeStatNoSca7TeV->Write();
	gaeSystNoSca7TeV->Write();
	fout->Write();
}


