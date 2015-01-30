#include <math.h>
#include <iostream>
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegendEntry.h"
#include "TPad.h"
#include "TF1.h"
#include "TColor.h"

#define NUM 5
//#define APb_ORI 208
#define APb 1


#define NUM_ATL 8

using namespace std;

#define BIN_NUM 5

void Comp_dsigma(int option){

	// 1(5.02 TeV pPb), 2(2.76 TeV), 3(5.02 TeV with 2.76 TeV pp data)
/*
	double APb;
	if (option==1 || option==2) APb=1.;
	else APb=APb_ORI;
*/
	std::string rmk;
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	//TFile* fin_FONLL=new TFile("../ResultsBplus/dSigmadpt_FONLL.root");
	//TFile* fin_ddFONLL=new TFile("../ResultsBplus/dSigmadpt_dataFONLL.root");

	TFile* fin_FONLL;
	fin_FONLL=new TFile("../ResultsBplus/dSigmadpt_Bplus.root");
	TFile* fin_ddFONLL;
	fin_ddFONLL=new TFile("../ResultsBplus/dSigmadpt_Bplus_ppF.root");

	if (option==1) {
		fin_FONLL=new TFile("../ResultsBplus/dSigmadpt_Bplus.root");
	fin_ddFONLL=new TFile("../ResultsBplus/dSigmadpt_Bplus_ppF.root");
	rmk="5TeV";
	}
	else if (option==2) {
		fin_FONLL=new TFile("../ResultsBplus_2760GeVpp/dSigmadpt_Bplus_2760GeVpp.root");
	fin_ddFONLL=new TFile("../ResultsBplus_2760GeVpp/dSigmadpt_Bplus_2760GeVppF.root");
		rmk="2760GeV";
		}
	else if (option==3) {
		fin_FONLL=new TFile("../ResultsBplus/dSigmadpt_Bplus.root");
	fin_ddFONLL=new TFile("../ResultsBplus/dSigmadpt_Bplus_ppFw2760.root");
		rmk="5TeV_w2760GeV";
		}

	TGraphAsymmErrors* gaeDatastat = (TGraphAsymmErrors*)fin_FONLL->Get("gSigmastat");
	/*
		 TGraphAsymmErrors* gaeDatasyst = (TGraphAsymmErrors*)fin_FONLL->Get("gSigmasyst");
		 TGraphAsymmErrors* gaeFONLL = (TGraphAsymmErrors*)fin_FONLL->Get("gaeBplusReference");
		 TGraphAsymmErrors* gaeddFONLL = (TGraphAsymmErrors*)fin_ddFONLL->Get("gaeBplusReference");
	 */
	TGraphAsymmErrors* href3 = (TGraphAsymmErrors*)fin_FONLL->Get("gSigmasyst");
	TGraphAsymmErrors* href2 = (TGraphAsymmErrors*)fin_FONLL->Get("gaeBplusReference");
	TGraphAsymmErrors* href1 = (TGraphAsymmErrors*)fin_ddFONLL->Get("gaeBplusReference");

	//ref 1 : pp+FONLL, ref 2: FONLL, ref 3 : pPb

	Double_t xbins_ref1[5] = {12.5-1.25,17.5-1.25,22.5-1.25,27.5-1.25,45.0-1.25};
	Double_t xbins_ref2[5] = {12.5+0.00,17.5+0.00,22.5+0.00,27.5+0.00,45.0+0.00};
	Double_t xbins_ref3[5] = {12.5+1.25,17.5+1.25,22.5+1.25,27.5+1.25,45.0+1.25};

	Double_t exl_ref[5] = {0.6,0.6,0.6,0.6,0.6};

	TGraphAsymmErrors* hR12 = (TGraphAsymmErrors*)href1->Clone();
	TGraphAsymmErrors* hR31 = (TGraphAsymmErrors*)href2->Clone();
	TGraphAsymmErrors* hR32 = (TGraphAsymmErrors*)href3->Clone();

	double BIN[BIN_NUM+1] = {10,15,20,25,30,60};

	for (int hh=0;hh<BIN_NUM;hh++) {
		href1->SetPoint(hh,xbins_ref1[hh],href1->GetY()[hh]*APb);
		href2->SetPoint(hh,xbins_ref2[hh],href2->GetY()[hh]*APb);
		href3->SetPoint(hh,xbins_ref3[hh],href3->GetY()[hh]);
		href1->SetPointEXlow(hh,exl_ref[hh]);
		href2->SetPointEXlow(hh,exl_ref[hh]);
		href3->SetPointEXlow(hh,exl_ref[hh]);
		href1->SetPointEXhigh(hh,exl_ref[hh]);
		href2->SetPointEXhigh(hh,exl_ref[hh]);
		href3->SetPointEXhigh(hh,exl_ref[hh]);

		hR12->SetPoint(hh,xbins_ref1[hh],(href1->GetY()[hh])/(href2->GetY()[hh]));
		hR31->SetPoint(hh,xbins_ref2[hh],(href3->GetY()[hh])/(href1->GetY()[hh])/APb);
		hR32->SetPoint(hh,xbins_ref3[hh],(href3->GetY()[hh])/(href2->GetY()[hh])/APb);

		double R12EYlow = sqrt(pow(href1->GetEYlow()[hh]/href1->GetY()[hh],2)+pow(href2->GetEYlow()[hh]/href2->GetY()[hh],2))*(href1->GetY()[hh]/href2->GetY()[hh]);
		double R12EYhigh = sqrt(pow(href1->GetEYhigh()[hh]/href1->GetY()[hh],2)+pow(href2->GetEYhigh()[hh]/href2->GetY()[hh],2))*(href1->GetY()[hh]/href2->GetY()[hh]);

		hR12->SetPointEYlow(hh,R12EYlow);
		hR31->SetPointEYlow(hh,(href1->GetEYhigh()[hh]/(href1->GetY()[hh]+href1->GetEYhigh()[hh]))*(href3->GetY()[hh]/href1->GetY()[hh])/APb);
		hR32->SetPointEYlow(hh,(href2->GetEYhigh()[hh]/(href2->GetY()[hh]+href2->GetEYhigh()[hh]))*(href3->GetY()[hh]/href2->GetY()[hh])/APb);

		hR12->SetPointEYhigh(hh,R12EYhigh);
		hR31->SetPointEYhigh(hh,(href1->GetEYlow()[hh]/(href1->GetY()[hh]-href1->GetEYlow()[hh]))*(href3->GetY()[hh]/href1->GetY()[hh])/APb);
		hR32->SetPointEYhigh(hh,(href2->GetEYlow()[hh]/(href2->GetY()[hh]-href2->GetEYlow()[hh]))*(href3->GetY()[hh]/href2->GetY()[hh])/APb);

		hR12->SetPointEXlow(hh,exl_ref[hh]);
		hR31->SetPointEXlow(hh,exl_ref[hh]);
		hR32->SetPointEXlow(hh,exl_ref[hh]);
		hR12->SetPointEXhigh(hh,exl_ref[hh]);
		hR31->SetPointEXhigh(hh,exl_ref[hh]);
		hR32->SetPointEXhigh(hh,exl_ref[hh]);

	}
	TCanvas* canvasRref = new TCanvas("canvasRref","",500,500);
	//canvasRref->Range(-1.989924,-0.2917772,25.49622,2.212202);
	canvasRref->SetFillColor(0);
	canvasRref->SetBorderMode(0);
	canvasRref->SetBorderSize(2);
	canvasRref->SetLeftMargin(0.16);
	canvasRref->SetRightMargin(0.02);
	//###canvasRref->SetTopMargin(0.01474576);//0.08474576
	canvasRref->SetTopMargin(0.08);//0.08474576
	canvasRref->SetBottomMargin(0.15);
	canvasRref->SetFrameBorderMode(0);
	canvasRref->SetFrameBorderMode(0);
	canvasRref->SetLogy();

	canvasRref->cd();
	canvasRref->SetLogy(1);
	//TH2F* hemptyCD=new TH2F("hemptyCD","", 10, 5., 65., 5., 0.0001, 50.0);    
	TH2F* hemptyCD;
	
	if (option==2) hemptyCD=new TH2F("hemptyCD","", 10, 5., 65., 5., 0.0001, 50.0);
		else hemptyCD=new TH2F("hemptyCD","", 10, 5., 65., 5., 0.1, 1000.0);    


	hemptyCD->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	//if(particle=="Bplus") hemptyCD->GetYaxis()->SetTitle("d#sigma / dp_{T} (B^{+}) (pb GeV^{-1}c)");
	//if(particle=="Bzero") hemptyCD->GetYaxis()->SetTitle("d#sigma / dp_{T} (B^{0}) (pb GeV^{-1}c)");
	//if(particle=="Bs") hemptyCD->GetYaxis()->SetTitle("d#sigma / dp_{T} (B_{s}) (pb GeV^{-1}c)");
	hemptyCD->GetXaxis()->CenterTitle();
	hemptyCD->GetYaxis()->CenterTitle();
	hemptyCD->GetYaxis()->SetTitle("d#sigma / dp_{T}( #mub GeV^{-1}c)");
	hemptyCD->GetXaxis()->SetTitleOffset(1.0);//###1.0
	hemptyCD->GetYaxis()->SetTitleOffset(1.0);//###1.3
	hemptyCD->GetXaxis()->SetTitleSize(0.070);//###0.055
	hemptyCD->GetYaxis()->SetTitleSize(0.070);//###0.055
	hemptyCD->GetXaxis()->SetTitleFont(42);
	hemptyCD->GetYaxis()->SetTitleFont(42);
	hemptyCD->GetXaxis()->SetLabelFont(42);
	hemptyCD->GetYaxis()->SetLabelFont(42);
	hemptyCD->GetXaxis()->SetLabelSize(0.060);//###0.055
	hemptyCD->GetYaxis()->SetLabelSize(0.060);//###0.055
	hemptyCD->SetMaximum(2);
	hemptyCD->SetMinimum(0.);
	hemptyCD->Draw();

	href1->SetFillColor(kGreen);
	href2->SetFillColor(5);
	href3->SetFillColor(0);
	href1->SetLineColor(kRed);
	href2->SetLineColor(kBlue);
	href3->SetLineColor(1);
	href3->SetFillStyle(0);
	href1->SetLineWidth(2);
	href2->SetLineWidth(2);
	href3->SetLineWidth(2);
	href1->SetMarkerStyle(24);
	href2->SetMarkerStyle(25);
	href3->SetMarkerStyle(26);
	href1->SetMarkerSize(1);
	href2->SetMarkerSize(1);
	href3->SetMarkerSize(1);

	href1->Draw("same2");
	href2->Draw("same2");
	href3->Draw("same2");
	href1->Draw("samep");
	href2->Draw("samep");
	href3->Draw("samep");

	TLegend *legendSigmaC=new TLegend(0.45,0.65,0.85,0.85,"");
	legendSigmaC->SetBorderSize(0);
	legendSigmaC->SetLineColor(0);
	legendSigmaC->SetFillColor(0);
	legendSigmaC->SetFillStyle(1001);
	legendSigmaC->SetTextFont(42);
	legendSigmaC->SetTextSize(0.055);//###0.045

	legendSigmaC->AddEntry(href1,"pp+FONLL","pf");
	legendSigmaC->AddEntry(href2,"FONLL","pf");
	legendSigmaC->AddEntry(href3,"pPb","pf");
	legendSigmaC->Draw("");
	//canvasRref->SaveAs("../ResultsBplus/hrefcomp_Comp5TeV.pdf");
	canvasRref->SaveAs(Form("../ResultsBplus/hrefcomp_Comp%s.pdf",rmk.c_str()));


	legendSigmaC->Clear();

	canvasRref->SetLogy(0);

	TH2F* hemptyCR;
	if (option==1 || option==2) hemptyCR=new TH2F("hemptyCR","", 10, 5., 65., 15., 0.5, 3.0);
	else if (option==3) hemptyCR=new TH2F("hemptyCR","", 10, 5., 65., 15., 0.0, 5.0);    
	//TH2F* hemptyCR=new TH2F("hemptyCR","", 10, 5., 65., 15., 0.5, 2.0);    


	hemptyCR->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	//if(particle=="Bplus") hemptyCR->GetYaxis()->SetTitle("d#sigma / dp_{T} (B^{+}) (pb GeV^{-1}c)");
	//if(particle=="Bzero") hemptyCR->GetYaxis()->SetTitle("d#sigma / dp_{T} (B^{0}) (pb GeV^{-1}c)");
	//if(particle=="Bs") hemptyCR->GetYaxis()->SetTitle("d#sigma / dp_{T} (B_{s}) (pb GeV^{-1}c)");
	hemptyCR->GetXaxis()->CenterTitle();
	hemptyCR->GetYaxis()->CenterTitle();
	hemptyCR->GetYaxis()->SetTitle("Ratios");
	hemptyCR->GetXaxis()->SetTitleOffset(1.0);//###1.0
	hemptyCR->GetYaxis()->SetTitleOffset(1.0);//###1.3
	hemptyCR->GetXaxis()->SetTitleSize(0.070);//###0.055
	hemptyCR->GetYaxis()->SetTitleSize(0.070);//###0.055
	hemptyCR->GetXaxis()->SetTitleFont(42);
	hemptyCR->GetYaxis()->SetTitleFont(42);
	hemptyCR->GetXaxis()->SetLabelFont(42);
	hemptyCR->GetYaxis()->SetLabelFont(42);
	hemptyCR->GetXaxis()->SetLabelSize(0.060);//###0.055
	hemptyCR->GetYaxis()->SetLabelSize(0.060);//###0.055
	hemptyCR->SetMaximum(2);
	hemptyCR->SetMinimum(0.);
	hemptyCR->Draw();

	hR12->SetFillColor(kGreen);
	hR31->SetFillColor(kPink+1);
	hR32->SetFillColor(5);
	hR32->SetFillStyle(1001);
	hR12->SetLineColor(kRed);
	hR31->SetLineColor(kBlue);
	hR32->SetLineColor(1);
	hR12->SetLineWidth(2);
	hR31->SetLineWidth(2);
	hR32->SetLineWidth(2);
	hR12->SetMarkerStyle(20);
	hR31->SetMarkerStyle(21);
	hR32->SetMarkerStyle(22);
	hR12->SetMarkerSize(1);
	hR31->SetMarkerSize(1);
	hR32->SetMarkerSize(1);



	TLine *l2 = new TLine(5.1,1.0,64.9,1.0);
	l2->SetLineStyle(2);
	l2->SetLineColor(kBlue);
	l2->SetLineStyle(2);  
	l2->SetLineWidth(2);

	hR12->Draw("samee2");
	hR31->Draw("samee2");
	hR32->Draw("samee2");
	hR12->Draw("samep");
	hR31->Draw("samep");
	hR32->Draw("samep");
	l2->Draw();
	legendSigmaC=new TLegend(0.40,0.65,0.80,0.90,"");
	legendSigmaC->SetBorderSize(0);
	legendSigmaC->SetLineColor(0);
	legendSigmaC->SetFillColor(0);
	legendSigmaC->SetFillStyle(1001);
	legendSigmaC->SetTextFont(42);
	legendSigmaC->SetTextSize(0.055);//###0.045
	legendSigmaC->AddEntry(hR12,"data+FONLL/FONLL","pf");
	legendSigmaC->AddEntry(hR31,"R_{pPb}^{data+FONLL}","pf");
	legendSigmaC->AddEntry(hR32,"R_{pPb}^{FONLL}","pf");
	legendSigmaC->Draw("");

	//canvasRref->SaveAs("../ResultsBplus/hRcomp_Comp5TeV.pdf");
	canvasRref->SaveAs(Form("../ResultsBplus/hRcomp_Comp%s.pdf",rmk.c_str()));


	for (int i=0;i<6;i++){
		std::cout << i << " : " << hR12->GetY()[i] << std::endl;
	}
}


