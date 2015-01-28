//### Code name :
//### Made by Hyunchul Kim at Jan. 16th. 2015
//### Purpose : Get the fitting function 

#include <iostream>
//#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraphAsymmErrors.h"

#define NUM 5	    // number of bins for CMS pp data
#define NUM_ATL 8   // number of bins for ATLAS pp data
//#define APb 208
#define APb 1

void iterpoints(int i, double oribin[], TF1* f, double stmpt[], double iterorder)
{
	//iterpoints(i, ptbin, fitft_opt2, stmpt, 2)
	double integ1,integ2,gapinteg,minpt;
	double mininteg=1000000000000000.0;
	double initval = pow(10,iterorder*(-1));
	//std::cout << "#### Inital value : " << initval << std::endl;
	if (iterorder==0) {
		for (double j=oribin[i]+1.0;j<oribin[i+1];j+=1.0){
			integ1 = f->Integral(oribin[i],j);
			integ2 = f->Integral(j,oribin[i+1]);
			gapinteg = fabs(integ1-integ2);
			//###std:: cout << j << " : " << integ1 << " , " << integ2 << std::endl;
			if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
		}
	}
	else {
		for (double j=stmpt[i]-initval*9;j<stmpt[i]+initval*10;j+=initval){
			integ1 = f->Integral(oribin[i],j);
			integ2 = f->Integral(j,oribin[i+1]);
			gapinteg = fabs(integ1-integ2);
			//###std:: cout << j << " : " << "----- " << integ1 << " , " << integ2 << std::endl;
			if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
		}
	}
	stmpt[i]=minpt;
	std::cout << i << " : " << std::string(iterorder*3+1, '-') << " gap : " << mininteg << ", minpt : " << minpt << std::endl;
}

void iterpoints_COM(int i, double lowv, double highv, TF1* f, double stmpt[], double iterorder)
{
	// iterpoints_COM(i, ptbin_COM_l, ptbin_COM_h, fitft_opt7, stminpt_COM2, 0);
	double integ1,integ2,gapinteg,minpt;
	double mininteg=1000000000000000.0;
	double initval = pow(10,iterorder*(-1));
	//std::cout << "#### Inital value : " << initval << std::endl;
	if (iterorder==0) {
		for (double j=lowv+1.0;j<highv;j+=1.0){
			integ1 = f->Integral(lowv,j);
			integ2 = f->Integral(j,highv);
			gapinteg = fabs(integ1-integ2);
			//###std:: cout << j << " : " << integ1 << " , " << integ2 << std::endl;
			if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
		}
	}
	else {
		for (double j=stmpt[i]-initval*9;j<stmpt[i]+initval*10;j+=initval){
			integ1 = f->Integral(lowv,j);
			integ2 = f->Integral(j,highv);
			gapinteg = fabs(integ1-integ2);
			//###std:: cout << j << " : " << "----- " << integ1 << " , " << integ2 << std::endl;
			if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
		}
	}
	stmpt[i]=minpt;
	std::cout << i << " : " << std::string(iterorder*3+1, '-') << " gap : " << mininteg << ", minpt : " << minpt << std::endl;
}

void ScaleCor_pt()
{

	//gROOT->SetStyle("plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	TFile* fout = new TFile("ScaleCor_reweighted7TeVdata.root","RECREATE");
	TH1F* hparafit = new TH1F("hparafit","",3,0.0,3.0);

	// CMS 7TeV, pp->(B+)X, unit : micro barn
	double ptbin[NUM+1]={5.,10.,13.,17.,24.,30.};
	double sigmapt[NUM]={4.07,1.47,0.412,0.181,0.042};
	double sigmapt_sta[NUM]={0.47,0.13,0.041,0.015,0.007};
	double sigmapt_sys[NUM]={0.31,0.09,0.026,0.012,0.004};

	// ATLAS 7TeV, dsigma/dp_T (B+) * BR(B+->J/psiK+) * BR (J/Psi->mu+mu-), unit : pico barn
	double ptbin_ATL[NUM_ATL+1]={9.,13.,16.,20.,25.,35.,50.,70.,120.};
	double sigmapt_ATL[NUM_ATL]={103.4,36.03,15.33,6.056,1.814,0.3477,0.06244,0.006099};
	double sigmapt_sta_ATL[NUM_ATL]={3.7,0.80,0.25,0.093,0.027,0.0084,0.00293,0.000561};
	double sigmapt_sys_ATL[NUM_ATL]={7.6,2.32,0.98,0.376,0.115,0.0280,0.00526,0.000666};

	// Branching fraction for ATLAS data 
	double BRs;
	BRs = (1.027*0.001)*(5.961*0.01);	//BR(B+->J/psiK+) * BR (J/Psi->mu+mu-)

	// call histograms for CMS and ATLAS data
	TH1D* hsigmapt_noerr = new TH1D("hsigmapt_noerr","",NUM,ptbin);
	TH1D* hsigmapt_staerr = new TH1D("hsigmapt_staerr","",NUM,ptbin);
	TH1D* hsigmapt_syserr = new TH1D("hsigmapt_syserr","",NUM,ptbin);
	TH1D* hsigmapt_noerr_ATL = new TH1D("hsigmapt_noerr_ATL","",NUM_ATL,ptbin_ATL);
	TH1D* hsigmapt_staerr_ATL = new TH1D("hsigmapt_staerr_ATL","",NUM_ATL,ptbin_ATL);
	TH1D* hsigmapt_syserr_ATL = new TH1D("hsigmapt_syserr_ATL","",NUM_ATL,ptbin_ATL);

	// store the corrected data in the histogram for combining CMS and ATLAS together
	for(int i=0;i<NUM;i++){
		hsigmapt_noerr->SetBinContent(i+1,sigmapt[i]*APb);
		hsigmapt_staerr->SetBinContent(i+1,sigmapt[i]*APb);
		hsigmapt_syserr->SetBinContent(i+1,sigmapt[i]*APb);
		hsigmapt_staerr->SetBinError(i+1,sigmapt_sta[i]*APb);
		hsigmapt_syserr->SetBinError(i+1,sigmapt_sys[i]*APb);
	}
	for(int i=0;i<NUM_ATL;i++){
		hsigmapt_noerr_ATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001);
		hsigmapt_staerr_ATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001);
		hsigmapt_syserr_ATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001);
		hsigmapt_staerr_ATL->SetBinError(i+1,sigmapt_sta_ATL[i]*APb/BRs*0.000001);
		hsigmapt_syserr_ATL->SetBinError(i+1,sigmapt_sys_ATL[i]*APb/BRs*0.000001);
	}

	// call for fitting function and set the related option
	TF1* fitft_opt1 = new TF1("fitft_opt1","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
	fitft_opt1->SetLineColor(kRed);
	TF1* fitft_opt2 = new TF1("fitft_opt2","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
	fitft_opt2->SetLineColor(kAzure-6);
	TF1* fitft_opt3 = new TF1("fitft_opt3","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
	fitft_opt3->SetLineColor(kOrange);
	fitft_opt3->SetLineStyle(2);
	TF1* fitft_opt4 = new TF1("fitft_opt4","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
	fitft_opt4->SetLineColor(kGreen);
	fitft_opt4->SetLineStyle(2);

	// parameter limits
	fitft_opt1->SetParLimits(1,-0.025,-0.015);
	fitft_opt2->SetParLimits(1,-0.025,-0.015);

	std::cout << "*** Fitting with 7TeV pp CMS " << std::endl;
	std::cout << "(1) Fit with CMS data" << std::endl;
	hsigmapt_staerr->Fit("fitft_opt1","M","",0.0,120.0);
	std::cout << "(2) Fit with ATLAS data" << std::endl;
	hsigmapt_staerr_ATL->Fit("fitft_opt2","M","",0.0,120.0);

	fitft_opt4->SetParameter(0,fitft_opt1->GetParameter(0));
	fitft_opt4->SetParameter(1,fitft_opt1->GetParameter(1));
	fitft_opt4->SetParameter(2,fitft_opt1->GetParameter(2));

	fitft_opt3->SetParameter(0,fitft_opt2->GetParameter(0));
	fitft_opt3->SetParameter(1,fitft_opt2->GetParameter(1));
	fitft_opt3->SetParameter(2,fitft_opt2->GetParameter(2));

	std::cout << "(3) Fit with CMS data based on (2)" << std::endl;
	hsigmapt_staerr->Fit("fitft_opt3","M","",0.0,120.0);
	std::cout << "(4) Fit with ATLAS data based on (1)" << std::endl;
	hsigmapt_staerr_ATL->Fit("fitft_opt4","M","",0.0,120.0);

	std::cout << std::endl;

	double stminpt[NUM];
	double sty[NUM],stexl[NUM],stexh[NUM];
	double steyl_sta[NUM],steyh_sta[NUM];
	double steyl_sys[NUM],steyh_sys[NUM];

	std::cout << std::string(50,'#') << std::endl; 
	for (int i=0;i<NUM;i++){
		std::cout << i << " : " << ", central point : " << (ptbin[i+1]+ptbin[i])/2 << std::endl;
		iterpoints(i, ptbin, fitft_opt4, stminpt, 0);
		iterpoints(i, ptbin, fitft_opt4, stminpt, 1);
		iterpoints(i, ptbin, fitft_opt4, stminpt, 2);
		sty[i]=hsigmapt_noerr->GetBinContent(i+1);
		stexl[i]=stminpt[i]-ptbin[i];
		stexh[i]=ptbin[i+1]-stminpt[i];
		steyl_sta[i]=hsigmapt_staerr->GetBinError(i+1);
		steyh_sta[i]=hsigmapt_staerr->GetBinError(i+1);
		steyl_sys[i]=hsigmapt_syserr->GetBinError(i+1);
		steyh_sys[i]=hsigmapt_syserr->GetBinError(i+1);
	}
	std::cout << std::string(50,'#') << std::endl; 

	// ALTAS iteration
	double stminpt_ATL[NUM_ATL];
	double sty_ATL[NUM_ATL],stexl_ATL[NUM_ATL],stexh_ATL[NUM_ATL];
	double steyl_sta_ATL[NUM_ATL],steyh_sta_ATL[NUM_ATL];
	double steyl_sys_ATL[NUM_ATL],steyh_sys_ATL[NUM_ATL];

	std::cout << std::endl;
	std::cout << std::string(50,'#') << std::endl; 
	for (int i=0;i<NUM_ATL;i++){
		std::cout << i << " : " << ", central point : " << (ptbin_ATL[i+1]+ptbin_ATL[i])/2 << std::endl;
		iterpoints(i, ptbin_ATL, fitft_opt4, stminpt_ATL, 0);
		iterpoints(i, ptbin_ATL, fitft_opt4, stminpt_ATL, 1);
		iterpoints(i, ptbin_ATL, fitft_opt4, stminpt_ATL, 2);
		sty_ATL[i]=hsigmapt_noerr_ATL->GetBinContent(i+1);
		stexl_ATL[i]=stminpt_ATL[i]-ptbin_ATL[i];
		stexh_ATL[i]=ptbin_ATL[i+1]-stminpt_ATL[i];
		steyl_sta_ATL[i]=hsigmapt_staerr_ATL->GetBinError(i+1);
		steyh_sta_ATL[i]=hsigmapt_staerr_ATL->GetBinError(i+1);
		steyl_sys_ATL[i]=hsigmapt_syserr_ATL->GetBinError(i+1);
		steyh_sys_ATL[i]=hsigmapt_syserr_ATL->GetBinError(i+1);
	}
	std::cout << std::string(50,'#') << std::endl; 

	// CMS, ATLAS merged points
	int NUM_COM=NUM+NUM_ATL;
	double stminpt_COM[NUM_COM];
	double sty_COM[NUM_COM],stexl_COM[NUM_COM],stexh_COM[NUM_COM];
	double steyl_sta_COM[NUM_COM],steyh_sta_COM[NUM_COM];
	double steyl_sys_COM[NUM_COM],steyh_sys_COM[NUM_COM];

	for (int i=0;i<NUM;i++) {
		stminpt_COM[i]=stminpt[i];
		sty_COM[i]=sty[i];
		stexl_COM[i]=stexl[i];
		stexh_COM[i]=stexh[i];
		steyl_sta_COM[i]=steyl_sta[i];
		steyh_sta_COM[i]=steyh_sta[i];
		steyl_sys_COM[i]=steyl_sys[i];
		steyh_sys_COM[i]=steyh_sys[i];
	}
	for (int i=0;i<NUM_ATL;i++) {
		stminpt_COM[NUM+i]=stminpt_ATL[i];
		sty_COM[NUM+i]=sty_ATL[i];
		stexl_COM[NUM+i]=stexl_ATL[i];
		stexh_COM[NUM+i]=stexh_ATL[i];
		steyl_sta_COM[NUM+i]=steyl_sta_ATL[i];
		steyh_sta_COM[NUM+i]=steyh_sta_ATL[i];
		steyl_sys_COM[NUM+i]=steyl_sys_ATL[i];
		steyh_sys_COM[NUM+i]=steyh_sys_ATL[i];
	}


	/////////////////////////////////////////////
	TGraphAsymmErrors* gaeCMSsta = new TGraphAsymmErrors(NUM,stminpt,sty,stexl,stexh,steyl_sta,steyh_sta);
	TGraphAsymmErrors* gaeCMSsys = new TGraphAsymmErrors(NUM,stminpt,sty,stexl,stexh,steyl_sys,steyh_sys);
	TGraphAsymmErrors* gaeATLsta = new TGraphAsymmErrors(NUM_ATL,stminpt_ATL,sty_ATL,stexl_ATL,stexh_ATL,steyl_sta_ATL,steyh_sta_ATL);
	TGraphAsymmErrors* gaeATLsys = new TGraphAsymmErrors(NUM_ATL,stminpt_ATL,sty_ATL,stexl_ATL,stexh_ATL,steyl_sys_ATL,steyh_sys_ATL);
	///////////////////
	TGraphAsymmErrors* gaeCOMsta = new TGraphAsymmErrors(NUM_COM,stminpt_COM,sty_COM,stexl_COM,stexh_COM,steyl_sta_COM,steyh_sta_COM);
	TGraphAsymmErrors* gaeCOMsys = new TGraphAsymmErrors(NUM_COM,stminpt_COM,sty_COM,stexl_COM,stexh_COM,steyl_sys_COM,steyh_sys_COM);
	gaeCMSsta->SetName("gaeCMSsta");
	gaeCMSsys->SetName("gaeCMSsys");
	gaeATLsta->SetName("gaeATLsta");
	gaeATLsys->SetName("gaeATLsys");
	gaeCOMsta->SetName("gaeCOMsta");
	gaeCOMsys->SetName("gaeCOMsys");





	////////////////
	TF1* fitft_opt5 = new TF1("fitft_opt5","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
	TF1* fitft_opt6 = new TF1("fitft_opt6","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
	TF1* fitft_opt7 = new TF1("fitft_opt7","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
	//###TF1* fitft_opt7 = new TF1("fitft_opt7","pow(10,[0]+[1]*x+[2]*x*x)",0.0,120.0);




	/*
		 fitft_opt5->SetParameter(0,fitft_opt1->GetParameter(0));
		 fitft_opt5->SetParameter(1,fitft_opt1->GetParameter(1));
		 fitft_opt5->SetParameter(2,fitft_opt1->GetParameter(2));
	//fitft_opt5->SetParameter(3,fitft_opt3->GetParameter(3));
	 */
	fitft_opt5->SetParLimits(0,5.0,7.0);
	fitft_opt5->SetParLimits(1,-0.025,-0.015);



	fitft_opt6->SetParameter(0,fitft_opt4->GetParameter(0));
	fitft_opt6->SetParameter(1,fitft_opt4->GetParameter(1));
	fitft_opt6->SetParameter(2,fitft_opt4->GetParameter(2));
	//fitft_opt6->SetParameter(3,fitft_opt4->GetParameter(3));



	fitft_opt5->SetLineColor(kOrange+7);
	fitft_opt5->SetLineStyle(1);

	fitft_opt6->SetLineColor(kGreen+4);
	fitft_opt6->SetLineStyle(1);

	fitft_opt7->SetLineColor(kViolet+5);
	fitft_opt7->SetLineStyle(1);


	std::cout << "(5) Fit with CMS weighted data, based on (6)" << std::endl;


	gaeCMSsta->Fit("fitft_opt5","M","",5.0,30.0);


	//std::cout << "(5) Fit with CMS weighted data based on (1)" << std::endl;
	std::cout << "(6) Fit with ATLAS weighted data based on (5)" << std::endl;
	gaeATLsta->Fit("fitft_opt6","M","",0.0,120.0);

	fitft_opt7->SetParameter(0,fitft_opt5->GetParameter(0));
	fitft_opt7->SetParameter(1,fitft_opt5->GetParameter(1));
	fitft_opt7->SetParameter(2,fitft_opt5->GetParameter(2));

	//fitft_opt7->SetParameter(0,6.30625e+00);
	//fitft_opt7->SetParameter(1,-2.41914e-02);

	fitft_opt7->SetParLimits(0,5.0,7.0);
	fitft_opt7->SetParLimits(1,-0.100,-0.005);



	std::cout << "(7) Fit with CMS+ATLAS weighted data" << std::endl;
	gaeCOMsta->Fit("fitft_opt7","M","",0.0,100.0);

	hparafit->SetBinContent(1,fitft_opt7->GetParameter(0));
	hparafit->SetBinContent(2,fitft_opt7->GetParameter(1));
	hparafit->SetBinContent(3,fitft_opt7->GetParameter(2));

	double stminpt_COM2[NUM_COM];
	double sty_COM2[NUM_COM],stexl_COM2[NUM_COM],stexh_COM2[NUM_COM];
	double steyl_sta_COM2[NUM_COM],steyh_sta_COM2[NUM_COM];
	double steyl_sys_COM2[NUM_COM],steyh_sys_COM2[NUM_COM];

	std::cout << std::endl;
	std::cout << std::string(50,'#') << std::endl; 
	for (int i=0;i<NUM_COM;i++){
		std::cout << i << " : merged, weighted point : " << stminpt_COM[i] << std::endl;
		double ptbin_COM_l=stminpt_COM[i]-stexl_COM[i];//low edge
		double ptbin_COM_h=stminpt_COM[i]+stexh_COM[i];//high edge
		iterpoints_COM(i, ptbin_COM_l, ptbin_COM_h, fitft_opt7, stminpt_COM2, 0);
		iterpoints_COM(i, ptbin_COM_l, ptbin_COM_h, fitft_opt7, stminpt_COM2, 1);
		iterpoints_COM(i, ptbin_COM_l, ptbin_COM_h, fitft_opt7, stminpt_COM2, 2);
		sty_COM2[i]=sty_COM[i];
		stexl_COM2[i]=stminpt_COM2[i]-stminpt_COM[i];
		stexh_COM2[i]=stminpt_COM[i]-stminpt_COM2[i];
		steyl_sta_COM2[i]=steyl_sta_COM[i];
		steyh_sta_COM2[i]=steyh_sta_COM[i];
		steyl_sys_COM2[i]=steyl_sys_COM[i];
		steyh_sys_COM2[i]=steyh_sys_COM[i];
	}
	std::cout << std::string(50,'#') << std::endl; 

	TGraphAsymmErrors* gaeCOMsta2 = new TGraphAsymmErrors(NUM_COM,stminpt_COM2,sty_COM2,stexl_COM2,stexh_COM2,steyl_sta_COM2,steyh_sta_COM2);
	TGraphAsymmErrors* gaeCOMsys2 = new TGraphAsymmErrors(NUM_COM,stminpt_COM2,sty_COM2,stexl_COM2,stexh_COM2,steyl_sys_COM2,steyh_sys_COM2);
	gaeCOMsta->SetName("gaeCOMsta2");
	gaeCOMsys->SetName("gaeCOMsys2");






	/*
		 TF1* fitft_opt3 = new TF1("fitft_opt3","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
		 fitft_opt3->SetLineColor(kOrange);
		 fitft_opt3->SetLineStyle(2);
		 TF1* fitft_opt4 = new TF1("fitft_opt4","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
		 std::cout << "(3) Fit with CMS data based on (2)" << std::endl;
	//hsigmapt_staerr->Fit("fitft_opt3","L","",0.0,120.0);
	hsigmapt_staerr->Fit("fitft_opt3","M","",0.0,60.0);

	std::cout << "(4) Fit with ATLAS data based on (1)" << std::endl;
	hsigmapt_staerr_ATL->Fit("fitft_opt4","M","",0.0,120.0);
	 */



	TCanvas *cSigma = new TCanvas("cSigma","",500,500);

	cSigma->SetFillColor(0);
	cSigma->SetBorderMode(0);
	cSigma->SetBorderSize(2);
	cSigma->SetFrameBorderMode(0);
	cSigma->SetLeftMargin(0.16);
	cSigma->SetRightMargin(0.02);
	cSigma->SetTopMargin(0.08);
	cSigma->SetBottomMargin(0.15);
	cSigma->SetLogy(1);

	//TH2D* hempty = new TH2D("hempty",";B^{+} p_{T} (GeV/c);d#sigma/dp_{T} A",13,0.0,130.0,7,0.01,10000.0);
	TH2D* hempty = new TH2D("hempty",";B^{+} p_{T} (GeV/c);d#sigma/dp_{T} A",13,0.0,130.0,10,0.00001,10000.0);

	//hempty->SetAxisRange(7.5, 42.5, "X");
	//hempty->SetAxisRange(0.1, 1000, "Y");
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();
	hempty->GetXaxis()->SetTitleOffset(1.0);//###1.0
	hempty->GetYaxis()->SetTitleOffset(1.0);//###1.3
	hempty->GetXaxis()->SetTitleSize(0.070);//###0.055
	hempty->GetYaxis()->SetTitleSize(0.070);//###0.055
	hempty->GetXaxis()->SetTitleFont(42);
	hempty->GetYaxis()->SetTitleFont(42);
	hempty->GetXaxis()->SetLabelFont(42);
	hempty->GetYaxis()->SetLabelFont(42);
	hempty->GetXaxis()->SetLabelSize(0.060);//###0.055
	hempty->GetYaxis()->SetLabelSize(0.060);//###0.055    

	hsigmapt_staerr->SetMarkerColor(kRed);
	hsigmapt_staerr->SetMarkerStyle(24);
	hsigmapt_staerr->SetMarkerSize(1);
	hsigmapt_staerr->SetLineWidth(2);
	hsigmapt_staerr->SetLineColor(kRed);

	hsigmapt_staerr_ATL->SetMarkerColor(kTeal+3);
	hsigmapt_staerr_ATL->SetMarkerStyle(25);
	hsigmapt_staerr_ATL->SetMarkerSize(1);
	hsigmapt_staerr_ATL->SetLineWidth(2);
	hsigmapt_staerr_ATL->SetLineColor(kTeal+3);

	hsigmapt_syserr->SetLineColor(kRed);
	hsigmapt_syserr->SetFillColor(5);
	hsigmapt_syserr->SetMarkerSize(0);

	hsigmapt_syserr_ATL->SetLineColor(kTeal+3);
	hsigmapt_syserr_ATL->SetFillColor(kTeal);//kGreen
	hsigmapt_syserr_ATL->SetMarkerSize(0);

	hempty->Draw("");
	hsigmapt_syserr->Draw("samee2");
	hsigmapt_syserr_ATL->Draw("samee2");
	hsigmapt_staerr->Draw("samee");
	hsigmapt_staerr_ATL->Draw("samee");
	fitft_opt1->Draw("same");
	fitft_opt2->Draw("same");
	fitft_opt3->Draw("same");
	fitft_opt4->Draw("same");


	TLegend* leg = new TLegend(0.30,0.65,0.60,0.90,"");

	leg->SetBorderSize(0);
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	leg->SetTextFont(42);
	leg->SetTextSize(0.040);

	hsigmapt_staerr->SetFillColor(5);
	hsigmapt_staerr_ATL->SetFillColor(kTeal);//kGreen

	TLegendEntry* ent_ppCMS7=leg->AddEntry(hsigmapt_staerr,"CMS 7 TeV, |y|<2.4","lpf");

	ent_ppCMS7->SetTextFont(42);
	ent_ppCMS7->SetFillStyle(1001);
	ent_ppCMS7->SetFillColor(5);
	ent_ppCMS7->SetMarkerColor(kRed);
	ent_ppCMS7->SetMarkerStyle(24);
	ent_ppCMS7->SetMarkerSize(15);
	ent_ppCMS7->SetLineWidth(2);
	ent_ppCMS7->SetLineColor(kRed);

	TLegendEntry* ent_ppATL7=leg->AddEntry(hsigmapt_staerr_ATL,"ATLAS 7 TeV, |y|<2.25","lpf");

	ent_ppATL7->SetTextFont(42);
	ent_ppATL7->SetFillColor(kTeal);
	ent_ppATL7->SetMarkerColor(kTeal+3);
	ent_ppATL7->SetMarkerStyle(25);
	ent_ppATL7->SetMarkerSize(15);
	ent_ppATL7->SetLineWidth(2);
	ent_ppATL7->SetLineColor(kTeal+3);

	leg->AddEntry(fitft_opt1,"(1) Fit with CMS data","l");
	leg->AddEntry(fitft_opt2,"(2) Fit with ATLAS data","l");
	leg->AddEntry(fitft_opt3,"(3) Fit with CMS data based on (2)","l");
	leg->AddEntry(fitft_opt4,"(4) Fit with ATLAS data based on (1)","l");


	leg->Draw("");

	cSigma->SaveAs("Plots/cSigma_CompFit.pdf");
	cSigma->Clear();
	//cSigma->SaveAs("cSigma_AddATLAS_wParLimit.pdf");//With ParLimit
	cSigma->SetFillColor(0);
	cSigma->SetBorderMode(0);
	cSigma->SetBorderSize(2);
	cSigma->SetFrameBorderMode(0);
	cSigma->SetLeftMargin(0.16);
	cSigma->SetRightMargin(0.02);
	cSigma->SetTopMargin(0.08);
	cSigma->SetBottomMargin(0.15);
	cSigma->SetLogy(1);

	//hempty->SetAxisRange(7.5, 42.5, "X");
	//hempty->SetAxisRange(0.1, 1000, "Y");
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();
	hempty->GetXaxis()->SetTitleOffset(1.0);//###1.0
	hempty->GetYaxis()->SetTitleOffset(1.0);//###1.3
	hempty->GetXaxis()->SetTitleSize(0.070);//###0.055
	hempty->GetYaxis()->SetTitleSize(0.070);//###0.055
	hempty->GetXaxis()->SetTitleFont(42);
	hempty->GetYaxis()->SetTitleFont(42);
	hempty->GetXaxis()->SetLabelFont(42);
	hempty->GetYaxis()->SetLabelFont(42);
	hempty->GetXaxis()->SetLabelSize(0.060);//###0.055
	hempty->GetYaxis()->SetLabelSize(0.060);//###0.055    

	hsigmapt_staerr->SetMarkerColor(kRed);
	hsigmapt_staerr->SetMarkerStyle(24);
	hsigmapt_staerr->SetMarkerSize(1);
	hsigmapt_staerr->SetLineWidth(2);
	hsigmapt_staerr->SetLineColor(kRed);

	hsigmapt_staerr_ATL->SetMarkerColor(kTeal+3);
	hsigmapt_staerr_ATL->SetMarkerStyle(25);
	hsigmapt_staerr_ATL->SetMarkerSize(1);
	hsigmapt_staerr_ATL->SetLineWidth(2);
	hsigmapt_staerr_ATL->SetLineColor(kTeal+3);

	hsigmapt_syserr->SetLineColor(kRed);
	hsigmapt_syserr->SetFillColor(5);
	hsigmapt_syserr->SetMarkerSize(0);

	hsigmapt_syserr_ATL->SetLineColor(kTeal+3);
	hsigmapt_syserr_ATL->SetFillColor(kTeal);//kGreen
	hsigmapt_syserr_ATL->SetMarkerSize(0);

	gaeCMSsta->SetLineColor(kRed);
	gaeCMSsta->SetMarkerColor(kRed);
	gaeCMSsta->SetMarkerStyle(20);
	gaeCMSsta->SetMarkerSize(1);
	gaeCMSsta->SetFillColor(5);

	gaeCMSsys->SetLineColor(kRed);
	gaeCMSsys->SetFillColor(5);
	gaeCMSsys->SetMarkerSize(0);

	gaeATLsta->SetLineColor(kTeal+3);
	gaeATLsta->SetMarkerColor(kTeal+3);
	gaeATLsta->SetMarkerStyle(21);
	gaeATLsta->SetMarkerSize(1);
	gaeATLsta->SetFillColor(kTeal);

	gaeATLsys->SetLineColor(kTeal+3);
	gaeATLsys->SetFillColor(kTeal);//kGreen
	gaeATLsys->SetMarkerSize(0);

	gaeCOMsta->SetLineColor(kViolet+5);
	gaeCOMsta->SetMarkerColor(kViolet+5);
	gaeCOMsta->SetMarkerStyle(33);
	gaeCOMsta->SetMarkerSize(2);
	gaeCOMsta->SetFillColor(kViolet-4);

	gaeCOMsys->SetLineColor(kViolet+5);
	gaeCOMsys->SetFillColor(kViolet-4);
	gaeCOMsys->SetMarkerSize(0);

	gaeCOMsta2->SetLineColor(kViolet+5);
	gaeCOMsta2->SetMarkerColor(kViolet+5);
	gaeCOMsta2->SetMarkerStyle(27);
	gaeCOMsta2->SetMarkerSize(2);
	gaeCOMsta2->SetFillColor(kViolet-4);


	gaeCOMsys2->SetLineColor(kViolet+5);
	gaeCOMsys2->SetFillColor(kViolet-4);
	gaeCOMsys2->SetMarkerSize(0);


	hempty->Draw("");
	/*
		 hsigmapt_syserr->Draw("samee2");
		 hsigmapt_syserr_ATL->Draw("samee2");
		 hsigmapt_staerr->Draw("samee");
		 hsigmapt_staerr_ATL->Draw("samee");
	 */
	/*
		 fitft_opt1->Draw("same");
		 fitft_opt2->Draw("same");
	 */
	fitft_opt3->Draw("same");
	fitft_opt4->Draw("same");
	gaeCMSsta->Draw("samee2");
	gaeATLsta->Draw("samee2");
	gaeCMSsta->Draw("samep");
	gaeATLsta->Draw("samep");
	fitft_opt3->Draw("same");
	fitft_opt4->Draw("same");
	fitft_opt5->Draw("same");
	fitft_opt6->Draw("same");

	leg->Clear();
	leg = new TLegend(0.30,0.55,0.60,0.90,"");

	leg->SetBorderSize(0);
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	leg->SetTextFont(42);
	leg->SetTextSize(0.040);

	leg->AddEntry(hsigmapt_staerr,"CMS 7 TeV, |y|<2.4","lpf");
	leg->AddEntry(hsigmapt_staerr_ATL,"ATLAS 7 TeV, |y|<2.25","lpf");
	leg->AddEntry(gaeCMSsta,"Reweighted CMS 7 TeV, |y|<2.4","lpf");
	leg->AddEntry(gaeATLsta,"Reweighted ATLAS 7 TeV, |y|<2.25","lpf");
	leg->AddEntry(fitft_opt3,"(3) Fit with CMS data based on (2)","l");
	leg->AddEntry(fitft_opt4,"(4) Fit with ATLAS data based on (1)","l");
	leg->AddEntry(fitft_opt5,"(5) Fit with reweighted CMS data","l");
	leg->AddEntry(fitft_opt6,"(6) Fit with reweighted ATLAS data","l");

	leg->Draw("");


	cSigma->SaveAs("Plots/cSigma_CompReweightedFit.pdf");

	cSigma->Clear();
	//cSigma->SaveAs("cSigma_AddATLAS_wParLimit.pdf");//With ParLimit
	cSigma->SetFillColor(0);
	cSigma->SetBorderMode(0);
	cSigma->SetBorderSize(2);
	cSigma->SetFrameBorderMode(0);
	cSigma->SetLeftMargin(0.16);
	cSigma->SetRightMargin(0.02);
	cSigma->SetTopMargin(0.08);
	cSigma->SetBottomMargin(0.15);
	cSigma->SetLogy(1);

	//hempty->SetAxisRange(7.5, 42.5, "X");
	//hempty->SetAxisRange(0.1, 1000, "Y");
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();
	hempty->GetXaxis()->SetTitleOffset(1.0);//###1.0
	hempty->GetYaxis()->SetTitleOffset(1.0);//###1.3
	hempty->GetXaxis()->SetTitleSize(0.070);//###0.055
	hempty->GetYaxis()->SetTitleSize(0.070);//###0.055
	hempty->GetXaxis()->SetTitleFont(42);
	hempty->GetYaxis()->SetTitleFont(42);
	hempty->GetXaxis()->SetLabelFont(42);
	hempty->GetYaxis()->SetLabelFont(42);
	hempty->GetXaxis()->SetLabelSize(0.060);//###0.055
	hempty->GetYaxis()->SetLabelSize(0.060);//###0.055    

	hsigmapt_staerr->SetMarkerColor(kRed);
	hsigmapt_staerr->SetMarkerStyle(24);
	hsigmapt_staerr->SetMarkerSize(1);
	hsigmapt_staerr->SetLineWidth(2);
	hsigmapt_staerr->SetLineColor(kRed);

	hsigmapt_staerr_ATL->SetMarkerColor(kTeal+3);
	hsigmapt_staerr_ATL->SetMarkerStyle(25);
	hsigmapt_staerr_ATL->SetMarkerSize(1);
	hsigmapt_staerr_ATL->SetLineWidth(2);
	hsigmapt_staerr_ATL->SetLineColor(kTeal+3);

	hsigmapt_syserr->SetLineColor(kRed);
	hsigmapt_syserr->SetFillColor(5);
	hsigmapt_syserr->SetMarkerSize(0);

	hsigmapt_syserr_ATL->SetLineColor(kTeal+3);
	hsigmapt_syserr_ATL->SetFillColor(kTeal);//kGreen
	hsigmapt_syserr_ATL->SetMarkerSize(0);


	hempty->Draw("");
	/*
		 hsigmapt_syserr->Draw("samee2");
		 hsigmapt_syserr_ATL->Draw("samee2");
		 hsigmapt_staerr->Draw("samee");
		 hsigmapt_staerr_ATL->Draw("samee");
	 */
	/*
		 fitft_opt1->Draw("same");
		 fitft_opt2->Draw("same");

		 fitft_opt3->Draw("same");
		 fitft_opt4->Draw("same");

		 gaeCMSsta->Draw("samep");
		 gaeATLsta->Draw("samep");
		 fitft_opt5->Draw("same");
		 fitft_opt6->Draw("same");
	 */
	gaeCOMsys->Draw("samee2");
	gaeCOMsta->Draw("samep");
	fitft_opt7->Draw("same");
	gaeCOMsys2->Draw("samee2");
	gaeCOMsta2->Draw("samep");
	leg->Clear();
	leg = new TLegend(0.20,0.70,0.50,0.90,"");

	leg->SetBorderSize(0);
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(1001);
	leg->SetTextFont(42);
	leg->SetTextSize(0.040);


	leg->AddEntry(gaeCOMsta,"Reweighted CMS+ATLAS 7 TeV data","lpf");
	leg->AddEntry(gaeCOMsta2,"Twice reweighted CMS+ATLAS 7 TeV data","lpf");
	leg->AddEntry(fitft_opt7,"(7) Fit with reweighted CMS+ATLAS data","l");
	leg->Draw("");

	cSigma->SaveAs("Plots/cSigma_CompCombinedFit.pdf");

	fout->cd();
	hparafit->Write();
	gaeCMSsta->Write();
	gaeATLsta->Write();
	gaeCOMsta->Write();
	gaeCMSsys->Write();
	gaeATLsys->Write();
	gaeCOMsys->Write();
	gaeCOMsta2->Write();
	gaeCOMsys2->Write();


	fout->Close();
}
