//Made by Gian Michele, modified by Hyunchul Kim

#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include <cmath>
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include <fstream>
#include <iostream>
#include "TMath.h"

#define BIN_NUM 460		// 5~120 GeV, interval : 0.25 GeV

#define HMIN 5   
#define HMAX 120    

#define REBIN_bin0 6	// CMS pPb pt bin
#define REBIN_bin1 5	// CMS pp pt bin
#define REBIN_bin2 8	// ATLAS pp pt bin
#define REBIN 55			// For 5~60, bin width : 1 GeV     


void Bplusdsigmadpt_all(int option=5, bool isBinned=true, int bincase=0)
{

	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	gStyle->SetPadTopMargin(0.095);
	gStyle->SetPadLeftMargin(0.130);
	gStyle->SetPadRightMargin(0.045);

	std::string ename;
	ifstream getdata;

	if (option==5) {getdata.open("../FONLLInputs/fo_pPb_dpt_pt5120_rap24_all_5TeV.dat");ename="5TeV";}
	if (option==7) {getdata.open("../FONLLInputs/fo_pp_dpt_pt5120_rap24_all_7TeV.dat");ename="7TeV";}
	if (option==2) {getdata.open("../FONLLInputs/fo_pp_dpt_pt5120_rap24_all_2760GeV.dat");ename="2760GeV";}

	if (option==71) {getdata.open("../FONLLInputs/fo_pp_dpt_pt5120_rap24_all_7TeV_MSTW2008nlo68cl.dat");ename="7TeV_MSTW2008nlo68cl";}
	if (option==72) {getdata.open("../FONLLInputs/fo_pp_dpt_pt5120_rap24_all_7TeV_NNPDF30nlo_as0118.dat");ename="7TeV_NNPDF30nlo_as0118";}

	std::string isBinnedrmk;
	int REBINn;

	if (isBinned) {
		if (bincase==0) {isBinnedrmk="BinnedpPb";REBINn=REBIN_bin0;}
		if (bincase==1) {isBinnedrmk="BinnedCMS";REBINn=REBIN_bin1;}
		if (bincase==2) {isBinnedrmk="BinnedATL";REBINn=REBIN_bin2;}
	}
	else {isBinnedrmk="Unbinned";REBINn=REBIN;}//5~60, 1 GeV interval

	std::string tlatexrem;
	double tlatexptmin, tlatexptmax;

	if (option==5) {tlatexrem="5.02 TeV";tlatexptmin=-2.865;tlatexptmax=1.935;}
	if (option==7) {tlatexrem="7 TeV";tlatexptmin=-2.4;tlatexptmax=2.4;}
	if (option==2) {tlatexrem="2.76 TeV";tlatexptmin=-2.4;tlatexptmax=2.4;}
	if (option==71) {tlatexrem="7 TeV";tlatexptmin=-2.4;tlatexptmax=2.4;}
	if (option==72) {tlatexrem="7 TeV";tlatexptmin=-2.4;tlatexptmax=2.4;}

	TFile*foutput=new TFile(Form("../ResultsBplus/outputBplus_%s_%s.root",isBinnedrmk.c_str(),ename.c_str()),"recreate");

	if(!getdata.is_open()) cout<<"Opening the file fails"<<endl;

	float central[BIN_NUM];
	float min_all[BIN_NUM],max_all[BIN_NUM],min_sc[BIN_NUM],max_sc[BIN_NUM],min_mass[BIN_NUM],max_mass[BIN_NUM],min_pdf[BIN_NUM],max_pdf[BIN_NUM];
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
		// if (option!=2) {	
		getdata>>min_pdf[i];
		getdata>>max_pdf[i];
		// }
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
		// if (option!=2) {
		hminpdf->SetBinContent(i+1,min_pdf[i]);
		hmaxpdf->SetBinContent(i+1,max_pdf[i]);
		// }
		hfr_05_05->SetBinContent(i+1,fr_05_05[i]);
		hfr_20_20->SetBinContent(i+1,fr_20_20[i]);
		hfr_20_10->SetBinContent(i+1,fr_20_10[i]);
		hfr_10_20->SetBinContent(i+1,fr_10_20[i]);
		hfr_05_10->SetBinContent(i+1,fr_05_10[i]);
		hfr_10_05->SetBinContent(i+1,fr_10_05[i]);
	}

	// REBIN here

	//double rebiny[6] = {10,15,20,25,30,60};//rebin edge
	double rebiny[REBIN_bin0+1]		  = {5.,10.,15.,20.,25.,30.,60.};//rebin edge for CMS pPb binning
	double rebiny_CMS[REBIN_bin1+1] = {5.,10.,13.,17.,24.,30.};//rebin edge for CMS pp binning
	double rebiny_ATL[REBIN_bin2+1] = {9.,13.,16.,20.,25.,35.,50.,70.,120.};//rebin edge for ATLAS pp binning

	Double_t rebin[REBINn+1];
	Float_t apt[REBINn];
	Float_t aptl[REBINn];
	Float_t asigma[REBINn],aminall[REBINn],amaxall[REBINn],aminsc[REBINn],amaxsc[REBINn],aminmass[REBINn],amaxmass[REBINn],aminpdf[REBINn],amaxpdf[REBINn],aerrorl[REBINn],aerrorh[REBINn];
	// number of every rebinned bin
	double bin_num[REBINn];

	if (isBinned) {
		if (bincase==0) {
			for (i=0;i<REBIN_bin0;i++) {
				apt[i]=(rebiny[i+1]+rebiny[i])/2;//bin middle
				aptl[i]=(rebiny[i+1]-rebiny[i])/2;//bin half width
				rebin[i]=rebiny[i];//Rebin Edge
			}
			rebin[REBIN_bin0]=60.0;
		}
		if (bincase==1) {
			for (i=0;i<REBIN_bin1;i++) {
				apt[i]=(rebiny_CMS[i+1]+rebiny_CMS[i])/2;//bin middle
				aptl[i]=(rebiny_CMS[i+1]-rebiny_CMS[i])/2;//bin half width
				rebin[i]=rebiny_CMS[i];//Rebin Edge
			}
			rebin[REBIN_bin1]=30.0;
		}
		if (bincase==2) {
			for (i=0;i<REBIN_bin2;i++) {
				apt[i]=(rebiny_ATL[i+1]+rebiny_ATL[i])/2;//bin middle
				aptl[i]=(rebiny_ATL[i+1]-rebiny_ATL[i])/2;//bin half width
				rebin[i]=rebiny_ATL[i];//Rebin Edge
			}
			rebin[REBIN_bin2]=120.0;
		}
	}
	else {
		for (i=0;i<REBIN;i++) {
			apt[i]=((5.0+i*1.0)+(5.0+(i+1)*1.0))/2;//bin middle
			aptl[i]=0.5;//bin half width
			rebin[i]=5.0+i*1.0;//Rebin Edge
			//std::cout << i << " - " << rebin[i] << std::endl;
		}
		rebin[REBIN]=60.0;
	}

	std::cout << "##### REBINn" << REBINn << std::endl;

	TH1F* hpt_rebin = (TH1F*)hpt->Rebin(REBINn,"hpt_rebin",rebin);
	TH1F* hminall_rebin = (TH1F*)hminsc->Rebin(REBINn,"hminall_rebin",rebin);
	TH1F* hmaxall_rebin = (TH1F*)hmaxsc->Rebin(REBINn,"hmaxall_rebin",rebin);
	TH1F* hminsc_rebin = (TH1F*)hminsc->Rebin(REBINn,"hminsc_rebin",rebin);
	TH1F* hmaxsc_rebin = (TH1F*)hmaxsc->Rebin(REBINn,"hmaxsc_rebin",rebin);
	TH1F* hminmass_rebin = (TH1F*)hminmass->Rebin(REBINn,"hminmass_rebin",rebin);
	TH1F* hmaxmass_rebin = (TH1F*)hmaxmass->Rebin(REBINn,"hmaxmass_rebin",rebin);
	TH1F* hminpdf_rebin;
	TH1F* hmaxpdf_rebin;
	//if (option!=2) {
	hminpdf_rebin = (TH1F*)hminpdf->Rebin(REBINn,"hminpdf_rebin",rebin);
	hmaxpdf_rebin = (TH1F*)hmaxpdf->Rebin(REBINn,"hmaxpdf_rebin",rebin);
	//}

	int j;
	double norm=1.;

	for(j=0;j<REBINn;j++)
	{
		bin_num[j]=(rebin[j+1]-rebin[j])/0.25;//number of each rebinned bin

		tem = hpt_rebin->GetBinContent(j+1);
		asigma[j] = tem*norm/bin_num[j];

		tem = hminall_rebin->GetBinContent(j+1);
		aminall[j] = tem*norm/bin_num[j];

		//tem = hmaxsc_rebin->GetBinContent(j+1);
		tem = hmaxall_rebin->GetBinContent(j+1);
		amaxall[j] = tem*norm/bin_num[j];

		tem = hminsc_rebin->GetBinContent(j+1);
		aminsc[j] = tem*norm/bin_num[j];

		tem = hmaxsc_rebin->GetBinContent(j+1);
		amaxsc[j] = tem*norm/bin_num[j];

		tem = hminmass_rebin->GetBinContent(j+1);
		aminmass[j] = tem*norm/bin_num[j];

		tem = hmaxmass_rebin->GetBinContent(j+1);
		amaxmass[j] = tem*norm/bin_num[j];

		//		if (option!=2) {
		tem = hminpdf_rebin->GetBinContent(j+1);
		aminpdf[j] = tem*norm/bin_num[j];

		tem = hmaxpdf_rebin->GetBinContent(j+1);
		amaxpdf[j] = tem*norm/bin_num[j];
		//		}
		aerrorl[j] = asigma[j]-aminall[j];//all,sc,mass,pdf
		aerrorh[j] = amaxall[j]-asigma[j];//all,sc,mass,pdf
	}

	if (option==5) std::cout << "------- pPb 5.02 TeV -------, " << tlatexptmin << " < y_{CM} < " << tlatexptmax << " -------" << std::endl;
	else std::cout << "------- pp " << tlatexrem.c_str() << " -------, " << tlatexptmin << " < y_{CM} < " << tlatexptmax << " -------" << std::endl;

	std::cout << std::endl;

	std::cout << "##### REBINn" << REBINn << std::endl;

	TGraphAsymmErrors* gae = new TGraphAsymmErrors(REBINn, apt, asigma, aptl, aptl, aerrorl, aerrorh);

	TGraphAsymmErrors* gaeSigmaDecay=(TGraphAsymmErrors*)gae->Clone();
	gaeSigmaDecay->SetName("gaeSigmaDecay");

	//double BRchain=6.09604e-5;
	double BRchain=1.;
	double Fraction=0.401;

	norm=208.;
	double BRFraction=BRchain*Fraction;

	for (i=0;i<gaeSigmaDecay->GetN();i++){
		gaeSigmaDecay->GetY()[i] *= BRFraction*norm;
		gaeSigmaDecay->SetPointEYhigh(i,gaeSigmaDecay->GetErrorYhigh(i)*BRFraction*norm);
		gaeSigmaDecay->SetPointEYlow(i,gaeSigmaDecay->GetErrorYlow(i)*BRFraction*norm); 
		// Note - gaeSigmaDecay   : *BRchain*Fraction*N(208) -> only consider B+->J/Psi+K+
	}
	std::cout << "GetY: " << gaeSigmaDecay->GetY()[0] << std::endl;

	if (option==5) std::cout << "####### Fit with FONLL pPb 5.02 TeV, " << tlatexptmin << " < p_{T} < " << tlatexptmax << " , -2.865<y_{CM}<1.935" << std::endl;
	else std::cout << "####### Fit with FONLL pp " << tlatexrem.c_str() << " -------, " << tlatexptmin << " < p_{T} < " << tlatexptmax << " , -2.4<y_{CM}<2.4" << std::endl;

	gaeSigmaDecay->SetFillColor(2);
	gaeSigmaDecay->SetFillStyle(3001); 
	gaeSigmaDecay->SetTitle(";p_{T}(GeV/c);d#sigma(B^{+} full chain)/dp_{T} #times A (pb GeV^{-1}c)");

	//#TH2F* hempty=new TH2F("hempty","",10,0,70.,10.,1.,500000);  
	TH2F* hempty=new TH2F("hempty","",10,0,70.,10.,0.01,100000000000.0);  
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
	hempty->GetYaxis()->SetLabelSize(0.04);  
	hempty->GetYaxis()->SetTitle("d#sigma/dp_{T}(B^{+}) #times A (pb c/GeV)");
	//hempty->GetYaxis()->SetRangeUser(1.0,1000000.0);
	double rumin, rumax;
	double reval;
	reval=modf(log10(TMath::MinElement(gaeSigmaDecay->GetN(),gaeSigmaDecay->GetY())),&rumin);
	reval=modf(log10(TMath::MaxElement(gaeSigmaDecay->GetN(),gaeSigmaDecay->GetY())),&rumax);
	std::cout << TMath::MinElement(gaeSigmaDecay->GetN(),gaeSigmaDecay->GetY()) << " " << TMath::MaxElement(gaeSigmaDecay->GetN(),gaeSigmaDecay->GetY()) << " " << rumin << " " << rumax << std::endl;

	hempty->GetYaxis()->SetLimits(pow(10,rumin-1),pow(10,rumax+2));
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();

	TCanvas*canvas=new TCanvas("canvas","canvas",600,500);
	canvas->cd();
	canvas->Clear();
	canvas->SetLogy();
	hempty->Draw();
	gaeSigmaDecay->Draw("psame");

	TLatex * tlatex=new TLatex(0.2,0.85,"B^{+}=40.1%, |y_{LAB}|<2.4, BR unc not shown");
	//tlatex=new TLatex(0.2,0.85,"B^{+}=40.1%, |y_{LAB}|<2.4, BR unc not shown");
	tlatex->SetNDC();
	tlatex->SetTextColor(1);
	tlatex->SetTextFont(42);
	tlatex->SetTextSize(0.04);
	tlatex->Draw();

	gae->SetName("gaeBplus");
	gaeSigmaDecay->SetName("gaeSigmaDecayBplus");
	//###gaeSigmaDecayv2->SetName("gaeSigmaDecayv2Bplus");

	//canvas->SaveAs("../ResultsBplus/canvasBplusFinePt_7TeV.eps");
	canvas->SaveAs(Form("../ResultsBplus/canvasBplus%s_%s.pdf",isBinnedrmk.c_str(),ename.c_str()));

	// Draw the unbinned FONLL values with different settings

	TCanvas* c1 = new TCanvas("c1","",500,500);
	c1->SetLogy(1);

	TH1F* hbase = new TH1F("hbase","",13,5.0,65.0);
	hbase->GetXaxis()->SetTitle("B^{+} p_{T}");
	hbase->GetYaxis()->SetTitle(Form("FONLL expectation at %s",tlatexrem.c_str()));
	//hbase->GetYaxis()->SetRangeUser(100.0,1000000000.0);
	reval=modf(log10(hpt->GetMinimum()),&rumin);
	reval=modf(log10(hpt->GetMaximum()),&rumax);
	std::cout << hpt->GetMinimum() << " " << hpt->GetMaximum() << " " << rumin << " " << rumax << std::endl;
	hbase->GetYaxis()->SetRangeUser(pow(10,rumin+0),pow(10,rumax+1));
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
	c1->SaveAs(Form("../ResultsBplus/CompFONLLdsigmadpt_%s_%s.pdf",isBinnedrmk.c_str(),ename.c_str()));

	gae->SetTitle(";p_{T}(GeV/c);d#sigma (B admix) / dp_{T} (pb c/GeV)");
	gae->SetFillColor(2);
	gae->SetFillStyle(3001);

	TCanvas* cr = new TCanvas("cr","cr",600,500);
	cr->SetLogy();
	//TH2F* hempty=new TH2F("hempty","",10,5,60.,10,0.1,100000000.0);  
	hempty=new TH2F("hempty","",10,5,60.,10,0.1,100000000.0);  
	hempty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hempty->GetYaxis()->SetTitle("d#sigma(B admix)/dp_{T} (pb c/GeV)");
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
	//hempty->GetYaxis()->SetLimits(0.00,2.6*hmaxall->GetMaximum());
	reval=modf(log10(hminall->GetMinimum()),&rumin);
	reval=modf(log10(hmaxall->GetMaximum()),&rumax);
	std::cout << hminall->GetMinimum() << " " << hmaxall->GetMaximum() << " " << rumin << " " << rumax << std::endl;
	hempty->GetYaxis()->SetRangeUser(pow(10,rumin-0),pow(10,rumax+3));
	hempty->Draw();
	hminall->SetLineColor(2);
	hmaxall->SetLineColor(2);
	hpt->SetLineColor(2);
	hminall->Draw("same");
	hmaxall->Draw("same");
	hpt->Draw("same");
	gae->SetLineWidth(1);//3
	/*
		 TF1* fitft = new TF1("fitft","pow(10,[0]*exp([1]+[2]*x)+[3])",5.0,60.0);
		 fitft->SetParameter(0,2.92920e+00);
		 fitft->SetParameter(1,7.06534e-01);
		 fitft->SetParameter(2,-2.47338e-02);
		 fitft->SetParameter(3,2.00309e+00);

	// from option (false,7)
	1  p0           2.92920e+00   1.52820e-01   2.38787e-06  -1.84407e-02
	2  p1           7.06534e-01   5.21733e-02   8.14391e-07  -5.40158e-02
	3  p2          -2.47338e-02   1.04804e-03   2.39403e-08  -1.42800e+00
	4  p3           2.00309e+00   1.09026e-01   1.98609e-06  -1.82817e-02
	 */

	TF1* fitft = new TF1("fitft","pow(10,[0]*exp([1]*x)+[2])",5.0,60.0);

	std::cout << "##################################################################" << std::endl;
	gae->Fit(fitft,"","",5.0,60.0);
	std::cout << "##################################################################" << std::endl;

	gae->Draw("psame");

	//TLatex * tlatex;
	if (option==5) tlatex=new TLatex(0.18,0.85,Form("pp collisions at %s from FONLL, 10<p_{T}<60, %1.3f<y_{CM}<%1.3f",tlatexrem.c_str(),tlatexptmin,tlatexptmax));
	else tlatex=new TLatex(0.18,0.85,Form("pp collisions at %s from FONLL, 10<p_{T}<60, %1.1f<y_{CM}<%1.1f",tlatexrem.c_str(),tlatexptmin,tlatexptmax));

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
	cr->SaveAs(Form("../ResultsBplus/cBmesonPredFONLLBplus%s_%s.pdf",isBinnedrmk.c_str(),ename.c_str()));


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
	//foutput->Write();
}
