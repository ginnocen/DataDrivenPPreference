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

#define BIN_NUM 600

#define REBIN 5
#define REBINp 6

#define HMINp -3.0
#define HMAXp 3.0

#define REBINmax 49

void Bplusdsigmady_all(int option, bool isBinned)
{

	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	gStyle->SetPadTopMargin(0.095);
	gStyle->SetPadLeftMargin(0.130);
	gStyle->SetPadRightMargin(0.045);

	std::string ename;
	ifstream getdata;
	if (option==5) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_10_60_all_5TeV.dat");ename="5TeV_rap30_pt1060";}
	if (option==7) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_10_60_all_7TeV.dat");ename="7TeV_rap30_pt1060";}
	if (option==71) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_5_120_all_7TeV.dat");ename="7TeV_rap30_pt5120";}
	if (option==72) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_9_120_all_7TeV.dat");ename="7TeV_rap30_pt9120";}
	if (option==2) {getdata.open("../FONLLInputs/fo_pPb_y_rap30_pt_10_60_all_2760GeV.dat");ename="2760GeV_rap30_pt1060";}

	std::string isBinnedrmk;

	int REBINn;

	if (isBinned) {isBinnedrmk="Binned";REBINn=REBIN;}
	else {isBinnedrmk="Unbinned";REBINn=REBINmax-1;}

	std::string tlatexrem;
	int tlatexptmin, tlatexptmax;

	if (option==5) {tlatexrem="5.02 TeV";tlatexptmin=10;tlatexptmax=60;}
	if (option==7) {tlatexrem="7 TeV";tlatexptmin=10;tlatexptmax=60;}
	if (option==71) {tlatexrem="7 TeV";tlatexptmin=5;tlatexptmax=120;}
	if (option==72) {tlatexrem="7 TeV";tlatexptmin=9;tlatexptmax=120;}
	if (option==2) {tlatexrem="2.76 TeV";tlatexptmin=10;tlatexptmax=60;}


	TFile*foutput=new TFile(Form("../outputBplusy_%s_%s.root",isBinnedrmk.c_str(),ename.c_str()),"recreate");

	if(!getdata.is_open()) cout<<"Opening the file fails"<<endl;

	float central[BIN_NUM],y[BIN_NUM];
	float min_all[BIN_NUM],max_all[BIN_NUM],min_sc[BIN_NUM],max_sc[BIN_NUM],min_mass[BIN_NUM],max_mass[BIN_NUM],min_pdf[BIN_NUM],max_pdf[BIN_NUM];
	float fr_05_05[BIN_NUM],fr_20_20[BIN_NUM],fr_20_10[BIN_NUM],fr_10_20[BIN_NUM],fr_10_05[BIN_NUM],fr_05_10[BIN_NUM];
	float tem;

	int i;
	for(i=0;i<BIN_NUM;i++)
	{
		getdata>>tem;
		// all the value is corrected to "y_{CM}"
		if (option==5) y[i]=tem-0.465; else y[i]=tem;
		getdata>>central[i];
		getdata>>min_all[i];
		getdata>>max_all[i];
		getdata>>min_sc[i];
		getdata>>max_sc[i];
		getdata>>min_mass[i];
		getdata>>max_mass[i];
		if (option!=2) {
			getdata>>min_pdf[i];
			getdata>>max_pdf[i];
		}
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
		if (option!=2) {
			hminpdf->SetBinContent(i+1,min_pdf[i]);
			hmaxpdf->SetBinContent(i+1,max_pdf[i]);
		}
		hfr_05_05->SetBinContent(i+1,fr_05_05[i]);
		hfr_20_20->SetBinContent(i+1,fr_20_20[i]);
		hfr_20_10->SetBinContent(i+1,fr_20_10[i]);
		hfr_10_20->SetBinContent(i+1,fr_10_20[i]);
		hfr_05_10->SetBinContent(i+1,fr_05_10[i]);
		hfr_10_05->SetBinContent(i+1,fr_10_05[i]);
	}

	TCanvas* c1 = new TCanvas("c1","",500,500);

	TH1F* hbase = new TH1F("hbase","",60,-3.001,3.001);
	hbase->GetXaxis()->SetTitle("B^{+} y_{CM}");
	hbase->GetYaxis()->SetTitle(Form("FONLL expectation at %s",tlatexrem.c_str()));
	hbase->GetYaxis()->SetRangeUser(0.00,2.6*hpt->GetMaximum());
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
	//hpt->GetYaxis()->SetRangeUser(0,7000000);
	//hpt->GetYaxis()->SetRangeUser(0,2.0*hpt->GetMaximum());
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
	c1->SaveAs(Form("Plots/CompFONLLdsigmady_%s_%s.pdf",isBinnedrmk.c_str(),ename.c_str()));


	// REBIN here

	double rebiny[REBIN+1] = {-2.865,-1.935,-1.0,0.0,1.0,1.935};//CM frame

	Double_t rebin[61];
	Float_t apt[60];
	Float_t aptl[60];
	Float_t asigma[60],aminall[60],amaxall[60],aminsc[60],amaxsc[60],aminmass[60],amaxmass[60],aminpdf[60],amaxpdf[60],aerrorl[60],aerrorh[60];

	// number of every rebinned bin
	double bin_num[60];

	if (isBinned) {
		for (i=0;i<REBIN;i++) {
			apt[i]=(rebiny[i+1]+rebiny[i])/2;//bin middle
			aptl[i]=(rebiny[i+1]-rebiny[i])/2;//bin half width
			rebin[i]=rebiny[i];//Rebin edge
		}
		rebin[REBIN]=1.935;
	}
	else {
		for (i=0;i<60;i++) {
			if (option==5) {
				apt[i]=((-3.0-0.465+i*0.1)+(-3.0-0.465+(i+1)*0.1))/2;
				rebin[i]=-3.0-0.465+i*0.1;
			}

			else {
				apt[i]=((-3.0+i*0.1)+(-3.0+(i+1)*0.1))/2;
				rebin[i]=-3.0+i*0.1;
			}
			aptl[i]=0.05;
			//std::cout << i << " - " << rebin[i] << std::endl;
		}
		//rebin[60] = 1.935;
		if (option==5) rebin[60]=3.0-0.465; else rebin[60]=3.0;
	}

	//rebin[REBINn] = 1.935;

	//std::cout << "##### REBINn" << REBINn << std::endl;


	TH1F* hpt_rebin = (TH1F*)hpt->Rebin(REBINn,"hpt_rebin",rebin);
	TH1F* hminall_rebin = (TH1F*)hminsc->Rebin(REBINn,"hminall_rebin",rebin);
	TH1F* hmaxall_rebin = (TH1F*)hmaxsc->Rebin(REBINn,"hmaxall_rebin",rebin);
	TH1F* hminsc_rebin = (TH1F*)hminsc->Rebin(REBINn,"hminsc_rebin",rebin);
	TH1F* hmaxsc_rebin = (TH1F*)hmaxsc->Rebin(REBINn,"hmaxsc_rebin",rebin);
	TH1F* hminmass_rebin = (TH1F*)hminmass->Rebin(REBINn,"hminmass_rebin",rebin);
	TH1F* hmaxmass_rebin = (TH1F*)hmaxmass->Rebin(REBINn,"hmaxmass_rebin",rebin);
	TH1F* hminpdf_rebin;
	TH1F* hmaxpdf_rebin;
	if (option!=2) {
		hminpdf_rebin = (TH1F*)hminpdf->Rebin(REBINn,"hminpdf_rebin",rebin);
		hmaxpdf_rebin = (TH1F*)hmaxpdf->Rebin(REBINn,"hmaxpdf_rebin",rebin);
	}
	int j;
	double norm=1.;

	for(j=0;j<REBINn;j++)
	{
		bin_num[j]=(rebin[j+1]-rebin[j])/0.01;//number of every rebinned bin
		tem = hpt_rebin->GetBinContent(j+1);
		asigma[j] = tem*norm/bin_num[j];
		//std::cout << tem << " " << asigma[j] << std::endl;

		tem = hminall_rebin->GetBinContent(j+1);
		aminall[j] = tem*norm/bin_num[j];

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

		if (option!=2) {
			tem = hminpdf_rebin->GetBinContent(j+1);
			aminpdf[j] = tem*norm/bin_num[j];

			tem = hmaxpdf_rebin->GetBinContent(j+1);
			amaxpdf[j] = tem*norm/bin_num[j];
		}
		aerrorl[j] = asigma[j]-aminall[j];//all,sc,mass,pdf
		aerrorh[j] = amaxall[j]-asigma[j];//all,sc,mass,pdf
	}


	if (option==5) std::cout << "------- pPb 5.02 TeV -------, " << tlatexptmin << " < p_{T} < " << tlatexptmax << " -------" << std::endl;
	else std::cout << "------- pp " << tlatexrem.c_str() << " -------, " << tlatexptmin << " < p_{T} < " << tlatexptmax << " -------" << std::endl;

	std::cout << std::endl;

	std::cout << "##### REBINn" << REBINn << std::endl;

	TGraphAsymmErrors* gae = new TGraphAsymmErrors(REBINn, apt, asigma, aptl, aptl, aerrorl, aerrorh);
	TGraphAsymmErrors* gaeSigmaDecay = new TGraphAsymmErrors(REBINn, apt, asigma, aptl, aptl, aerrorl, aerrorh);
	TGraphAsymmErrors* gaeSigmaDecayv2 = new TGraphAsymmErrors(REBINn, apt, asigma, aptl, aptl, aerrorl, aerrorh);

	gae->SetTitle(";y_{CM};d#sigma (B admix) / dy (pb)");
	gae->SetFillColor(2);
	gae->SetFillStyle(3001);

	TCanvas* cr = new TCanvas("cr","cr",600,500);

	TH2F* hempty;
	hempty=new TH2F("hempty","",50,-2.965,2.035,100.,0,10000000);
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
	hempty->GetYaxis()->SetLimits(0.00,2.6*hmaxall->GetMaximum());
	hempty->Draw();
	hminall->SetLineColor(2);
	hmaxall->SetLineColor(2);
	hpt->SetLineColor(2);
	hminall->Draw("same");
	hmaxall->Draw("same");
	hpt->Draw("same");
	gae->SetLineWidth(1);

	//TF1* fitft = new TF1("fitft","[0]+[1]*pow((x-0.465),2)",-3.0,3.0);
	TF1* fitft = new TF1("fitft","[0]+[1]*pow(x,2)",-3.0,3.0);

	std::cout << "##################################################################" << std::endl;
	gae->Fit(fitft,"","",-2.865,1.935);
	std::cout << "##################################################################" << std::endl;

	gae->Draw("psame");

	TLatex * tlatex;
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
	//cr->SaveAs("Plots/cBmesonPredFONLLBplusyBinning.eps");
	cr->SaveAs(Form("Plots/cBmesonPredFONLLBplusy%s_%s.pdf",isBinnedrmk.c_str(),ename.c_str()));

	//TGraphAsymmErrors* gaeSigmaDecay=(TGraphAsymmErrors*)gae->Clone();
	gaeSigmaDecay->SetName("gaeSigmaDecay");

	//TGraphAsymmErrors* gaeSigmaDecayv2=(TGraphAsymmErrors*)gae->Clone();
	gaeSigmaDecayv2->SetName("gaeSigmaDecayv2");

	double BRchain=6.09604e-5;
	double Fraction=0.401;

	//double norm=208.;
	norm=208.;
	double BRFraction=BRchain*Fraction;

	for (i=0;i<gaeSigmaDecay->GetN();i++){
		gaeSigmaDecay->GetY()[i] *= BRFraction*norm;
		//gaeSigmaDecay->SetPoint(i,gaeSigmaDecay->GetX()[i],gaeSigmaDecay->GetY()[i]*BRFraction*norm);
		gaeSigmaDecay->SetPointEYhigh(i,gaeSigmaDecay->GetErrorYhigh(i)*BRFraction*norm);
		gaeSigmaDecay->SetPointEYlow(i,gaeSigmaDecay->GetErrorYlow(i)*BRFraction*norm); 

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
		hpt_rebinv2->SetBinContent(i+1,hpt_rebin->GetBinContent(i+1)*Fraction*norm*(1.0e-6)/bin_num[i]);
		hminall_rebinv2->SetBinContent(i+1,hminall_rebin->GetBinContent(i+1)*Fraction*norm*(1.0e-6)/bin_num[i]);
		hmaxall_rebinv2->SetBinContent(i+1,hmaxall_rebin->GetBinContent(i+1)*Fraction*norm*(1.0e-6)/bin_num[i]);
	}

	std::cout << "GetY: " << gaeSigmaDecay->GetY()[0] << std::endl; 
	//###TF1* fitft2 = new TF1("fitft2","[0]+[1]*pow(x,2)",-3.0,3.0);


	if (option==5) std::cout << "####### Fit with FONLL pPb 5.02 TeV, " << tlatexptmin << " < p_{T} < " << tlatexptmax << " , -2.865<y_{CM}<1.935" << std::endl;
	else std::cout << "####### Fit with FONLL pp " << tlatexrem.c_str() << " -------, " << tlatexptmin << " < p_{T} < " << tlatexptmax << " , -2.865<y_{CM}<1.935" << std::endl;

	gaeSigmaDecay->SetFillColor(2);
	gaeSigmaDecay->SetFillStyle(3001); 
	//  if (option==5 || option==7 || option==71 || option==72) gaeSigmaDecay->SetTitle(";y_{LAB};d#sigma(B^{+} full chain)/dy #times A");
	//  else gaeSigmaDecay->SetTitle(";y_{CM};d#sigma(B^{+} full chain)/dy #times A");
	gaeSigmaDecay->SetTitle(";y_{CM};d#sigma(B^{+} full chain)/dy #times A");

	//if (option==5 || option==7 || option==71 || option==72) hempty=new TH2F("hempty","",60,-3.0,3.0,150.,0,150000);
	//else hempty=new TH2F("hempty","",50,-2.965,2.035,15.,0,15000);
	hempty=new TH2F("hempty","",60,-3.0,3.0,150.,0,150000);

	//if (option==5 || option==7 || option==71 || option==72) hempty->GetXaxis()->SetTitle("y_{LAB}");
	//else hempty->GetXaxis()->SetTitle("y_{CM}");
	hempty->GetXaxis()->SetTitle("y_{CM}");
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
	hempty->GetYaxis()->SetTitle("d#sigma/dy(B^{+}) #times A (pb)");
	TCanvas*canvas=new TCanvas("canvas","canvas",600,500);
	hempty->Draw();
	gaeSigmaDecay->Draw("psame");

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
	canvas->SaveAs(Form("Plots/canvasBplusy%s_%s.pdf",isBinnedrmk.c_str(),ename.c_str()));

	hmaxall_rebinv2->GetYaxis()->SetRangeUser(0.00,2.6*hmaxall_rebinv2->GetMaximum());
	//hmaxall_rebinv2->GetYaxis()->SetRangeUser(0,500);
	hmaxall_rebinv2->Draw();
	hpt_rebinv2->Draw("same");
	hminall_rebinv2->Draw("same");
	canvas->SaveAs(Form("Plots/canvasBplusyv2%s_%s.pdf",isBinnedrmk.c_str(),ename.c_str()));

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
/*
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
*/
	foutput->Write();

}
