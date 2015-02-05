#include <iostream>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>

#define APb 208

#define REBIN_bin0 5
#define REBIN_bin1 5
#define REBIN_bin2 8
#define REBIN_Fine 60 

void CheckScaleFac_v2() {

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetPadTopMargin(0.095);
	gStyle->SetPadLeftMargin(0.130);
	gStyle->SetPadRightMargin(0.045);

	TFile* fin_F2 = new TFile("../ResultsBplus/outputBplus_Binned_2760GeV.root");
	TFile* fin_F7 = new TFile("../ResultsBplus/outputBplus_Binned_7TeV.root");
	TFile* fin_D2 = new TFile("../ResultsBplus/Estimatedpp5TeV_with2760GeV.root");
	TFile* fin_D7 = new TFile("../ResultsBplus/Estimatedpp2760GeV_with7TeV.root");

// CMS 7TeV, pp->(B+)X, unit : micro barn
	double ptbin[REBIN_bin1+1]={5.,10.,13.,17.,24.,30.};
	double sigmapt[REBIN_bin1]={4.07,1.47,0.412,0.181,0.042};
	double sigmapt_sta[REBIN_bin1]={0.47,0.13,0.041,0.015,0.007};
	double sigmapt_sys[REBIN_bin1]={0.31,0.09,0.026,0.012,0.004};

	// ATLAS 7TeV, dsigma/dp_T (B+) * BR(B+->J/psiK+) * BR (J/Psi->mu+mu-), unit : pico barn
	double ptbin_ATL[REBIN_bin2+1]={9.,13.,16.,20.,25.,35.,50.,70.,120.};
	double sigmapt_ATL[REBIN_bin2]={103.4,36.03,15.33,6.056,1.814,0.3477,0.06244,0.006099};
	double sigmapt_sta_ATL[REBIN_bin2]={3.7,0.80,0.25,0.093,0.027,0.0084,0.00293,0.000561};
	double sigmapt_sys_ATL[REBIN_bin2]={7.6,2.32,0.98,0.376,0.115,0.0280,0.00526,0.000666};

	TFile* fin_F7_CMS = new TFile("../ResultsBplus/outputBplus_BinnedCMS_7TeV.root");
	//TFile* fin_F7_ATL = new TFile("../ResultsBplus/outputBplus_BinnedATL_7TeV.root");
	TFile* fin_F7_ATL = new TFile("../ResultsBplus/outputBplus_BinnedATL_7TeVATL.root");
	//TFile* fin_F7_ATL = new TFile("../ResultsBplus/outputBplus_BinnedATL_7TeVATLv2.root");



	TGraphAsymmErrors* gaeFONL7_CMS = (TGraphAsymmErrors*)fin_F7_CMS->Get("gaeSigmaDecayv2Bplus");//already times A(208)
	TGraphAsymmErrors* gaeFONL7_ATL = (TGraphAsymmErrors*)fin_F7_ATL->Get("gaeSigmaDecayv2Bplus");//already times A(208)
	gaeFONL7_CMS->SetName("gaeFONL7_CMS");
	gaeFONL7_ATL->SetName("gaeFONL7_ATL");

	TF1* fitft = new TF1 ("fitft","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);

	TFile* fin_fitpar = new TFile("../ResultsBplus/ScaleCor_reweighted7TeVdata.root");
	TH1D* hparafit = (TH1D*)fin_fitpar->Get("hparafit");
	fitft->SetParameter(0,hparafit->GetBinContent(1));
	fitft->SetParameter(1,hparafit->GetBinContent(2));
	fitft->SetParameter(2,hparafit->GetBinContent(3));

	TGraphAsymmErrors* gaeStat7_CMS = (TGraphAsymmErrors*)gaeFONL7_CMS->Clone("gaeStat7_CMS");
	TGraphAsymmErrors* gaeSyst7_CMS = (TGraphAsymmErrors*)gaeFONL7_CMS->Clone("gaeSyst7_CMS");
	TGraphAsymmErrors* gaeStat7_ATL = (TGraphAsymmErrors*)gaeFONL7_ATL->Clone("gaeStat7_ATL");
	TGraphAsymmErrors* gaeSyst7_ATL = (TGraphAsymmErrors*)gaeFONL7_ATL->Clone("gaeSyst7_ATL");

	TGraphAsymmErrors* gaeFitT7_CMS = (TGraphAsymmErrors*)gaeFONL7_CMS->Clone("gaeFitT7_CMS");
	TGraphAsymmErrors* gaeFitT7_ATL = (TGraphAsymmErrors*)gaeFONL7_ATL->Clone("gaeFitT7_ATL");

	TGraphAsymmErrors* gae_StavsFON_CMS = (TGraphAsymmErrors*)gaeFONL7_CMS->Clone("gae_StavsFON_CMS");
	TGraphAsymmErrors* gae_StavsFON_ATL = (TGraphAsymmErrors*)gaeFONL7_ATL->Clone("gae_StavsFON_ATL");
	TGraphAsymmErrors* gae_SysvsFON_CMS = (TGraphAsymmErrors*)gaeFONL7_CMS->Clone("gae_SysvsFON_CMS");
	TGraphAsymmErrors* gae_SysvsFON_ATL = (TGraphAsymmErrors*)gaeFONL7_ATL->Clone("gae_SysvsFON_ATL");

	TGraphAsymmErrors* gae_StavsFit_CMS = (TGraphAsymmErrors*)gaeFONL7_CMS->Clone("gae_StavsFit_CMS");
	TGraphAsymmErrors* gae_StavsFit_ATL = (TGraphAsymmErrors*)gaeFONL7_ATL->Clone("gae_StavsFit_ATL");
	TGraphAsymmErrors* gae_SysvsFit_CMS = (TGraphAsymmErrors*)gaeFONL7_CMS->Clone("gae_SysvsFit_CMS");
	TGraphAsymmErrors* gae_SysvsFit_ATL = (TGraphAsymmErrors*)gaeFONL7_ATL->Clone("gae_SysvsFit_ATL");

	TGraphAsymmErrors* gae_FitvsFON_CMS = (TGraphAsymmErrors*)gaeFONL7_CMS->Clone("gae_FitvsFON_CMS");
	TGraphAsymmErrors* gae_FitvsFON_ATL = (TGraphAsymmErrors*)gaeFONL7_ATL->Clone("gae_FitvsFON_ATL");

	for (int i=0;i<REBIN_bin1;i++) {
		gaeStat7_CMS->SetPoint(i,gaeStat7_CMS->GetX()[i],sigmapt[i]*1.0);
		gaeStat7_CMS->SetPointEYhigh(i,sigmapt_sta[i]*1.0);
		gaeStat7_CMS->SetPointEYlow(i,sigmapt_sta[i]*1.0);
		gaeSyst7_CMS->SetPoint(i,gaeSyst7_CMS->GetX()[i],sigmapt[i]*1.0);
		gaeSyst7_CMS->SetPointEYhigh(i,sigmapt_sys[i]*1.0);
		gaeSyst7_CMS->SetPointEYlow(i,sigmapt_sys[i]*1.0);
		gaeFitT7_CMS->SetPoint(i,gaeFitT7_CMS->GetX()[i],(fitft->Integral(ptbin[i],ptbin[i+1])));
		gaeFitT7_CMS->SetPointEYhigh(i,0.0);
		gaeFitT7_CMS->SetPointEYlow(i,0.0);

		gae_StavsFON_CMS->SetPoint(i,gae_StavsFON_CMS->GetX()[i],gaeStat7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb);
		gae_StavsFON_CMS->SetPointEYhigh(i,(sqrt(pow(gaeStat7_CMS->GetErrorYhigh(i)/gaeStat7_CMS->GetY()[i],2)+pow(gaeFONL7_CMS->GetErrorYhigh(i)/gaeFONL7_CMS->GetY()[i],2))*(gaeStat7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb)));
		gae_StavsFON_CMS->SetPointEYlow(i,(sqrt(pow(gaeStat7_CMS->GetErrorYlow(i)/gaeStat7_CMS->GetY()[i],2)+pow(gaeFONL7_CMS->GetErrorYlow(i)/gaeFONL7_CMS->GetY()[i],2))*(gaeStat7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb)));

		gae_SysvsFON_CMS->SetPoint(i,gae_SysvsFON_CMS->GetX()[i],gaeSyst7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb);
		gae_SysvsFON_CMS->SetPointEYhigh(i,(sqrt(pow(gaeSyst7_CMS->GetErrorYhigh(i)/gaeSyst7_CMS->GetY()[i],2)+pow(gaeFONL7_CMS->GetErrorYhigh(i)/gaeFONL7_CMS->GetY()[i],2))*(gaeSyst7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb)));
		gae_SysvsFON_CMS->SetPointEYlow(i,(sqrt(pow(gaeSyst7_CMS->GetErrorYlow(i)/gaeSyst7_CMS->GetY()[i],2)+pow(gaeFONL7_CMS->GetErrorYlow(i)/gaeFONL7_CMS->GetY()[i],2))*(gaeSyst7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb)));

		gae_StavsFit_CMS->SetPoint(i,gae_StavsFit_CMS->GetX()[i],gaeStat7_CMS->GetY()[i]/gaeFitT7_CMS->GetY()[i]);
		gae_StavsFit_CMS->SetPointEYhigh(i,gaeStat7_CMS->GetErrorYhigh(i)*(gaeStat7_CMS->GetY()[i]/gaeFitT7_CMS->GetY()[i]));
		gae_StavsFit_CMS->SetPointEYlow(i,gaeStat7_CMS->GetErrorYlow(i)*(gaeStat7_CMS->GetY()[i]/gaeFitT7_CMS->GetY()[i]));

		gae_SysvsFit_CMS->SetPoint(i,gae_SysvsFit_CMS->GetX()[i],gaeSyst7_CMS->GetY()[i]/gaeFitT7_CMS->GetY()[i]);
		gae_SysvsFit_CMS->SetPointEYhigh(i,gaeSyst7_CMS->GetErrorYhigh(i)*(gaeSyst7_CMS->GetY()[i]/gaeFitT7_CMS->GetY()[i]));
		gae_SysvsFit_CMS->SetPointEYlow(i,gaeSyst7_CMS->GetErrorYlow(i)*(gaeSyst7_CMS->GetY()[i]/gaeFitT7_CMS->GetY()[i]));

		gae_FitvsFON_CMS->SetPoint(i,gae_FitvsFON_CMS->GetX()[i],gaeFitT7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb);
		gae_FitvsFON_CMS->SetPointEYhigh(i,gaeFitT7_CMS->GetErrorYlow(i)/(gaeFONL7_CMS->GetY()[i]-gaeFitT7_CMS->GetErrorYlow(i))*(gaeFitT7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb));
		gae_FitvsFON_CMS->SetPointEYlow(i,gaeFitT7_CMS->GetErrorYhigh(i)/(gaeFONL7_CMS->GetY()[i]+gaeFitT7_CMS->GetErrorYhigh(i))*(gaeFitT7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]*APb));
/*
		std::cout << gaeStat7_CMS->GetX()[i] << std::endl;
		std::cout << gaeStat7_CMS->GetY()[i] << std::endl;
		std::cout << gaeStat7_CMS->GetErrorYhigh(i) << std::endl;
		std::cout << gaeStat7_CMS->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
		std::cout << gaeSyst7_CMS->GetX()[i] << std::endl;
		std::cout << gaeSyst7_CMS->GetY()[i] << std::endl;
		std::cout << gaeSyst7_CMS->GetErrorYhigh(i) << std::endl;
		std::cout << gaeSyst7_CMS->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
		std::cout << gaeFONL7_CMS->GetX()[i] << std::endl;
		std::cout << gaeFONL7_CMS->GetY()[i] << std::endl;
		std::cout << gaeFONL7_CMS->GetErrorYhigh(i) << std::endl;
		std::cout << gaeFONL7_CMS->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
		std::cout << gaeFitT7_CMS->GetX()[i] << std::endl;
		std::cout << gaeFitT7_CMS->GetY()[i] << std::endl;
		std::cout << gaeFitT7_CMS->GetErrorYhigh(i) << std::endl;
		std::cout << gaeFitT7_CMS->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
*/
	}

	for (int i=0;i<REBIN_bin2;i++) {

		double BRs;
		BRs = (1.027*0.001)*(5.961*0.01);

		gaeStat7_ATL->SetPoint(i,gaeStat7_ATL->GetX()[i],sigmapt_ATL[i]*1.0/BRs*0.000001);
		gaeStat7_ATL->SetPointEYhigh(i,sigmapt_sta_ATL[i]*1.0/BRs*0.000001);
		gaeStat7_ATL->SetPointEYlow(i,sigmapt_sta_ATL[i]*1.0/BRs*0.000001);
		gaeSyst7_ATL->SetPoint(i,gaeSyst7_ATL->GetX()[i],sigmapt_ATL[i]/BRs*0.000001);
		gaeSyst7_ATL->SetPointEYhigh(i,sigmapt_sys_ATL[i]*1.0/BRs*0.000001);
		gaeSyst7_ATL->SetPointEYlow(i,sigmapt_sys_ATL[i]*1.0/BRs*0.000001);
		gaeFitT7_ATL->SetPoint(i,gaeFitT7_ATL->GetX()[i],(fitft->Integral(ptbin_ATL[i],ptbin_ATL[i+1])));
		gaeFitT7_ATL->SetPointEYhigh(i,0.0);
		gaeFitT7_ATL->SetPointEYlow(i,0.0);
/*
		std::cout << gaeStat7_ATL->GetX()[i] << std::endl;
		std::cout << gaeStat7_ATL->GetY()[i] << std::endl;
		std::cout << gaeStat7_ATL->GetErrorYhigh(i) << std::endl;
		std::cout << gaeStat7_ATL->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
		std::cout << gaeSyst7_ATL->GetX()[i] << std::endl;
		std::cout << gaeSyst7_ATL->GetY()[i] << std::endl;
		std::cout << gaeSyst7_ATL->GetErrorYhigh(i) << std::endl;
		std::cout << gaeSyst7_ATL->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
		std::cout << gaeFONL7_ATL->GetX()[i] << std::endl;
		std::cout << gaeFONL7_ATL->GetY()[i] << std::endl;
		std::cout << gaeFONL7_ATL->GetErrorYhigh(i) << std::endl;
		std::cout << gaeFONL7_ATL->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
		std::cout << gaeFitT7_ATL->GetX()[i] << std::endl;
		std::cout << gaeFitT7_ATL->GetY()[i] << std::endl;
		std::cout << gaeFitT7_ATL->GetErrorYhigh(i) << std::endl;
		std::cout << gaeFitT7_ATL->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
*/


		gae_StavsFON_ATL->SetPoint(i,gae_StavsFON_ATL->GetX()[i],gaeStat7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb);
		gae_StavsFON_ATL->SetPointEYhigh(i,(sqrt(pow(gaeStat7_ATL->GetErrorYhigh(i)/gaeStat7_ATL->GetY()[i],2)+pow(gaeFONL7_ATL->GetErrorYhigh(i)/gaeFONL7_ATL->GetY()[i],2))*(gaeStat7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb)));
		gae_StavsFON_ATL->SetPointEYlow(i,(sqrt(pow(gaeStat7_ATL->GetErrorYlow(i)/gaeStat7_ATL->GetY()[i],2)+pow(gaeFONL7_ATL->GetErrorYlow(i)/gaeFONL7_ATL->GetY()[i],2))*(gaeStat7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb)));
/*
		std::cout << gae_StavsFON_ATL->GetX()[i] << std::endl;
		std::cout << gae_StavsFON_ATL->GetY()[i] << std::endl;
		std::cout << gae_StavsFON_ATL->GetErrorYhigh(i) << std::endl;
		std::cout << gae_StavsFON_ATL->GetErrorYlow(i) << std::endl;
		std::cout << std::endl;
*/

		gae_SysvsFON_ATL->SetPoint(i,gae_SysvsFON_ATL->GetX()[i],gaeSyst7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb);
		gae_SysvsFON_ATL->SetPointEYhigh(i,(sqrt(pow(gaeSyst7_ATL->GetErrorYhigh(i)/gaeSyst7_ATL->GetY()[i],2)+pow(gaeFONL7_ATL->GetErrorYhigh(i)/gaeFONL7_ATL->GetY()[i],2))*(gaeSyst7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb)));
		gae_SysvsFON_ATL->SetPointEYlow(i,(sqrt(pow(gaeSyst7_ATL->GetErrorYlow(i)/gaeSyst7_ATL->GetY()[i],2)+pow(gaeFONL7_ATL->GetErrorYlow(i)/gaeFONL7_ATL->GetY()[i],2))*(gaeSyst7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb)));

/*
		gae_StavsFON_CMS->SetPoint(i,gae_StavsFON_CMS->GetX()[i],gaeStat7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]);
		gae_StavsFON_CMS->SetPointEYhigh(i,(sqrt(pow(gaeStat7_CMS->GetErrorYhigh(i)/gaeStat7_CMS->GetY()[i],2)+pow(gaeFONL7_CMS->GetErrorYhigh(i)/gaeFONL7_CMS->GetY()[i],2))*(gaeStat7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i])));
		gae_StavsFON_CMS->SetPointEYlow(i,(sqrt(pow(gaeStat7_CMS->GetErrorYlow(i)/gaeStat7_CMS->GetY()[i],2)+pow(gaeFONL7_CMS->GetErrorYlow(i)/gaeFONL7_CMS->GetY()[i],2))*(gaeStat7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i])));

		gae_SysvsFON_CMS->SetPoint(i,gae_SysvsFON_CMS->GetX()[i],gaeSyst7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i]);
		gae_SysvsFON_CMS->SetPointEYhigh(i,(sqrt(pow(gaeSyst7_CMS->GetErrorYhigh(i)/gaeSyst7_CMS->GetY()[i],2)+pow(gaeFONL7_CMS->GetErrorYhigh(i)/gaeFONL7_CMS->GetY()[i],2))*(gaeSyst7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i])));
		gae_SysvsFON_CMS->SetPointEYlow(i,(sqrt(pow(gaeSyst7_CMS->GetErrorYlow(i)/gaeSyst7_CMS->GetY()[i],2)+pow(gaeFONL7_CMS->GetErrorYlow(i)/gaeFONL7_CMS->GetY()[i],2))*(gaeSyst7_CMS->GetY()[i]/gaeFONL7_CMS->GetY()[i])));


*/
		gae_StavsFit_ATL->SetPoint(i,gae_StavsFit_ATL->GetX()[i],gaeStat7_ATL->GetY()[i]/gaeFitT7_ATL->GetY()[i]);
		gae_StavsFit_ATL->SetPointEYhigh(i,gaeStat7_ATL->GetErrorYhigh(i)*(gaeStat7_ATL->GetY()[i]/gaeFitT7_ATL->GetY()[i]));
		gae_StavsFit_ATL->SetPointEYlow(i,gaeStat7_ATL->GetErrorYlow(i)*(gaeStat7_ATL->GetY()[i]/gaeFitT7_ATL->GetY()[i]));

		gae_SysvsFit_ATL->SetPoint(i,gae_SysvsFit_ATL->GetX()[i],gaeSyst7_ATL->GetY()[i]/gaeFitT7_ATL->GetY()[i]);
		gae_SysvsFit_ATL->SetPointEYhigh(i,gaeSyst7_ATL->GetErrorYhigh(i)*(gaeSyst7_ATL->GetY()[i]/gaeFitT7_ATL->GetY()[i]));
		gae_SysvsFit_ATL->SetPointEYlow(i,gaeSyst7_ATL->GetErrorYlow(i)*(gaeSyst7_ATL->GetY()[i]/gaeFitT7_ATL->GetY()[i]));

		gae_FitvsFON_ATL->SetPoint(i,gae_FitvsFON_ATL->GetX()[i],gaeFitT7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb);
		gae_FitvsFON_ATL->SetPointEYhigh(i,gaeFitT7_ATL->GetErrorYlow(i)/(gaeFONL7_ATL->GetY()[i]-gaeFitT7_ATL->GetErrorYlow(i))*(gaeFitT7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb));
		gae_FitvsFON_ATL->SetPointEYlow(i,gaeFitT7_ATL->GetErrorYhigh(i)/(gaeFONL7_ATL->GetY()[i]+gaeFitT7_ATL->GetErrorYhigh(i))*(gaeFitT7_ATL->GetY()[i]/gaeFONL7_ATL->GetY()[i]*APb));


	}

//////////////////////////////////////////////////////////////////////////////////////////////
//1) ratio of ATLAS/ pure FONLL, ratio of CMS/ pure FONLL (all at 7 TeV)
//2) ratio of CMS+ATLAS data / fit function ATLAS+CMS (at 7 TeV)
//3) ratio of fit function ATLAS+CMS/FONLL (at 7 TeV)
/////////////////////////////////////////////////////////////////////////////////////////////
	TFile* fout = new TFile("test.root","RECREATE");

	TGraphAsymmErrors* gaeStat2 = (TGraphAsymmErrors*)fin_D2->Get("gaeStatNoSca2760GeV");
	TGraphAsymmErrors* gaeSyst2 = (TGraphAsymmErrors*)fin_D2->Get("gaeSystNoSca2760GeV");
	TGraphAsymmErrors* gaeStat7 = (TGraphAsymmErrors*)fin_D7->Get("gaeStatNoSca7TeV");
	TGraphAsymmErrors* gaeSyst7 = (TGraphAsymmErrors*)fin_D7->Get("gaeSystNoSca7TeV");
	TGraphAsymmErrors* gaeFONL2 = (TGraphAsymmErrors*)fin_F2->Get("gaeSigmaDecayv2Bplus");//already times A(208)
	TGraphAsymmErrors* gaeFONL7 = (TGraphAsymmErrors*)fin_F7->Get("gaeSigmaDecayv2Bplus");//already times A(208)

	gaeStat2->SetName("gaeStat2");
	gaeSyst2->SetName("gaeSyst2");
	gaeStat7->SetName("gaeStat7");
	gaeSyst7->SetName("gaeSyst7");
	gaeFONL2->SetName("gaeFONL2");
	gaeFONL7->SetName("gaeFONL7");

//////////////////////////////////////////////////////////////////

	//TGraphAsymmErrors* gae_2Fvs7F = (TGraphAsymmErrors*)gaeStat2->Clone();
	TGraphAsymmErrors* gae_2Dvs7D = (TGraphAsymmErrors*)gaeStat2->Clone();
	TGraphAsymmErrors* gae_2Dvs2F = (TGraphAsymmErrors*)gaeStat2->Clone();
	TGraphAsymmErrors* gae_7Dvs7F = (TGraphAsymmErrors*)gaeStat2->Clone();

	//gae_2Fvs7F->SetName("gae_2Fvs7F");
	gae_2Dvs7D->SetName("gae_2Dvs7D");
	gae_2Dvs2F->SetName("gae_2Dvs2F");
	gae_7Dvs7F->SetName("gae_7Dvs7F");

////////////////////////////////////////////////////////////

	TGraphAsymmErrors* gae_2Fvs7F			 = (TGraphAsymmErrors*)gaeStat2->Clone();

	TGraphAsymmErrors* gae_2Dvs7D_SystErr = (TGraphAsymmErrors*)gaeStat2->Clone();
	TGraphAsymmErrors* gae_2Dvs7D_StatErr = (TGraphAsymmErrors*)gaeStat2->Clone();

	TGraphAsymmErrors* gae_2Dvs2F_StatErr = (TGraphAsymmErrors*)gaeStat2->Clone();
	TGraphAsymmErrors* gae_2Dvs2F_SystErr = (TGraphAsymmErrors*)gaeStat2->Clone();
	TGraphAsymmErrors* gae_2Dvs2F_FONLErr = (TGraphAsymmErrors*)gaeStat2->Clone();

	TGraphAsymmErrors* gae_7Dvs7F_StatErr = (TGraphAsymmErrors*)gaeStat2->Clone();
	TGraphAsymmErrors* gae_7Dvs7F_SystErr = (TGraphAsymmErrors*)gaeStat2->Clone();
	TGraphAsymmErrors* gae_7Dvs7F_FONLErr = (TGraphAsymmErrors*)gaeStat2->Clone();

	gae_2Fvs7F->SetName("gae_2Fvs7F");

	gae_2Dvs7D_SystErr->SetName("gae_2Dvs7D_SystErr");
	gae_2Dvs7D_StatErr->SetName("gae_2Dvs7D_StatErr");

	gae_2Dvs2F_SystErr->SetName("gae_2Dvs2F_SystErr");
	gae_2Dvs2F_StatErr->SetName("gae_2Dvs2F_StatErr");
	gae_2Dvs2F_FONLErr->SetName("gae_2Dvs2F_FONLErr");

	gae_7Dvs7F_SystErr->SetName("gae_7Dvs7F_SystErr");
	gae_7Dvs7F_StatErr->SetName("gae_7Dvs7F_StatErr");
	gae_7Dvs7F_FONLErr->SetName("gae_7Dvs7F_FONLErr");

/*
	for(int i=0;i<6;i++){
		gaeFONL2->SetPoint(i,gaeFONL2->GetX()[i],gaeFONL2->GetY()[i]/APb);
		gaeFONL7->SetPoint(i,gaeFONL7->GetX()[i],gaeFONL7->GetY()[i]/APb);
	}
*/
	//TFile* fin_Val_FONLL2vs7=new TFile("../ResultsBplus/CompFONLL_Bplus_Binned_Val_2760GeVvs7TeV.root");
	TFile* fin_Norm_FONLL2vs7=new TFile("../ResultsBplus/CompFONLL_Bplus_Binned_Norm_2760GeVvs7TeV.root");

	//TH1D* hfrval=(TH1D*)fin_Val_FONLL2vs7->Get("hfrval");
	TH1D* hfrmaxerr=(TH1D*)fin_Norm_FONLL2vs7->Get("hfrmaxerr");
	TH1D* hfrminerr=(TH1D*)fin_Norm_FONLL2vs7->Get("hfrminerr");


	for(int i=0;i<5;i++){
		//gae_2Fvs7F : error - 2F/7F FONLL error
		gae_2Fvs7F->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeFONL2->GetY()[i+1]/gaeFONL7->GetY()[i+1]);
		//gae_2Fvs7F->SetPointEYhigh(i,sqrt(pow(gaeFONL2->GetErrorYhigh(i+1)/gaeFONL2->GetY()[i+1],2)+pow(gaeFONL2->GetErrorYhigh(i+1)/gaeFONL2->GetY()[i+1],2))*(gaeFONL2->GetY()[i+1]/gaeFONL7->GetY()[i+1]));
		//gae_2Fvs7F->SetPointEYlow(i,sqrt(pow(gaeFONL2->GetErrorYlow(i+1)/gaeFONL2->GetY()[i+1],2)+pow(gaeFONL2->GetErrorYlow(i+1)/gaeFONL2->GetY()[i+1],2))*(gaeFONL2->GetY()[i+1]/gaeFONL7->GetY()[i+1]));
		gae_2Fvs7F->SetPointEYhigh(i,hfrmaxerr->GetBinContent(i+2)*(gaeFONL2->GetY()[i+1]/gaeFONL7->GetY()[i+1]));
		gae_2Fvs7F->SetPointEYlow(i,hfrminerr->GetBinContent(i+2)*(gaeFONL2->GetY()[i+1]/gaeFONL7->GetY()[i+1]));

//	gae_2Dvs7D_SystErr->SetName("gae_2Dvs7D_SystErr");
//	gae_2Dvs7D_StatErr->SetName("gae_2Dvs7D_StatErr");

		//gae_2Dvs7D : error - 2D/7D sys error
		gae_2Dvs7D_SystErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeSyst2->GetY()[i]/gaeSyst7->GetY()[i+1]);
		gae_2Dvs7D_StatErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeStat2->GetY()[i]/gaeStat7->GetY()[i+1]);


		gae_2Dvs7D_SystErr->SetPointEYhigh(i,sqrt(pow(gaeSyst2->GetErrorYhigh(i)/gaeSyst2->GetY()[i],2)+pow(gaeSyst7->GetErrorYhigh(i+1)/gaeSyst7->GetY()[i+1],2))*(gaeSyst2->GetY()[i]/gaeSyst7->GetY()[i+1]));
		gae_2Dvs7D_SystErr->SetPointEYlow(i,sqrt(pow(gaeSyst2->GetErrorYlow(i)/gaeSyst2->GetY()[i],2)+pow(gaeSyst7->GetErrorYlow(i+1)/gaeSyst7->GetY()[i+1],2))*(gaeSyst2->GetY()[i]/gaeSyst7->GetY()[i+1]));
		gae_2Dvs7D_StatErr->SetPointEYhigh(i,sqrt(pow(gaeStat2->GetErrorYhigh(i)/gaeStat2->GetY()[i],2)+pow(gaeStat7->GetErrorYhigh(i+1)/gaeStat7->GetY()[i+1],2))*(gaeStat2->GetY()[i]/gaeStat7->GetY()[i+1]));
		gae_2Dvs7D_StatErr->SetPointEYlow(i,sqrt(pow(gaeStat2->GetErrorYlow(i)/gaeStat2->GetY()[i],2)+pow(gaeStat7->GetErrorYlow(i+1)/gaeStat7->GetY()[i+1],2))*(gaeStat2->GetY()[i]/gaeStat7->GetY()[i+1]));


		//gae_2Dvs2F : error - 2D sys, 2F FONLL error

		//gae_2Dvs2F_FONLErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeSyst2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb);
		//gae_2Dvs2F_FONLErr->SetPointEYhigh(i,sqrt(pow(gaeSyst2->GetErrorYhigh(i)/gaeSyst2->GetY()[i],2)+pow(gaeFONL2->GetErrorYhigh(i+1)/gaeFONL2->GetY()[i+1],2))*(gaeSyst2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb));
		//gae_2Dvs2F_FONLErr->SetPointEYlow(i,sqrt(pow(gaeSyst2->GetErrorYlow(i)/gaeSyst2->GetY()[i],2)+pow(gaeFONL2->GetErrorYlow(i+1)/gaeFONL2->GetY()[i+1],2))*(gaeSyst2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb));

		gae_2Dvs2F_SystErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeSyst2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb);
		gae_2Dvs2F_SystErr->SetPointEYhigh(i,gaeSyst2->GetErrorYhigh(i)*(gaeSyst2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb));
		gae_2Dvs2F_SystErr->SetPointEYlow(i,gaeSyst2->GetErrorYlow(i)*(gaeSyst2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb));

		gae_2Dvs2F_StatErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeStat2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb);
		gae_2Dvs2F_StatErr->SetPointEYhigh(i,gaeStat2->GetErrorYhigh(i)*(gaeStat2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb));
		gae_2Dvs2F_StatErr->SetPointEYlow(i,gaeStat2->GetErrorYlow(i)*(gaeStat2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb));

		gae_2Dvs2F_FONLErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeStat2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb);
		gae_2Dvs2F_FONLErr->SetPointEYhigh(i,gaeFONL2->GetErrorYlow(i+1)/(gaeFONL2->GetY()[i+1]-gaeFONL2->GetErrorYlow(i+1))*(gaeStat2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb));
		gae_2Dvs2F_FONLErr->SetPointEYlow(i,gaeFONL2->GetErrorYhigh(i+1)/(gaeFONL2->GetY()[i+1]+gaeFONL2->GetErrorYhigh(i+1))*(gaeStat2->GetY()[i]/gaeFONL2->GetY()[i+1]*APb));

		//gae_7Dvs7F : error - 7D sys, 7F FONLL error

		//gae_7Dvs7F->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeSyst7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb);
		//gae_7Dvs7F->SetPointEYhigh(i,sqrt(pow(gaeSyst7->GetErrorYhigh(i+1)/gaeSyst7->GetY()[i+1],2)+pow(gaeFONL7->GetErrorYhigh(i+1)/gaeFONL7->GetY()[i+1],2))*(gaeSyst7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb));
		//gae_7Dvs7F->SetPointEYlow(i,sqrt(pow(gaeSyst7->GetErrorYlow(i+1)/gaeSyst7->GetY()[i+1],2)+pow(gaeFONL7->GetErrorYlow(i+1)/gaeFONL7->GetY()[i+1],2))*(gaeSyst7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb));

		gae_7Dvs7F_SystErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeSyst7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb);
		gae_7Dvs7F_SystErr->SetPointEYhigh(i,gaeSyst7->GetErrorYhigh(i+1)*(gaeSyst7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb));
		gae_7Dvs7F_SystErr->SetPointEYlow(i,gaeSyst7->GetErrorYlow(i+1)*(gaeSyst7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb));

		gae_7Dvs7F_StatErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeStat7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb);
		gae_7Dvs7F_StatErr->SetPointEYhigh(i,gaeStat7->GetErrorYhigh(i+1)*(gaeStat7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb));
		gae_7Dvs7F_StatErr->SetPointEYlow(i,gaeStat7->GetErrorYlow(i+1)*(gaeStat7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb));

		gae_7Dvs7F_FONLErr->SetPoint(i,gae_2Fvs7F->GetX()[i],gaeStat7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb);
		gae_7Dvs7F_FONLErr->SetPointEYhigh(i,gaeFONL7->GetErrorYlow(i+1)/(gaeFONL7->GetY()[i+1]-gaeFONL7->GetErrorYlow(i+1))*(gaeStat7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb));
		gae_7Dvs7F_FONLErr->SetPointEYlow(i,gaeFONL7->GetErrorYhigh(i+1)/(gaeFONL7->GetY()[i+1]+gaeFONL7->GetErrorYhigh(i+1))*(gaeStat7->GetY()[i+1]/gaeFONL7->GetY()[i+1]*APb));

	}

	gae_2Fvs7F->SetMarkerColor(kRed);
	gae_2Fvs7F->SetMarkerStyle(20);
	gae_2Fvs7F->SetMarkerSize(1);
	gae_2Fvs7F->SetFillColor(5);
	gae_2Fvs7F->SetFillStyle(1001);
	gae_2Fvs7F->SetLineColor(1);
	gae_2Fvs7F->SetLineWidth(2);

	gae_2Dvs7D_SystErr->SetLineColor(kAzure);
	gae_2Dvs7D_SystErr->SetMarkerColor(kAzure);
	gae_2Dvs7D_SystErr->SetMarkerStyle(21);
	gae_2Dvs7D_SystErr->SetMarkerSize(1);
	gae_2Dvs7D_SystErr->SetFillColor(0);
	gae_2Dvs7D_SystErr->SetFillStyle(0);

	gae_2Dvs7D_StatErr->SetLineColor(kAzure);
	gae_2Dvs7D_StatErr->SetMarkerColor(kAzure);
	gae_2Dvs7D_StatErr->SetMarkerStyle(21);
	gae_2Dvs7D_StatErr->SetMarkerSize(1);

	TCanvas* c1=new TCanvas("c1","",500,500);

	TLegend* leg=new TLegend(0.15,0.70,0.35,0.85);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.045);
	leg->AddEntry(gae_2Fvs7F,"2.76 TeV / 7 TeV of FONLL","pf");
	leg->AddEntry(gae_2Dvs7D_StatErr,"2.76 TeV / 7 TeV of pp data with Stat. Err.","lp");
	leg->AddEntry(gae_2Dvs7D_SystErr,"2.76 TeV / 7 TeV of pp data with Syst. Err.","f");



	TH2F* hempty = new TH2F("hempty",";B^{+} p_{T} (GeV/c);Ratio",12,5,70,10,0.0,1.0);
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();
	hempty->GetXaxis()->SetTitleOffset(1.0);
	hempty->GetYaxis()->SetTitleOffset(1.4);
	hempty->SetLineColor(0);
	hempty->GetXaxis()->SetTitleSize(0.045);
	hempty->GetYaxis()->SetTitleSize(0.045);
	hempty->GetXaxis()->SetTitleFont(42);
	hempty->GetYaxis()->SetTitleFont(42);
	hempty->GetXaxis()->SetLabelSize(0.040);
	hempty->GetYaxis()->SetLabelSize(0.040);
	hempty->GetXaxis()->SetLabelFont(42);
	hempty->GetYaxis()->SetLabelFont(42);

	hempty->Draw("");
	gae_2Fvs7F->Draw("same2p");
	gae_2Dvs7D_SystErr->Draw("same2e");

	gae_2Fvs7F->SetLineColor(kRed);
	gae_2Fvs7F->SetLineWidth(1);

	TGraphAsymmErrors* gae_2Fvs7F_v2 = (TGraphAsymmErrors*)gae_2Fvs7F->Clone("gae_2Fvs7F_v2");
	gae_2Fvs7F_v2->SetMarkerColor(kRed);
	gae_2Fvs7F_v2->SetMarkerStyle(20);
	gae_2Fvs7F_v2->SetMarkerSize(1);
	gae_2Fvs7F_v2->SetFillColor(0);
	gae_2Fvs7F_v2->SetFillStyle(0);
	gae_2Fvs7F_v2->SetLineColor(1);
	gae_2Fvs7F_v2->SetLineStyle(2);

	gae_2Fvs7F_v2->Draw("same2");

	gae_2Dvs7D_StatErr->Draw("samep");

	leg->Draw("");

	c1->SaveAs("CheckAccSys_20150203_re1.pdf");


	gae_2Dvs2F_StatErr->SetLineColor(kViolet);
	gae_2Dvs2F_StatErr->SetMarkerColor(kViolet);
	gae_2Dvs2F_StatErr->SetMarkerStyle(22);
	gae_2Dvs2F_StatErr->SetMarkerSize(1);
	gae_2Dvs2F_StatErr->SetLineWidth(1);

	gae_2Dvs2F_SystErr->SetLineColor(kViolet);
	gae_2Dvs2F_SystErr->SetMarkerColor(kViolet);
	gae_2Dvs2F_SystErr->SetMarkerStyle(22);
	gae_2Dvs2F_SystErr->SetMarkerSize(1);
	gae_2Dvs2F_SystErr->SetFillColor(0);
	gae_2Dvs2F_SystErr->SetFillStyle(0);
	gae_2Dvs2F_SystErr->SetLineWidth(2);


	gae_2Dvs2F_FONLErr->SetLineColor(kRed-9);
	gae_2Dvs2F_FONLErr->SetMarkerColor(kViolet);
	gae_2Dvs2F_FONLErr->SetMarkerStyle(22);
	gae_2Dvs2F_FONLErr->SetMarkerSize(1);
	gae_2Dvs2F_FONLErr->SetFillColor(kRed-9);
	gae_2Dvs2F_FONLErr->SetFillStyle(3244);//3444

	gae_7Dvs7F_StatErr->SetLineColor(kTeal+3);
	gae_7Dvs7F_StatErr->SetMarkerColor(kTeal+3);
	gae_7Dvs7F_StatErr->SetMarkerStyle(23);
	gae_7Dvs7F_StatErr->SetMarkerSize(1);
	gae_7Dvs7F_StatErr->SetLineWidth(1);


	gae_7Dvs7F_SystErr->SetLineColor(kTeal+3);
	gae_7Dvs7F_SystErr->SetMarkerColor(kTeal+3);
	gae_7Dvs7F_SystErr->SetMarkerStyle(22);
	gae_7Dvs7F_SystErr->SetMarkerSize(1);
	gae_7Dvs7F_SystErr->SetFillColor(0);
	gae_7Dvs7F_SystErr->SetFillStyle(0);
	gae_7Dvs7F_SystErr->SetLineWidth(2);

	gae_7Dvs7F_FONLErr->SetLineColor(kTeal-4);
	gae_7Dvs7F_FONLErr->SetMarkerColor(kTeal+3);
	gae_7Dvs7F_FONLErr->SetMarkerStyle(22);
	gae_7Dvs7F_FONLErr->SetMarkerSize(1);
	gae_7Dvs7F_FONLErr->SetFillColor(kTeal-4);
	gae_7Dvs7F_FONLErr->SetFillStyle(3209);//3409


//////////////////////////////////////////////////

	c1->Clear();
/*
	leg->DeleteEntry();
	leg->DeleteEntry();
	leg->DeleteEntry();
*/
	leg=new TLegend(0.15,0.65,0.35,0.90);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.045);

	leg->AddEntry(gae_2Dvs2F_StatErr,"Data / FONLL at 2.76 TeV with Stat. Err.","lp");
	leg->AddEntry(gae_2Dvs2F_SystErr,"Data / FONLL at 2.76 TeV with Syst. Err.","f");
	leg->AddEntry(gae_2Dvs2F_FONLErr,"Data / FONLL at 2.76 TeV from FONLL","f");
	leg->AddEntry(gae_7Dvs7F_StatErr,"Data / FONLL at 7 TeV with Stat. Err.","lp");
	leg->AddEntry(gae_7Dvs7F_SystErr,"Data / FONLL at 7 TeV with Syst. Err.","f");
	leg->AddEntry(gae_7Dvs7F_FONLErr,"Data / FONLL at 7 TeV from FONLL","f");

	hempty = new TH2F("hempty",";B^{+} p_{T} (GeV/c);Ratio",12,5,70,25,0.6,2.6);
	hempty->GetXaxis()->CenterTitle();
	hempty->GetYaxis()->CenterTitle();
	hempty->GetXaxis()->SetTitleOffset(1.0);
	hempty->GetYaxis()->SetTitleOffset(1.4);
	hempty->SetLineColor(0);
	hempty->GetXaxis()->SetTitleSize(0.045);
	hempty->GetYaxis()->SetTitleSize(0.045);
	hempty->GetXaxis()->SetTitleFont(42);
	hempty->GetYaxis()->SetTitleFont(42);
	hempty->GetXaxis()->SetLabelSize(0.040);
	hempty->GetYaxis()->SetLabelSize(0.040);
	hempty->GetXaxis()->SetLabelFont(42);
	hempty->GetYaxis()->SetLabelFont(42);

	hempty->Draw("");
	gae_2Dvs2F_FONLErr->Draw("same2");
	gae_7Dvs7F_FONLErr->Draw("same2");

	TGraphAsymmErrors* gae_2Dvs2F_FONLErr_v2 = (TGraphAsymmErrors*)gae_2Dvs2F_FONLErr->Clone("gae_2Dvs2F_FONLErr_v2");
	TGraphAsymmErrors* gae_7Dvs7F_FONLErr_v2 = (TGraphAsymmErrors*)gae_7Dvs7F_FONLErr->Clone("gae_7Dvs7F_FONLErr_v2");
	gae_2Dvs2F_FONLErr_v2->SetFillColor(0);
	gae_2Dvs2F_FONLErr_v2->SetFillStyle(0);
	gae_2Dvs2F_FONLErr_v2->SetLineWidth(2);
	gae_7Dvs7F_FONLErr_v2->SetFillColor(0);
	gae_7Dvs7F_FONLErr_v2->SetFillStyle(0);
	gae_7Dvs7F_FONLErr_v2->SetLineWidth(2);
	gae_2Dvs2F_FONLErr_v2->Draw("same2");
	gae_7Dvs7F_FONLErr_v2->Draw("same2");
	gae_2Dvs2F_SystErr->Draw("same2");
	gae_7Dvs7F_SystErr->Draw("same2");
	gae_2Dvs2F_StatErr->Draw("samep");
	gae_7Dvs7F_StatErr->Draw("samep");

	TLine* l1 = new TLine(5.0,1.0,70.0,1.0);
	l1->SetLineColor(kBlue);
	l1->SetLineStyle(2);
	l1->SetLineWidth(2);
	l1->Draw("");
	leg->Draw("");
	c1->SaveAs("CheckAccSys_20150203_re2.pdf");

	c1->cd();c1->Clear();
	gae_StavsFON_CMS->Draw("ape");c1->SaveAs("gae_StavsFON_CMS.pdf");
	c1->cd();c1->Clear();
	gae_SysvsFON_CMS->Draw("ape");c1->SaveAs("gae_SysvsFON_CMS.pdf");
	c1->cd();c1->Clear();
	gae_StavsFit_CMS->Draw("ape");c1->SaveAs("gae_StavsFit_CMS.pdf");
	c1->cd();c1->Clear();
	gae_SysvsFit_CMS->Draw("ape");c1->SaveAs("gae_SysvsFit_CMS.pdf");
	c1->cd();c1->Clear();
	gae_FitvsFON_CMS->Draw("ape");c1->SaveAs("gae_FitvsFON_CMS.pdf");

	c1->cd();c1->Clear();
	gae_StavsFON_ATL->Draw("ape");c1->SaveAs("gae_StavsFON_ATL.pdf");
	c1->cd();c1->Clear();
	gae_SysvsFON_ATL->Draw("ape");c1->SaveAs("gae_SysvsFON_ATL.pdf");
	c1->cd();c1->Clear();
	gae_StavsFit_ATL->Draw("ape");c1->SaveAs("gae_StavsFit_ATL.pdf");
	c1->cd();c1->Clear();
	gae_SysvsFit_ATL->Draw("ape");c1->SaveAs("gae_SysvsFit_ATL.pdf");
	c1->cd();c1->Clear();
	gae_FitvsFON_ATL->Draw("ape");c1->SaveAs("gae_FitvsFON_ATL.pdf");


	fout->cd();
	gaeStat2->Write();
	gaeSyst2->Write(); 
	gaeStat7->Write(); 
	gaeSyst7->Write(); 
	gaeFONL2->Write(); 
	gaeFONL7->Write(); 

	gae_2Fvs7F->Write(); 
	gae_2Dvs7D->Write(); 
	gae_2Dvs2F->Write(); 
	gae_7Dvs7F->Write(); 

	gaeFONL7_CMS->Write();
	gaeFONL7_ATL->Write();

	gaeStat7_CMS->Write();
	gaeSyst7_CMS->Write();
	gaeStat7_ATL->Write();
	gaeSyst7_ATL->Write();

	gaeFitT7_CMS->Write();
	gaeFitT7_ATL->Write();

	gae_StavsFON_CMS->Write();
	gae_StavsFON_ATL->Write();
	gae_SysvsFON_CMS->Write();
	gae_SysvsFON_ATL->Write();

	gae_StavsFit_CMS->Write();
	gae_StavsFit_ATL->Write();
	gae_SysvsFit_CMS->Write();
	gae_SysvsFit_ATL->Write();

	gae_FitvsFON_CMS->Write();
	gae_FitvsFON_ATL->Write();


	}
