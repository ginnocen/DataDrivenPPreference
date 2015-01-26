//### Code name :
//### Made by Hyunchul Kim at Jan. 16th. 2015
//### Purpose : Get the fitting function 

#include <iostream>
//#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

#define NUM 9	    // number of bins for CMS pp data
#define NUM_ATL 7   // number of bins for ATLAS pp data
#define NUM_BANA 5
#define APb 208

void ScaleCor_BplusY_v3()
{

    //gROOT->SetStyle("plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // CMS 7TeV, pp->(B+)X, unit : micro barn
    double ptbin[NUM+1]={-2.4,-1.8,-1.45,-1.10,-0.60,0.60,1.10,1.45,1.80,2.40};
    double sigmapt[NUM]={3.31,5.01,7.11,6.11,7.39,6.11,7.11,5.01,3.31};
    double sigmapt_sta[NUM]={0.42,0.55,0.69,0.64,0.65,0.64,0.69,0.55,0.42};
    double sigmapt_sys[NUM]={0.28,0.42,0.59,0.47,0.53,0.47,0.59,0.42,0.28};

    // ATLAS 7TeV, dsigma/dp_T (B+) * BR(B+->J/psiK+) * BR (J/Psi->mu+mu-), unit : pico barn
    double ptbin_ATL[NUM_ATL+1]={-2.25,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.25};
    double sigmapt_ATL[NUM_ATL]={132.0,143.9,142.7,153.7,142.7,143.9,132.0};
    double sigmapt_sta_ATL[NUM_ATL]={7.4,7.5,5.5,4.6,5.5,7.5,7.4};
    double sigmapt_sys_ATL[NUM_ATL]={10.8,9.8,8.5,10.1,8.5,9.8,10.8};

    // scale factor : 7TeV (10,60) / 7TeV (5,120)
    double scafac1[NUM]={0.240936,0.251162,0.256777,0.261614,0.264938,0.261707,0.256914,0.251345,0.241179};
    //double scafac1_higherrPerc[NUM]={0.0415639,0.0432057,0.0441331,0.0449673,0.0454538,0.0449755,0.044155,0.0432352,0.0416145};
    double scafac1_higherrPerc[NUM]={0.104572,0.106311,0.107527,0.108475,0.10918,0.108489,0.107558,0.106333,0.104632};
    double scafac1_lowerrPerc[NUM]={0.0467926,0.04699,0.04717,0.0474607,0.0475526,0.047478,0.0471655,0.0470032,0.0467807};

    // scale factor : 7TeV (9,120) / 7TeV (5,120)
    double scafac2[NUM]={0.315041,0.325914,0.331878,0.337008,0.340543,0.337106,0.332023,0.326108,0.315300};
    double scafac2_higherrPerc[NUM]={0.0382577,0.0396277,0.0403873,0.0411094,0.0415108,0.0411201,0.040406,0.03965,0.0382986};

/*
    // scale factor : 5TeV (10,60) / 7TeV (5,120)
    double scafac1[NUM]={0.115581,0.133784,0.145557,0.157723,0.174997,0.185269,0.187314,0.187551,0.185627};
    double scafac1_higherrPerc[NUM]={0.0774542,0.0767504,0.0757453,0.0744357,0.0709480,0.0664896,0.0640618,0.0618008,0.0584849};

    // scale factor : 5TeV (10,60) / 7TeV (9,120)
    double scafac2[NUM]={0.366874,0.410489,0.438586,0.468009,0.513877,0.549585,0.564160,0.575119,0.588730};
    double scafac2_higherrPerc[NUM]={0.0377523,0.0357077,0.0339855,0.0320103,0.0282639,0.0243676,0.0227370,0.0213062,0.0194416};
*/

    // scale factor : 7TeV (10,60) / 7TeV (9,120)
    double scafacATL[NUM_ATL]={0.767797,0.773887,0.776727,0.778075,0.776765,0.773962,0.767920};
    //double scafacATL_higherrPerc[NUM_ATL]={0.00331509,0.00358474,0.00373018,0.00379133,0.00373411,0.00359035,0.00331652};
    double scafacATL_higherrPerc[NUM_ATL]={0.0142932,0.0144848,0.0146981,0.0148175,0.0147026,0.0144899,0.0142972};
    double scafacATL_lowerrPerc[NUM_ATL]={0.00905752,0.00900733,0.00896227,0.0089069,0.00895989,0.00900924,0.00905448};


    TFile* fin = new TFile("../outputBplusyUnbinned_7TeV_rap30_pt1060.root");
    TGraphAsymmErrors* fonllref = (TGraphAsymmErrors*)fin->Get("gaeSigmaDecayv2Bplus");
    TFile* fin2 = new TFile("../outputBplusyUnbinned_5TeV_rap30_pt1060_Binned.root");
    TGraphAsymmErrors* fonllref2 = (TGraphAsymmErrors*)fin2->Get("gaeSigmaDecayv2Bplus");


/*
    // scale factor : 5TeV (10,60) / 7TeV (9,120)
    double scafacATL[NUM_ATL]={0.388793,0.440473,0.474273,0.514172,0.545919,0.563254,0.582096};
    double scafacATL_higherrPerc[NUM_ATL]={0.036585,0.0338129,0.0315913,0.0282863,0.0248507,0.0228102,0.0203317};
*/

    // TF1
    //TF1* f1fit_7TeV_pt1060_y30 = new TF1("f1fit_7TeV_pt1060_y30","((2.11951e+04)+(-1.44821e+03)*x*x)*(1.0e-06)/(6.09604e-5)",-3.0,3.0);
    TF1* f1fit_7TeV_pt1060_y30 = new TF1("f1fit_7TeV_pt1060_y30","((2.11913e+04)+(-1.44672e+03)*x*x)*(1.0e-06)/(6.09604e-5)",-3.0,3.0);

    TF1* fonllf = new TF1("fonllf","[0]+[1]*x*x",-3.0,3.0);
    fonllf->SetLineWidth(2);
    fonllf->SetLineStyle(2);
    fonllf->SetLineColor(kRed);


    fonllref->Fit(fonllf,"L","",-3.0,3.0);
    	  
    // Branching fraction for ATLAS data 
    double BRs;
    BRs = (1.027*0.001)*(5.961*0.01);	//BR(B+->J/psiK+) * BR (J/Psi->mu+mu-)

    // call histograms for CMS and ATLAS data
    TH1D* hsigmapt_noerr = new TH1D("hsigmapt_noerr","",NUM,ptbin);
    TH1D* hsigmapt_staerr = new TH1D("hsigmapt_staerr","",NUM,ptbin);
    TH1D* hsigmapt_syserr = new TH1D("hsigmapt_syserr","",NUM,ptbin);

    TH1D* hsigmapt_noerr_sf1 = new TH1D("hsigmapt_noerr_sf1","",NUM,ptbin);
    TH1D* hsigmapt_staerr_sf1 = new TH1D("hsigmapt_staerr_sf1","",NUM,ptbin);
    TH1D* hsigmapt_syserr_sf1 = new TH1D("hsigmapt_syserr_sf1","",NUM,ptbin);

    TH1D* hsigmapt_noerr_sf2 = new TH1D("hsigmapt_noerr_sf2","",NUM,ptbin);
    TH1D* hsigmapt_staerr_sf2 = new TH1D("hsigmapt_staerr_sf2","",NUM,ptbin);
    TH1D* hsigmapt_syserr_sf2 = new TH1D("hsigmapt_syserr_sf2","",NUM,ptbin);

    TH1D* hsigmapt_noerr_sfATL = new TH1D("hsigmapt_noerr_sfATL","",NUM_ATL,ptbin_ATL);
    TH1D* hsigmapt_staerr_sfATL = new TH1D("hsigmapt_staerr_sfATL","",NUM_ATL,ptbin_ATL);
    TH1D* hsigmapt_syserr_sfATL = new TH1D("hsigmapt_syserr_sfATL","",NUM_ATL,ptbin_ATL);

    TH1D* hsigmapt_noerr_ATL = new TH1D("hsigmapt_noerr_ATL","",NUM_ATL,ptbin_ATL);
    TH1D* hsigmapt_staerr_ATL = new TH1D("hsigmapt_staerr_ATL","",NUM_ATL,ptbin_ATL);
    TH1D* hsigmapt_syserr_ATL = new TH1D("hsigmapt_syserr_ATL","",NUM_ATL,ptbin_ATL);

    // store the corrected data in the histogram for combining CMS and ATLAS together
	std::cout << "################ SCALED CMS DATA ########################" << std::endl;
    for(int i=0;i<NUM;i++){
	hsigmapt_noerr->SetBinContent(i+1,sigmapt[i]*APb);
	hsigmapt_staerr->SetBinContent(i+1,sigmapt[i]*APb);
	hsigmapt_syserr->SetBinContent(i+1,sigmapt[i]*APb);
	hsigmapt_staerr->SetBinError(i+1,sigmapt_sta[i]*APb);
	hsigmapt_syserr->SetBinError(i+1,sigmapt_sys[i]*APb);

	hsigmapt_noerr_sf1->SetBinContent(i+1,sigmapt[i]*APb*scafac1[i]);
	hsigmapt_staerr_sf1->SetBinContent(i+1,sigmapt[i]*APb*scafac1[i]);
	hsigmapt_syserr_sf1->SetBinContent(i+1,sigmapt[i]*APb*scafac1[i]);
	hsigmapt_staerr_sf1->SetBinError(i+1,sqrt(pow(sigmapt_sta[i]/sigmapt[i],2)+pow(scafac1_higherrPerc[i],2))*(sigmapt[i]*APb*scafac1[i]));
	hsigmapt_syserr_sf1->SetBinError(i+1,sqrt(pow(sigmapt_sys[i]/sigmapt[i],2)+pow(scafac1_higherrPerc[i],2))*(sigmapt[i]*APb*scafac1[i]));
   
	std::cout << i << " : " << hsigmapt_noerr_sf1->GetBinContent(i+1) << " + " << hsigmapt_staerr_sf1->GetBinError(i+1) << " (sta.) + " << hsigmapt_syserr_sf1->GetBinError(i+1) << " (sys.)" << std::endl;
    	hsigmapt_noerr_sf2->SetBinContent(i+1,sigmapt[i]*APb*scafac2[i]);
	hsigmapt_staerr_sf2->SetBinContent(i+1,sigmapt[i]*APb*scafac2[i]);
	hsigmapt_syserr_sf2->SetBinContent(i+1,sigmapt[i]*APb*scafac2[i]);
	hsigmapt_staerr_sf2->SetBinError(i+1,sqrt(pow(sigmapt_sta[i]/sigmapt[i],2)+pow(scafac2_higherrPerc[i],2))*(sigmapt[i]*APb*scafac2[i]));
	hsigmapt_syserr_sf2->SetBinError(i+1,sqrt(pow(sigmapt_sys[i]/sigmapt[i],2)+pow(scafac2_higherrPerc[i],2))*(sigmapt[i]*APb*scafac2[i]));
    }

	std::cout << "################ SCALED ATLAS DATA ########################" << std::endl;

    for(int i=0;i<NUM_ATL;i++){
	hsigmapt_noerr_ATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001);
	hsigmapt_staerr_ATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001);
	hsigmapt_syserr_ATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001);
	hsigmapt_staerr_ATL->SetBinError(i+1,sigmapt_sta_ATL[i]*APb/BRs*0.000001);
	hsigmapt_syserr_ATL->SetBinError(i+1,sigmapt_sys_ATL[i]*APb/BRs*0.000001);

    	hsigmapt_noerr_sfATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001*scafacATL[i]);
	hsigmapt_staerr_sfATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001*scafacATL[i]);
	hsigmapt_syserr_sfATL->SetBinContent(i+1,sigmapt_ATL[i]*APb/BRs*0.000001*scafacATL[i]);
	hsigmapt_staerr_sfATL->SetBinError(i+1,sqrt(pow(sigmapt_sta_ATL[i]/sigmapt_ATL[i],2)+pow(scafac2_higherrPerc[i],2))*(sigmapt_ATL[i]*APb/BRs*0.000001*scafacATL[i]));
	hsigmapt_syserr_sfATL->SetBinError(i+1,sqrt(pow(sigmapt_sys_ATL[i]/sigmapt_ATL[i],2)+pow(scafac2_higherrPerc[i],2))*(sigmapt_ATL[i]*APb/BRs*0.000001*scafacATL[i]));
    	std::cout << i << " : " << hsigmapt_noerr_sfATL->GetBinContent(i+1) << " + " << hsigmapt_staerr_sfATL->GetBinError(i+1) << " (sta.) + " << hsigmapt_syserr_sfATL->GetBinError(i+1) << " (sys.)" << std::endl;
    	
    }

    // call for fitting function and set the related option
    //###TF1* fitft = new TF1("fitft","pow(10,[0]*exp([1]+[2]*x)+[3])",0.0,120.0);
    TF1* fitft = new TF1("fitft","pol2",-3.0,3.0);
    //TF1* fitft = new TF1("fitft","expo",0.0,120.0);
    //fitft->SetParLimits(2,-0.033,0.00);
    fitft->SetLineColor(kAzure-6);

    //###TF1* fitft_old = new TF1("fitft_old","pow(10,(8.00042e-01)*exp((2.07267e+00)+(-2.42747e-02)*x)+(-2.36680e+00))",0.0,120.0);//60.0
    TF1* fitft_old = new TF1("fitft_old","pol2",-3.0,3.0);
    fitft_old->SetLineColor(kRed);

    //###TF1* fitft_c = new TF1("fitft_c","pow(10,(4.00284e+00)*exp((4.62641e-01)+(-2.42728e-02)*x)+(-2.36722e+00))",0.0,120.0);//60.0
    //###TF1* fitft_c = new TF1("fitft_c","pol2",-2.5,2.5);
    //TF1* fitft_c = new TF1("fitft_c","pow(10,(5.33672e+00)*exp((1.74899e-01)+(-2.42769e-02)*x)+(-2.36635e+00))",0.0,120.0);//with SetParLimit
/*
    fitft_c->SetLineColor(kViolet-4);
    fitft_c->SetLineStyle(2);
*/
    hsigmapt_staerr->Fit("fitft_old","L","",-3.0,3.0);

    std::cout << "*** Fitting with 7TeV pp ATLAS " << std::endl;
    //###hsigmapt_staerr_ATL->Fit("fitft","L","",0.0,120.0);
    hsigmapt_staerr_ATL->Fit("fitft","L","",-3.0,3.0);
    std::cout << std::endl;
    std::cout << "*** Fitting with 7TeV pp CMS " << std::endl;
    //###hsigmapt_staerr->Fit("fitft","L","",0.0,120.0);
/*
    hsigmapt_staerr_ATL->Fit("fitft_c","L","",-2.4,2.4);
    std::cout << std::endl;
    std::cout << "*** Fitting with 7TeV pp CMS " << std::endl;
    hsigmapt_staerr->Fit("fitft_c","L","",-2.4,2.4);
*/

    TCanvas *cSigma = new TCanvas("cSigma","",500,500);

    cSigma->SetFillColor(0);
    cSigma->SetBorderMode(0);
    cSigma->SetBorderSize(2);
    cSigma->SetFrameBorderMode(0);
    cSigma->SetLeftMargin(0.16);
    cSigma->SetRightMargin(0.04);//0.02
    cSigma->SetTopMargin(0.08);
    cSigma->SetBottomMargin(0.15);
    //cSigma->SetLogy(1);

    //###TH2D* hempty = new TH2D("hempty",";B^{+} p_{T} (GeV/c);d#sigma/dp_{T} A",13,0.0,130.0,7,0.01,10000.0);
    TH2D* hempty = new TH2D("hempty",";B^{+} y_{lab};d#sigma/dy A (#mub)",60,-3.0,3.0,35,0,3500.0);

    //hempty->SetAxisRange(7.5, 42.5, "X");
    //hempty->SetAxisRange(0.1, 1000, "Y");
    hempty->GetXaxis()->CenterTitle();
    hempty->GetYaxis()->CenterTitle();
    hempty->GetXaxis()->SetTitleOffset(1.0);//###1.0
    hempty->GetYaxis()->SetTitleOffset(1.2);//###1.0
    hempty->GetXaxis()->SetTitleSize(0.065);//###0.070
    hempty->GetYaxis()->SetTitleSize(0.065);//###0.070
    hempty->GetXaxis()->SetTitleFont(42);
    hempty->GetYaxis()->SetTitleFont(42);
    hempty->GetXaxis()->SetLabelFont(42);
    hempty->GetYaxis()->SetLabelFont(42);
    hempty->GetXaxis()->SetLabelSize(0.055);//###0.060
    hempty->GetYaxis()->SetLabelSize(0.050);//###0.060    

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

    hsigmapt_staerr_sf1->SetMarkerColor(kViolet);
    hsigmapt_staerr_sf1->SetMarkerStyle(20);
    hsigmapt_staerr_sf1->SetMarkerSize(1);
    hsigmapt_staerr_sf1->SetLineWidth(2);
    hsigmapt_staerr_sf1->SetLineColor(kViolet);

    hsigmapt_syserr_sf1->SetLineColor(kViolet);
    hsigmapt_syserr_sf1->SetFillColor(kViolet-4);
    hsigmapt_syserr_sf1->SetMarkerSize(0);
    hsigmapt_syserr_sf1->SetFillStyle(3004);

    hsigmapt_staerr_sf2->SetMarkerColor(kGreen+4);
    hsigmapt_staerr_sf2->SetMarkerStyle(26);
    hsigmapt_staerr_sf2->SetMarkerSize(1);
    hsigmapt_staerr_sf2->SetLineWidth(2);
    hsigmapt_staerr_sf2->SetLineColor(kGreen+4);

    hsigmapt_syserr_sf2->SetLineColor(kGreen+4);
    hsigmapt_syserr_sf2->SetFillColor(kGreen-8);
    hsigmapt_syserr_sf2->SetMarkerSize(0);

    hsigmapt_staerr_sfATL->SetMarkerColor(kGreen+4);
    hsigmapt_staerr_sfATL->SetMarkerStyle(21);
    hsigmapt_staerr_sfATL->SetMarkerSize(1);
    hsigmapt_staerr_sfATL->SetLineWidth(2);
    hsigmapt_staerr_sfATL->SetLineColor(kGreen+4);

    hsigmapt_syserr_sfATL->SetLineColor(kGreen+4);
    hsigmapt_syserr_sfATL->SetFillColor(kGreen-8);
    hsigmapt_syserr_sfATL->SetMarkerSize(0);
    hsigmapt_syserr_sfATL->SetFillStyle(3005);

    //TF1* fitft_sf1 = new TF1("fitft_sf1","pol2",-3.0,3.0);
    TF1* fitft_sf1 = new TF1("fitft_sf1","[0]+[1]*x*x",-3.0,3.0);


    fitft_sf1->SetLineColor(kViolet);
    std::cout << std::endl << std::endl;
    std::cout << "### Fit parameters with CMS data #########################" << std::endl;
    hsigmapt_staerr_sf1->Fit("fitft_sf1","L","",-3.0,3.0);
/*
    TF1* fitft_sf2 = new TF1("fitft_sf2","pol2",-2.5,2.5);
    fitft_sf2->SetLineColor(kGreen+4);
    hsigmapt_staerr_sf2->Fit("fitft_sf2","L","",-2.4,2.4);
*/
    //TF1* fitft_sfATL = new TF1("fitft_sfATL","pol2",-3.0,3.0);
    TF1* fitft_sfATL = new TF1("fitft_sfATL","[0]+[1]*x*x",-3.0,3.0);

    fitft_sfATL->SetLineColor(kGreen+4);
    std::cout << std::endl << std::endl;
    std::cout << "### Fit parameters with ALTAS data #########################" << std::endl;
    hsigmapt_staerr_sfATL->Fit("fitft_sfATL","L","",-3.0,3.0);

    TH1F* hsigmapt_syserr_sf1_v2 = (TH1F*)hsigmapt_syserr_sf1->Clone();
    TH1F* hsigmapt_syserr_sfATL_v2 = (TH1F*)hsigmapt_syserr_sfATL->Clone();

    hsigmapt_syserr_sf1_v2->SetFillStyle(0);
    hsigmapt_syserr_sfATL_v2->SetFillStyle(0);

    hempty->Draw("");
    hsigmapt_syserr->Draw("samee2");
    hsigmapt_syserr_ATL->Draw("samee2");
    hsigmapt_staerr->Draw("samee");
    hsigmapt_staerr_ATL->Draw("samee");
    fitft_old->Draw("same");
    //fitft_c->Draw("same");
    hsigmapt_syserr_sf1->Draw("samee2");
    hsigmapt_syserr_sf1_v2->Draw("samee2");
    hsigmapt_staerr_sf1->Draw("samee");
    //hsigmapt_syserr_sf2->Draw("samee2");
    //hsigmapt_staerr_sf2->Draw("samee");
    hsigmapt_syserr_sfATL->Draw("samee2");
    hsigmapt_syserr_sfATL_v2->Draw("samee2");
    hsigmapt_staerr_sfATL->Draw("samee");
    //###f1fit_7TeV_pt1060_y30->Draw("same");
     fonllf->Draw("same");
   
    //###TLegend* leg = new TLegend(0.30,0.65,0.60,0.90,"");
    TLegend* leg = new TLegend(0.20,0.55,0.50,0.90,"");


    leg->SetBorderSize(0);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    hsigmapt_staerr->SetFillColor(5);
    hsigmapt_staerr_ATL->SetFillColor(kTeal);//kGreen
 
    TLegendEntry* ent_ppCMS7=leg->AddEntry(hsigmapt_staerr,"CMS 7 TeV, p_{T} > 5 GeV/c","lpf");

    ent_ppCMS7->SetTextFont(42);
    ent_ppCMS7->SetFillStyle(1001);
    ent_ppCMS7->SetFillColor(5);
    ent_ppCMS7->SetMarkerColor(kRed);
    ent_ppCMS7->SetMarkerStyle(24);
    ent_ppCMS7->SetMarkerSize(15);
    ent_ppCMS7->SetLineWidth(2);
    ent_ppCMS7->SetLineColor(kRed);
   
    TLegendEntry* ent_ppATL7=leg->AddEntry(hsigmapt_staerr_ATL,"ATLAS 7 TeV, 9 <  p_{T} <  120 GeV/c","lpf");

    ent_ppATL7->SetTextFont(42);
    ent_ppATL7->SetFillColor(kTeal);
    ent_ppATL7->SetMarkerColor(kTeal+3);
    ent_ppATL7->SetMarkerStyle(25);
    ent_ppATL7->SetMarkerSize(15);
    ent_ppATL7->SetLineWidth(2);
    ent_ppATL7->SetLineColor(kTeal+3);

    TLegendEntry* ent_ppCMS7sf1=leg->AddEntry(hsigmapt_staerr_sf1,"CMS 7 TeV, Scaled","pf");
/*
    ent_ppCMS7sf1->SetTextFont(42);
    ent_ppCMS7sf1->SetFillStyle(3004);
    ent_ppCMS7sf1->SetFillColor(kViolet-4);
    ent_ppCMS7sf1->SetMarkerColor(kViolet);
    ent_ppCMS7sf1->SetMarkerStyle(21);
    ent_ppCMS7sf1->SetMarkerSize(15);
    ent_ppCMS7sf1->SetLineWidth(2);
    ent_ppCMS7sf1->SetLineColor(kViolet);
*/   
    TLegendEntry* ent_ppATL7sf=leg->AddEntry(hsigmapt_staerr_sfATL,"ATLAS 7 TeV, scaled","pf");
/*
    ent_ppATL7sf->SetTextFont(42);
    ent_ppATL7sf->SetFillColor(kTeal);
    ent_ppATL7sf->SetMarkerColor(kTeal+3);
    ent_ppATL7sf->SetMarkerStyle(25);
    ent_ppATL7sf->SetMarkerSize(15);
    ent_ppATL7sf->SetLineWidth(2);
    ent_ppATL7sf->SetLineColor(kTeal+3);
*/

    TLegendEntry* ent_fitold=leg->AddEntry(fitft_old,"(1) Fit with CMS data","l");
    TLegendEntry* ent_fitnew=leg->AddEntry(fitft,"(2) Fit with ATLAS data","l");
    //TLegendEntry* ent_fitc=leg->AddEntry(fitft_c,"(3) Fit with CMS data based on (2)","l");
    TLegendEntry* ent_fitsf1=leg->AddEntry(fitft_sf1,"(3) Fit with scaled CMS data, p_{T} (10,60)","l");
    TLegendEntry* ent_fitsfA=leg->AddEntry(fitft_sfATL,"(4) Fit with scaled ATLAS data, p_{T} (10,60)","l");
    TLegendEntry* ent_fitfon=leg->AddEntry(fonllf,"(5) Fit with 7TeV FONLL, p_{T} (10,60)","l");




    leg->Draw("");

    //###cSigma->SaveAs("cSigma_AddATLAS.pdf");
    //cSigma->SaveAs("cSigma_AddATLAS_wParLimit.pdf");//With ParLimit
    //#####cSigma->SaveAs("cSigma_BplusY_v2_5vs7.pdf");
    cSigma->SaveAs("cSigma_BplusY_v2.pdf");


    hempty = new TH2D("hempty",";B^{+} y_{lab};d#sigma/dy A (#mub)",60,-3.0,3.0,70,0,700.0);
    hempty->GetXaxis()->CenterTitle();
    hempty->GetYaxis()->CenterTitle();
    hempty->GetXaxis()->SetTitleOffset(1.0);//###1.0
    hempty->GetYaxis()->SetTitleOffset(1.2);//###1.0
    hempty->GetXaxis()->SetTitleSize(0.065);//###0.070
    hempty->GetYaxis()->SetTitleSize(0.065);//###0.070
    hempty->GetXaxis()->SetTitleFont(42);
    hempty->GetYaxis()->SetTitleFont(42);
    hempty->GetXaxis()->SetLabelFont(42);
    hempty->GetYaxis()->SetLabelFont(42);
    hempty->GetXaxis()->SetLabelSize(0.055);//###0.060
    hempty->GetYaxis()->SetLabelSize(0.050);//###0.060    

    hempty->Draw("");

    //fonllref->SetLineColor(kAzure+4);
    fonllref->SetFillColor(kYellow-9);
    fonllref->SetLineColor(kYellow-9);
    fonllref->SetLineWidth(0);
    //fonllref->SetFillColor(0);
    fonllref->SetFillStyle(1001);
    fonllref->Draw("samee3");
    fonllf->Draw("same");

    fitft_sf1->SetLineWidth(2);
    fitft_sfATL->SetLineWidth(2);

    hsigmapt_syserr_sf1->Draw("samee2");
    hsigmapt_syserr_sf1_v2->Draw("samee2");
    hsigmapt_staerr_sf1->Draw("samee");
    hsigmapt_syserr_sfATL->Draw("samee2");
    hsigmapt_syserr_sfATL_v2->Draw("samee2");
    hsigmapt_staerr_sfATL->Draw("samee");
    fonllf->Draw("same");

    //f1fit_7TeV_pt1060_y30->Draw("same");
    //#####cSigma->SaveAs("cSigma_BplusY_v2_fitfts_5vs7.pdf");
    
    leg->Clear();
    leg = new TLegend(0.20,0.67,0.50,0.91,"");
    leg->SetBorderSize(0);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    leg->AddEntry(hsigmapt_staerr_sf1,"CMS 7 TeV, scaled","pf");
    leg->AddEntry(hsigmapt_staerr_sfATL,"ATLAS 7 TeV, scaled","pf");
    leg->AddEntry(fonllref,"FONLL assumption, p_{T} (10,60)","f");
    leg->AddEntry(fitft_sf1,"(3) Fit with scaled CMS data, p_{T} (10,60)","l");
    leg->AddEntry(fitft_sfATL,"(4) Fit with scaled ATLAS data, p_{T} (10,60)","l");
    leg->AddEntry(fonllf,"(5) Fit with 7TeV FONLL, p_{T} (10,60)","l");

    leg->Draw(); 
    cSigma->SaveAs("cSigma_BplusY_v2re_fitfts.pdf");

    //double ptbinpPb[NUM_BANA+1]={-2.865,-1.935,-1.0,0.0,1.0,1.935};
    double ptbinpPb[NUM_BANA+1]={-2.4,-1.470,-0.535,0.465,1.465,2.4};



/*
(-2.865 , -1.935) : 
----------- 93.2602 312.341 209.217
(-1.935 , -1) : 
----------- 291.061 363.435 294.787
(-1 , 0) : 
----------- 394.912 390.261 339.714
(0 , 1) : 
----------- 394.912 390.261 339.714
(1 , 1.935) : 
----------- 291.061 363.435 294.787
*/

    double rat5vs7[NUM_BANA]={0.614124,0.661836,0.673195,0.672551,0.647923};
    double rat5vs7_higherr[NUM_BANA]={0.0267497,0.0263131,0.0261087,0.0261546,0.0263078};
    double rat5vs7_lowerr[NUM_BANA]={0.0192216,0.0182598,0.0176511,0.0176945,0.0182942};

    double CMS_sc_sysPerc[NUM_BANA]={0.216967456,0.202607258,0.195379494,0.195385517,0.202904093};
    double ATL_sc_sysPerc[NUM_BANA]={0.120591615,0.107314368,0.091820064,0.092316837,0.108632875};

    TH1D* hbin_CMS7=new TH1D("hbin_CMS7","",NUM_BANA,ptbinpPb);
    TH1D* hbin_ATL7=new TH1D("hbin_ATL7","",NUM_BANA,ptbinpPb);
    TH1D* hbin_fonl=new TH1D("hbin_fonl","",NUM_BANA,ptbinpPb);

    double xfac=0.465;
    for (int i=0;i<NUM_BANA;i++){
	double integ3,integ4,integ5;
	integ3=fitft_sf1->Integral(ptbinpPb[i]-xfac,ptbinpPb[i+1]-xfac)/(ptbinpPb[i+1]-ptbinpPb[i]);
	integ4=fitft_sfATL->Integral(ptbinpPb[i]-xfac,ptbinpPb[i+1]-xfac)/(ptbinpPb[i+1]-ptbinpPb[i]);
	integ5=fonllf->Integral(ptbinpPb[i]-xfac,ptbinpPb[i+1]-xfac)/(ptbinpPb[i+1]-ptbinpPb[i]);

	hbin_CMS7->SetBinContent(i+1,integ3*rat5vs7[i]);
	hbin_ATL7->SetBinContent(i+1,integ4*rat5vs7[i]);
	hbin_fonl->SetBinContent(i+1,integ5*rat5vs7[i]);

	hbin_CMS7->SetBinError(i+1,sqrt(pow(CMS_sc_sysPerc[i],2)+pow(rat5vs7_higherr[i],2))*(integ3*rat5vs7[i]));
	hbin_ATL7->SetBinError(i+1,sqrt(pow(ATL_sc_sysPerc[i],2)+pow(rat5vs7_higherr[i],2))*(integ4*rat5vs7[i]));

	std::cout << "(" << ptbinpPb[i] << " , " << ptbinpPb[i+1] << ") : " << std::endl;
	std::cout << "----------- " << integ3 << " " << integ4 << " " << integ5 << std::endl;
	std::cout << "----------------- " << hbin_CMS7->GetBinContent(i+1) << " " << hbin_ATL7->GetBinContent(i+1) << " " << hbin_fonl->GetBinContent(i+1) << std::endl;


	
    }

///////////////////////////////////////////////////////////////////////////
    hempty = new TH2D("hempty",";B^{+} y_{lab};d#sigma/dy A (#mub)",48,-2.4,2.4,50,0,500.0);
    hempty->GetXaxis()->CenterTitle();
    hempty->GetYaxis()->CenterTitle();
    hempty->GetXaxis()->SetTitleOffset(1.0);//###1.0
    hempty->GetYaxis()->SetTitleOffset(1.2);//###1.0
    hempty->GetXaxis()->SetTitleSize(0.065);//###0.070
    hempty->GetYaxis()->SetTitleSize(0.065);//###0.070
    hempty->GetXaxis()->SetTitleFont(42);
    hempty->GetYaxis()->SetTitleFont(42);
    hempty->GetXaxis()->SetLabelFont(42);
    hempty->GetYaxis()->SetLabelFont(42);
    hempty->GetXaxis()->SetLabelSize(0.055);//###0.060
    hempty->GetYaxis()->SetLabelSize(0.050);//###0.060    

    gStyle->SetHatchesLineWidth(2.0);

    hbin_CMS7->SetLineColor(kViolet);
    hbin_ATL7->SetLineColor(kGreen+3);
    hbin_fonl->SetLineColor(kRed);

    hbin_CMS7->SetLineWidth(4);
    hbin_ATL7->SetLineWidth(4);
    hbin_fonl->SetLineWidth(2);
    fonllref2->SetLineWidth(4.5);

    hbin_CMS7->SetLineStyle(2);

    hbin_CMS7->SetFillColor(kViolet-9);
    hbin_ATL7->SetFillColor(kGreen-10);
    fonllref2->SetFillColor(5);

    //hbin_CMS7->SetFillStyle(4104);
    //hbin_ATL7->SetFillStyle(4105);
    fonllref2->SetFillStyle(1001);

    hbin_CMS7->SetMarkerStyle(20);
    hbin_ATL7->SetMarkerStyle(21);
    fonllref2->SetMarkerStyle(26);//22

    hbin_CMS7->SetMarkerColor(kViolet);
    hbin_ATL7->SetMarkerColor(kGreen+3);
    fonllref2->SetMarkerColor(kBlack);

    hbin_CMS7->SetMarkerSize(1.4);
    hbin_ATL7->SetMarkerSize(1.5);
    fonllref2->SetMarkerSize(2);

    hempty->Draw("");
    fonllref2->Draw("samee2");
    hbin_ATL7->Draw("samee2");
    hbin_CMS7->Draw("samee2");
    fonllref2->Draw("sameep");
    hbin_ATL7->Draw("samee");
    hbin_CMS7->Draw("samee");

    //hbin_fonl->Draw("same");
    leg->Clear();
    leg = new TLegend(0.20,0.67,0.50,0.91,"");
    leg->SetBorderSize(0);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    leg->AddEntry(hbin_CMS7,"scaled 5 TeV pp, boosted with CMS data","pf");
    leg->AddEntry(hbin_ATL7,"scaled 5 TeV pp, boosted with ATLAS data","pf");
    //leg->AddEntry(fonllref2,"5TeV pp, boosted FONLL assumption, p_{T} (10,60)","pf");
    leg->AddEntry(fonllref2,"5TeV pp, boosted FONLL assumption","pf");
    leg->Draw(); 

    TLatex* latexpt=new TLatex(0.55,0.23,"10 < p_{T} < 60 GeV/c");//for only check
    latexpt->SetNDC();
    latexpt->SetTextColor(1);
    latexpt->SetTextFont(42);//42
    latexpt->SetTextSize(0.050);//0.045
    latexpt->Draw();

    cSigma->SaveAs("cSigma_BplusY_v2re_fitfts2.pdf");
    cSigma->Clear();
    fonllref2->Draw("");

    cSigma->SaveAs("cSigma_BplusY_v2re_fitfts3.pdf");

    TFile* fout = new TFile("Compdsigmady_ScaleCor_v3.root","RECREATE");
    fout->cd();
    fonllref2->Write();
    hbin_ATL7->Write();
    hbin_CMS7->Write();

}
