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

#define NUM 5	    // number of bins for CMS pp data
#define NUM_ATL 8   // number of bins for ATLAS pp data
#define APb 208

void ScaleCor_pt()
{

    //gROOT->SetStyle("plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

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
    TF1* fitft = new TF1("fitft","pow(10,[0]*exp([1]+[2]*x)+[3])",0.0,120.0);
    //TF1* fitft = new TF1("fitft","expo",0.0,120.0);
    //fitft->SetParLimits(2,-0.033,0.00);
    fitft->SetLineColor(kAzure-6);

    TF1* fitft_old = new TF1("fitft_old","pow(10,(8.00042e-01)*exp((2.07267e+00)+(-2.42747e-02)*x)+(-2.36680e+00))",0.0,120.0);//60.0
    fitft_old->SetLineColor(kRed);

    TF1* fitft_c = new TF1("fitft_c","pow(10,(4.00284e+00)*exp((4.62641e-01)+(-2.42728e-02)*x)+(-2.36722e+00))",0.0,120.0);//60.0
    //TF1* fitft_c = new TF1("fitft_c","pow(10,(5.33672e+00)*exp((1.74899e-01)+(-2.42769e-02)*x)+(-2.36635e+00))",0.0,120.0);//with SetParLimit


    fitft_c->SetLineColor(kOrange);
    fitft_c->SetLineStyle(2);

    std::cout << "*** Fitting with 7TeV pp CMS " << std::endl;
    hsigmapt_staerr->Fit("fitft","L","",0.0,120.0);

    //#####TF1* fitft = new TF1("fitft","pow(10,[0]*exp([1]+[2]*x)+[3])",0.0,120.0);


    std::cout << std::endl;
// CMS iteration
    //double ptbin[NUM+1]={5.,10.,13.,17.,24.,30.};

    double integ1,integ2,gapinteg,mininteg,minpt;
    double stminpt[NUM];
    for (int i=0;i<NUM;i++){
    mininteg=100000000000.0;
	for (double j=ptbin[i]+1.0;j<ptbin[i+1];j+=1.0){
	integ1 = fitft->Integral(ptbin[i],j);
	integ2 = fitft->Integral(j,ptbin[i+1]);
	gapinteg = fabs(integ1-integ2);
	std:: cout << j << " : " << integ1 << " , " << integ2 << std::endl;
	if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
	}
    stminpt[i]=minpt;
    std::cout << i << " : " << "gap : " << mininteg << ", minpt : " << minpt << std::endl;
	for (double j=stminpt[i]-0.9;j<stminpt[i]+1.0;j+=0.1){
	integ1 = fitft->Integral(ptbin[i],j);
	integ2 = fitft->Integral(j,ptbin[i+1]);
	gapinteg = fabs(integ1-integ2);
	std:: cout << j << " : " << "----- " << integ1 << " , " << integ2 << std::endl;

	if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
	}
    stminpt[i]=minpt;

    std::cout << i << " : " << "--------- gap : " << mininteg << ", minpt : " << minpt << std::endl;
    }



    std::cout << "*** Fitting with 7TeV pp ATLAS " << std::endl;
    hsigmapt_staerr_ATL->Fit("fitft","L","",0.0,120.0);

// ALTAS iteration
    //double ptbin_ATL[NUM_ATL+1]={9.,13.,16.,20.,25.,35.,50.,70.,120.};
    double stminpt_ATL[NUM_ATL];
    for (int i=0;i<NUM_ATL;i++){
    mininteg=100000000000.0;
	for (double j=ptbin_ATL[i]+1.0;j<ptbin_ATL[i+1];j+=1.0){
	integ1 = fitft->Integral(ptbin_ATL[i],j);
	integ2 = fitft->Integral(j,ptbin_ATL[i+1]);
	gapinteg = fabs(integ1-integ2);
	std:: cout << j << " : " << integ1 << " , " << integ2 << std::endl;
	if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
	}
    stminpt_ATL[i]=minpt;
    std::cout << i << " : " << "gap : " << mininteg << ", minpt : " << minpt << std::endl;
	for (double j=stminpt_ATL[i]-0.9;j<stminpt_ATL[i]+1.0;j+=0.1){
	integ1 = fitft->Integral(ptbin_ATL[i],j);
	integ2 = fitft->Integral(j,ptbin_ATL[i+1]);
	gapinteg = fabs(integ1-integ2);
	std:: cout << j << " : " << "----- " << integ1 << " , " << integ2 << std::endl;

	if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
	}
    stminpt_ATL[i]=minpt;
    std::cout << i << " : " << "gap : " << mininteg << ", minpt : " << minpt << std::endl;
	for (double j=stminpt_ATL[i]-0.09;j<stminpt_ATL[i]+0.1;j+=0.01){
	integ1 = fitft->Integral(ptbin_ATL[i],j);
	integ2 = fitft->Integral(j,ptbin_ATL[i+1]);
	gapinteg = fabs(integ1-integ2);
	std:: cout << j << " : " << "-------- " << integ1 << " , " << integ2 << std::endl;

	if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
	}
    stminpt_ATL[i]=minpt;

    std::cout << i << " : " << "--------- gap : " << mininteg << ", minpt : " << minpt << std::endl;
    }










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

    TH2D* hempty = new TH2D("hempty",";B^{+} p_{T} (GeV/c);d#sigma/dp_{T} A",13,0.0,130.0,7,0.01,10000.0);
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
    fitft_old->Draw("same");
    fitft_c->Draw("same");
   
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

    TLegendEntry* ent_fitold=leg->AddEntry(fitft_old,"(1) Fit with CMS data","l");
    TLegendEntry* ent_fitnew=leg->AddEntry(fitft,"(2) Fit with ATLAS data","l");
    TLegendEntry* ent_fitc=leg->AddEntry(fitft_c,"(3) Fit with CMS data based on (2)","l");

    leg->Draw("");

    cSigma->SaveAs("cSigma_AddATLAS.pdf");
    //cSigma->SaveAs("cSigma_AddATLAS_wParLimit.pdf");//With ParLimit


}
