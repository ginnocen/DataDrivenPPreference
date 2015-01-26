#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

#define BIN_NUM 220
#define HMIN 5
#define HMAX 60
#define REBIN 5
#define REBIN7 9

void CompFONLL_BplusdsigmadyUnbinned_v6(int isBinned=0,int isNorm=1,int codeNum=1,int codeDen=2){

    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    gStyle->SetPadTopMargin(0.095);
    gStyle->SetPadLeftMargin(0.130);
    gStyle->SetPadRightMargin(0.045);

    //gStyle->SetTitleOffset(0.9,"X");
    //gStyle->SetTitleOffset(5.0,"Y");

    TFile* fin_name[5];
    fin_name[0]=new TFile("../outputBplusyUnbinned_5TeV_rap30_pt1060.root");//already boosted
    fin_name[1]=new TFile("../outputBplusyUnbinned_7TeV_rap30_pt1060.root");
    fin_name[2]=new TFile("../outputBplusyUnbinned_7TeV_rap30_pt5120.root");
    fin_name[3]=new TFile("../outputBplusyUnbinned_7TeVboosted.root");
    fin_name[4]=new TFile("../outputBplusyUnbinned_7TeV_rap30_pt9120.root");

    //fin_name[0]=new TFile("../outputBplusyUnbinned_5TeV.root");
    //fin_name[1]=new TFile("../outputBplusyUnbinned_7TeV.root");
    //fin_name[2]=new TFile("../outputBplusyUnbinned_7TeV_pt5120.root");
    //fin_name[4]=new TFile("../outputBplusyUnbinned_7TeV_pt9120.root");


    bool isBoosted[5]={true,false,false,true,false};
    //std::string codeName[5]={"5TeV","7TeV","7TeVpt5120","7TeVboosted","7TeVpt9120"};
    //###std::string codeName[5]={"5TeV","7TeVrap30pt1060","7TeVrap30pt5120","7TeVboosted","7TeVrap30pt9120"};
    std::string codeName[5]={"5TeVrap30pt1060","7TeVrap30pt1060","7TeVrap30pt5120","7TeVboosted","7TeVrap30pt9120"};

    TFile* fin_5TeV = fin_name[codeNum];
    TFile* fin_7TeV = fin_name[codeDen];

    TH1F* hfr_10_10_5TeV=(TH1F*)fin_5TeV->Get("hpt");
    TH1F* hfr_05_05_5TeV=(TH1F*)fin_5TeV->Get("hfr_05_05");
    TH1F* hfr_20_20_5TeV=(TH1F*)fin_5TeV->Get("hfr_20_20");
    TH1F* hfr_20_10_5TeV=(TH1F*)fin_5TeV->Get("hfr_20_10");
    TH1F* hfr_10_20_5TeV=(TH1F*)fin_5TeV->Get("hfr_10_20");
    TH1F* hfr_10_05_5TeV=(TH1F*)fin_5TeV->Get("hfr_10_05");
    TH1F* hfr_05_10_5TeV=(TH1F*)fin_5TeV->Get("hfr_05_10");

    TH1F* hfr_10_10_7TeV=(TH1F*)fin_7TeV->Get("hpt");
    TH1F* hfr_05_05_7TeV=(TH1F*)fin_7TeV->Get("hfr_05_05");
    TH1F* hfr_20_20_7TeV=(TH1F*)fin_7TeV->Get("hfr_20_20");
    TH1F* hfr_20_10_7TeV=(TH1F*)fin_7TeV->Get("hfr_20_10");
    TH1F* hfr_10_20_7TeV=(TH1F*)fin_7TeV->Get("hfr_10_20");
    TH1F* hfr_10_05_7TeV=(TH1F*)fin_7TeV->Get("hfr_10_05");
    TH1F* hfr_05_10_7TeV=(TH1F*)fin_7TeV->Get("hfr_05_10");

// v4
//    double rebin[REBIN+1] = {-2.4,-1.470,-0.535,0.465,1.465,2.4};
//    double rebin5[REBIN+1] = {-2.865,-1.935,-1.0,0.0,1.0,1.935};
// v5
//    double rebin7[REBIN7+1] = {-2.4,-1.8,-1.45,-1.10,-0.60,0.60,1.10,1.45,1.80,2.40};
// v6
    //double rebin7[REBIN+1] = {-2.865,-1.935,-1.0,0.0,1.0,1.935};



    double rebin5[REBIN+1];
    double rebin7[REBIN+1] = {-2.4,-1.470,-0.535,0.465,1.465,2.4};

    int REBINn;
    REBINn=REBIN;
    REBINn=REBIN7;

    REBINn=REBIN;

    double rebinnum[REBINn+1];
    double rebinden[REBINn+1];

    for (int i=0;i<REBINn+1;i++) {
    rebin5[i]=rebin7[i]-0.465;
    if (isBoosted[codeNum]) rebinnum[i]=rebin7[i]; else rebinnum[i]=rebin5[i];
    if (isBoosted[codeDen]) rebinden[i]=rebin7[i]; else rebinden[i]=rebin5[i];
    }
    if (isBinned==1) {
    hfr_10_10_5TeV=(TH1F*)hfr_10_10_5TeV->Rebin(REBINn,"hfr_10_10_5TeV",rebinnum);
    hfr_05_05_5TeV=(TH1F*)hfr_05_05_5TeV->Rebin(REBINn,"hfr_05_05_5TeV",rebinnum);
    hfr_20_20_5TeV=(TH1F*)hfr_20_20_5TeV->Rebin(REBINn,"hfr_20_20_5TeV",rebinnum);
    hfr_20_10_5TeV=(TH1F*)hfr_20_10_5TeV->Rebin(REBINn,"hfr_20_10_5TeV",rebinnum);
    hfr_10_20_5TeV=(TH1F*)hfr_10_20_5TeV->Rebin(REBINn,"hfr_10_20_5TeV",rebinnum);
    hfr_10_05_5TeV=(TH1F*)hfr_10_05_5TeV->Rebin(REBINn,"hfr_10_05_5TeV",rebinnum);
    hfr_05_10_5TeV=(TH1F*)hfr_05_10_5TeV->Rebin(REBINn,"hfr_05_10_5TeV",rebinnum);

    hfr_10_10_7TeV=(TH1F*)hfr_10_10_7TeV->Rebin(REBINn,"hfr_10_10_7TeV",rebinden);
    hfr_05_05_7TeV=(TH1F*)hfr_05_05_7TeV->Rebin(REBINn,"hfr_05_05_7TeV",rebinden);
    hfr_20_20_7TeV=(TH1F*)hfr_20_20_7TeV->Rebin(REBINn,"hfr_20_20_7TeV",rebinden);
    hfr_20_10_7TeV=(TH1F*)hfr_20_10_7TeV->Rebin(REBINn,"hfr_20_10_7TeV",rebinden);
    hfr_10_20_7TeV=(TH1F*)hfr_10_20_7TeV->Rebin(REBINn,"hfr_10_20_7TeV",rebinden);
    hfr_10_05_7TeV=(TH1F*)hfr_10_05_7TeV->Rebin(REBINn,"hfr_10_05_7TeV",rebinden);
    hfr_05_10_7TeV=(TH1F*)hfr_05_10_7TeV->Rebin(REBINn,"hfr_05_10_7TeV",rebinden);
    }


    TH1F* hfr_10_10_Comp=(TH1F*)hfr_10_10_7TeV->Clone();
    TH1F* hfr_05_05_Comp=(TH1F*)hfr_05_05_7TeV->Clone();
    TH1F* hfr_20_20_Comp=(TH1F*)hfr_20_20_7TeV->Clone();
    TH1F* hfr_20_10_Comp=(TH1F*)hfr_20_10_7TeV->Clone();
    TH1F* hfr_10_20_Comp=(TH1F*)hfr_10_20_7TeV->Clone();
    TH1F* hfr_10_05_Comp=(TH1F*)hfr_10_05_7TeV->Clone();
    TH1F* hfr_05_10_Comp=(TH1F*)hfr_05_10_7TeV->Clone();

    TH1F* hfr_05_05_5TeV_C;
    TH1F* hfr_20_20_5TeV_C;
    TH1F* hfr_20_10_5TeV_C;
    TH1F* hfr_10_20_5TeV_C;
    TH1F* hfr_10_05_5TeV_C;
    TH1F* hfr_05_10_5TeV_C;
    TH1F* hfr_10_10_5TeV_C;
 
    if (isNorm==1) {
    hfr_05_05_5TeV->Divide(hfr_05_05_5TeV,hfr_10_10_5TeV,1,1,"B");
    hfr_20_20_5TeV->Divide(hfr_20_20_5TeV,hfr_10_10_5TeV,1,1,"B");
    hfr_20_10_5TeV->Divide(hfr_20_10_5TeV,hfr_10_10_5TeV,1,1,"B");
    hfr_10_20_5TeV->Divide(hfr_10_20_5TeV,hfr_10_10_5TeV,1,1,"B");
    hfr_10_05_5TeV->Divide(hfr_10_05_5TeV,hfr_10_10_5TeV,1,1,"B");
    hfr_05_10_5TeV->Divide(hfr_05_10_5TeV,hfr_10_10_5TeV,1,1,"B");
    hfr_10_10_5TeV->Divide(hfr_10_10_5TeV,hfr_10_10_5TeV,1,1,"B");

    hfr_05_05_7TeV->Divide(hfr_05_05_7TeV,hfr_10_10_7TeV,1,1,"B");
    hfr_20_20_7TeV->Divide(hfr_20_20_7TeV,hfr_10_10_7TeV,1,1,"B");
    hfr_20_10_7TeV->Divide(hfr_20_10_7TeV,hfr_10_10_7TeV,1,1,"B");
    hfr_10_20_7TeV->Divide(hfr_10_20_7TeV,hfr_10_10_7TeV,1,1,"B");
    hfr_10_05_7TeV->Divide(hfr_10_05_7TeV,hfr_10_10_7TeV,1,1,"B");
    hfr_05_10_7TeV->Divide(hfr_05_10_7TeV,hfr_10_10_7TeV,1,1,"B");
    hfr_10_10_7TeV->Divide(hfr_10_10_7TeV,hfr_10_10_7TeV,1,1,"B");

    }

    hfr_05_05_5TeV_C=(TH1F*)hfr_05_05_7TeV->Clone();
    hfr_20_20_5TeV_C=(TH1F*)hfr_20_20_7TeV->Clone();
    hfr_20_10_5TeV_C=(TH1F*)hfr_20_10_7TeV->Clone();
    hfr_10_20_5TeV_C=(TH1F*)hfr_10_20_7TeV->Clone();
    hfr_10_05_5TeV_C=(TH1F*)hfr_10_05_7TeV->Clone();
    hfr_05_10_5TeV_C=(TH1F*)hfr_05_10_7TeV->Clone();
    hfr_10_10_5TeV_C=(TH1F*)hfr_10_10_7TeV->Clone();

    for (int i=0;i<5;i++){
    hfr_05_05_5TeV_C->SetBinContent(i+1,hfr_05_05_5TeV->GetBinContent(i+1));
    hfr_20_20_5TeV_C->SetBinContent(i+1,hfr_20_20_5TeV->GetBinContent(i+1));
    hfr_20_10_5TeV_C->SetBinContent(i+1,hfr_20_10_5TeV->GetBinContent(i+1));
    hfr_10_20_5TeV_C->SetBinContent(i+1,hfr_10_20_5TeV->GetBinContent(i+1));
    hfr_10_05_5TeV_C->SetBinContent(i+1,hfr_10_05_5TeV->GetBinContent(i+1));
    hfr_05_10_5TeV_C->SetBinContent(i+1,hfr_05_10_5TeV->GetBinContent(i+1));
    hfr_10_10_5TeV_C->SetBinContent(i+1,hfr_10_10_5TeV->GetBinContent(i+1));
    }

    hfr_10_10_Comp->Divide(hfr_10_10_5TeV_C,hfr_10_10_7TeV,1,1,"B");
    hfr_05_05_Comp->Divide(hfr_05_05_5TeV_C,hfr_05_05_7TeV,1,1,"B");
    hfr_20_20_Comp->Divide(hfr_20_20_5TeV_C,hfr_20_20_7TeV,1,1,"B");
    hfr_20_10_Comp->Divide(hfr_20_10_5TeV_C,hfr_20_10_7TeV,1,1,"B");
    hfr_10_20_Comp->Divide(hfr_10_20_5TeV_C,hfr_10_20_7TeV,1,1,"B");
    hfr_10_05_Comp->Divide(hfr_10_05_5TeV_C,hfr_10_05_7TeV,1,1,"B");
    hfr_05_10_Comp->Divide(hfr_05_10_5TeV_C,hfr_05_10_7TeV,1,1,"B");

    for (int i=1;i<REBINn+1;i++) {
	double maxbinv,minbinv;
	maxbinv = max(hfr_10_10_Comp->GetBinContent(i),max(hfr_05_05_Comp->GetBinContent(i),max(hfr_20_20_Comp->GetBinContent(i),max(hfr_20_10_Comp->GetBinContent(i),max(hfr_10_20_Comp->GetBinContent(i),max(hfr_10_05_Comp->GetBinContent(i),hfr_05_10_Comp->GetBinContent(i)))))));
	minbinv = min(hfr_10_10_Comp->GetBinContent(i),min(hfr_05_05_Comp->GetBinContent(i),min(hfr_20_20_Comp->GetBinContent(i),min(hfr_20_10_Comp->GetBinContent(i),min(hfr_10_20_Comp->GetBinContent(i),min(hfr_10_05_Comp->GetBinContent(i),hfr_05_10_Comp->GetBinContent(i)))))));
	std::cout << hfr_10_10_Comp->GetBinContent(i) << "+" << maxbinv-hfr_10_10_Comp->GetBinContent(i) << "-" << hfr_10_10_Comp->GetBinContent(i)-minbinv << std::endl;
	//std::cout << hfr_10_10_Comp->GetBinContent(i) << "+" << hfr_05_05_Comp->GetBinContent(i)-hfr_10_10_Comp->GetBinContent(i) << "-" << hfr_10_10_Comp->GetBinContent(i)-hfr_20_20_Comp->GetBinContent(i) << std::endl;
    }
    TCanvas*c1 = new TCanvas("c1","",500,500);

    //TH1F* hbase = new TH1F("hbase",";B^{+} p_{t} (gev/c);(variation/central)5TeV / (variation/central)7TeV",50,10.0,60.0);
    TH1F* hbase;
    if (isBoosted[codeNum]) hbase = new TH1F("hbase","",50,-2.965,2.035); else hbase = new TH1F("hbase","",50,-2.5,2.5);
 
    //hfr_10_10_comp->getyaxis()->settitle("(variation/central)5tev / (variation/central)7tev");
    hbase->GetXaxis()->SetTitle("B^{+} y_{LAB}");
    if (isNorm==1) {
	hbase->GetYaxis()->SetTitle("variation/central)5TeV / (variation/central)7TeV");
        hbase->GetYaxis()->SetRangeUser(0.95,1.08);
   }
    else {
	hbase->GetYaxis()->SetTitle("variation,5TeV / variation,7TeV");
        //hbase->GetYaxis()->SetRangeUser(0.52,0.72);
        hbase->GetYaxis()->SetRangeUser(0.00,1.40);
    }
    hbase->GetXaxis()->CenterTitle();
    hbase->GetYaxis()->CenterTitle();
    hbase->GetXaxis()->SetTitleOffset(1.0);
    hbase->GetYaxis()->SetTitleOffset(1.5);
    hbase->SetLineColor(0);

    hfr_05_05_Comp->SetLineWidth(2);
    hfr_20_20_Comp->SetLineWidth(2);

    hbase->Draw("");
    hfr_10_10_Comp->Draw("same");
    hfr_05_05_Comp->Draw("same");
    hfr_20_20_Comp->Draw("same");
    hfr_20_10_Comp->Draw("same");
    hfr_10_20_Comp->Draw("same");
    hfr_10_05_Comp->Draw("same");
    hfr_05_10_Comp->Draw("same");

    TLegend* leg = new TLegend(0.75,0.58,0.93,0.90);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(hbase,"#mu_{F}/#mu_{0} , #mu_{R}/#mu_{0}","l");
    leg->AddEntry(hfr_10_10_Comp,"   1.0 ,   1.0","l");
    leg->AddEntry(hfr_05_05_Comp,"   0.5 ,   0.5","l");
    leg->AddEntry(hfr_20_20_Comp,"   2.0 ,   2.0","l");
    leg->AddEntry(hfr_20_10_Comp,"   2.0 ,   1.0","l");
    leg->AddEntry(hfr_10_20_Comp,"   1.0 ,   2.0","l");
    leg->AddEntry(hfr_10_05_Comp,"   1.0 ,   0.5","l");
    leg->AddEntry(hfr_05_10_Comp,"   0.5 ,   1.0","l");
    leg->Draw("");
/*
    if (isBinned==1) {
	if (isNorm==1) c1->SaveAs("Comp_Binnedy_Norm_FONLL.pdf");
	else c1->SaveAs("Comp_Binnedy_Val_FONLL.pdf");
    }
    else {
	if (isNorm==1) c1->SaveAs("Comp_Unbinnedy_Norm_FONLL.pdf");
	else c1->SaveAs("Comp_Unbinnedy_Val_FONLL.pdf");
    }
*/
/*
    if (isBinned==1) {
	if (isNorm==1) c1->SaveAs(Form("Comp_Binnedy_Norm_FONLL_v3_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
	else c1->SaveAs(Form("Comp_Binnedy_Val_FONLL_v3_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
    }
    else {
	if (isNorm==1) c1->SaveAs(Form("Comp_Unbinnedy_Norm_FONLL_v3_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
	else c1->SaveAs(Form("Comp_Unbinnedy_Val_FONLL_v3_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
    }
*/
    if (isBinned==1) {
	if (isNorm==1) c1->SaveAs(Form("Comp_Binnedy_Norm_FONLL_v6_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
	else c1->SaveAs(Form("Comp_Binnedy_Val_FONLL_v6_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
    }
    else {
	if (isNorm==1) c1->SaveAs(Form("Comp_Unbinnedy_Norm_FONLL_v6_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
	else c1->SaveAs(Form("Comp_Unbinnedy_Val_FONLL_v6_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
    }

/*
    c1->Clear();
    //hfr_10_10_Comp_rebin->GetYaxis()->SetRangeUser(0.95,1.05);
    hfr_10_10_Comp_rebin->GetYaxis()->SetTitle("(variation/central)5TeV / (variation/central)7TeV");
    hfr_10_10_Comp_rebin->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
    hfr_10_10_Comp_rebin->GetXaxis()->CenterTitle();
    hfr_10_10_Comp_rebin->GetYaxis()->CenterTitle();
    hfr_10_10_Comp_rebin->GetXaxis()->SetTitleOffset(1.0);
    hfr_10_10_Comp_rebin->GetYaxis()->SetTitleOffset(1.5);

    hfr_05_05_Comp_rebin->SetLineWidth(2);
    hfr_20_20_Comp_rebin->SetLineWidth(2);

    hfr_10_10_Comp_rebin->Draw("");
    hfr_05_05_Comp_rebin->Draw("same");
    hfr_20_20_Comp_rebin->Draw("same");
    hfr_20_10_Comp_rebin->Draw("same");
    hfr_10_20_Comp_rebin->Draw("same");
    hfr_10_05_Comp_rebin->Draw("same");
    hfr_10_10_Comp_rebin->Draw("same");

    c1->SaveAs("ComprebinFONLL.pdf");
*/
}
