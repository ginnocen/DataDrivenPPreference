#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

#define BIN_NUM 460   // 5~120 GeV, interval : 0.25 GeV

#define HMIN 5   
#define HMAX 120    

#define REBIN_bin0 6  // CMS pPb pt bin
#define REBIN_bin1 5  // CMS pp pt bin
#define REBIN_bin2 8  // ATLAS pp pt bin
#define REBIN 55      // For 5~60, bin width : 1 GeV     

void CompFONLL_Bplusdsigmadpt(int isBinned=0,int isNorm=1, int codeNum=1, int codeDen=2){

	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	gStyle->SetPadTopMargin(0.095);
	gStyle->SetPadLeftMargin(0.130);
	gStyle->SetPadRightMargin(0.045);

	//gStyle->SetTitleOffset(0.9,"X");
	//gStyle->SetTitleOffset(5.0,"Y");

	TFile* fin_name[5];
	/*
	fin_name[0]=new TFile("../ResultsBplus/outputBplus_Unbinned_5TeV.root");
	fin_name[1]=new TFile("../ResultsBplus/outputBplus_Unbinned_7TeV.root");
	fin_name[2]=new TFile("../ResultsBplus/outputBplus_Unbinned_2760GeV.root");
	*/
	fin_name[0]=new TFile("../CodeOrg/Rootf/outputBplus_pp_pt_rap24_5p02TeV_PPbBin.root");
	fin_name[1]=new TFile("../CodeOrg/Rootf/outputBplus_pp_pt_rap24_7TeV_PPbBin.root");
	fin_name[2]=new TFile("../CodeOrg/Rootf/outputBplus_pp_pt_rap24_2p76TeV_PPbBin.root");
	fin_name[3]=new TFile("../ResultsBplus/outputBplus_Unbinned_7TeV_MSTW2008nlo68cl.root");
	fin_name[4]=new TFile("../ResultsBplus/outputBplus_Unbinned_7TeV_NNPDF30nlo_as0118.root");

	std::string codeName[5]={"5TeV","7TeV","2760GeV","7TeV_MSTW2008nlo68cl","7TeV_NNPDF30nlo_as0118"};

	int REBINn;

	std::string tlatexrem;

	TFile* fin_5TeV = fin_name[codeNum];
	TFile* fin_7TeV = fin_name[codeDen];

	std::string isBinnedrmk, isNormst;
		if (isBinned==0) {isBinnedrmk="BinnedpPb";REBINn=REBIN_bin0;}
		else if (isBinned==1) {isBinnedrmk="BinnedCMS";REBINn=REBIN_bin1;}
		else if (isBinned==2) {isBinnedrmk="BinnedATL";REBINn=REBIN_bin2;}
		else {isBinnedrmk="Fine";REBINn=BIN_NUM;isBinned=99;}//5~60, 1 GeV interval
	if (isNorm) isNormst="Norm"; else isNormst="Val";

	std::cout << "isBinned : " << isBinned << std::endl;
	TFile* fout = new TFile(Form("../ResultsBplus/CompFONLL_Bplus_%s_%s_%svs%s.root",isBinnedrmk.c_str(),isNormst.c_str(),codeName[codeNum].c_str(),codeName[codeDen].c_str()),"recreate");

	TH1D* hfr_10_10_5TeV=(TH1D*)fin_5TeV->Get("hpt");
	TH1D* hfr_05_05_5TeV=(TH1D*)fin_5TeV->Get("hfr_05_05");
	TH1D* hfr_20_20_5TeV=(TH1D*)fin_5TeV->Get("hfr_20_20");
	TH1D* hfr_20_10_5TeV=(TH1D*)fin_5TeV->Get("hfr_20_10");
	TH1D* hfr_10_20_5TeV=(TH1D*)fin_5TeV->Get("hfr_10_20");
	TH1D* hfr_10_05_5TeV=(TH1D*)fin_5TeV->Get("hfr_10_05");
	TH1D* hfr_05_10_5TeV=(TH1D*)fin_5TeV->Get("hfr_05_10");

	TH1D* hfr_10_10_7TeV=(TH1D*)fin_7TeV->Get("hpt");
	TH1D* hfr_05_05_7TeV=(TH1D*)fin_7TeV->Get("hfr_05_05");
	TH1D* hfr_20_20_7TeV=(TH1D*)fin_7TeV->Get("hfr_20_20");
	TH1D* hfr_20_10_7TeV=(TH1D*)fin_7TeV->Get("hfr_20_10");
	TH1D* hfr_10_20_7TeV=(TH1D*)fin_7TeV->Get("hfr_10_20");
	TH1D* hfr_10_05_7TeV=(TH1D*)fin_7TeV->Get("hfr_10_05");
	TH1D* hfr_05_10_7TeV=(TH1D*)fin_7TeV->Get("hfr_05_10");

	hfr_10_10_5TeV->SetName("hfr_10_10_5TeV");
	hfr_05_05_5TeV->SetName("hfr_05_05_5TeV");
	hfr_20_20_5TeV->SetName("hfr_20_20_5TeV");
	hfr_20_10_5TeV->SetName("hfr_20_10_5TeV");
	hfr_10_20_5TeV->SetName("hfr_10_20_5TeV");
	hfr_10_05_5TeV->SetName("hfr_10_05_5TeV");
	hfr_05_10_5TeV->SetName("hfr_05_10_5TeV");

	hfr_10_10_7TeV->SetName("hfr_10_10_7TeV");
	hfr_05_05_7TeV->SetName("hfr_05_05_7TeV");
	hfr_20_20_7TeV->SetName("hfr_20_20_7TeV");
	hfr_20_10_7TeV->SetName("hfr_20_10_7TeV");
	hfr_10_20_7TeV->SetName("hfr_10_20_7TeV");
	hfr_10_05_7TeV->SetName("hfr_10_05_7TeV");
	hfr_05_10_7TeV->SetName("hfr_05_10_7TeV");

  // REBIN here

  //double rebiny[6] = {10,15,20,25,30,60};//rebin edge
  double rebiny[REBIN_bin0+1]     = {5.,10.,15.,20.,25.,30.,60.};//rebin edge for CMS pPb binning
  double rebiny_CMS[REBIN_bin1+1] = {5.,10.,13.,17.,24.,30.};//rebin edge for CMS pp binning
  double rebiny_ATL[REBIN_bin2+1] = {9.,13.,16.,20.,25.,35.,50.,70.,120.};//rebin edge for ATLAS pp binning

	double rebin[REBINn+1];

		if (isBinned==0) {for (int i=0;i<REBIN_bin0+1;i++) {rebin[i]=rebiny[i];}}
		else if (isBinned==1) {for (int i=0;i<REBIN_bin1+1;i++) {rebin[i]=rebiny_CMS[i];}}
		else if (isBinned==2) {for (int i=0;i<REBIN_bin2+1;i++) {rebin[i]=rebiny_ATL[i];}}
	  //else {for (int i=0;i<REBIN+1;i++) {rebin[i]=5.0+i*1.0;}}

	std::cout << "##### REBINn" << REBINn << std::endl;

	if (isBinned!=99) {
		hfr_10_10_5TeV=(TH1D*)hfr_10_10_5TeV->Rebin(REBINn,"hfr_10_10_5TeV",rebin);
		hfr_05_05_5TeV=(TH1D*)hfr_05_05_5TeV->Rebin(REBINn,"hfr_05_05_5TeV",rebin);
		hfr_20_20_5TeV=(TH1D*)hfr_20_20_5TeV->Rebin(REBINn,"hfr_20_20_5TeV",rebin);
		hfr_20_10_5TeV=(TH1D*)hfr_20_10_5TeV->Rebin(REBINn,"hfr_20_10_5TeV",rebin);
		hfr_10_20_5TeV=(TH1D*)hfr_10_20_5TeV->Rebin(REBINn,"hfr_10_20_5TeV",rebin);
		hfr_10_05_5TeV=(TH1D*)hfr_10_05_5TeV->Rebin(REBINn,"hfr_10_05_5TeV",rebin);
		hfr_05_10_5TeV=(TH1D*)hfr_05_10_5TeV->Rebin(REBINn,"hfr_05_10_5TeV",rebin);

		hfr_10_10_7TeV=(TH1D*)hfr_10_10_7TeV->Rebin(REBINn,"hfr_10_10_7TeV",rebin);
		hfr_05_05_7TeV=(TH1D*)hfr_05_05_7TeV->Rebin(REBINn,"hfr_05_05_7TeV",rebin);
		hfr_20_20_7TeV=(TH1D*)hfr_20_20_7TeV->Rebin(REBINn,"hfr_20_20_7TeV",rebin);
		hfr_20_10_7TeV=(TH1D*)hfr_20_10_7TeV->Rebin(REBINn,"hfr_20_10_7TeV",rebin);
		hfr_10_20_7TeV=(TH1D*)hfr_10_20_7TeV->Rebin(REBINn,"hfr_10_20_7TeV",rebin);
		hfr_10_05_7TeV=(TH1D*)hfr_10_05_7TeV->Rebin(REBINn,"hfr_10_05_7TeV",rebin);
		hfr_05_10_7TeV=(TH1D*)hfr_05_10_7TeV->Rebin(REBINn,"hfr_05_10_7TeV",rebin);
	}


	TH1D* hfr_10_10_Comp=(TH1D*)hfr_10_10_5TeV->Clone();
	TH1D* hfr_05_05_Comp=(TH1D*)hfr_05_05_5TeV->Clone();
	TH1D* hfr_20_20_Comp=(TH1D*)hfr_20_20_5TeV->Clone();
	TH1D* hfr_20_10_Comp=(TH1D*)hfr_20_10_5TeV->Clone();
	TH1D* hfr_10_20_Comp=(TH1D*)hfr_10_20_5TeV->Clone();
	TH1D* hfr_10_05_Comp=(TH1D*)hfr_10_05_5TeV->Clone();
	TH1D* hfr_05_10_Comp=(TH1D*)hfr_05_10_5TeV->Clone();

	TH1D* hfr_05_05_5TeV_C;
	TH1D* hfr_20_20_5TeV_C;
	TH1D* hfr_20_10_5TeV_C;
	TH1D* hfr_10_20_5TeV_C;
	TH1D* hfr_10_05_5TeV_C;
	TH1D* hfr_05_10_5TeV_C;
	TH1D* hfr_10_10_5TeV_C;

	hfr_10_10_Comp->SetName("hfr_10_10_Comp");
	hfr_05_05_Comp->SetName("hfr_05_05_Comp");
	hfr_20_20_Comp->SetName("hfr_20_20_Comp");
	hfr_20_10_Comp->SetName("hfr_20_10_Comp");
	hfr_10_20_Comp->SetName("hfr_10_20_Comp");
	hfr_10_05_Comp->SetName("hfr_10_05_Comp");
	hfr_05_10_Comp->SetName("hfr_05_10_Comp");


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

	hfr_05_05_5TeV_C=(TH1D*)hfr_05_05_5TeV->Clone();
	hfr_20_20_5TeV_C=(TH1D*)hfr_20_20_5TeV->Clone();
	hfr_20_10_5TeV_C=(TH1D*)hfr_20_10_5TeV->Clone();
	hfr_10_20_5TeV_C=(TH1D*)hfr_10_20_5TeV->Clone();
	hfr_10_05_5TeV_C=(TH1D*)hfr_10_05_5TeV->Clone();
	hfr_05_10_5TeV_C=(TH1D*)hfr_05_10_5TeV->Clone();
	hfr_10_10_5TeV_C=(TH1D*)hfr_10_10_5TeV->Clone();

	hfr_10_10_5TeV_C->SetName("hfr_10_10_5TeV_C");
	hfr_05_05_5TeV_C->SetName("hfr_05_05_5TeV_C");
	hfr_20_20_5TeV_C->SetName("hfr_20_20_5TeV_C");
	hfr_20_10_5TeV_C->SetName("hfr_20_10_5TeV_C");
	hfr_10_20_5TeV_C->SetName("hfr_10_20_5TeV_C");
	hfr_10_05_5TeV_C->SetName("hfr_10_05_5TeV_C");
	hfr_05_10_5TeV_C->SetName("hfr_05_10_5TeV_C");


	for (int i=0;i<REBIN;i++){
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

	TH1D* hfrval=(TH1D*)hfr_10_10_Comp->Clone();
	TH1D* hfrmaxerr=(TH1D*)hfr_10_10_Comp->Clone();
	TH1D* hfrminerr=(TH1D*)hfr_10_10_Comp->Clone();

	hfrval->SetName("hfrval");
	hfrmaxerr->SetName("hfrmaxerr");
	hfrminerr->SetName("hfrminerr");

	for (int i=1;i<REBINn+1;i++) {
		double maxbinv,minbinv;
		maxbinv = max(hfr_10_10_Comp->GetBinContent(i),max(hfr_05_05_Comp->GetBinContent(i),max(hfr_20_20_Comp->GetBinContent(i),max(hfr_20_10_Comp->GetBinContent(i),max(hfr_10_20_Comp->GetBinContent(i),max(hfr_10_05_Comp->GetBinContent(i),hfr_05_10_Comp->GetBinContent(i)))))));
		minbinv = min(hfr_10_10_Comp->GetBinContent(i),min(hfr_05_05_Comp->GetBinContent(i),min(hfr_20_20_Comp->GetBinContent(i),min(hfr_20_10_Comp->GetBinContent(i),min(hfr_10_20_Comp->GetBinContent(i),min(hfr_10_05_Comp->GetBinContent(i),hfr_05_10_Comp->GetBinContent(i)))))));
		std::cout << hfr_10_10_Comp->GetBinContent(i) << "+" << maxbinv-hfr_10_10_Comp->GetBinContent(i) << "-" << hfr_10_10_Comp->GetBinContent(i)-minbinv << std::endl;
	hfrval->SetBinContent(i,hfr_10_10_Comp->GetBinContent(i));
	hfrmaxerr->SetBinContent(i,maxbinv-hfr_10_10_Comp->GetBinContent(i));
	hfrminerr->SetBinContent(i,hfr_10_10_Comp->GetBinContent(i)-minbinv);
	}
	TCanvas*c1 = new TCanvas("c1","",500,500);

	//TH1D* hbase = new TH1D("hbase",";B^{+} p_{t} (gev/c);(variation/central)5TeV / (variation/central)7TeV",50,10.0,60.0);
	TH1D* hbase = new TH1D("hbase","",55,5.0,60.0);

	//hfr_10_10_comp->getyaxis()->settitle("(variation/central)5tev / (variation/central)7tev");
	hbase->GetXaxis()->SetTitle("B^{+} p_{T} (Gev/c)");
	if (isNorm==1) {
		hbase->GetYaxis()->SetTitle(Form("(variation/central)%s / (variation/central)%s",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
		//hbase->GetYaxis()->SetRangeUser(0.95,1.06);
		hbase->GetYaxis()->SetRangeUser(0.90,1.20);
	}
	else {
		hbase->GetYaxis()->SetTitle(Form("variation,%s / variation,%s",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
		//hbase->GetYaxis()->SetRangeUser(0.52,0.72);
		//hbase->GetYaxis()->SetRangeUser(0.00,1.2);
		hbase->GetYaxis()->SetRangeUser(0.00,1.2*hfr_10_10_Comp->GetMaximum());
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
	hfr_10_10_Comp->Draw("same");
/*
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
	leg->AddEntry(hfr_10_10_Comp,"   1.0 ,   1.0","l");
*/

  TLegend* leg;
	//if (isNorm==1) leg = new TLegend(0.68,0.51,0.86,0.90); else leg = new TLegend(0.68,0.13,0.86,0.52);
	if (isNorm==1) leg = new TLegend(0.68,0.51,0.86,0.90); else leg = new TLegend(0.68,0.11,0.86,0.46);


  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.040);
  leg->AddEntry(hbase,"#mu_{F}/#mu_{0} , #mu_{R}/#mu_{0}","l");
  leg->AddEntry(hfr_10_10_Comp,"  1.0   ,   1.0","l");
  leg->AddEntry(hfr_05_05_Comp,"  0.5   ,   0.5","l");
  leg->AddEntry(hfr_20_20_Comp,"  2.0   ,   2.0","l");
  leg->AddEntry(hfr_20_10_Comp,"  2.0   ,   1.0","l");
  leg->AddEntry(hfr_10_20_Comp,"  1.0   ,   2.0","l");
  leg->AddEntry(hfr_10_05_Comp,"  1.0   ,   0.5","l");
  leg->AddEntry(hfr_05_10_Comp,"  0.5   ,   1.0","l");

	leg->Draw("");

  if (isBinned!=99) {
    if (isNorm==1) c1->SaveAs(Form("../ResultsBplus/CompFONLL_Bplus_%s_Norm_%svs%s.pdf",isBinnedrmk.c_str(),codeName[codeNum].c_str(),codeName[codeDen].c_str()));
    else c1->SaveAs(Form("../ResultsBplus/CompFONLL_Bplus_%s_Val_%svs%s.pdf",isBinnedrmk.c_str(),codeName[codeNum].c_str(),codeName[codeDen].c_str()));
  }
  else {
    if (isNorm==1) c1->SaveAs(Form("../ResultsBplus/CompFONLL_Bplus_Unbinned_Norm_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
    else c1->SaveAs(Form("../ResultsBplus/CompFONLL_Bplus_Unbinned_Val_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
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

	fout->cd();
	hfrval->Write();
	hfrmaxerr->Write();
	hfrminerr->Write();

	hfr_10_10_Comp->Write();
	hfr_05_05_Comp->Write();
	hfr_20_20_Comp->Write();
	hfr_20_10_Comp->Write();
	hfr_10_20_Comp->Write();
	hfr_05_10_Comp->Write();
	hfr_10_05_Comp->Write();
	hfr_10_10_5TeV_C->Write();
	hfr_05_05_5TeV_C->Write();
	hfr_20_20_5TeV_C->Write();
	hfr_20_10_5TeV_C->Write();
	hfr_10_20_5TeV_C->Write();
	hfr_05_10_5TeV_C->Write();
	hfr_10_05_5TeV_C->Write();
	hfr_10_10_7TeV->Write();
	hfr_05_05_7TeV->Write();
	hfr_20_20_7TeV->Write();
	hfr_20_10_7TeV->Write();
	hfr_10_20_7TeV->Write();
	hfr_05_10_7TeV->Write();
	hfr_10_05_7TeV->Write();
	hfr_10_10_5TeV->Write();
	hfr_05_05_5TeV->Write();
	hfr_20_20_5TeV->Write();
	hfr_20_10_5TeV->Write();
	hfr_10_20_5TeV->Write();
	hfr_05_10_5TeV->Write();
	hfr_10_05_5TeV->Write();


	//fout->Write();
}
