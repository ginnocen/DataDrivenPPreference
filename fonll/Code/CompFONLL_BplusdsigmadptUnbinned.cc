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

void CompFONLL_BplusdsigmadptUnbinned(int isBinned=0,int isNorm=1, int codeNum=1, int codeDen=2){

	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	gStyle->SetPadTopMargin(0.095);
	gStyle->SetPadLeftMargin(0.130);
	gStyle->SetPadRightMargin(0.045);

	//gStyle->SetTitleOffset(0.9,"X");
	//gStyle->SetTitleOffset(5.0,"Y");

	TFile* fin_name[3];
	fin_name[0]=new TFile("../outputBplus_Unbinned_5TeV.root");
	fin_name[1]=new TFile("../outputBplus_Unbinned_7TeV.root");
	fin_name[2]=new TFile("../outputBplus_Unbinned_2760GeV.root");

	std::string codeName[3]={"5TeV","7TeV","2760GeV"};

	TFile* fin_5TeV = fin_name[codeNum];
	TFile* fin_7TeV = fin_name[codeDen];

	std::string isBinnedst, isNormst;
	if (isBinned) isBinnedst="Binned"; else isBinnedst="Fine";
	if (isNorm) isNormst="Norm"; else isNormst="Val";

	TFile* fout = new TFile(Form("CompFONLL_Bplus_%s_%s_%svs%s.root",isBinnedst.c_str(),isNormst.c_str(),codeName[codeNum].c_str(),codeName[codeDen].c_str()),"recreate");

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

	double rebin[REBIN+1] = {10.0,15.0,20.0,25.0,30.0,60.0};

	if (isBinned==1) {
		hfr_10_10_5TeV=(TH1F*)hfr_10_10_5TeV->Rebin(REBIN,"hfr_10_10_5TeV",rebin);
		hfr_05_05_5TeV=(TH1F*)hfr_05_05_5TeV->Rebin(REBIN,"hfr_05_05_5TeV",rebin);
		hfr_20_20_5TeV=(TH1F*)hfr_20_20_5TeV->Rebin(REBIN,"hfr_20_20_5TeV",rebin);
		hfr_20_10_5TeV=(TH1F*)hfr_20_10_5TeV->Rebin(REBIN,"hfr_20_10_5TeV",rebin);
		hfr_10_20_5TeV=(TH1F*)hfr_10_20_5TeV->Rebin(REBIN,"hfr_10_20_5TeV",rebin);
		hfr_10_05_5TeV=(TH1F*)hfr_10_05_5TeV->Rebin(REBIN,"hfr_10_05_5TeV",rebin);
		hfr_05_10_5TeV=(TH1F*)hfr_05_10_5TeV->Rebin(REBIN,"hfr_05_10_5TeV",rebin);

		hfr_10_10_7TeV=(TH1F*)hfr_10_10_7TeV->Rebin(REBIN,"hfr_10_10_7TeV",rebin);
		hfr_05_05_7TeV=(TH1F*)hfr_05_05_7TeV->Rebin(REBIN,"hfr_05_05_7TeV",rebin);
		hfr_20_20_7TeV=(TH1F*)hfr_20_20_7TeV->Rebin(REBIN,"hfr_20_20_7TeV",rebin);
		hfr_20_10_7TeV=(TH1F*)hfr_20_10_7TeV->Rebin(REBIN,"hfr_20_10_7TeV",rebin);
		hfr_10_20_7TeV=(TH1F*)hfr_10_20_7TeV->Rebin(REBIN,"hfr_10_20_7TeV",rebin);
		hfr_10_05_7TeV=(TH1F*)hfr_10_05_7TeV->Rebin(REBIN,"hfr_10_05_7TeV",rebin);
		hfr_05_10_7TeV=(TH1F*)hfr_05_10_7TeV->Rebin(REBIN,"hfr_05_10_7TeV",rebin);
	}


	TH1F* hfr_10_10_Comp=(TH1F*)hfr_10_10_5TeV->Clone();
	TH1F* hfr_05_05_Comp=(TH1F*)hfr_05_05_5TeV->Clone();
	TH1F* hfr_20_20_Comp=(TH1F*)hfr_20_20_5TeV->Clone();
	TH1F* hfr_20_10_Comp=(TH1F*)hfr_20_10_5TeV->Clone();
	TH1F* hfr_10_20_Comp=(TH1F*)hfr_10_20_5TeV->Clone();
	TH1F* hfr_10_05_Comp=(TH1F*)hfr_10_05_5TeV->Clone();
	TH1F* hfr_05_10_Comp=(TH1F*)hfr_05_10_5TeV->Clone();

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

	TH1F* hfrval=(TH1F*)hfr_10_10_Comp->Clone();
	TH1F* hfrmaxerr=(TH1F*)hfr_10_10_Comp->Clone();
	TH1F* hfrminerr=(TH1F*)hfr_10_10_Comp->Clone();

	hfrval->SetName("hfrval");
	hfrmaxerr->SetName("hfrmaxerr");
	hfrminerr->SetName("hfrminerr");

	for (int i=1;i<REBIN+1;i++) {
		double maxbinv,minbinv;
		maxbinv = max(hfr_10_10_Comp->GetBinContent(i),max(hfr_05_05_Comp->GetBinContent(i),max(hfr_20_20_Comp->GetBinContent(i),max(hfr_20_10_Comp->GetBinContent(i),max(hfr_10_20_Comp->GetBinContent(i),max(hfr_10_05_Comp->GetBinContent(i),hfr_05_10_Comp->GetBinContent(i)))))));
		minbinv = min(hfr_10_10_Comp->GetBinContent(i),min(hfr_05_05_Comp->GetBinContent(i),min(hfr_20_20_Comp->GetBinContent(i),min(hfr_20_10_Comp->GetBinContent(i),min(hfr_10_20_Comp->GetBinContent(i),min(hfr_10_05_Comp->GetBinContent(i),hfr_05_10_Comp->GetBinContent(i)))))));
		std::cout << hfr_10_10_Comp->GetBinContent(i) << "+" << maxbinv-hfr_10_10_Comp->GetBinContent(i) << "-" << hfr_10_10_Comp->GetBinContent(i)-minbinv << std::endl;
	hfrval->SetBinContent(i,hfr_10_10_Comp->GetBinContent(i));
	hfrmaxerr->SetBinContent(i,maxbinv-hfr_10_10_Comp->GetBinContent(i));
	hfrminerr->SetBinContent(i,hfr_10_10_Comp->GetBinContent(i)-minbinv);
	}
	TCanvas*c1 = new TCanvas("c1","",500,500);

	//TH1F* hbase = new TH1F("hbase",";B^{+} p_{t} (gev/c);(variation/central)5TeV / (variation/central)7TeV",50,10.0,60.0);
	TH1F* hbase = new TH1F("hbase","",50,10.0,60.0);

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
	if (isNorm==1) leg = new TLegend(0.68,0.51,0.86,0.90); else leg = new TLegend(0.68,0.15,0.86,0.54);

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

  if (isBinned==1) {
    if (isNorm==1) c1->SaveAs(Form("CompFONLL_Bplus_Binned_Norm_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
    else c1->SaveAs(Form("CompFONLL_Bplus_Binned_Val_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
  }
  else {
    if (isNorm==1) c1->SaveAs(Form("CompFONLL_Bplus_Unbinned_Norm_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
    else c1->SaveAs(Form("CompFONLL_Bplus_Unbinned_Val_%svs%s.pdf",codeName[codeNum].c_str(),codeName[codeDen].c_str()));
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
}
