TString particle="Bplus";
const int nbins=5;
Double_t xbins[nbins]={12.5,17.5,22.5,27.5,45.};
Double_t exl[nbins]={2.5,2.5,2.5,2.5,15.};
Double_t exl0[nbins]={0.,0.,0.,0.,0.};
Double_t yPercSigmapPbSystTotHigh[nbins]={0.163,0.150,0.146,0.142,0.140};
Double_t yPercSigmapPbSystTotLow[nbins]={0.163,0.150,0.146,0.142,0.140};
Double_t commonErrorP = TMath::Sqrt(0.0445*0.0445);
Double_t commonErrorN = TMath::Sqrt(0.0445*0.0445);
Double_t FFsysterror=0.7/40.1;
Double_t tagandprobcorrection[nbins]={1.049,1.030,1.019,1.012,1.006};


TString fofrom;

void ReferencePPOverFO(int option=2){

  if(option==2)  fofrom= "2760GeV";
  if(option==7)  fofrom= "7TeV";

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  TFile*filePPReference=new TFile(Form("../ResultsBplus/Estimatedpp5TeV_with%s.root",fofrom.Data()));  
  TGraphAsymmErrors*gaeBplusPPReference=(TGraphAsymmErrors*)filePPReference->Get("gaeEstimatedpp5TeV");
  gaeBplusPPReference->SetName(Form("gae%sPPReference",particle.Data()));
  TFile*fileFOReference=new TFile(Form("../ResultsBplus/output%s.root",particle.Data()));  
  TGraphAsymmErrors*gaeBplusFOReference=(TGraphAsymmErrors*)fileFOReference->Get(Form("gaeSigmaDecay%s",particle.Data()));
  gaeBplusFOReference->SetName(Form("gae%sFOReference",particle.Data()));
  
  double rebin[nbins+1] = {10,15,20,25,30,60};
  TH1F* hratio_rebin = new TH1F("hratio_rebin","hratio_rebin",5,rebin);

  double yPP,yFO,xvalue;

  for (int i=0;i<nbins;i++)
    {
      yPP=-1.;
      yFO=-1.;
      xvalue=-1.;
      if(fofrom=="7TeV") gaeBplusPPReference->GetPoint(i+1,xvalue,yPP);
      else if(fofrom=="2760GeV") gaeBplusPPReference->GetPoint(i,xvalue,yPP);
      xvalue=-1.;
      gaeBplusFOReference->GetPoint(i,xvalue,yFO);
      hratio_rebin->SetBinContent(i+1,yPP*(1.e+6)*208/yFO);
    }
  
  TCanvas*cratio=new TCanvas("cratio","cratio",500,500);
  hratio_rebin->SetMaximum(2.);
  hratio_rebin->SetMinimum(0.5);
  hratio_rebin->SetXTitle("p_{T}(GeV/c)");
  hratio_rebin->SetYTitle("Reference_{Data Driven} / Reference_{FONLL}");
  hratio_rebin->SetTitleOffset(1.2,"Y");
  hratio_rebin->SetLineColor(kRed);
  hratio_rebin->SetFillStyle(3004);
  hratio_rebin->SetFillColor(kRed);
  hratio_rebin->SetLineWidth(3);
  hratio_rebin->Draw();

  TLegend *leg = new TLegend(0.5,0.75,0.9,0.9);
  leg->AddEntry((TObject*)0,"Reference 5.02TeV","");
  leg->AddEntry((TObject*)0,"|y_{LAB}|<2.4","");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("same");

  TLine* lin0=new TLine(9,1,120,1);
  lin0->SetLineStyle(2);
  lin0->SetLineColor(1);
  lin0->SetLineWidth(3);
  lin0->Draw("same");

  cratio->SaveAs(Form("Plots/cratio_ReferencePPOverFO_%s.pdf",fofrom.Data()));

}
