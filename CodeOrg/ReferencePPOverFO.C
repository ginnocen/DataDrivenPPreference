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
  
  double yPP[nbins],yFO[nbins],yPPerrorup[nbins],yPPerrorlow[nbins],yFOerrorup[nbins],yFOerrorlow[nbins],ratio[nbins];
  
  double yPPval,yFOval,xvalue;

  for (int i=0;i<nbins;i++)
    {      
      if(fofrom=="7TeV") {
        gaeBplusPPReference->GetPoint(i+1,xvalue,yPPval);
        yPP[i]=yPPval;
        yPPerrorup[i] =gaeBplusPPReference->GetErrorYhigh(i+1)/yPPval;
        yPPerrorlow[i]=gaeBplusPPReference->GetErrorYlow(i+1)/yPPval;
      }
      else if(fofrom=="2760GeV") {
        gaeBplusPPReference->GetPoint(i,xvalue,yPPval);
        yPP[i]=yPPval;
        yPPerrorup[i] =gaeBplusPPReference->GetErrorYhigh(i)/yPPval;
        yPPerrorlow[i]=gaeBplusPPReference->GetErrorYlow(i)/yPPval;

      }
      xvalue=-1.;
      gaeBplusFOReference->GetPoint(i,xvalue,yFOval);
      yFO[i]=yFOval;
      yFOerrorup[i]=gaeBplusFOReference->GetErrorYhigh(i)/yFOval;
      yFOerrorlow[i]=gaeBplusFOReference->GetErrorYlow(i)/yFOval;
      
    } 
    
  for (int i=0;i<nbins;i++)
    {      
    ratio[i]=yPP[i]*(1.e+6)*208/yFO[i];
    yPPerrorup[i]=yPPerrorup[i] *ratio[i];
    yPPerrorlow[i]=yPPerrorlow[i] *ratio[i];
    yFOerrorup[i]=yFOerrorup[i] *ratio[i];
    yFOerrorlow[i]=yFOerrorlow[i] *ratio[i];
    hratio_rebin->SetBinContent(i+1,yPP[i]*(1.e+6)*208/yFO[i]);

    }   
    
  TGraphAsymmErrors *gaeDD = new TGraphAsymmErrors(nbins,xbins,ratio,exl,exl,yPPerrorlow,yPPerrorup);  
  TGraphAsymmErrors *gaeFONLL = new TGraphAsymmErrors(nbins,xbins,ratio,exl,exl,yFOerrorlow,yFOerrorup);        


  
  TCanvas*cratio=new TCanvas("cratio","cratio",500,500);
  hratio_rebin->SetMaximum(3);
  hratio_rebin->SetMinimum(0);
  hratio_rebin->SetXTitle("p_{T}(GeV/c)");
  hratio_rebin->SetYTitle("Reference_{Data Driven} / Reference_{FONLL}");
  hratio_rebin->SetTitleOffset(1.2,"Y");
  //hratio_rebin->SetLineColor(kRed);
  //hratio_rebin->SetFillStyle(3004);
  //hratio_rebin->SetFillColor(kRed);
  hratio_rebin->SetLineWidth(3);
  hratio_rebin->Draw("p");


  gaeFONLL->SetMarkerColor(1);
  gaeFONLL->SetMarkerStyle(21);  
  gaeFONLL->SetFillColor(5);
  gaeFONLL->SetFillStyle(1001);
  gaeFONLL->SetLineColor(1);
  gaeFONLL->SetLineWidth(5);
  gaeFONLL->Draw("2");
  
  gaeDD->SetLineColor(1);
  gaeDD->SetMarkerColor(1);
  gaeDD->Draw("epsame");
  gaeDD->SetLineWidth(3);
  

  TLegend *legend=new TLegend(0.1975806,0.6109937,0.4959677,0.8012685,"");
  legend->SetBorderSize(0);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(1001);
  legend->SetTextFont(42);
  legend->SetTextSize(0.045);

  TLegendEntry *ent_FONLL=legend->AddEntry(gaeFONLL,"syst fonll uncertainty","f");
  ent_FONLL->SetTextFont(42);
  ent_FONLL->SetLineColor(5);
  ent_FONLL->SetMarkerColor(5);

  TLegendEntry *ent_DD=legend->AddEntry(gaeDD,"Tot data driven reference unc","ple");
  ent_DD->SetTextFont(42);
  ent_DD->SetLineColor(1);
  ent_DD->SetMarkerColor(1);

  TLine* lin0=new TLine(9,1,120,1);
  lin0->SetLineStyle(2);
  lin0->SetLineColor(1);
  lin0->SetLineWidth(3);
  lin0->Draw("same");
  legend->Draw("same");

  cratio->SaveAs(Form("Plots/cratio_ReferencePPOverFO_%s.pdf",fofrom.Data()));

}
