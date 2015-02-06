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


TString fofrom = "7TeV";
//TString fofrom = "2760GeV";

void NuclearModificationDD(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  TFile*filePPReference=new TFile(Form("../ResultsBplus/Estimatedpp5TeV_with%s.root",fofrom.Data()));  
  TGraphAsymmErrors*gaeBplusReferenceOrigin=(TGraphAsymmErrors*)filePPReference->Get("gaeEstimatedpp5TeV");
  gaeBplusReferenceOrigin->SetName(Form("gae%sReferenceOrigin",particle.Data()));
  
  TFile*filepPb=new TFile(Form("../Results%s/Sigma%s.root",particle.Data(),particle.Data()));
  TH1F*hSigmapPbStat=(TH1F*)filepPb->Get("hPtSigma");  
  TH1F*hPt=(TH1F*)filepPb->Get("hPt");
  TH1F*hEff=(TH1F*)filepPb->Get("hEff");
  
  double scalingfactor=1.;
  //double scalingfactor=1e-6;
  
  double yvalue,xvalue,yerrorhigh,yerrorlow;

  double ytem[nbins],ylow[nbins],yhigh[nbins];
  for (int i=0;i<nbins;i++)
    {
      hSigmapPbStat->SetBinContent(i+1,(1.e-6)*hSigmapPbStat->GetBinContent(i+1));
      hSigmapPbStat->SetBinError(i+1,(1.e-6)*hSigmapPbStat->GetBinError(i+1));
      
      yvalue=-1.;
      xvalue=-1.;
      yerrorhigh=-1.;
      yerrorlow=-1.;

      if(fofrom=="7TeV")
	{
	  gaeBplusReferenceOrigin->GetPoint(i+1,xvalue,yvalue);
	  yerrorhigh=gaeBplusReferenceOrigin->GetEYhigh()[i+1];
	  yerrorlow=gaeBplusReferenceOrigin->GetEYlow()[i+1];
	}
      if(fofrom=="2760GeV")
	{
	  gaeBplusReferenceOrigin->GetPoint(i,xvalue,yvalue);
	  yerrorhigh=gaeBplusReferenceOrigin->GetEYhigh()[i];
	  yerrorlow=gaeBplusReferenceOrigin->GetEYlow()[i];
	}
      ytem[i] = yvalue*scalingfactor*208;
      ylow[i] = yerrorhigh*scalingfactor*208;
      yhigh[i] = yerrorlow*scalingfactor*208;
    } 
  TGraphAsymmErrors *gaeBplusReference = new TGraphAsymmErrors(nbins,xbins,ytem,exl,exl,ylow,yhigh);  
  
  for (int i=0;i<nbins;i++)
    {
      hSigmapPbStat->SetBinContent(i+1,(1./tagandprobcorrection[i])*hSigmapPbStat->GetBinContent(i+1));
      hSigmapPbStat->SetBinError(i+1,(1./tagandprobcorrection[i])*hSigmapPbStat->GetBinError(i+1));
    } 

  Double_t yRefPP[nbins];                        //value y reference
  Double_t xRefPP[nbins];                        //value x reference
  Double_t yPPsystFONLLhigh[nbins];              //y err syst FONLL high
  Double_t yPPsystFONLLlow[nbins];               //y err syst FONLL low
  Double_t yPercPPsystFONLLhigh[nbins];          //y percentuale err syst FONLL high
  Double_t yPercPPsystFONLLlow[nbins];           //y percentuale err syst FONLL low
    
  Double_t ySigmapPb[nbins];                     //value y pPb 
  Double_t xSigmapPb[nbins];                     //value x pPb
  Double_t ySigmapPbStat[nbins];                 //y err stat pPb
  Double_t yPercSigmapPbStat[nbins];             //y err stat pPb
  
  Double_t yFONLL[nbins];                        //1
  Double_t yRpA[nbins];                          //value y RpA 
  Double_t yRpAStat[nbins];                      //y err stat RpA 
  Double_t yRpAsystFONLLhigh[nbins];             //y err syst FONLL RpA high
  Double_t yRpAsystFONLLlow[nbins];              //y err syst FONLL RpA lzow
  Double_t yPercRpAsystFONLLhigh[nbins];         //y percentuale err syst FONLL RpA high
  Double_t yPercRpAsystFONLLlow[nbins];          //y percentuale err syst FONLL RpA low
    
  Double_t ySigmapPbSystTotHigh[nbins];              //y percentuale err syst pPb TOT
  Double_t ySigmapPbSystTotLow[nbins];              //y percentuale err syst pPb TOT

  Double_t yPercRpPbSystTotHigh[nbins];          //y percentuale err syst RpPb TOT
  Double_t yPercRpPbSystTotLow[nbins];          //y percentuale err syst RpPb TOT
  
  Double_t yRpPbSystTotHigh[nbins];              //y percentuale err syst RpPb TOT
  Double_t yRpPbSystTotLow[nbins];              //y percentuale err syst RpPb TOT

  
  double x,y;
  for (Int_t i=0;i<nbins;i++) {
    gaeBplusReference->GetPoint(i,xRefPP[i],yRefPP[i]);
    yPPsystFONLLhigh[i]=gaeBplusReference->GetEYhigh()[i];
    yPPsystFONLLlow[i]=gaeBplusReference->GetEYlow()[i];
    yPercPPsystFONLLhigh[i]=yPPsystFONLLhigh[i]/yRefPP[i];
    yPercPPsystFONLLlow[i]=yPPsystFONLLlow[i]/yRefPP[i];
    yPercPPsystFONLLhigh[i]=TMath::Sqrt(yPercPPsystFONLLhigh[i]*yPercPPsystFONLLhigh[i]+FFsysterror*FFsysterror);
    yPercPPsystFONLLlow[i]=TMath::Sqrt(yPercPPsystFONLLlow[i]*yPercPPsystFONLLlow[i]+FFsysterror*FFsysterror);
    
  }
  
  for(Int_t i=0;i<nbins;i++) {
    ySigmapPb[i]=hSigmapPbStat->GetBinContent(i+1);
    ySigmapPbStat[i]=hSigmapPbStat->GetBinError(i+1);
    yPercSigmapPbStat[i]=ySigmapPbStat[i]/ySigmapPb[i];
    ySigmapPbSystTotHigh[i]=yPercSigmapPbSystTotHigh[i]*ySigmapPb[i];
    ySigmapPbSystTotLow[i]=yPercSigmapPbSystTotLow[i]*ySigmapPb[i];
  }
  
  for(Int_t i=0;i<nbins;i++) {
    yRpA[i]=ySigmapPb[i]/yRefPP[i];
    yRpAStat[i]=ySigmapPbStat[i]/yRefPP[i];
    yFONLL[i]=yRpA[i];
    yPercRpAsystFONLLhigh[i]=(yPercPPsystFONLLlow[i]/(1-yPercPPsystFONLLlow[i]));
    yPercRpAsystFONLLlow[i]=(yPercPPsystFONLLhigh[i]/(1+yPercPPsystFONLLhigh[i]));
    yRpAsystFONLLhigh[i]=yPercRpAsystFONLLhigh[i]*yRpA[i];
    yRpAsystFONLLlow[i]=yPercRpAsystFONLLlow[i]*yRpA[i];

    yRpPbSystTotHigh[i]=yPercSigmapPbSystTotHigh[i]*yRpA[i];
    yRpPbSystTotLow[i]=yPercSigmapPbSystTotLow[i]*yRpA[i];
    //cout<<yRpPbSystTot[i]<<endl;
    
  }
  
  TGraphAsymmErrors *gSigmasyst = new TGraphAsymmErrors(nbins,xbins,ySigmapPb,exl,exl,ySigmapPbSystTotLow,ySigmapPbSystTotHigh);
  gSigmasyst->SetTitle("Sigma syst uncertainty from pPb");
  gSigmasyst->SetMarkerColor(1);
  gSigmasyst->SetLineColor(1);
  gSigmasyst->SetLineWidth(2);   
  gSigmasyst->SetMarkerStyle(21);
  gSigmasyst->SetMarkerColor(1);

  TGraphAsymmErrors *gSigmastat = new TGraphAsymmErrors(nbins,xbins,ySigmapPb,exl,exl,ySigmapPbStat,ySigmapPbStat);
  gSigmastat->SetTitle("Sigma stat uncertainty from pPb");
  gSigmastat->SetMarkerColor(1);
  gSigmastat->SetLineColor(1);
  gSigmastat->SetLineWidth(1);   
  gSigmastat->SetMarkerStyle(21);
  gSigmastat->SetMarkerColor(1);
  
  gSigmastat->SetFillColor(0);
  gSigmastat->SetFillStyle(0);


  TCanvas *canvasSigma=new TCanvas("canvasSigma","canvasSigma",500,500);   
  canvasSigma->cd();
  canvasSigma->Range(-1.989924,-0.2917772,25.49622,2.212202);
  canvasSigma->SetFillColor(0);
  canvasSigma->SetBorderMode(0);
  canvasSigma->SetBorderSize(2);
  canvasSigma->SetLeftMargin(0.1451613);
  canvasSigma->SetRightMargin(0.05443548);
  canvasSigma->SetTopMargin(0.08474576);
  canvasSigma->SetBottomMargin(0.1165254);
  canvasSigma->SetFrameBorderMode(0);
  canvasSigma->SetFrameBorderMode(0);
  canvasSigma->SetLogy();
  
  TH2F* hempty=new TH2F("hempty","",10,0.,70,10.,0.1,1e3.);  
  hempty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //if(particle=="Bplus") hempty->GetYaxis()->SetTitle("d#sigma / dp_{T} (B^{+}) (pb GeV^{-1}c)");
  //if(particle=="Bzero") hempty->GetYaxis()->SetTitle("d#sigma / dp_{T} (B^{0}) (pb GeV^{-1}c)");
  //if(particle=="Bs") hempty->GetYaxis()->SetTitle("d#sigma / dp_{T} (B_{s}) (pb GeV^{-1}c)");
  
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();

  
  hempty->GetYaxis()->SetTitle("d#sigma / dp_{T}( #mub GeV^{-1}c)");
  

  hempty->GetXaxis()->SetTitleOffset(1.);
  hempty->GetYaxis()->SetTitleOffset(1.3);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.04);  
  hempty->SetMaximum(2);
  hempty->SetMinimum(0.);
  hempty->Draw();
    
  
  gaeBplusReference->SetMarkerColor(1);
  gaeBplusReference->SetMarkerStyle(21);  
  gaeBplusReference->SetFillColor(5);
  gaeBplusReference->SetFillStyle(1001);
  gaeBplusReference->SetLineColor(1);
  gaeBplusReference->SetLineWidth(5);
  
  
  gSigmastat->SetMarkerColor(1);
  gSigmastat->SetLineColor(1);
  gSigmastat->SetLineWidth(2);   
  gSigmastat->SetMarkerStyle(21);
  gSigmastat->SetMarkerColor(1);

  gaeBplusReference->Draw("2same");
  gSigmastat->SetFillColor(0);
  gSigmastat->Draw("epsame");  

  //coord. for B+ and B0
  //TLegend *legendSigma=new TLegend(0.5745968,0.4756871,0.8729839,0.6490486,"");
  //coord.  B0
  TLegend *legendSigma=new TLegend(0.5100806,0.5868644,0.8084677,0.7605932,"");
  //coord. for B_s
  //TLegend *legendSigma=new TLegend(0.2580645,0.6122881,0.5564516,0.7860169,"");
  legendSigma->SetBorderSize(0);
  legendSigma->SetLineColor(0);
  legendSigma->SetFillColor(0);
  legendSigma->SetFillStyle(1001);
  legendSigma->SetTextFont(42);
  legendSigma->SetTextSize(0.045);
  
  TBox *c = new TBox(3,1-commonErrorN,7,1+commonErrorP);
  c->SetLineColor(5);
  c->SetFillColor(5);
  c->Draw();

  
  //TLegendEntry *ent_Sigmapp=legendSigma->AddEntry(gaeBplusReference,"pp reference","PLF");
  
  TLegendEntry *ent_SigmapPb=legendSigma->AddEntry(gSigmastat,"pPb","pf");
  ent_SigmapPb->SetTextFont(42);
  ent_SigmapPb->SetLineColor(1);
  //ent_SigmapPb->SetFillColor(0);
  ent_SigmapPb->SetMarkerColor(1);
  
  //TLegendEntry *ent_SigmapPbSyst=legendSigma->AddEntry(gSigmasyst,"pPb","f");
  //ent_SigmapPbSyst->SetTextFont(42);
  //ent_SigmapPbSyst->SetLineColor(1);
  //ent_SigmapPbSyst->SetMarkerColor(1);
  
  TLegendEntry *ent_Sigmapp=legendSigma->AddEntry(c,"FONLL pp ref.","f");
  ent_Sigmapp->SetTextFont(42);
  ent_Sigmapp->SetLineColor(5);
  ent_Sigmapp->SetMarkerColor(1);

  legendSigma->Draw("same");
  
  gSigmasyst->SetFillColor(0);
  gSigmasyst->SetFillStyle(0);
  
  //hSigmapPbStat->Draw("same");
  gSigmasyst->SetFillColor(0);
  gSigmasyst->SetFillStyle(0);
  gSigmasyst->Draw("2same");

  TBox *d = new TBox(3,1-commonErrorN,7,1+commonErrorP);
  d->SetLineColor(1);
  d->SetFillColor(0);
  d->Draw();
    
  TLatex * tlatex1=new TLatex(0.1612903,0.8625793,"CMS Preliminary     pPb #sqrt{s_{NN}}= 5.02 TeV");
  tlatex1->SetNDC();
  tlatex1->SetTextColor(1);
  tlatex1->SetTextFont(42);
  tlatex1->SetTextSize(0.045);
  tlatex1->Draw();

  TString mypar;
  if(particle=="Bplus") mypar="B^{+}";
  if(particle=="Bzero") mypar="B^{0}";
  if(particle=="Bs") mypar="B_{s}";
  
  TLatex * tlatexlumi=new TLatex(0.671371,0.7801268,"L = 34.8 nb^{-1}");
  tlatexlumi->SetNDC();
  tlatexlumi->SetTextColor(1);
  tlatexlumi->SetTextFont(42);
  tlatexlumi->SetTextSize(0.045);
  tlatexlumi->Draw();

 double xpos=0.8528226;
 double ypos=0.6849894;
  
  TLatex * tlatex3=new TLatex(xpos,ypos,mypar.Data());
   tlatex3->SetNDC();
  tlatex3->SetTextColor(1);
  tlatex3->SetTextFont(42);
  tlatex3->SetTextSize(0.06);
  tlatex3->Draw();
  
  canvasSigma->SaveAs(Form("Plots/canvasSigma%s_with%s.pdf",particle.Data(),particle.Data(),fofrom.Data()));  
  
  //TGraphAsymmErrors *gRpAstat = new TGraphAsymmErrors(nbins,xbins,yRpA,exl,exl,yRpAStat,yRpAStat);
  TGraphAsymmErrors *gRpAstat = new TGraphAsymmErrors(nbins,xbins,yRpA,exl0,exl0,yRpAStat,yRpAStat);
  gRpAstat->SetTitle("RpA stat uncertainty from pPb");
  gRpAstat->SetMarkerColor(1);
  gRpAstat->SetLineColor(1);
  gRpAstat->SetLineWidth(2);   
  gRpAstat->SetMarkerStyle(21);
  gRpAstat->SetMarkerColor(1);
  
  TGraphAsymmErrors *gRpAsyst = new TGraphAsymmErrors(nbins,xbins,yRpA,exl,exl,yRpPbSystTotLow,yRpPbSystTotHigh);
  gRpAsyst->SetTitle("RpA syst uncertainty from pPb");
  gRpAsyst->SetMarkerColor(1);
  gRpAsyst->SetLineColor(1);
  gRpAsyst->SetLineWidth(2);   
  gRpAsyst->SetMarkerStyle(21);
  gRpAsyst->SetMarkerColor(1);

   
  TGraphAsymmErrors *gRpAsystFONLL = new TGraphAsymmErrors(nbins,xbins,yFONLL,exl,exl,yRpAsystFONLLlow,yRpAsystFONLLhigh);
  gRpAsystFONLL->SetTitle("RpA syst uncertainty from FONLL reference");
  gRpAsystFONLL->SetMarkerColor(2);
  gRpAsystFONLL->SetLineColor(2);
  gRpAsystFONLL->SetLineWidth(2); 
  gRpAsystFONLL->SetMarkerStyle(21);
  gRpAsystFONLL->SetMarkerColor(2);

  TCanvas *canvasRpA=new TCanvas("canvasRpA","canvasRpA",500,500);   
  
  canvasRpA->Range(-1.989924,-0.2917772,25.49622,2.212202);
  canvasRpA->SetFillColor(0);
  canvasRpA->SetBorderMode(0);
  canvasRpA->SetBorderSize(2);
  canvasRpA->SetLeftMargin(0.1451613);
  canvasRpA->SetRightMargin(0.05443548);
  canvasRpA->SetTopMargin(0.08474576);
  canvasRpA->SetBottomMargin(0.1165254);
  canvasRpA->SetFrameBorderMode(0);
  canvasRpA->SetFrameBorderMode(0);
  
  TLegend *legendRpA=new TLegend(0.1975806,0.6109937,0.4959677,0.8012685,"");
  legendRpA->SetBorderSize(0);
  legendRpA->SetLineColor(0);
  legendRpA->SetFillColor(0);
  legendRpA->SetFillStyle(1001);
  legendRpA->SetTextFont(42);
  legendRpA->SetTextSize(0.045);

  TH2F* hempty=new TH2F("hempty","",10,0.,70,10.,0.,2.5);  
  
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  
  hempty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hempty->GetYaxis()->SetTitle("R^{FONLL}_{pA}");
  
  hempty->GetXaxis()->SetTitleOffset(1.1);
  hempty->GetYaxis()->SetTitleOffset(1.1);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.045);
  hempty->GetYaxis()->SetLabelSize(0.045);  
  hempty->SetMaximum(2);
  hempty->SetMinimum(0.);
  hempty->Draw();
  
  TLine *l = new TLine(0,1,70,1);
  l->SetLineStyle(2);
  
  legendRpA->Draw();
  gRpAstat->SetMarkerStyle(21);
  gRpAstat->SetLineColor(1);
  gRpAstat->SetMarkerColor(1);

  gRpAsystFONLL->SetLineStyle(3);
  gRpAsystFONLL->SetLineWidth(3);
  gRpAsystFONLL->SetFillColor(5);
  gRpAsystFONLL->SetLineColor(5);//kAzure-3);
  gRpAsystFONLL->SetMarkerColor(5);//kAzure-3);
  gRpAsystFONLL->Draw("2same");
  gRpAstat->Draw("psame");
  gRpAsyst->SetFillColor(0);
  gRpAsyst->SetFillStyle(0);
  
  gRpAstat->SetFillColor(0);
  gRpAstat->SetFillStyle(0);
  gRpAsyst->Draw("2same");
  
  TBox *a = new TBox(3,1-commonErrorN,7,1+commonErrorP);
  a->SetLineColor(1);
  a->SetFillColor(0);
  a->Draw();


  TBox *b = new TBox(3,1-commonErrorN,7,1+commonErrorP);
  b->SetLineColor(1);
  b->SetFillColor(kGray);
  b->Draw();


  TLegendEntry *ent_RpAstat=legendRpA->AddEntry(gRpAstat,"R^{FONLL}_{pA}","pf");
  ent_RpAstat->SetTextFont(42);
  ent_RpAstat->SetLineColor(2);
  ent_RpAstat->SetMarkerColor(2);
  //TLegendEntry *ent_RpAsystData=legendRpA->AddEntry(gRpAsyst,"            syst. unc.","f");
  //ent_RpAsystData->SetTextFont(42);
  //ent_RpAsystData->SetLineColor(1);
  //ent_RpAsystData->SetMarkerColor(1);
  TLegendEntry *ent_RpAsystData=legendRpA->AddEntry(b,"Syst. L+BR","f");
  ent_RpAsystData->SetTextFont(42);
  ent_RpAsystData->SetLineColor(2);
  ent_RpAsystData->SetMarkerColor(2);


  
  TLegendEntry *ent_RpAsystFONLL=legendRpA->AddEntry(gRpAsystFONLL,"Syst. err. from FONLL pp ref.","f");
  ent_RpAsystFONLL->SetTextFont(42);
  ent_RpAsystFONLL->SetLineColor(5);
  ent_RpAsystFONLL->SetLineStyle(1);
  ent_RpAsystFONLL->SetMarkerColor(5);
  
  legendRpA->Draw();
  
  TLatex * tlatex1=new TLatex(0.1612903,0.8625793,"CMS Preliminary     pPb #sqrt{s_{NN}}= 5.02 TeV");
  tlatex1->SetNDC();
  tlatex1->SetTextColor(1);
  tlatex1->SetTextFont(42);
  tlatex1->SetTextSize(0.045);
  tlatex1->Draw();
  
  TLatex * tlatex2=new TLatex(0.671371,0.7801268,"L = 34.8 nb^{-1}");
  tlatex2->SetNDC();
  tlatex2->SetTextColor(1);
  tlatex2->SetTextFont(42);
  tlatex2->SetTextSize(0.045);
  tlatex2->Draw();
  
  tlatex3->Draw();

//  l->Draw();  
  canvasRpA->SaveAs(Form("Plots/canvasRpA%s_%s.pdf",particle.Data(),particle.Data(),fofrom.Data()));  
  
  TFile *fout=new TFile(Form("Rootf/fileRpA%s_%s.root",particle.Data(),particle.Data(),fofrom.Data()),"recreate");  
  fout->cd();
  gSigmasyst->SetName("gSigmasyst");
  gSigmasyst->Write();
  gaeBplusReference->Write();
  hSigmapPbStat->Write();
}
