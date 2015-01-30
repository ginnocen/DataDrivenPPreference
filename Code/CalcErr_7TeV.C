#include <iostream>
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"

#define NUM 5
#define NUM_ATL 8
#define NUM_pPb 6

void iterpoints(int i, double oribin[], TF1* f, double stmpt[], double iterorder)
{
  //iterpoints(i, ptbin, fitft_opt2, stmpt, 2)
  double integ1,integ2,gapinteg,minpt;
  double mininteg=1000000000000000.0;
  double initval = pow(10,iterorder*(-1));
  //std::cout << "#### Inital value : " << initval << std::endl;
  if (iterorder==0) {
    for (double j=oribin[i]+1.0;j<oribin[i+1];j+=1.0){
      integ1 = f->Integral(oribin[i],j);
      integ2 = f->Integral(j,oribin[i+1]);
      gapinteg = fabs(integ1-integ2);
      //###std:: cout << j << " : " << integ1 << " , " << integ2 << std::endl;
      if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
    }
  }
  else {
    for (double j=stmpt[i]-initval*9;j<stmpt[i]+initval*10;j+=initval){
      integ1 = f->Integral(oribin[i],j);
      integ2 = f->Integral(j,oribin[i+1]);
      gapinteg = fabs(integ1-integ2);
      //###std:: cout << j << " : " << "----- " << integ1 << " , " << integ2 << std::endl;
      if (gapinteg<mininteg) {mininteg=gapinteg;minpt=j;}
    }
  }
  stmpt[i]=minpt;
  std::cout << i << " : " << std::string(iterorder*3+1, '-') << " gap : " << mininteg << ", minpt : " << minpt << std::endl;
}

void CalcErr_7TeV() {

	// CMS 7TeV, pp->(B+)X, unit : micro barn
	//double ptbin[NUM+1]={5.,10.,13.,17.,24.,30.};
	double sigmapt[NUM]={4.07,1.47,0.412,0.181,0.042};
	double sigmapt_sta[NUM]={0.47,0.13,0.041,0.015,0.007};
	double sigmapt_sys[NUM]={0.31,0.09,0.026,0.012,0.004};

	double sigmapt_sta_pPb[NUM_pPb];
	double sigmapt_sys_pPb[NUM_pPb];

	// ATLAS 7TeV, dsigma/dp_T (B+) * BR(B+->J/psiK+) * BR (J/Psi->mu+mu-), unit : pico barn
/*
	double ptbin_ATL[NUM_ATL+1]={9.,13.,16.,20.,25.,35.,50.,70.,120.};
	double sigmapt_ATL[NUM_ATL]={103.4,36.03,15.33,6.056,1.814,0.3477,0.06244,0.006099};
	double sigmapt_sta_ATL[NUM_ATL]={3.7,0.80,0.25,0.093,0.027,0.0084,0.00293,0.000561};
	double sigmapt_sys_ATL[NUM_ATL]={7.6,2.32,0.98,0.376,0.115,0.0280,0.00526,0.000666};

	double sigmapt_sta_ATL_pPb[NUM_pPb];
	double sigmapt_sys_ATL_pPb[NUM_pPb];
*/
	// Our pPb analysis binning
	double ptbin_pPb[NUM_pPb+1]={5.,10.,15.,20.,25.,30.,60.};

	TF1* fitft = new TF1 ("fitft","pow(10,[0]*exp([1]*x)+[2])",0.0,120.0);
	TFile* fin_fitpar = new TFile("../ResultsBplus/ScaleCor_reweighted7TeVdata.root");
	TH1D* hparafit = (TH1D*)fin_fitpar->Get("hparafit");

	TFile* fout = new TFile("../ResultsBplus/CalcErrPerc_7TeV.root","RECREATE");

	fitft->SetParameter(0,hparafit->GetBinContent(1));
	fitft->SetParameter(1,hparafit->GetBinContent(2));
	fitft->SetParameter(2,hparafit->GetBinContent(3));

  double stminpt_pPb[NUM_pPb];
  //double sty_pPb[NUM_pPb],stexl_pPb[NUM_pPb],stexh_pPb[NUM_pPb];
  //double steyl_sta_pPb[NUM_pPb],steyh_sta_pPb[NUM_pPb];
  //double steyl_sys_pPb[NUM_pPb],steyh_sys_pPb[NUM_pPb];

	//double ptbin_pPb[NUM_pPb+1]={5.,10.,15.,20.,25.,30.,60.};

  std::cout << std::endl;
  std::cout << std::string(50,'#') << std::endl;
  for (int i=0;i<NUM_pPb;i++){
    std::cout << i << " : " << ", central point : " << (ptbin_pPb[i+1]+ptbin_pPb[i])/2 << std::endl;
    iterpoints(i, ptbin_pPb, fitft, stminpt_pPb, 0);
    iterpoints(i, ptbin_pPb, fitft, stminpt_pPb, 1);
    iterpoints(i, ptbin_pPb, fitft, stminpt_pPb, 2);
/*
    sty_pPb[i]=hsigmapt_noerr_pPb->GetBinContent(i+1);
    stexl_pPb[i]=stminpt_pPb[i]-ptbin_pPb[i];
    stexh_pPb[i]=ptbin_pPb[i+1]-stminpt_pPb[i];
    steyl_sta_pPb[i]=hsigmapt_staerr_pPb->GetBinError(i+1);
    steyh_sta_pPb[i]=hsigmapt_staerr_pPb->GetBinError(i+1);
    steyl_sys_pPb[i]=hsigmapt_syserr_pPb->GetBinError(i+1);
    steyh_sys_pPb[i]=hsigmapt_syserr_pPb->GetBinError(i+1);
*/  
}
  std::cout << std::string(50,'#') << std::endl;

	sigmapt_sta_pPb[0]=(sigmapt_sta[0]/sigmapt[0]);//5,10
	sigmapt_sta_pPb[1]=((sigmapt_sta[1]/sigmapt[1])*3+(sigmapt_sta[2]/sigmapt[2])*2)/(ptbin_pPb[2]-ptbin_pPb[1]);//10,15
	sigmapt_sta_pPb[2]=((sigmapt_sta[2]/sigmapt[2])*2+(sigmapt_sta[3]/sigmapt[3])*3)/(ptbin_pPb[3]-ptbin_pPb[2]);//15,20
	sigmapt_sta_pPb[3]=((sigmapt_sta[3]/sigmapt[3])*4+(sigmapt_sta[4]/sigmapt[4])*1)/(ptbin_pPb[4]-ptbin_pPb[3]);//20,25
	sigmapt_sta_pPb[4]=((sigmapt_sta[4]/sigmapt[4])*5)/(ptbin_pPb[5]-ptbin_pPb[4]);//25,30
	sigmapt_sta_pPb[5]=((sigmapt_sta[4]/sigmapt[4])-(sigmapt_sta[3]/sigmapt[3]))/(stminpt_pPb[4]-stminpt_pPb[3])*(stminpt_pPb[5]-stminpt_pPb[3])+(sigmapt_sta[3]/sigmapt[3]);//30,60 rough estimation

	sigmapt_sys_pPb[0]=(sigmapt_sys[0]/sigmapt[0]);//5,10
	sigmapt_sys_pPb[1]=((sigmapt_sys[1]/sigmapt[1])*3+(sigmapt_sys[2]/sigmapt[2])*2)/(ptbin_pPb[2]-ptbin_pPb[1]);//10,15
	sigmapt_sys_pPb[2]=((sigmapt_sys[2]/sigmapt[2])*2+(sigmapt_sys[3]/sigmapt[3])*3)/(ptbin_pPb[3]-ptbin_pPb[2]);//15,20
	sigmapt_sys_pPb[3]=((sigmapt_sys[3]/sigmapt[3])*4+(sigmapt_sys[4]/sigmapt[4])*1)/(ptbin_pPb[4]-ptbin_pPb[3]);//20,25
	sigmapt_sys_pPb[4]=((sigmapt_sys[4]/sigmapt[4])*5)/(ptbin_pPb[5]-ptbin_pPb[4]);//25,30
	sigmapt_sys_pPb[5]=((sigmapt_sys[4]/sigmapt[4])-(sigmapt_sys[3]/sigmapt[3]))/(stminpt_pPb[4]-stminpt_pPb[3])*(stminpt_pPb[5]-stminpt_pPb[3])+(sigmapt_sys[3]/sigmapt[3]);//30,60 rough estimation

  TH1D* hstaerrPerc = new TH1D("hstaerrPrec","",NUM_pPb,ptbin_pPb);
  TH1D* hsyserrPerc = new TH1D("hsyserrPrec","",NUM_pPb,ptbin_pPb);
  TH1D* hWgtcen = new TH1D("hWgtcen","",NUM_pPb,ptbin_pPb);

	for (int i=0;i<NUM_pPb;i++) {
		hstaerrPerc->SetBinContent(i+1,sigmapt_sta_pPb[i]);
		hsyserrPerc->SetBinContent(i+1,sigmapt_sys_pPb[i]);
		hWgtcen->SetBinContent(i+1,stminpt_pPb[i]);
		std::cout << i << " --- sta.Err.Perc. : " << sigmapt_sta_pPb[i] << " , sys.Err.Perc. : " << sigmapt_sys_pPb[i] << std::endl;
	}
	fout->cd();
	hstaerrPerc->Write();
	hsyserrPerc->Write();
	hWgtcen->Write();
}
