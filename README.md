# DataDrivenPPreference

## **Sequence of running codes for pT dependence**

### **QUICK SERVICE** : Just run ./auto_refstudy_Bplus.sh 
After running, you can get the results for the case of 

(1). 5TeV pPb results (with pure FONLL and data-driven FONLL with 7 TeV pp data)

(2). 2.76TeV pp results (with pure FONLL and data-driven FONLL)
 
(3). 5TeV pPb results (with pure FONLL and data-driven FONLL with 2.76 TeV pp data) 
 
### **DETAILED INSTRUCTION**

#### 0. Overall notice 
* All the executive codes(.C or .cc) are in /Code 
* All the FONLL input file(.dat) are in /FONLLInputs 
* Except the step to get RpA and dsigma/dpt like public results, all the resulted root and pdf files are stored in /ResultsBplus 

1. Make the root file from dat files including FONLL expectation

* usage : root -l -b -q 'Bplusdsigmadpt_all.cc+(5)' 
* optional parameters 
  1. beam energy : 5(5TeV), 7(7TeV), 2(2.76TeV) 
  2. is it with analysis binning or not(bin width : 1GeV) : true(analysis binned), false(fine binned)
* resulted file 
	* ../ResultsBplus/outputBplus_Unbinned(Binned)_5TeV(7TeV, 2760GeV).root

2. Calculate the ratio of FONLL expectation (A vs. B) with our binning and related systematics

* usage : root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,0,0,1)';root -l -b -q 'CompFONLL_BplusdsigmadptUnbinned.cc+(1,1,0,1)'
* optional parameters 
  1. isBinned(default:0) : binned with our analysis binning(1) or with 1GeV fine binning(0)
  2. isNorm(default:1) : normalized by central value(1) or central value itself(0)
  3. What is numerator? : 5TeV(0), 7TeV(1), 2.76TeV(2) 
  4. What is denominator? : 5TeV(0), 7TeV(1), 2.76TeV(2) 
* resulted file
	* ../ResultsBplus/CompFONLL_Bplus_Binned(Fine)_Val(Norm)_5TeV(7TeV,2760GeV)vs7TeV(5TeV,2760GeV).root
	* ../ResultsBplus/CompFONLL_Bplus_Binned(Unbinned)_Val(Norm)_5TeV(7TeV,2760GeV)vs7TeV(5TeV,2760GeV).pdf - Ratio itself or relative errors from FONLL comparison
	 and valuesvalues to check! - save as histogram in root file 
* example in stored information in root file
  1. CompFONLL_Bplus_Binned_Val_7TeVvs2760GeV.root - with (1,0) : central value 
	**0.664924**+0.0171475-0.0126994 

	**0.631362**+0.013487-0.0108966 

  2. CompFONLL_Bplus_Binned_Norm_7TeVvs2760GeV.root - with (1,1) : check plus minus error 
	1+**0.0257887**-**0.019099** 

	1+**0.0213618**-**0.0172589** 

3. Fit on pp data with their binning with power law function reweighting center of bin and combining CMS+ATLAS data

* usage : root -l -b -q ScaleCor_pt.C+ 
* fitting parameters are stored in histogram 
* resulted file 
	* ../ResultsBplus/ScaleCor_reweighted7TeVdata.root 
	* ../ResultsBplus/cSigma_CompFit.pdf 
	* ../ResultsBplus/cSigma_CompReweightedFit.pdf 
	* ../ResultsBplus/cSigma_CompCombinedFit.pdf 

4. With fitting function, get the pp data with our binning

* usage : root -l -b -q 'CalcErr_7TeV.C+' 
* calculate statistical and systematical errors for our binning with considering trend of published 7 TeV data 
* resulted file 
	* ../ResultsBplus/CalcErrPerc_7TeV.root

5. Get the pp data-driven reference (pp+FONLL)

* usage 
	* For 5.02 TeV pPb, with 7 TeV pp data  : root -l -b -q 'GetBinned.C+(0,1)' 
	* For 2.76 TeV pp, with 7 TeV pp data  : root -l -b -q 'GetBinned.C+(2,1)'  
	* For 5.02 TeV pPb, with 2.76 TeV pp data  : root -l -b -q 'GetBinned_w2760GeVpp.C+(0,2)' 
* resulted file 
	* common results 
		* ../ResultsBplus/Estimatedpp5TeV(2760GeV)_with7TeV(2760GeV).root 
		* ../ResultsBplus/gaeEstimatedpp5TeV(2760GeV)_from7TeV(2760GeV).pdf  
	* only from GetBinned.C 
		* ../ResultsBplus/CompRecal5TeV(2760GeV)vs7TeV(2760GeV).pdf 
		* ../ResultsBplus/CompRecalRatio5TeV(2760GeV)vs7TeV(2760GeV).pdf  
	
6. Final results

* usage : root -l -b -q 'CrossSection_forPPref.C+(1~5)' 
* resulted file
	* ../ResultsBplus(Bplus_pp2760GeV)/dSigmadpt_Bplus(Bplus_ppF,Bplus_ppFw2760,Bplus_2760GeVpp,Bplus_2760GeVppF).root 
	* ../ResultsBplus(Bplus_pp2760GeV)/canvasSigmaBplus(Bplus_ppF,Bplus_ppFw2760,Bplus_2760GeVpp,Bplus_2760GeVppF).pdf 
	* ../ResultsBplus(Bplus_pp2760GeV)/canvasPPOverFONLLBplus(Bplus_ppF,Bplus_ppFw2760,Bplus_2760GeVpp,Bplus_2760GeVppF).pdf 

7. Comparison between with pure FONLL and data-driven FONLL

* usage : root -l -b -q 'Comp_dsigma.C+(1~3)'
* resulted file 
	* ../ResultsBplus/hrefcomp_Comp5TeV(2760GeV,5TeV_w2760GeV).pdf 
	* ../ResultsBplus/hRcomp_Comp5TeV(2760GeV,5TeV_w2760GeV).pdf 

8. Check some ratios between FONLL, data and fitting

* usage : root -l -b -q 'CheckScaleFac_v2.C+'

For ATLAS FONLL, root -l -b -q 'Bplusdsigmadpt_ATL.cc+(72,true,2)'



---------------------------------------------------------------------------------

## **Sequence of running codes for y dependence**

1. Make the root file from dat files including FONLL expectation

* usage : root -l -b -q 'Bplusdsigmady_all.cc+(5)' 
* optional parameters 
  1. beam energy and pT range: 5(5TeV), 7(7TeV), 2(2.76TeV), 71(7TeV,(5,120)-CMS pp binning), 72(7TeV, (9,12)-ATLAS pp binning)
  2. is it with analysis binning or not(bin width : 1GeV) : true(analysis binned), false(fine binned)
* resulted file 
	* ../ResultsBplus_y/outputBplusy_Unbinned(Binned)_5TeV(7TeV, 2760GeV)_rap30_pt1060(5120,9120,1060).root

2. Calculate the ratio of FONLL expectation (A vs. B) with our binning and related systematics

* usage : root -l -b -q 'CompFONLL_Bplusdsigmady.cc+(1,0,0,1)';root -l -b -q 'CompFONLL_BplusdsigmadptUnbinned.cc+(1,1,0,1)'
* optional parameters 
  1. isBinned(default:0) : binned with our analysis binning(1) or with 1GeV fine binning(0)
  2. isNorm(default:1) : normalized by central value(1) or central value itself(0)
  3. What is numerator? : 5TeV_rap30_pt1060(0), 7TeV_rap30_pt1060(1), 7TeV_rap30_pt5120(2), 7TeV_rap30_pt9120(3), 2760GeV_rap30_pt1060(4) 
  4. What is denominator? : 5TeV_rap30_pt1060(0), 7TeV_rap30_pt1060(1), 7TeV_rap30_pt5120(2), 7TeV_rap30_pt9120(3), 2760GeV_rap30_pt1060(4) 
* resulted file
	* ../ResultsBplus/CompFONLL_Bplus_Binned(Fine)_Val(Norm)_5TeV(7TeV,2760GeV)vs7TeV(5TeV,2760GeV).root
	* ../ResultsBplus/CompFONLL_Bplus_Binned(Unbinned)_Val(Norm)_5TeV(7TeV,2760GeV)vs7TeV(5TeV,2760GeV).pdf - Ratio itself or relative errors from FONLL comparison
	 and valuesvalues to check! - save as histogram in root file 
* example in stored information in root file
  1. CompFONLL_Bplus_Binned_Val_7TeVvs2760GeV.root - with (1,0) : central value 
	**0.664924**+0.0171475-0.0126994 

	**0.631362**+0.013487-0.0108966 

  2. CompFONLL_Bplus_Binned_Norm_7TeVvs2760GeV.root - with (1,1) : check plus minus error 
	1+**0.0257887**-**0.019099** 

	1+**0.0213618**-**0.0172589** 



