# DataDrivenPPreference

**Sequence of running codes for pT dependence**

1. Make the root file from dat files including FONLL expectation

* file location : /fonll/Code
 
* usage : root -l -b -q 'Bplusdsigmadpt_all.cc+(5)'

* optional parameters

  1. beam energy : 5(5TeV), 7(7TeV), 2(2.76TeV)

  2. is it with analysis binning or not(bin width : 1GeV) : true(analysis binned), false(fine binned)

* resulted file and values

2. Fit on pp data with their binning with power law function

* file location : /fonll/Code 
* usage : root -l -b -q ScaleCor_pt.C+ 
* parameters to check!

*** Fitting with 7TeV pp ATLAS  
  NO.   NAME      VALUE             
   1  p0           3.68067e+00  
   2  p1           3.80198e-01   
   3  p2          -3.76518e-02  
   4  p3          -1.01557e+00 

*** Fitting with 7TeV pp CMS 
  NO.   NAME      VALUE         
   1  p0           4.00284e+00  
   2  p1           4.62641e-01  
   3  p2          -2.42728e-02  
   4  p3          -2.36722e+00  

3. With fitting function, get the pp data with our binning

* file location : /fonll/Code/.xml  
* With some thought, calculated centrl value and errors for our binning 

4. Calculate the ratio of FONLL expectation (A vs. B) with our binning and related systematics

* file location : /fonll/Code 
* usage : root -l -b -q 'CompFONLL_BplusdsigmadptUnbinned.cc+(1,0,1,2)';root -l -b -q 'CompFONLL_BplusdsigmadptUnbinned.cc+(1,1,1,2)'
* optional parameters 
  1. isBinned(default:0) : binned with our analysis binning(1) or with 1GeV fine binning(0)
  2. isNorm(default:1) : normalized by central value(1) or central value itself(0)
  3. What is numerator? : 5TeV(0), 7TeV(1), 2.76TeV(2) 
  4. What is denominator? : 5TeV(0), 7TeV(1), 2.76TeV(2) 
* values to check! - save as histogram in root file 

  1. CompFONLL_Bplus_Binned_Val_7TeVvs2760GeV.root - with (1,0) : central value 
	**0.664924**+0.0171475-0.0126994 
	**0.631362**+0.013487-0.0108966 
	**0.60587**+0.0110387-0.00962627 
	**0.585093**+0.00921929-0.00869071 
	**0.55274**+0.00729316-0.00759625 

  2. CompFONLL_Bplus_Binned_Norm_7TeVvs2760GeV.root - with (1,1) : check plus minus error
	1+**0.0257887**-**0.019099**
	1+**0.0213618**-**0.0172589**
	1+**0.0182195**-**0.0158883**
	1+**0.015757**-**0.0148537**
	1+**0.0131946**-**0.0137429** 

5. Get the pp data-driven reference (pp+FONLL)

* file location : Bntuple/CrossSection/Analysis/
* usage : root -l GetBinned.C+ 
* result 
# 0 : 245.456 #pm 23.7691 
# 1 : 63.1326 #pm 6.42757 
# 2 : 20.576 #pm 3.52816 
# 3 : 8.15009 #pm 1.67695 
# 4 : 1.45249 #pm 0.276782 

6. Final results

* file location : Bntuple/CrossSection/Analysis/ 
* usage : root -l NuclearModification_DrawOnSamePad_v3.C+






**Sequence of running codes for y dependence**

1. Make the root file from dat files including FONLL expectation

* file location : /fonll/Code

* usage : root -l -b -q 'Bplusdsigmady_all.cc+(5)'

* optional parameters

  1. beam energy and pT range: 5(5TeV), 7(7TeV), 2(2.76TeV), 71(7TeV,(5,120)-CMS pp binning), 72(7TeV, (9,12)-ATLAS pp binning)

  2. is it with analysis binning or not(bin width : 1GeV) : true(analysis binned), false(fine binned)

* resulted file and values


