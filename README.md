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

* file location : /fonll/code

* usage : root -l -b -q ScaleCor.C+

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



**Sequence of running codes for y dependence**

1. Make the root file from dat files including FONLL expectation

* file location : /fonll/Code

* usage : root -l -b -q 'Bplusdsigmady_all.cc+(5)'

* optional parameters

  1. beam energy and pT range: 5(5TeV), 7(7TeV), 2(2.76TeV), 71(7TeV,(5,120)-CMS pp binning), 72(7TeV, (9,12)-ATLAS pp binning)

  2. is it with analysis binning or not(bin width : 1GeV) : true(analysis binned), false(fine binned)

* resulted file and values


