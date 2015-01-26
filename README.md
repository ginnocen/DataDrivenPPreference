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

 FCN=0.287191 FROM HESSE     STATUS=OK             23 CALLS         443 TOTAL

                     EDM=6.76274e-10    STRATEGY= 1      ERROR MATRIX ACCURATE 

  EXT PARAMETER                                   STEP         FIRST   

  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 

   1  p0           3.68067e+00   6.49876e+00   2.49686e-06   2.94775e-04

   2  p1           3.80198e-01   1.76055e+00   6.78578e-07   1.08482e-03

   3  p2          -3.76518e-02   1.07197e-02   5.35272e-08   2.76702e-03

   4  p3          -1.01557e+00   8.01984e-01   2.26098e-06   1.98043e-04

                               ERR DEF= 0.5


*** Fitting with 7TeV pp CMS 

 FCN=5.52214 FROM HESSE     STATUS=NOT POSDEF     23 CALLS         421 TOTAL

                     EDM=2.73951e-07    STRATEGY= 1      ERR MATRIX NOT POS-DEF

  EXT PARAMETER                APPROXIMATE        STEP         FIRST   

  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 

   1  p0           4.00284e+00   1.19329e-01   3.23700e-06   5.21426e-02

   2  p1           4.62641e-01   2.97749e-02   8.09033e-07   2.08719e-01

   3  p2          -2.42728e-02   1.19544e-03   8.57120e-08   1.67634e+00

   4  p3          -2.36722e+00   1.65639e-01   4.10697e-06   3.98973e-02

                               ERR DEF= 0.5



**Sequence of running codes for y dependence**

1. Make the root file from dat files including FONLL expectation

* file location : /fonll/Code

* usage : root -l -b -q 'Bplusdsigmady_all.cc+(5)'

* optional parameters

  1. beam energy and pT range: 5(5TeV), 7(7TeV), 2(2.76TeV), 71(7TeV,(5,120)-CMS pp binning), 72(7TeV, (9,12)-ATLAS pp binning)

  2. is it with analysis binning or not(bin width : 1GeV) : true(analysis binned), false(fine binned)

* resulted file and values


