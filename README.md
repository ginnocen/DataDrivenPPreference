# DataDrivenPPreference

- Sequence of running codes for pT dependence

1. Make the rrot file from dat files including FONLL expectation
 
* usage : root -l -b -q 'Bplusdsigmadpt_all.cc+(5)'

* optional parameters

  1. beam energy : 5(5TeV), 7(7TeV), 2(2.76TeV)

  2. is it with analysis binning or not(bin width : 1GeV) : true(analysis binned), false(fine binned)

resulted file and values


- Sequence of running codes for y dependence

1. Make the root file from dat files including FONLL expectation

* usage : root -l -b -q 'Bplusdsigmady_all.cc+(5)'

* optional parameters

  1. beam energy and pT range: 5(5TeV), 7(7TeV), 2(2.76TeV), 71(7TeV,(5,120)-CMS pp binning), 72(7TeV, (9,12)-ATLAS pp binning)

  2. is it with analysis binning or not(bin width : 1GeV) : true(analysis binned), false(fine binned)

* resulted file and values


