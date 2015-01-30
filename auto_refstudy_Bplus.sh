#!/bin/bash

### For 5.02 TeV pPb Comparison
cd Code;root -l -b -q 'Bplusdsigmadpt_all.cc+(5,false)'
cd Code;root -l -b -q 'Bplusdsigmadpt_all.cc+(7,false)'
cd Code;root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,0,0,1)'
cd Code;root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,1,0,1)'
cd Code;root -l -b -q 'ScaleCor_pt.C+'
cd Code;root -l -b -q 'CalcErr_7TeV.C+'
cd Code;root -l -b -q 'GetBinned.C+(0,1)'
cd Code;root -l -b -q 'CrossSection_forPPref.C+(1)'
cd Code;root -l -b -q 'CrossSection_forPPref.C+(2)'
cd Code;root -l -b -q 'Comp_dsigma.C+(1)'

### For 2.76 TeV Comparison
cd -
cd Code;root -l -b -q 'Bplusdsigmadpt_all.cc+(7,false)'
cd Code;root -l -b -q 'Bplusdsigmadpt_all.cc+(2,false)'
cd Code;root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,0,2,1)'
cd Code;root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,1,2,1)'
cd Code;root -l -b -q 'ScaleCor_pt.C+'
cd Code;root -l -b -q 'CalcErr_7TeV.C+'
cd Code;root -l -b -q 'GetBinned.C+(2,1)'
cd Code;root -l -b -q 'CrossSection_forPPref.C+(3)'
cd Code;root -l -b -q 'CrossSection_forPPref.C+(4)'
cd Code;root -l -b -q 'Comp_dsigma.C+(2)'

### For 5.02 TeV pPb with 2.76 TeV pp
cd -
cd Code;root -l -b -q 'Bplusdsigmadpt_all.cc+(5,false)'
cd Code;root -l -b -q 'Bplusdsigmadpt_all.cc+(2,false)'
cd Code;root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,0,0,2)'
cd Code;root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,1,0,2)'
cd Code;root -l -b -q 'GetBinned_w2760GeVpp.C+(0,2)'
cd Code;root -l -b -q 'CrossSection_forPPref.C+(1)'
cd Code;root -l -b -q 'CrossSection_forPPref.C+(5)'
cd Code;root -l -b -q 'Comp_dsigma.C+(3)'
