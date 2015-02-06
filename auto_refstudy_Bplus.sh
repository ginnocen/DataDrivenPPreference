#!/bin/bash

### For 5.02 TeV pPb Comparison
cd CodeOrg
root -l -b -q 'Bplusdsigmadpt_pp_7TeV_PPbBin.cc'
root -l -b -q 'Bplusdsigmadpt_pp_5p02TeV_PPbBin.cc'
root -l -b -q 'Bplusdsigmadpt_pp_2p76TeV_PPbBin.cc'
cd ../Code
root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,0,0,1)'
root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,1,0,1)'
root -l -b -q 'ScaleCor_pt.C+'
root -l -b -q 'CalcErr_7TeV.C+'
root -l -b -q 'GetBinned.C+(0,1)'
root -l -b -q 'CrossSection_forPPref.C+(1)'
root -l -b -q 'CrossSection_forPPref.C+(2)'
root -l -b -q 'Comp_dsigma.C+(1)'

### For 2.76 TeV Comparison
root -l -b -q 'Bplusdsigmadpt_all.cc+(7,false)'
root -l -b -q 'Bplusdsigmadpt_all.cc+(2,false)'
root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,0,2,1)'
root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,1,2,1)'
root -l -b -q 'ScaleCor_pt.C+'
root -l -b -q 'CalcErr_7TeV.C+'
root -l -b -q 'GetBinned.C+(2,1)'
root -l -b -q 'CrossSection_forPPref.C+(3)'
root -l -b -q 'CrossSection_forPPref.C+(4)'
root -l -b -q 'Comp_dsigma.C+(2)'

### For 5.02 TeV pPb with 2.76 TeV pp
root -l -b -q 'Bplusdsigmadpt_all.cc+(5,false)'
root -l -b -q 'Bplusdsigmadpt_all.cc+(2,false)'
root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,0,0,2)'
root -l -b -q 'CompFONLL_Bplusdsigmadpt.cc+(1,1,0,2)'
root -l -b -q 'GetBinned_w2760GeVpp.C+(0,2)'
root -l -b -q 'CrossSection_forPPref.C+(1)'
root -l -b -q 'CrossSection_forPPref.C+(5)'
root -l -b -q 'Comp_dsigma.C+(3)'
