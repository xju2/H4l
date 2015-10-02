#!/bin/bash
base_dir="$GROUPEOSDIR/monoH4l/minitrees/mc15_13TeV_v1"
file_tag="mono2"
signal="${base_dir}/mc_341748_zphxx_mzp200_mx1_fullsim_${file_tag}.root"
ggF="${base_dir}/mc_PowhegPythia8EvtGen_CT10_AZNLOCTEQ6L1_ggH125_ZZ4lep_noTau_combined_${file_tag}.root"
VBF="${base_dir}/mc_PowhegPythia8EvtGen_CT10_AZNLOCTEQ6L1_VBFH125_ZZ4lep_noTau_combined_${file_tag}.root"
WH="${base_dir}/mc_Pythia8EvtGen_A14NNPDF23LO_WH125_ZZ4l_combined_${file_tag}.root"
ZH="${base_dir}/mc_Pythia8EvtGen_A14NNPDF23LO_ZH125_ZZ4l_combined_${file_tag}.root"
Zee="zee.list"
Zmumu="zmumu.list"
ZHllvv="${base_dir}/mc_Pythia8EvtGen_A14NNPDF23LO_ZH125_ZZllvv_combined.root"

#h4l_add_weight zz.list combined_zz.root
#h4l_add_weight $signal combined_zpxx_mzp200_mx1.root
#h4l_add_weight $ggF combined_ggH125.root
#h4l_add_weight $VBF combined_VBFH125.root

#h4l_add_weight $Zee combined_zee.root
#h4l_add_weight $Zmumu combined_zmumu.root

#h4l_add_weight $WH combined_WH125.root
#h4l_add_weight ${ZHllvv} combined_ZH125_llvv.root
#h4l_add_weight $ZH combined_ZH125.root

h4l_add_weight zh_lvlv.list combined_ZH125_lvlv.root
