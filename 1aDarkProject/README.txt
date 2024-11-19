This is the directory that contains programs related to Nick and Max's dark matter analysis for PHY-0082, Fall 2024

Nick's Dark n-Tuple File Path
-d /cluster/tufts/wongjiradlabnu/ndahle01/ubphoton/1aDarkProject/FinalDarkMatterInput.root

Taritree's Dark n-tuple Path
-d /cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/ntuple_epem_DarkNu_BenchmarkD_reco_v3lmshowerkp_gen2ntuple_evisvertex.root

Current Standard n-Tuple File Path
-s /cluster/tufts/wongjiradlabnu/nutufts/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge_reco_v2me06_gen2ntuple_photon_edep_fix_preview_v2.root

Current Cosmic n-Tuple File Path
-c /cluster/tufts/wongjiradlabnu/ndahle01/ubphoton/CosmicInput.root

PASTE THIS INPUT STACK: 
-d /cluster/tufts/wongjiradlabnu/ndahle01/ubphoton/1aDarkProject/FinalDarkMatterInput.root -s /cluster/tufts/wongjiradlabnu/nutufts/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge_reco_v2me06_gen2ntuple_photon_edep_fix_preview_v2.root -c /cluster/tufts/wongjiradlabnu/ndahle01/ubphoton/CosmicInput.root



This file contains simulated darkNu events made with a different version of GENIE than the regular n-tuples use. 

Our coding goal is to design a series of cuts that look over the darkNu file (signal) and the regular n-tuple file (background)
and is able to parse the events correctly. 