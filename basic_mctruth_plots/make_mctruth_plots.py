import os,sys
sys.path.append("../")
import ROOT as rt

import tools
from tools.plotmaker import PlotMaker

plotter = PlotMaker(targetpot=6.6e20)
plotter.load_ntuples()

#outroot = rt.TFile("out_temp_edepfix.root","recreate")
outroot = rt.TFile("out_temp_run3_CV_500k.root","recreate")
outroot.cd()

# true Enu
# define selection function
def select_fv_true_numu( tree, is_mc ):
    if not is_mc:
        return True
    
    if abs(tree.trueNuPDG)!=14:
        return False

    fv_dist = 17.0
    
    if tree.trueVtxX<fv_dist or tree.trueVtxX>256.0-fv_dist:
        return False
    if tree.trueVtxY<-116.5+fv_dist or tree.trueVtxY>116.5-fv_dist:
        return False
    if tree.trueVtxZ<fv_dist or tree.trueVtxZ>1036.0-fv_dist:
        return False
    
    return True

# define fill variable
def fill_true_Enu( tree, is_mc ):
    if not is_mc:
        return 0.0
    return tree.trueNuE*1.0e3

# define weights
def get_default_weight( tree, is_mc ):
    if not is_mc:
        return 1.0
    
    w = tree.xsecWeight
    return w

# above threshold gamma deposition



# define histogram
hEnu = rt.TH1D("hEnu", ";true E_{#nu_{#mu}};", 200, 0, 20000.0 )
plotter.addHist( "trueEnumu", hEnu, fill_true_Enu, select_fv_true_numu, get_default_weight )
plotter.runloop()

outroot.Close()
