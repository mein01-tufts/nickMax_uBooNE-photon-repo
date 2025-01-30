import os,sys
import numpy as np
import ROOT as rt
sys.path.append("../")
import tools
from tools.plotmaker import PlotMaker


plotter = PlotMaker(targetpot=6.6e20)
plotter.load_ntuples()

outroot = rt.TFile("out_gamma_truth_temp.root","recreate")
outroot.cd()

# true photon deposition
def event_has_truephoton( tree, is_mc ):
    if not is_mc:
        return True
    nuvtx = np.array( (tree.trueVtxX, tree.trueVtxY, tree.trueVtxZ) )

    n_primary_photon = 0

    for i in range(tree.nTrueSimParts):
        if tree.trueSimPartPDG[i]!=22:
            continue

        partvtx = np.array( (tree.trueSimPartX[i], tree.trueSimPartY[i], tree.trueSimPartZ[i]) )
        dvtx = nuvtx-partvtx
        dist = np.sqrt( (dvtx*dvtx).sum() )

        # check its a primary vertex
        if dist>0.5:
            # distance too far from vertex
            continue

        n_primary_photon += 1

    if n_primary_photon>0:
        return True

    return False

def select_gv_true( tree, is_mc ):
    """ defines a gamma volume based on distance outside of TPC -- but not used as interaction volume is the cryostat"""
    if not is_mc:
        return True

    if tree.trueVtxX<-60.0 or tree.trueVtxX>320.0:
        return False
    if tree.trueVtxY<-180.0 or tree.trueVtxY>160.0:
        return False
    if tree.trueVtxZ<-100.0 or tree.trueVtxZ>1100.0:
        return False
    
    return True

def event_has_single_truephoton_deposition_infv( tree, is_mc ):
    passes_hastrue_photon = event_has_truephoton(tree, is_mc)  
    if not passes_hastrue_photon:
        return False
    else:
        return True
    nuvtx = np.array( (tree.trueVtxX, tree.trueVtxY, tree.trueVtxZ) )
    ingv = select_gv_true( tree, is_mc )
    return ingv

    
# define fill variable
def fill_true_X( tree, is_mc ):
    if not is_mc:
        return 0.0
    return tree.trueVtxX

def fill_true_Y( tree, is_mc ):
    if not is_mc:
        return 0.0
    return tree.trueVtxY

def fill_true_Z( tree, is_mc ):
    if not is_mc:
        return 0.0
    return tree.trueVtxZ

# define weights
def get_default_weight( tree, is_mc ):
    if not is_mc:
        return 1.0
    
    w = tree.xsecWeight
    return w

# define histogram
hX = rt.TH1D("hX_truephoton", ";true vertex X (cm);", 75,  -250.0, 500.0 )
hY = rt.TH1D("hY_truephoton", ";true vertex Y (cm);", 50,  -250.0, 250.0 )
hZ = rt.TH1D("hZ_truephoton", ";true vertex Z (cm);", 150, -250.0, 1250.0 )
plotter.addHist( "hX", hX, fill_true_X, event_has_truephoton, get_default_weight )
plotter.addHist( "hY", hY, fill_true_Y, event_has_truephoton, get_default_weight )
plotter.addHist( "hZ", hZ, fill_true_Z, event_has_truephoton, get_default_weight )

plotter.runloop()

outroot.Close()
