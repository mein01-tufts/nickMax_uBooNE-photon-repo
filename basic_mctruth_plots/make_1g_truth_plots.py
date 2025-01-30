import os,sys
import numpy as np
import ROOT as rt
sys.path.append("../")
import tools
from tools.plotmaker import PlotMaker


plotter = PlotMaker(targetpot=6.6e20)
plotter.load_ntuples()

outroot = rt.TFile("out_1gtruth_temp.root","recreate")
outroot.cd()

# true photon deposition
def event_has_truephoton_deposition( tree, is_mc ):
    if not is_mc:
        return True
    nuvtx = np.array( (tree.trueVtxX, tree.trueVtxY, tree.trueVtxZ) )

    n_vis_photon = 0

    for i in range(tree.nTrueSimParts):
        if tree.trueSimPartPDG[i]!=22:
            continue

        partvtx = np.array( (tree.trueSimPartX[i], tree.trueSimPartY[i], tree.trueSimPartZ[i]) )
        dvtx = nuvtx-partvtx
        dist = np.sqrt( (dvtx*dvtx).sum() )

        # check its a primary vertex
        if dist>0.5:
            continue

        # get the energy deposit
        pixsum = np.zeros(3) 
        pixsum[0] = tree.trueSimPartPixelSumUplane[i]*0.016
        pixsum[1] = tree.trueSimPartPixelSumVplane[i]*0.016
        pixsum[2] = tree.trueSimPartPixelSumYplane[i]*0.016
        #print(pixsum)
        pixsum_max = np.max(pixsum)

        if pixsum_max<10:
            continue

        n_vis_photon += 1

    if n_vis_photon>0:
        return True

    return False

# true photon deposition
def event_has_single_truephoton_deposition( tree, is_mc ):
    if not is_mc:
        return True
    nuvtx = np.array( (tree.trueVtxX, tree.trueVtxY, tree.trueVtxZ) )

    n_vis_photon = 0

    for i in range(tree.nTrueSimParts):
        if tree.trueSimPartPDG[i]!=22:
            continue

        partvtx = np.array( (tree.trueSimPartX[i], tree.trueSimPartY[i], tree.trueSimPartZ[i]) )
        dvtx = nuvtx-partvtx
        dist = np.sqrt( (dvtx*dvtx).sum() )

        # check its a primary vertex
        if dist>0.5:
            continue

        # get the energy deposit
        pixsum = np.zeros(3) 
        pixsum[0] = tree.trueSimPartPixelSumUplane[i]*0.016
        pixsum[1] = tree.trueSimPartPixelSumVplane[i]*0.016
        pixsum[2] = tree.trueSimPartPixelSumYplane[i]*0.016
        #print(pixsum)
        pixsum_max = np.max(pixsum)

        if pixsum_max<10:
            continue

        n_vis_photon += 1

    if n_vis_photon==1:
        return True

    return False

def select_fv_true( tree, is_mc ):
    if not is_mc:
        return True

    fv_dist = 10.0
    
    if tree.trueVtxX<fv_dist or tree.trueVtxX>256.0-fv_dist:
        return False
    if tree.trueVtxY<-116.5+fv_dist or tree.trueVtxY>116.5-fv_dist:
        return False
    if tree.trueVtxZ<fv_dist or tree.trueVtxZ>1036.0-fv_dist:
        return False
    
    return True

def event_has_single_truephoton_deposition_infv( tree, is_mc ):
    passes_single_photon = event_has_single_truephoton_deposition(tree, is_mc)  
    if not passes_single_photon:
        return False
    nuvtx = np.array( (tree.trueVtxX, tree.trueVtxY, tree.trueVtxZ) )
    infv = select_fv_true( tree, is_mc )
    return infv

    
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
hX = rt.TH1D("hX_visphoton", ";true vertex X (cm);", 75, -250.0, 500.0 )
hY = rt.TH1D("hY_visphoton", ";true vertex Y (cm);", 50, -250.0, 250.0 )
hZ = rt.TH1D("hZ_visphoton", ";true vertex Z (cm);", 150, -250.0, 1250.0 )
plotter.addHist( "hX", hX, fill_true_X, event_has_truephoton_deposition, get_default_weight )
plotter.addHist( "hY", hY, fill_true_Y, event_has_truephoton_deposition, get_default_weight )
plotter.addHist( "hZ", hZ, fill_true_Z, event_has_truephoton_deposition, get_default_weight )

hX1g = rt.TH1D("hX_onevisphoton", ";true vertex X (cm);", 75, -250.0, 500.0 )
hY1g = rt.TH1D("hY_onevisphoton", ";true vertex Y (cm);", 50, -250.0, 250.0 )
hZ1g = rt.TH1D("hZ_onevisphoton", ";true vertex Z (cm);", 50, -250.0, 1250.0 )
plotter.addHist( "hX1g", hX1g, fill_true_X, event_has_single_truephoton_deposition, get_default_weight )
plotter.addHist( "hY1g", hY1g, fill_true_Y, event_has_single_truephoton_deposition, get_default_weight )
plotter.addHist( "hZ1g", hZ1g, fill_true_Z, event_has_single_truephoton_deposition, get_default_weight )

hX1g_fv = rt.TH1D("hX_onevisphoton_fv", ";true vertex X (cm);", 75, -250.0, 500.0 )
hY1g_fv = rt.TH1D("hY_onevisphoton_fv", ";true vertex Y (cm);", 50, -250.0, 250.0 )
hZ1g_fv = rt.TH1D("hZ_onevisphoton_fv", ";true vertex Z (cm);", 50, -250.0, 1250.0 )
plotter.addHist( "hX1g_fv", hX1g_fv, fill_true_X, event_has_single_truephoton_deposition_infv, get_default_weight )
plotter.addHist( "hY1g_fv", hY1g_fv, fill_true_Y, event_has_single_truephoton_deposition_infv, get_default_weight )
plotter.addHist( "hZ1g_fv", hZ1g_fv, fill_true_Z, event_has_single_truephoton_deposition_infv, get_default_weight )

plotter.runloop()

outroot.Close()
