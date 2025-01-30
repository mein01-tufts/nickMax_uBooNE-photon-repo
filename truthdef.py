import os,sys
from cuts import *

def truthdef_1gamma_cuts( eventTree, photonEDepThreshold, fiducialData, return_on_fail=True ):

    passes = False
    truthcuts = {"vertexFV":False,
                 "muonBelowThreshold":False,
                 "electronBelowThreshold":False,
                 "pionsBelowThreshold":False,
                 "protonsBelowThreshold":False,
                 "only1Photon":False,
                 "NEdepPhotons":0,
                 "truePhotonIDs":[],
                 "trueOpeningAngle":0.0,
                 "pionCount":0,
                 "protonCount":0}

    # count the number of true photons with energy deposits in the TPC
    truePhotonIDs = truePhotonList(eventTree, fiducialData, threshold=photonEDepThreshold)
    nEdepPhotons = len(truePhotonIDs)
    trueOpeningAngle = 0.0
    if nEdepPhotons==2:
        # WC inclusive opening angle acceptance
        trueOpeningAngle = trueTwoPhotonOpeningAngle( eventTree, truePhotonIDs[0], truePhotonIDs[1] )
        if trueOpeningAngle<20.0:
            nEdepPhotons=1
    truthcuts["truePhotonIDs"]=truePhotonIDs
    truthcuts["nEdepPhotons"]=nEdepPhotons
    truthcuts["trueOpeningAngle"]=trueOpeningAngle
    

    #Selecting events using truth
    if trueCutFiducials(eventTree, fiducialData) == False:
        if return_on_fail:
            return passes, truthcuts
    else:
        truthcuts["vertexFV"] = True
  
    if trueCutMuons(eventTree) == False:
        if return_on_fail:
            return passes, truthcuts
    else:
        truthcuts["muonBelowThreshold"] = True

    if trueCutElectrons(eventTree) == False:
        if return_on_fail:
            return passes, truthcuts
    else:
        truthcuts["electronBelowThreshold"] = True

    pionCount, protonCount = trueCutPionProton(eventTree)
    truthcuts["pionCount"] = pionCount
    truthcuts["protonCount"] = protonCount
    if pionCount > 0:
        if return_on_fail:
            return passes,truthcuts
    else:
        truthcuts["pionsBelowThreshold"]=True
        
    if protonCount>1:
        if return_on_fail:
            return passes,truthcuts
    else:
        truthcuts["protonsBelowThreshold"]=True

    if nEdepPhotons != 1:
        if return_on_fail:
            return passes, truthcuts
    else:
        truthcuts["only1Photon"]=True
    
    passes = True
    return passes, truthcuts
