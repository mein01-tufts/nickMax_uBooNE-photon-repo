import os,sys
from cuts import *

def run_1g1p_reco_selection_cuts( eventTree, classificationThreshold, fiducialData, return_on_fail=True ):
    passes = True
    recoList = []
    recoTrackList = []
    recoProtonCount = 0
    cuts_passed = {"novertex":True,
                   "noPrimaryMuon":True,
                   "noPrimaryElectron":True,
                   "cutCosmicPixels":True,
                   "vertexInFV":True,
                   "noChargedPion":True,
                   "noMoreThan1Proton":True,
                   "only1Photon":True,
                   "showerFromCharge":True,
                   "primaryScore":True,
                   "cutLongTracks":True,
                   "cutMuonCompleteness":True,
                   "maxVertexIntimePixels":True,
                   "cutShowerCompleteness":True,
                   "AllCuts":True}

    #Selecting events with reco
    #See if the event has a vertex
    if recoNoVertex(eventTree) == False:
        cuts_passed["noVertex"] = False
        passes = False
        if return_on_fail:
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    #See if the event is neutral current
    if recoCutMuons(eventTree, classificationThreshold) == False:
        cuts_passed["noPrimaryMuon"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    if recoCutElectrons(eventTree, classificationThreshold) == False:
        cuts_passed["noPrimaryElectron"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    #Use Matt's Cosmic Cut
    if trueCutCosmic(eventTree) == False:
        cuts_passed["cutCosmicPixels"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    #Make sure the event is within the fiducial volume
    if recoFiducials(eventTree, fiducialData) == False:
        cuts_passed["vertexInFV"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    #Cut events with suitably energetic charged pions 
    if recoPion(eventTree, classificationThreshold) == False:
        cuts_passed["noChargedPion"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        
        
    #Cut events with far-travelling protons
    recoProtonCount = recoProton(eventTree, classificationThreshold)
    if recoProtonCount > 1:
        cuts_passed["noMoreThan1Proton"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    #See if there are any photons in the event - if so, list them
    recoList = recoPhotonListFiducial(fiducialData, eventTree, 10.0)

    #List all photons classified as tracks
    recoTrackList = recoPhotonListTracks(fiducialData, eventTree, classificationThreshold)
    #recoTrackList = []

    if len(recoList) + len(recoTrackList) != 1:
        cuts_passed["only1Photon"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    #Try cutting based on data for Shower from Charged
    if recoCutShowerFromChargeScore(eventTree, recoList, recoTrackList) == False:
        cuts_passed["showerFromCharge"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        
  
    #Cut based on primary score
    if recoCutPrimary(eventTree, recoList, recoTrackList) == False:
        cuts_passed["primaryScore"] = False
        passes = False
        if return_on_fail:
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount

    #Cut based on the presence of tracks over 20 cm
    if recoCutLongTracks(eventTree, fiducialData) == False:
        cuts_passed["cutLongTracks"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    #Cut based on the completeness of known Muons
    if recoCutMuonCompleteness(eventTree) == False:
        cuts_passed["cutMuonCompleteness"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount        

    if recoCutMaxInTime(eventTree, recoProtonCount) == False:
        cuts_passed["maxVertexIntimePixels"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount

    # Cut based on completeness
    if recoCutCompleteness(eventTree, recoList, recoTrackList)==False:
        cuts_passed["cutShowerCompleteness"] = False
        passes = False
        if return_on_fail:        
            return passes, cuts_passed, recoList, recoTrackList, recoProtonCount
    
    cuts_passed["AllCuts"] = passes
    return passes, cuts_passed, recoList, recoTrackList, recoProtonCount

