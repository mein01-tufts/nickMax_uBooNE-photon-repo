# Welcome to darkCuts: It's like cuts.py, but like..... dark

#      _.-'''-._
#     .'   .-'``|'.
#    /    /    -*- \
#   ;   <{      |   ;
#   |    _\ |       | 
#   ;   _\ -*- |    ;
#    \   \  | -*-  /
#     '._ '.__ |_.'
#        '-----'

# Some imported libraries we may need
import sys, argparse
import numpy as np
import ROOT as rt
import math

# Placeholder
def cut(eventTree):
    for i in range(len(eventTree)):
        print("hello world" + str(i))
    return True

#Goes through an event, iterates over all tracks and showers, and creates two lists for each particle - one for any entries in tracks with that particle, and one for any entries in showers
def particleTallies(eventTree):
    electronShowerList, photonShowerList, muonShowerList, pionShowerList, protonShowerList, electronTrackList, photonTrackList, muonTrackList, pionTrackList, protonTrackList = []

    for x in range(eventTree.nshowers):
        if abs(eventTree.showerRecoPID[x]) == 11:
            electronShowerList.append(x)
        elif eventTree.showerRecoPID[x] == 22:
            photonShowerList.append(x)
        elif eventTree.showerRecoPID[x] == 13:
            muonShowerList.append(x)
        elif abs(eventTree.showerRecoPID[x]) == 211:
            pionShowerList.append(x)
        elif eventTree.showerRecoPID[x] == 2212:
            protonShowerList.append(x)

    for x in range(eventTree.nTracks):
        if abs(eventTree.trackRecoPID[x]) == 11:
            electronTrackList.append(x)
        elif eventTree.trackRecoPID[x] == 22:
            photonTrackList.append(X)
        elif eventTree.trackRecoPID[x] == 13:
            muonTrackList.append(x)
        elif abs(eventTree.trackRecoPID[x]) == 211:
            pionTrackList.append(x)
        elif eventTree.trackRecoPID[x] == 2212:
            protonTrackList.append(x)

return electronShowerList, photonShowerList, muonShowerList, pionShowerList, protonShowerList, electronTrackList, photonTrackList, muonTrackList, pionTrackList, protonTrackList

def protonCut(protonShowerList, protonTrackList, tolerance = 0):
    keep = True
    if len(protonShowerList) + len(protonTrackList) > tolerance:
        keep = False
    return keep

def muonCut(muonShowerList, muonTrackList, tolerance = 0):
    keep = True
    if len(muonShowerList) + len(muonTrackList) > tolerance:
        keep = False
    return keep

def pionCut(pionShowerList, pionTrackList, tolerance = 0):
    keep = True
    if len(pionShowerList) + len(pionTrackList) > tolerance:
        keep = False
    return keep

def muonCut(photonShowerList, photonTrackList, tolerance = 0):
    keep = True
    if len(photonShowerList) + len(photonTrackList) > tolerance:
        keep = False
    return keep

def electronCut(muonShowerList, muonTrackList, rangeLow = 1, rangeHigh = 2):
    keep = True
    if len(muonShowerList) + len(muonTrackList) < rangeLow or len(muonShowerList) + len(muonTrackList) > rangeHigh:
        keep = False
    return keep