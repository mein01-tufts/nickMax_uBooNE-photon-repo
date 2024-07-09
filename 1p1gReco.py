import sys, argparse
import numpy as np
import ROOT as rt

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="NpNgRecoLogged.root", help="output root file name")
args = parser.parse_args()

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")

#set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

total_1G = 0
total_2G = 0
total_3G = 0
trueFiducial_1G = 0
trueFiducial_2G = 0
trueFiducial_3G = 0
trueCC_1G = 0
trueCC_2G = 0
trueCC_3G = 0

avgVtxDist = 0
totalEvents = 0
avgVtxDist_1G = 0
avgVtxDist_2G = 0
avgVtxDist_3G = 0
#begin loop over events in ntuple file
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
    if eventTree.foundVertex != 1:
        continue

#reco fiducial cut
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
        eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
        eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
        #if abs(eventTree.vtxDistToTrue) <= 3:
        continue

       
#reco cosmic cut - same as true
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        continue

#reco charged-current cut:        
#iterate through all tracks in event,
#look for non-secondary tracks identified as muons
    primaryMuonFound = False
    primaryElectronTrackFound = False
    for i in range(eventTree.nTracks):
        if eventTree.trackIsSecondary[i] == 0:
            if abs(eventTree.trackPID[i]) == 13:
                primaryMuonFound = True
            if abs(eventTree.trackPID[i]) == 11:
                primaryElectronTrackFound = True
#cut events w/ primary muons
    if primaryMuonFound or primaryElectronTrackFound:    
        continue
    
#reco charged-current cut phase 2:
#look for non-secondary showers identified as electrions
    primaryElectronFound = False
    for i in range(eventTree.nShowers):
        if eventTree.showerIsSecondary[i] == 0:
            if abs(eventTree.showerPID[i]) == 11:
                primaryElectronFound == True
#cut events w/ primary electrons
    if primaryElectronFound:
        continue

#reco pi+ finder:
    piPlusPresent = False
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 211:
            piPlusPresent = True
            break
    if piPlusPresent:
        continue

#reco protons: find and tally
    nProtonsPresent = 0
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 2212: 
            nProtonsPresent += 1
    if nProtonsPresent == 0:
        continue

#reco photons: find and tally
    nPhotonsPresent = 0
    for i in range(eventTree.nShowers):
        if eventTree.showerPID[i] == 22:
            nPhotonsPresent += 1
    if nPhotonsPresent == 0:
        continue

    avgVtxDist += abs(eventTree.vtxDistToTrue)
    totalEvents += 1

#------ Now we truth-match all reco signal ---------------------------#
    if nPhotonsPresent != 1:
        continue


#----- end of event loop ---------------------------------------------#