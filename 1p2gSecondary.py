import sys, argparse
import numpy as np
import ROOT as rt

from cuts import histStackFill, kineticEnergyCalculator, sStackFillS

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gSecondaries1.root", help="output root file name")
args = parser.parse_args()

# open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

# scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles

# calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

# define histograms to be filled

muPiScoreCC = rt.TH1F("muPiCC", "Muon+Pion Score of Reco NC, True CC Tracks",60,0,50)
muPiScoreNC = rt.TH1F("muPiNC", "Muon+Pion Score of Reco NC, True NC Tracks",60,0,50)


xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

totalTracks = 0
recoSignalEvents = 0
ccTracks = 0
ncTracks = 0
ccEvents, ncEvents = 0, 0
ncMuons, ccMuons = 0, 0
ncSecondaryMuons, ccSecondaryMuons = 0, 0
ncMuonsUnclassified, ccMuonsUnclassified = 0, 0
trueNCMuons, trueCCMuons = 0, 0

# begin loop over events in ntuple file:
# start with culling for reco-defined 1p2g signal, then reiterate using true to
# determine what reco actually saw
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

# event vertex reconstructed cut
    if eventTree.foundVertex != 1:
        continue

# fiducial volume cut
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
        eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
        eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
        continue

# reco "cosmic cut"
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        continue

# reco charged-current cut:        
# iterate through all tracks/showers in event,
# look for non-secondary tracks identified as muons or electrons
# and for non-secondary showers identified as electrons
    recoPrimaryMuonFound = False
    recoPrimaryElectronTrackFound = False
    recoPrimaryElectronShowerFound = False
    for i in range(eventTree.nTracks):
        if eventTree.trackIsSecondary[i] == 0:
            if abs(eventTree.trackPID[i]) == 13:
                recoPrimaryMuonFound = True
            if abs(eventTree.trackPID[i]) == 11:
                recoPrimaryElectronTrackFound = True
    for i in range(eventTree.nShowers):
        if eventTree.showerIsSecondary[i] == 0:
            if abs(eventTree.showerPID[i]) == 11:
                recoPrimaryElectronShowerFound == True
# cut all events with primary electrons/muons
    if recoPrimaryMuonFound or recoPrimaryElectronTrackFound or recoPrimaryElectronShowerFound:    
        continue

# reco pi+ cut: cut all events w/ pi+ of KE >= 30MeV
    recoPiPlusPresent = False
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 211:
                if eventTree.trackRecoE[i] >= 30:
                    recoPiPlusPresent = True
                    break
    if recoPiPlusPresent:
        continue

# reco protons cut: all protons of KE >= 60MeV 
    nRecoProtons = 0
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 2212:
            if eventTree.trackRecoE[i] >= 60.0:
                nRecoProtons += 1
                recoProtonTID = eventTree.trackTrueTID[i]
    if nRecoProtons != 1:
        continue

# reco photon cut: find all photons, then append truth-matched TID to a list
    recoPhotonTIDList = []
    recoPhotonIndexList = []
    for i in range(eventTree.nShowers):
        if eventTree.showerPID[i] == 22:
            recoPhotonTIDList.append(eventTree.showerTrueTID[i])
            recoPhotonIndexList.append(i)
# cut events w/o two reco photons    
    if len(recoPhotonTIDList) != 2:
        continue

    recoSignalEvents += 1
    totalTracks += eventTree.nTracks

    recoTrackScoreList = []
    recoTrackPIDList = []
    recoTrackIsSecondaryList = []
    recoTrackClassifiedList = []
    for i in range(eventTree.nTracks):
        muPiScore = abs(eventTree.trackMuScore[i])
        recoTrackScoreList.append(muPiScore)
        recoTrackPIDList.append(eventTree.trackTruePID[i])
        recoTrackIsSecondaryList.append(eventTree.trackIsSecondary[i])
        recoTrackClassifiedList.append(eventTree.trackClassified[i])

    #-------------- end reco selection --------------#
    #------------- begin truth matching -------------#

# true charged-current selection
    if eventTree.trueNuCCNC != 0:
        ncTracks += len(recoTrackScoreList)
        ncEvents += 1
        for i in range(len(recoTrackScoreList)):
            muPiScoreNC.Fill(recoTrackScoreList[i], eventTree.xsecWeight)
        for i in range(len(recoTrackPIDList)):
            if recoTrackPIDList[i] == 13:
                if recoTrackClassifiedList[i] == 0:
                    ncMuonsUnclassified += 1
                if recoTrackIsSecondaryList[i] == 1:
                    ncSecondaryMuons += 1
                ncMuons += 1
        for i in range(eventTree.nTrueSimParts):
            if eventTree.trueSimPartProcess[i] == 0 and eventTree.trueSimPartPDG[i] == 13:
                trueNCMuons += 1

    else:
        ccTracks += len(recoTrackScoreList)
        ccEvents += 1
        for i in range(len(recoTrackScoreList)):
            muPiScoreCC.Fill(recoTrackScoreList[i], eventTree.xsecWeight)
        for i in range(len(recoTrackPIDList)):
            if recoTrackPIDList[i] == 13:
                ccMuons += 1
                if recoTrackClassifiedList[i] == 0:
                    ccMuonsUnclassified += 1
                if recoTrackIsSecondaryList[i] == 1:
                    ccSecondaryMuons += 1
        for i in range(eventTree.nTrueSimParts):
            if eventTree.trueSimPartProcess[i] == 0 and eventTree.trueSimPartPDG[i] == 13:
                trueCCMuons += 1

#----- end of event loop ---------------------------------------------#

muPiCanvas, muPiStack, muPiLegend, muPiInt = histStackFill("Muon+PionScores", [muPiScoreNC, muPiScoreCC], "muon + pion score: (", "muon + pion score", "tracks per 6.67e+20POT")

outFile = rt.TFile(args.outfile, "RECREATE")
muPiCanvas.Write()

print(str(recoSignalEvents) + " reco signal events")
print(str(totalTracks) + " total reco signal tracks")
print(str(ccEvents) + " true cc events, with " + str(ccTracks) + " tracks")
print(str(ncEvents) + " true nc events, with " + str(ncTracks) + " tracks")
print(str(ncMuons)+" nc muon tracks, " + str(ncMuonsUnclassified) + " are unclassified")
print(str(trueNCMuons) + " true nc primary muons")
print(str(ccMuons)+" cc muon tracks, " + str(ccMuonsUnclassified) + " are unclassified")
print(str(trueCCMuons) + " true cc primary muons")