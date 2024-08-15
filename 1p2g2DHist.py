import sys, argparse
import numpy as np
import ROOT as rt

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, histStackFill, trueCCCutLoose, recoCCCutLoose

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gTracks.root", help="output root file name")
args = parser.parse_args()

# Grab ntuple file, eventTree, and potTree for reference
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

# Scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles

# Calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

# Define hists to be filled:
trueNCTrackHist = rt.TH1F("nctracks", "nc tracks per event",10,0,10)
trueCCTrackHist = rt.TH1F("cctracks", "cc tracks per event",10,0,10)
# Define fiducial width
fiducialWidth = 10

recoNCtrueNC = 0
recoNCtrueCC = 0
trueNCunclassified = 0
trueCCunclassified = 0
trueNCCclassified = 0
trueCCclassified = 0
trueNCtotalTracks = 0
trueCCtotalTracks = 0

## Event loop: Select reco signal, then of those events which are CC, 
## Plot each non-proton track in terms of its energy and length(stop length at det edge)
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
# Reco CC cut
    if recoCCCut(eventTree):
        continue
# Reco loose CC cut
#    if recoCCCutLoose(eventTree):
#        continue
# Reco fiducial cut
    if recoFiducialCut(eventTree, fiducialWidth):
        continue
# Reco cosmic cut
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        continue
# Reco pi+ cut
    if recoPiPlusCut(eventTree):
        continue
# Reco proton selection
    nRecoProtons, recoProtonTID = recoProtonSelection(eventTree)
    if nRecoProtons != 1:
        continue
# Reco photon selection
    recoPhotonTIDList, recoLeadingPhotonEnergy = recoPhotonSelection(eventTree, fiducialWidth) 
    if len(recoPhotonTIDList) != 2:
        continue

# true CC: cut all non-cc events
    if trueCCCut(eventTree) == False:
        recoNCtrueNC += 1
        eventTracks = 0
        for i in range(eventTree.nTracks):
            eventTracks += 1
            if eventTree.trackClassified[i] == 0:
                trueNCunclassified += 1
            elif eventTree.trackClassified[i] == 1:
                trueNCCclassified += 1
        trueNCtotalTracks += eventTree.nTracks
        trueNCTrackHist.Fill(eventTracks, eventTree.xsecWeight)
        continue

    recoNCtrueCC += 1

# Calculate track lengths, store track energies for all 
    eventTracks = 0
    for i in range(eventTree.nTracks):
        eventTracks += 1
        unclassifiedInEvent = 0
        if eventTree.trackPID[i] != 2212 and eventTree.trackClassified[i] == 1:
            x=2
        if eventTree.trackClassified[i] == 0:
            trueCCunclassified += 1
            unclassifiedInEvent += 1
        if eventTree.trackClassified[i] == 1:
            trueCCclassified += 1
    trueCCTrackHist.Fill(eventTracks, eventTree.xsecWeight)
        
    trueCCtotalTracks += eventTree.nTracks
        


trackCanvas, trackStack, trackLegend, trackInt = histStackFill("num of tracks per event", [trueNCTrackHist, trueCCTrackHist], "number of tracks: (", "num of tracks", "tracks per 6.67e+20POT", ntuplePOTsum)


print("of " + str(recoNCtrueCC) + " total reco nc true CC events, " + str(trueCCunclassified/trueCCtotalTracks*100) + " percent of tracks are unclassified" + str(trueCCtotalTracks) + " total tracks")
print("of " + str(recoNCtrueNC) + " total reco nc true NC events, " + str(trueNCunclassified/trueNCtotalTracks*100) + " percent of tracks are unclassified" + str(trueNCtotalTracks) + " total tracks")

outFile = rt.TFile(args.outfile, "RECREATE")
trackCanvas.Write("Num of Tracks Hist")