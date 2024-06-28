import sys, argparse
import numpy as np
import ROOT as rt

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="NpNgOutput.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
args = parser.parse_args()

#needed for proper scaling of error bars:
#rt.TH1.SetDefaultSumw2(rt.kTRUE)

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

#scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles

#calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT


#define histograms to fill
#we will write histograms to output file for:
onePhotonHist = rt.TH1F("1Photon_nProtonHist", "Energy of NC events with N proton(s) and 1 photon",60,0,6)
twoPhotonHist = rt.TH1F("2Photon_nProtonHist", "Energy of NC events with N proton(s) and 2 photons",60,0,6)
threePhotonHist = rt.TH1F("3Photon_nProtonHist", "Energy of NC events with N proton(s) and 3+ photons",60,0,6)

#set histogram axis titles and increase line width
def configureHist(h):
    h.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
    h.GetXaxis().SetTitle("energy (GeV)")
    h.SetLineWidth(2)
    return h

#scale the histograms based on total good POT
onePhotonHist = configureHist(onePhotonHist)
twoPhotonHist = configureHist(twoPhotonHist)
threePhotonHist = configureHist(threePhotonHist)

#set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10


recoFiducialCut = 0
recoCosmicCut = 0

#begin loop over events in ntuple file
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

    #reco fiducial cut
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
        eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
        eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
        recoFiducialCut += 1
        continue

    #reco cosmic cut
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        recoCosmicCut += 1
        continue

    #iterate through all tracks in event
    #look for non-secondary tracks identified as muons
    primaryMuonFound = False
    for i in range(eventTree.nTracks):
        if eventTree.trackIsSecondary[i] == 0:
            if abs(eventTree.trackPID[i]) == 13:
                primaryMuonFound = True
                break
    #cut events w/ primary muons
    if primaryMuonFound:    
        continue
    
    #iterate through showers in event
    #look for non-secondary showers identified as electrions
    primaryElectronFound = False
    for i in range(eventTree.nShowers):
        if eventTree.showerIsSecondary == 0:
            if abs(eventTree.showerPID[i]) == 11:
                primaryElectronFound == True
                break
    #cut events w/ primary electrons
    if primaryElectronFound:
        continue

    