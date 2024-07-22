import sys, argparse
import numpy as np
import ROOT as rt

from cuts import histStackFill, kineticEnergyCalculator, sStackFillNS

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gCC.root", help="output root file name")
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
muonEnergyHist = rt.TH1F("muonEnergyHist", "Muon Energy of True CC Reco NC Events",60,0,2200)
trueChargedCurrentHist = rt.TH1F("trueCCHist", "True event is charged-current",60,0,1500)
trueOutFiducialHist = rt.TH1F("trueOutFiducialHist", "True vertex outside fiducial volume",60,0,1500)
truePiPlusHist = rt.TH1F("truePiPlusHist", "True event contains a pi+",60,0,1500)
trueNoProtonHist = rt.TH1F("trueNoProtonHist", "True event contains no protons",60,0,1500)
trueManyProtonHist = rt.TH1F("trueManyProtonHist", "True event contains multiple protons",60,0,1500)
trueWrongProtonHist = rt.TH1F("trueWrongProtonHist", "True event proton was not TID-matched by reco",60,0,1500)
trueNoPhotonHist = rt.TH1F("trueNoPhotonHist", "True event contains no photons",60,0,1500)
trueOnePhotonHist = rt.TH1F("trueOnePhotonHist", "True event contains one photon",60,0,1500)
trueThreePhotonHist = rt.TH1F("trueThreePhotonHist", "True event contains 3+ photons",60,0,1500)
trueNoTIDMatchedPhotonHist = rt.TH1F("trueNoTIDMatchedPhotonHist", "Neither reco photon was TID-matched to true",60,0,1500)
trueOneTIDMatchedPhotonHist = rt.TH1F("trueOneTIDMatchedPhotonHist", "One reco photon was TID-matched to true",60,0,1500)
trueSignalHist = rt.TH1F("trueSignalHist", "Event was fully truth-matched",60,0,1500)

xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

chargeCurrentEvents = 0
muonsObserved = 0
electronsObserved = 0
muonKEList = []
muonEList = []

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

# find leading photon energy in reco:
    recoPhotonEnergyList = []
    for i in range(len(recoPhotonIndexList)):
        recoPhotonEnergyMeV = eventTree.showerRecoE[recoPhotonIndexList[i]] 
        recoPhotonEnergyList.append(recoPhotonEnergyMeV)
    leadingPhotonEnergy = max(recoPhotonEnergyList)



    #-------------- end reco selection --------------#
    #------------- begin truth matching -------------#

# true charged-current selection
    
    if eventTree.trueNuCCNC != 1:
        trueChargedCurrentHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        chargeCurrentEvents += 1
        for i in range(eventTree.nTruePrimParts):
            if abs(eventTree.truePrimPartPDG[i]) == 13:
                muonsObserved += 1
                muonKE_MeV = (kineticEnergyCalculator(eventTree, i))*1000
                muonEnergyHist.Fill(muonKE_MeV, eventTree.xsecWeight)

            if abs(eventTree.truePrimPartPDG[i]) == 11:
                electronsObserved += 1


#----- end of event loop ---------------------------------------------#

#recoSignalCanvas, recoSignalStack, recoSignalLegend, recoSignalInt = \
#    sStackFillS("Reco-Identified NC 1 proton 2 gamma Events", recoSignalHist, rt.kBlue, "recoSignalCanvas")

#histList = [trueSignalHist, trueChargedCurrentHist, trueOutFiducialHist, truePiPlusHist, trueNoProtonHist, \
#            trueManyProtonHist, trueWrongProtonHist, trueNoPhotonHist, trueOnePhotonHist, trueThreePhotonHist, \
#                trueNoTIDMatchedPhotonHist, trueOneTIDMatchedPhotonHist]

#purityCanvas, purityStack, purityLegend, purityInt = histStackFill("True IDs of Reco NC 1 Proton 2 Gamma Events", histList, \
#                                                                   "Truth Identification: (", "Leading Reco Photon Energy (MeV)", \
#                                                                    "Events per 6.67e+20 POT")

muonEnergyCanvas, muonEnergyStack, muonEnergyLegend, muonEnergyInt = sStackFillNS("Muon Energy of Reco NC, True CC Events", muonEnergyHist, rt.kBlue, "muon energy")

outFile = rt.TFile(args.outfile, "RECREATE")
muonEnergyCanvas.Write()

print(str(chargeCurrentEvents) + " total cc mis-ID events")
print(str(muonsObserved) + " total primary muons found")
print(str(electronsObserved) + " total primary electrons found")
