import sys, argparse
import numpy as np
import ROOT as rt

from cuts import kineticEnergyCalculator, efficiencyPlot

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gRecoTIDMatchEnergy.root", help="output root file name")
args = parser.parse_args()

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
trueSignalHist = rt.TH1F("trueSignalHist", "True NC 1 proton 2 gamma Events",60,0,1500)

noVtxFoundHist = rt.TH1F("No Vertex Found", "Reco couldn't find a vertex",60,0,1500)
outFiducialHist = rt.TH1F("Futside Fiducial", "Reco placed vertex outside fiducial volume",60,0,1500)
chargedCurrentHist = rt.TH1F("Charged Current", "Reco identified as charged-current",60,0,1500)
piPlusHist = rt.TH1F("pi+", "Reco found a pi+",60,0,1500)
noProtonHist = rt.TH1F("No Protons", "Reco found no protons",60,0,1500)
pluralProtonHist = rt.TH1F("Multiple Protons", "Reco found multiple protons",60,0,1500)
wrongProtonHist = rt.TH1F("1 Wrong Proton", "Reco found 1 proton but failed TID matching",60,0,1500)
noPhotonHist = rt.TH1F("No Photons", "Reco found no photons",60,0,1500)
onePhotonHist = rt.TH1F("One Photon", "Reco found 1 photon",60,0,1500)
manyPhotonHist = rt.TH1F("Many Photons", "Reco found 3+ photons",60,0,1500)
recoSignalHist = rt.TH1F("Reco = True", "Correctly identified by Reco",60,0,1500)
truePhotonHist = rt.TH1F("True EDep Dist of all Photons", "EDep-Vtx Distance of all True Photons",60,0,1500)
photonReconstructedHist = rt.TH1F("reco reconstructed the photon", "Photon was TID-matched in Reco",60,0,1500)
photonNotReconstructedHist = rt.TH1F("reco didn't reconstruct the photon", "Photon not TID-matched in Reco",60,0,1500)

#set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

totaltrueevents = 0

totalrecophotons = 0
truthMatchedPhotons = 0
nonTruthMatchedPhotons = 0
exhaustive = 0
alltimeEDepList = []
#begin loop over events in ntuple file:
# if we start with truth then reco-match, we can ID reco's false negatives
# if we start with reco then truth-match, we can ID reco's false positives

#start by using truth to cut all background, then do the same checks in reco
#after reco-matching, sort events by how reco identifies the true signal, 
#then create a stacked histogram w/ x-axis of photon energy
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

#true charged-current cut
    if eventTree.trueNuCCNC != 1:
        continue

#true fiducial volume cut
    if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or \
        eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or \
        eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
        continue

#true cosmic background cut
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        continue

#true proton & pi+ check: include all photons above 60MeV, exclude all pi+ above 30MeV
    nProtons = 0
    piPlusPresent = False
    for i in range(len(eventTree.truePrimPartPDG)):
        #proton checker
        if abs(eventTree.truePrimPartPDG[i]) == 2212:
            kE = kineticEnergyCalculator(eventTree, i)
            if kE >= 0.06:
                nProtons += 1
        #pi+ checker        
        elif abs(eventTree.truePrimPartPDG[i]) == 211:
            kE = kineticEnergyCalculator(eventTree, i)
            if kE >= 0.03:
                piPlusPresent = True
                break
    if nProtons != 1 or piPlusPresent:
        continue

#true photon check: look through primary particles to count the photon daughters
    photonInSecondary = False
    primList = []
    photonIndexList = []
    truePhotonTIDList = []
    truePhotonEnergyList = []
    for i in range(len(eventTree.trueSimPartTID)):
        if eventTree.trueSimPartTID[i] == eventTree.trueSimPartMID[i]:
            primList.append(eventTree.trueSimPartTID[i])
    for i in range(len(eventTree.trueSimPartPDG)):
        if eventTree.trueSimPartPDG[i] == 22:
            if eventTree.trueSimPartMID[i] in primList:
                photonIndexList.append(i)
                truePhotonTIDList.append(eventTree.trueSimPartTID[i])
                truePhotonEnergyList.append(eventTree.trueSimPartE[i])
                photonInSecondary = True
            elif abs(eventTree.trueSimPartX[i] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[i] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[i] -eventTree.trueVtxZ) <= 0.15:
                photonIndexList.append(i)
                truePhotonTIDList.append(eventTree.trueSimPartTID[i])
                truePhotonEnergyList.append(eventTree.trueSimPartE[i])
                photonInSecondary = True

    if photonInSecondary == False:
        continue
    if len(photonIndexList) != 2:
        continue

    photonEDepOutsideFiducial = False
    for i in range(len(photonIndexList)):
        if eventTree.trueSimPartEDepX[photonIndexList[i]] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[photonIndexList[i]] >= (xMax - fiducialWidth)\
            or eventTree.trueSimPartEDepY[photonIndexList[i]] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[photonIndexList[i]] >= (yMax - fiducialWidth)\
            or eventTree.trueSimPartEDepZ[photonIndexList[i]] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[photonIndexList[i]] >= (zMax - fiducialWidth):
            photonEDepOutsideFiducial = True
    if photonEDepOutsideFiducial == True:
        continue
    
#find leading photon energy
    photonEnergyList = []
    for i in range(len(photonIndexList)):
        photonEnergyMeV = eventTree.trueSimPartE[photonIndexList[i]] 
        photonEnergyList.append(photonEnergyMeV)
    leadingPhotonEnergy = max(photonEnergyList)

    trueSignalHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)

#find and store the true trackID of the proton
    for i in range(len(eventTree.trueSimPartTID)):
        if eventTree.trueSimPartProcess[i] == 0:
            if abs(eventTree.trueSimPartPDG[i]) == 2212:
                momentumVector = np.square(eventTree.trueSimPartPx[i]) + \
                    np.square(eventTree.trueSimPartPy[i]) + np.square(eventTree.trueSimPartPz[i])
                kineticMeV = eventTree.trueSimPartE[i] - np.sqrt((np.square(eventTree.trueSimPartE[i])) - momentumVector)
                if kineticMeV >= 60:
                    trueProtonTID = eventTree.trueSimPartTID[i]

    #-------- end of truth selection --------#
    #-------- start of reco matching --------#

#reco vertex check: does reco even find an event
    if eventTree.foundVertex != 1:
        noVtxFoundHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue

#reco fiducial histogram fill
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
        eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
        eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
        outFiducialHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue

# don't need reco cosmic cut - it's the same as the true

#reco charged-current cut:        
#iterate through all tracks/showers in event,
#look for non-secondary tracks identified as muons or electrons
#and for non-secondary showers identified as electrons
    recoPrimaryMuonFound = False
    recoPrimaryElectronTrackFound = False
    recoPrimaryElectronFound = False
    for i in range(eventTree.nTracks):
        if eventTree.trackIsSecondary[i] == 0:
            if abs(eventTree.trackPID[i]) == 13:
                recoPrimaryMuonFound = True
            if abs(eventTree.trackPID[i]) == 11:
                recoPrimaryElectronTrackFound = True
    for i in range(eventTree.nShowers):
        if eventTree.showerIsSecondary[i] == 0:
            if abs(eventTree.showerPID[i]) == 11:
                recoPrimaryElectronFound == True
#fill hist w/ events that have primary electrons/muons (charged-current)
    if recoPrimaryMuonFound or recoPrimaryElectronTrackFound or recoPrimaryElectronFound:    
        chargedCurrentHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue

#reco pi+: fill and continue all events w/ pi+ of KE >= 30MeV
    recoPiPlusPresent = False
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 211:
                if eventTree.trackRecoE[i] >= 30:
                    recoPiPlusPresent = True
                    break
    if recoPiPlusPresent:
        piPlusHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue

#reco protons: find all protons of KE >= 60MeV and truthmatch TID
#then fill hists based on proton number
    recoNumProtons = 0
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 2212:
            if eventTree.trackRecoE[i] >= 60.0:
                recoNumProtons += 1
                recoProtonTID = eventTree.trackTrueTID[i]
    if recoNumProtons == 0:
        noProtonHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if recoNumProtons >= 2:
        pluralProtonHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if recoNumProtons == 1:
        if recoProtonTID != trueProtonTID:
            wrongProtonHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
            continue
    
#reco photons: find all photons, then append truth-matched TID to a list
    recoPhotonTIDList = []
    for i in range(eventTree.nShowers):
        if eventTree.showerPID[i] == 22:
            recoPhotonTIDList.append(eventTree.showerTrueTID[i])


#starting dictionary builder: all i do right now is fill the keys
    trueToRecoMap = {}
    for i in range (len(truePhotonTIDList)):
        trueToRecoMap[truePhotonTIDList[i]] = []

#calculate each photon's edep distance from vtx
    trueEDepDistanceList = []
    for i in range(len(truePhotonTIDList)):
        for j in range(len(eventTree.trueSimPartTID)):
            if eventTree.trueSimPartTID[j] == truePhotonTIDList[i]:
                deltaX = eventTree.trueSimPartEDepX[j] - eventTree.trueSimPartX[j]
                deltaY = eventTree.trueSimPartEDepY[j] - eventTree.trueSimPartY[j]
                deltaZ = eventTree.trueSimPartEDepZ[j] - eventTree.trueSimPartZ[j]
                trueEDepDistanceList.append(np.sqrt(np.square(deltaX) + np.square(deltaY) + np.square(deltaZ)))
                
#truth match checker: currently, i'm just checking that each true tid appears in the reco tid
    for i in range(len(truePhotonTIDList)):
        truePhotonHist.Fill(truePhotonEnergyList[i], eventTree.xsecWeight)
        if truePhotonTIDList[i] in recoPhotonTIDList:
            photonReconstructedHist.Fill(truePhotonEnergyList[i], eventTree.xsecWeight)
        else:
            photonNotReconstructedHist.Fill(truePhotonEnergyList[i], eventTree.xsecWeight)

#----- end of event loop ---------------------------------------------#

truthMatchEnergyHist = rt.TH1F("Photon Reconstruction Efficiency", "Efficiency of photon reconstruction as a function of photon energy",60,0,1500)

truthMatchedEfficiencyCanvas, truthMatchedEfficiencyStack = efficiencyPlot(truePhotonHist, \
    photonReconstructedHist, truthMatchEnergyHist, "Efficiency of photon reconstruction as a function of photon energy", \
        "Photon Energy (MeV)")

outFile = rt.TFile(args.outfile, "RECREATE")
truthMatchedEfficiencyCanvas.Write()
