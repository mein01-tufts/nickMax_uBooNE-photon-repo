import sys, argparse
import numpy as np
import ROOT as rt

from cuts import histStackFill, sStackFillS, sStackFillNS

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gMissingPhotonsTest.root", help="output root file name")
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
trueSignalHist = rt.TH1F("trueSignalHist", "True NC 1 proton 2 gamma Events",60,0,600)

noPhotonHist = rt.TH1F("No Photons", "Reco saw neither true photon",60,0,600)
onePhotonHist = rt.TH1F("One Photon", "Reco saw one true photon",60,0,600)
manyPhotonHist = rt.TH1F("Many Photons", "Reco saw extra photons",60,0,600)
recoSignalHist = rt.TH1F("Reco = True", "Reco saw both true photons",60,0,600)

#set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

totalEvents = 0
recoEvents = 0
recoNoPhotons = 0
reco1Photon = 0
reco2Photon = 0
reco3Photon = 0

recoTotalPhotons = 0
recoTruthMatchedPhotons = 0
recoTruthMatchedPhotons2 = 0
recoTruthMatchedPhotons3 = 0
trueTotalPhotons = 0
totalPhotonsRecoSees = 0

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

#true proton check: does its kinetic energy exceed 60mev? - count it if so
    nProtons = 0
    for i in range(len(eventTree.truePrimPartPDG)):
        if eventTree.truePrimPartPDG[i] == 2212:
            pVector = np.square(eventTree.truePrimPartPx[i]) + np.square(eventTree.truePrimPartPy[i]) + np.square(eventTree.truePrimPartPz[i])
            kE = eventTree.truePrimPartE[i] - np.sqrt((np.square(eventTree.truePrimPartE[i])) - pVector)
            if kE >= 0.06:
                nProtons += 1
    if nProtons != 1:
        continue

#true pi+ check: - we want to exclude any above 30MeV
    piPlusPresent = False
    for i in range(len(eventTree.truePrimPartPDG)):
        if abs(eventTree.truePrimPartPDG[i] == 211):
            pVector = np.square(eventTree.truePrimPartPx[i]) + np.square(eventTree.truePrimPartPy[i]) + np.square(eventTree.truePrimPartPz[i])
            kE = eventTree.truePrimPartE[i] - np.sqrt((np.square(eventTree.truePrimPartE[i])) - pVector)
            if kE >= 0.03:
                piPlusPresent = True
                break
    if piPlusPresent == True:
        continue

#true photon check: look through primary particles to count the photon daughters
    photonInSecondary = False
    primList = []
    photonIndexList = []
    photonTIDlist = []
    for i in range(len(eventTree.trueSimPartTID)):
        if eventTree.trueSimPartTID[i] == eventTree.trueSimPartMID[i]:
            primList.append(eventTree.trueSimPartTID[i])
    for i in range(len(eventTree.trueSimPartPDG)):
        if eventTree.trueSimPartPDG[i] == 22:
            if eventTree.trueSimPartMID[i] in primList:
                photonIndexList.append(i)
                photonTIDlist.append(eventTree.trueSimPartTID[i])
                photonInSecondary = True
            elif abs(eventTree.trueSimPartX[i] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[i] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[i] -eventTree.trueVtxZ) <= 0.15:
                photonIndexList.append(i)
                photonTIDlist.append(eventTree.trueSimPartTID[i])
                photonInSecondary =True
    if photonInSecondary == False:
        continue
    if len(photonIndexList) != 2:
        continue

    trueTotalPhotons += 2
#find leading photon energy
    photonEnergyList = []
    for i in range(len(photonIndexList)):
        photonEnergyMeV = eventTree.trueSimPartE[photonIndexList[i]] 
        photonEnergyList.append(photonEnergyMeV)

    avgEnergy = np.mean(photonEnergyList)

    trueSignalHist.Fill(avgEnergy, eventTree.xsecWeight)

    #-------- end of truth selection --------#
    #-------- start of reco matching --------#

#reco vertex cut:
    if eventTree.foundVertex != 1:
        continue

#reco fiducial cut
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
        eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
        eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
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
#cut events w/ primary muons and/or electrons
    if recoPrimaryMuonFound or recoPrimaryElectronTrackFound or recoPrimaryElectronFound:    
        continue

#reco protons: find all protons of KE >= 60MeV, cut if != 1
    recoNumProtons = 0
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 2212:
            if eventTree.trackRecoE[i] >= 60:
                recoNumProtons += 1
    if recoNumProtons != 1:
        continue
    
#reco pi+: cut all events w/ pi+ of KE >= 30MeV
    recoPiPlusPresent = False
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 211:
                if eventTree.trackRecoE[i] >= 30:
                    recoPiPlusPresent = True
                    break
    if recoPiPlusPresent:
        continue

#fill to truth all events where reco 
#fill each photon energy into true


    totalEvents += 1
#reco photons: find and tally all photons
    recoNumPhotons = 0
    recoEvents += 1
    totalPhotonsRecoSees += 2
    for i in range(eventTree.nShowers):
        if eventTree.showerPID[i] == 22:
            recoTotalPhotons += 1
            recoNumPhotons += 1
            if eventTree.showerTruePID[i] == 22:
                recoTruthMatchedPhotons += 1
                if eventTree.showerTruePID[i] == 22 and eventTree.showerTrueTID[i] in photonIndexList:
                    recoTruthMatchedPhotons3 += 1

        if eventTree.showerTrueTID[i] in photonIndexList:
            recoTruthMatchedPhotons2 += 1

    if recoNumPhotons == 0:
        recoNoPhotons += 1
        noPhotonHist.Fill(avgEnergy, eventTree.xsecWeight)
    if recoNumPhotons == 1:
        reco1Photon += 1
        onePhotonHist.Fill(avgEnergy, eventTree.xsecWeight)
    if recoNumPhotons == 2:
        reco2Photon += 1
        recoSignalHist.Fill(avgEnergy, eventTree.xsecWeight)
    if recoNumPhotons >= 3:
        reco3Photon += 1
        manyPhotonHist.Fill(avgEnergy, eventTree.xsecWeight)

    true2recoMap = {}
    for i in photonIndexList:
        true2recoMap[i] = []
    


#----- end of event loop ---------------------------------------------#
histList = [recoSignalHist, noPhotonHist, onePhotonHist, manyPhotonHist,]
signalCanvas, signalStack, signalLegend, signalInt = histStackFill("Reco ID of Photons in True NC 1 Proton, 2 Gamma Events", histList, "Reco Identification: (", "(True) Average Photon Energy (MeV)", "Events per 6.67e+20 POT")

trueCanvas, trueStack, trueLegend, trueInt = sStackFillS("True NC 1 Proton, 2 Gamma Events", trueSignalHist, rt.kBlue, "trueCombined")

outFile = rt.TFile(args.outfile, "RECREATE")
signalCanvas.Write()
trueCanvas.Write()

#print("Reco's truth-matching photon efficiency: "+ str(recoTruthMatchedPhotons/recoTotalPhotons * 100) + " %")
#print("Reco's efficiency at finding and truth-matching a true photon: " + str(recoTruthMatchedPhotons/totalPhotonsRecoSees*100) + "%")
#print("Reco's efficiency at finding the same number of photons it should see: " + str(recoTotalPhotons/totalPhotonsRecoSees*100) + "%")
#print("reco IDs a shower whose true ID is photon: " + str(recoTruthMatchedPhotons2))
#print("reco IDs a photon shower whose true ID is also photon: " + str(recoTruthMatchedPhotons))
#print("reco IDs a photon shower whose true ID is also photon, and whose TID is in the true photon list: " + str(recoTruthMatchedPhotons3))

#print("\nreco IDs a shower as non-photon, but that shower's true ID is photon: " + str(recoTruthMatchedPhotons2-recoTruthMatchedPhotons))
#print("reco IDs a photon shower whose true ID is also photon but is not a secondary particle: " + str(recoTruthMatchedPhotons-recoTruthMatchedPhotons3) + "\n")