import sys, argparse
import numpy as np
import ROOT as rt

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, histStackFill, trueCCCutLoose, recoCCCutLoose, recoInvariantMassCalculations

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-c", "--cosmicFile", type=str, required=True, help="input cosmic ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="2025-01-21_CombinedOutput4.4.root", help="output root file name")
args = parser.parse_args()

# Grab ntuple file, eventTree, and potTree for reference
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

# Scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 4.4e+19
targetPOTstring = "4.4e+19" #for plot axis titles

# Get relevant info from cosmic ntuple file
cosmic_ntuple_file = rt.TFile(args.cosmicFile)
cosmicTree = cosmic_ntuple_file.Get("EventTree")
cosmicPOTsum = 5.28e+19

# Calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

# Define histograms to be filled through loop

#Efficiency histograms: named reco bc reco called truth events x
trueTotalHist = rt.TH1F("trueTotalHist", "All True Signal Events", 60, 0, 1400)
recoNoVertexHist = rt.TH1F("NoVertexFound", "Vertex not reconstructed", 60, 0, 1400)
recoOutFiducialHist = rt.TH1F("outsideFiducial", "Vertex reconstructed outside fiducial volume", 60, 0, 1400)
recoCCHist = rt.TH1F("recoCC", "Event reconstructed as charged-current", 60, 0, 1400)
recoPiPlusHist = rt.TH1F("recoPiPlus", "Charged pion reconstructed", 60, 0, 1400)
recoNoProtonHist = rt.TH1F("recoNoProton", "No proton reconstructed", 60, 0, 1400)
recoPluralProtonHist = rt.TH1F("recoPluralProton", "2+ protons reconstructed", 60, 0, 1400)
recoWrongProtonHist = rt.TH1F("recoWrongProton", "Reconstructed proton failed TID-matching", 60, 0, 1400)
recoNoPhotonHist = rt.TH1F("recoNoPhoton", "No photon reconstructed", 60, 0, 1400)
recoOnePhotonHist = rt.TH1F("recoOnePhoton", "One photon reconstructed", 60, 0, 1400)
recoManyPhotonHist = rt.TH1F("recoManyPhoton", "3+ photons reconstructed", 60, 0, 1400)
recoSignalHist = rt.TH1F("recoSignal", "1p2g successfully reconstructed", 60, 0, 1400)


#Purity histograms: named true bc true called the reco events x
purityxMin = 0
purityxMax = 1200
purityBins = 60
recoTotalHist = rt.TH1F("recoTotalHist", "All Reco Signal Events", purityBins, purityxMin, purityxMax)
trueOutFiducialHist = rt.TH1F("trueOutsideFiducial", "True vertex outside fiducial volume", purityBins, purityxMin, purityxMax)
trueCCHist = rt.TH1F("trueCC", "True event charged-current", purityBins, purityxMin, purityxMax)
truePiPlusHist = rt.TH1F("truePiPlus", "True charged pion", purityBins, purityxMin, purityxMax)
trueNoProtonHist = rt.TH1F("trueNoProton", "No true proton", purityBins, purityxMin, purityxMax)
truePluralProtonHist = rt.TH1F("truePluralProton", "2+ true protons", purityBins, purityxMin, purityxMax)
trueWrongProtonHist = rt.TH1F("trueWrongProton", "Reconstructed proton failed TID-matching", purityBins, purityxMin, purityxMax)
trueNoPhotonHist = rt.TH1F("trueNoPhoton", "No true photon", purityBins, purityxMin, purityxMax)
trueOnePhotonHist = rt.TH1F("trueOnePhoton", "One true photon", purityBins, purityxMin, purityxMax)
trueManyPhotonHist = rt.TH1F("trueManyPhoton", "3+ true photons", purityBins, purityxMin, purityxMax)
trueCosmicsHist = rt.TH1F("trueCosmics", "Cosmic background reconstructed as signal", purityBins, purityxMin, purityxMax)
trueSignalHist = rt.TH1F("trueSignal", "1p2g successfully truth-matched", purityBins, purityxMin, purityxMax)

# Define min, max, fiducial and data function
fiducialWidth = 10

recoNoVertex = 0
trueTotalTally = 0
#---------------- Efficiency Loop ----------------#
#---------------- Efficiency Loop ----------------#
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
# True charged-current cut
    if trueCCCut(eventTree):
        continue
# True fiducial cut
    if trueFiducialCut(eventTree, fiducialWidth):
        continue
# True cosmic cut
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        continue
# True pi+ cut
    if truePiPlusCut(eventTree):
        continue
# True proton selection
    nTrueProtons, trueProtonTID, trueProtonIndex = trueProtonSelection(eventTree)
    if nTrueProtons != 1:
        continue
# True photon selection
    truePhotonTIDList, trueLeadingPhotonEnergy, photonIndexList = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) != 2:
        continue

    trueTotalHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
    trueTotalTally += 1
#------------ Reco Matching ------------#

# start with whether vertex was found
    if eventTree.foundVertex == 0:
        recoNoVertexHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        recoNoVertex += 1
        continue
# reco cc check
    if recoCCCut(eventTree):
        recoCCHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# reco photon selection up top
    recophotonTIDList, recoLeadingPhotonEnergy, recoPhotonIndexList = recoPhotonSelection(eventTree, fiducialWidth)
    if len(recophotonTIDList) == 0:
        recoNoPhotonHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if len(recophotonTIDList) == 1:
        recoOnePhotonHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if len(recophotonTIDList) >= 3:
        recoManyPhotonHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# reco fiducial check
    if recoFiducialCut(eventTree, fiducialWidth):
        recoOutFiducialHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# reco pi+ check
    if recoPiPlusCut(eventTree):
        recoPiPlusHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# reco proton selection and filling
    nRecoProtons, recoProtonTID, recoProtonIndex = recoProtonSelection(eventTree)
    if nRecoProtons == 0:
        recoNoProtonHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if nRecoProtons >= 2:
        recoPluralProtonHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# uncomment below for proton tid-matching
#    if recoProtonTID != trueProtonTID:
#        recoWrongProtonHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue
    recoSignalHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)

#------------------ Purity Loop ------------------#
#------------------ Purity Loop ------------------#

for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
# start with whether vertex was found
    if eventTree.foundVertex == 0:
        continue
# Reco CC cut
    if recoCCCut(eventTree):
        continue
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
    nRecoProtons, recoProtonTID, recoProtonIndex = recoProtonSelection(eventTree)
    if nRecoProtons != 1:
        continue
# Reco photon selection
    recoPhotonTIDList, recoLeadingPhotonEnergy, recoPhotonIndexList = recoPhotonSelection(eventTree, fiducialWidth) 
    if len(recoPhotonTIDList) != 2:
        continue
# Reco invariant mass calc
    # recoProtonInv, recoPiZeroInv, recoDeltaInv = recoInvariantMassCalculations(eventTree, recoProtonIndex, recoPhotonIndexList)

    recoTotalHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)

#------------ TruthMatching ------------#
# true cc check
    if trueCCCut(eventTree):
        trueCCHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
#true fiducial check
    if trueFiducialCut(eventTree, fiducialWidth):
        trueOutFiducialHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# true pi+ check
    if truePiPlusCut(eventTree):
        truePiPlusHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# true proton selection
    nTrueProtons, trueProtonTID, trueProtonIndex = trueProtonSelection(eventTree)
    if nTrueProtons == 0:
        trueNoProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if nTrueProtons >= 2:
        truePluralProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# uncomment below for proton tid-matching
#    if trueProtonTID != trueProtonTID:
#        trueWrongProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue    
# true photon selection up top
    truePhotonTIDList, trueLeadingPhotonEnergy, truePhotonIndexList = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) == 0:
        trueNoPhotonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if len(truePhotonTIDList) == 1:
        trueOnePhotonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if len(truePhotonTIDList) >= 3:
        trueManyPhotonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue

    trueSignalHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)

#------------------ Cosmic Loop ------------------#
#------------------ Cosmic Loop ------------------#
# Perform reco cuts on cosmic data to determine number of cosmics to be tagged as signal
cosmicEventsTally = 0
for i in range(cosmicTree.GetEntries()):
    cosmicTree.GetEntry(i)
# Keep a tally of all cosmic events    
    cosmicEventsTally += 1
# start with whether vertex was found
    if cosmicTree.foundVertex == 0:
        continue
# Cosmic CC Cut:
    recoPrimaryMuonTrackFound = False
    recoPrimaryMuonShowerFound = False
    recoPrimaryElectronTrackFound = False
    recoPrimaryElectronShowerFound = False
    cc = False
    unclassifiedTracks = 0
    for i in range(eventTree.nTracks):
        if eventTree.trackClassified[i] == 0:
            unclassifiedTracks += 1
        if eventTree.trackIsSecondary[i] == 0:
            if abs(eventTree.trackPID[i]) == 13:
                recoPrimaryMuonTrackFound = True
            elif abs(eventTree.trackPID[i]) == 11:
                recoPrimaryElectronTrackFound = True
    for i in range(eventTree.nShowers):
        if eventTree.showerIsSecondary[i] == 0:
            if abs(eventTree.showerPID[i]) == 11:
                recoPrimaryElectronShowerFound == True
            elif abs(eventTree.showerPID[i]) == 13:
                recoPrimaryMuonShowerFound == True
    if eventTree.nTracks != 0:
        if unclassifiedTracks/eventTree.nTracks >= 0.55:
            cc = True
    if recoPrimaryMuonTrackFound or recoPrimaryMuonShowerFound or recoPrimaryElectronTrackFound or recoPrimaryElectronShowerFound:   
        cc = True
    if eventTree.nTracks >= 4:
        cc = True
    if cc == True:
        continue
# Cosmic Fiducial Cut:
    if recoFiducialCut(cosmicTree, fiducialWidth):
        continue
# Reco cosmic cut
    if cosmicTree.vtxFracHitsOnCosmic >= 1.:
        continue
# Cosmic pi+ cut
    if recoPiPlusCut(cosmicTree):
        continue
# Proton selection
    nCosmicProtons = 0
    CosmicProtonIndex = 0
    for i in range(cosmicTree.nTracks):
        if cosmicTree.trackPID[i] == 2212 and cosmicTree.trackRecoE[i] >= 60:
            nCosmicProtons += 1
            CosmicProtonIndex = i
    if nCosmicProtons != 1:
        continue
# Photon Selection
    cosmic = 0
    cosmicPhotonIndexList = []
    xMin, xMax = 0, 256
    yMin, yMax = -116.5, 116.5
    zMin, zMax = 0, 1036
    for i in range(cosmicTree.nShowers):
        if cosmicTree.showerPID[i] == 22:
            if cosmicTree.showerStartPosX[i] <= (xMin + fiducialWidth) or cosmicTree.showerStartPosX[i] >= (xMax - fiducialWidth) or \
                cosmicTree.showerStartPosY[i] <= (yMin + fiducialWidth) or cosmicTree.showerStartPosY[i] >= (yMax - fiducialWidth) or \
                    cosmicTree.showerStartPosZ[i] <= (zMin + fiducialWidth) or cosmicTree.showerStartPosZ[i] >= (zMax - fiducialWidth):
                cosmic += 1
            else:
                cosmicPhotonIndexList.append(i)

    recoPhotonEnergyList = [0]
    recoLeadingPhotonEnergy = 0.
  
    for i in range(len(cosmicPhotonIndexList)):
        recoPhotonEnergyMeV = cosmicTree.showerRecoE[cosmicPhotonIndexList[i]] 
        recoPhotonEnergyList.append(recoPhotonEnergyMeV)
        recoLeadingPhotonEnergy = max(recoPhotonEnergyList)
 
    if len(cosmicPhotonIndexList) != 2:
        continue

    # recoProtonInv, recoPiZeroInv, recoDeltaInv = recoInvariantMassCalculations(cosmicTree, CosmicProtonIndex, cosmicPhotonIndexList)

    trueCosmicsHist.Fill(recoLeadingPhotonEnergy, 0.4)
 
#------------------ End of Loops ------------------#
trueCosmicsHist.Scale(ntuplePOTsum/cosmicPOTsum)

print(str(trueTotalTally) + "truth events")
print(str(recoNoVertex) + "truth where reco can't find vertex")

trueTotalHistInt = trueTotalHist.Integral(1, 60)
recoTotalHistInt = recoTotalHist.Integral(1, purityBins)

recoSignalHistInt = recoSignalHist.Integral(1, 60)
trueSignalHistInt = trueSignalHist.Integral(1, purityBins)

print("Efficiency: " + str(recoSignalHistInt / trueTotalHistInt * 100) + "%")
print("Purity: " + str(trueSignalHistInt / recoTotalHistInt * 100) + "%")

efficiencyHistList = [recoSignalHist, recoNoVertexHist, recoCCHist, recoNoPhotonHist, recoOnePhotonHist, recoManyPhotonHist, recoOutFiducialHist, recoPiPlusHist, recoNoProtonHist, recoPluralProtonHist]
purityHistList = [trueSignalHist, trueCCHist, trueOutFiducialHist, truePiPlusHist, trueNoProtonHist, truePluralProtonHist, trueNoPhotonHist, trueOnePhotonHist, trueManyPhotonHist, trueCosmicsHist]

efficiencyCanvas, efficiencyStack, efficiencyLegend, efficiencyInt = \
    histStackFill(" NC 2g + 1p Event Reconstruction Efficiency", efficiencyHistList, "All Truth Events: (", "True Leading Photon Energy (MeV)", "Events per 4.4e+19 POT", ntuplePOTsum)

purityCanvas, purityStack, purityLegend, purityInt = \
    histStackFill("NC 2g + 1p Event Reconstruction Purity", purityHistList, "All Reconstructed Events: (", "Reconstructed Leading Photon Energy (MeV)", "Events per 4.4e+19 POT", ntuplePOTsum)

# Write Efficiency and Purity Canvases to Outfile
outFile = rt.TFile(args.outfile, "RECREATE")
efficiencyCanvas.Write("EfficiencyHist")
purityCanvas.Write("PurityHist")