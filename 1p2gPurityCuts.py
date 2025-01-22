import sys, argparse
import numpy as np
import ROOT as rt

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, histStackFill, trueCCCutLoose, recoCCCutLoose, recoInvariantMassCalculations

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
# parser.add_argument("-c", "--cosmicFile", type=str, required=True, help="input cosmic ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="2025-01-20_vtxDist.root", help="output root file name")
args = parser.parse_args()

# Grab ntuple file, eventTree, and potTree for reference
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

# Scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67+20
targetPOTstring = "6.67e+20" #for plot axis titles

# Get relevant info from cosmic ntuple file
# cosmic_ntuple_file = rt.TFile(args.cosmicFile)
# cosmicTree = cosmic_ntuple_file.Get("EventTree")
# cosmicPOTsum = 5.28e+19

# Calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

#Purity histograms: named true bc true called the reco events x
purityxMin = -105
purityxMax = 95
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
# trueCosmicsHist = rt.TH1F("trueCosmics", "Cosmic background reconstructed as signal", purityBins, purityxMin, purityxMax)
trueSignalHist = rt.TH1F("trueSignal", "1p2g successfully truth-matched", purityBins, purityxMin, purityxMax)

# Define min, max, fiducial
fiducialWidth = 10

recoEvents = 0
#------------------ Purity Loop ------------------#
#------------------ Purity Loop ------------------#
# We want to obtain our 1p2g events using reconstruction, then make signal events green and background events red, sorting by diff variables

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
    recoEvents += 1
    recoTotalHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
    
#------------ TruthMatching ------------#
# true cc check
    if trueCCCut(eventTree):
        trueCCHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
        continue
#true fiducial check
    if trueFiducialCut(eventTree, fiducialWidth):
        trueOutFiducialHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
        continue
# true pi+ check
    if truePiPlusCut(eventTree):
        truePiPlusHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
        continue
# true proton selection
    nTrueProtons, trueProtonTID, trueProtonIndex = trueProtonSelection(eventTree)
    if nTrueProtons == 0:
        trueNoProtonHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
        continue
    if nTrueProtons >= 2:
        truePluralProtonHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
        continue
# uncomment below for proton tid-matching
#    if trueProtonTID != trueProtonTID:
#        trueWrongProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue    
# true photon selection up top
    truePhotonTIDList, trueLeadingPhotonEnergy, truePhotonIndexList = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) == 0:
        trueNoPhotonHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
        continue
    if len(truePhotonTIDList) == 1:
        trueOnePhotonHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
        continue
    if len(truePhotonTIDList) >= 3:
        trueManyPhotonHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)
        continue

    trueSignalHist.Fill(eventTree.vtxDistToTrue, eventTree.xsecWeight)

#------------------ Cosmic Loop ------------------#
#------------------ Cosmic Loop ------------------#
# Perform reco cuts on cosmic data to determine number of cosmics to be tagged as signal
# cosmicEventsTally = 0
# for i in range(cosmicTree.GetEntries()):
#     cosmicTree.GetEntry(i)
# # Keep a tally of all cosmic events    
#     cosmicEventsTally += 1
# # start with whether vertex was found
#     if cosmicTree.foundVertex == 0:
#         continue
# # Cosmic CC Cut:
#     if recoCCCut(cosmicTree):
#         continue
# # Cosmic Fiducial Cut:
#     if recoFiducialCut(cosmicTree, fiducialWidth):
#         continue
# # Reco cosmic cut
#     if cosmicTree.vtxFracHitsOnCosmic >= 1.:
#         continue
# # Cosmic pi+ cut
#     if recoPiPlusCut(cosmicTree):
#         continue
# # Proton selection
#     nCosmicProtons = 0
#     CosmicProtonIndex = 0
#     for i in range(cosmicTree.nTracks):
#         if cosmicTree.trackPID[i] == 2212 and cosmicTree.trackRecoE[i] >= 60:
#             nCosmicProtons += 1
#             CosmicProtonIndex = i
#     if nCosmicProtons != 1:
#         continue
# # Photon Selection
#     cosmic = 0
#     cosmicPhotonIndexList = []
#     xMin, xMax = 0, 256
#     yMin, yMax = -116.5, 116.5
#     zMin, zMax = 0, 1036
#     for i in range(cosmicTree.nShowers):
#         if cosmicTree.showerPID[i] == 22:
#             if cosmicTree.showerStartPosX[i] <= (xMin + fiducialWidth) or cosmicTree.showerStartPosX[i] >= (xMax - fiducialWidth) or \
#                 cosmicTree.showerStartPosY[i] <= (yMin + fiducialWidth) or cosmicTree.showerStartPosY[i] >= (yMax - fiducialWidth) or \
#                     cosmicTree.showerStartPosZ[i] <= (zMin + fiducialWidth) or cosmicTree.showerStartPosZ[i] >= (zMax - fiducialWidth):
#                 cosmic += 1
#             else:
#                 cosmicPhotonIndexList.append(i)

#     recoPhotonEnergyList = [0]
#     recoLeadingPhotonEnergy = 0.
  
#     for i in range(len(cosmicPhotonIndexList)):
#         recoPhotonEnergyMeV = cosmicTree.showerRecoE[cosmicPhotonIndexList[i]] 
#         recoPhotonEnergyList.append(recoPhotonEnergyMeV)
#         recoLeadingPhotonEnergy = max(recoPhotonEnergyList)
 
#     if len(cosmicPhotonIndexList) != 2:
#         continue

#     trueCosmicsHist.Fill(cosmicTree.vtxMaxIntimePixelSum, 0.4)
 
#------------------ End of Loops ------------------#

recoTotalHistInt = recoTotalHist.Integral(1, purityBins)
trueSignalHistInt = trueSignalHist.Integral(1, purityBins)

print("Purity: " + str(trueSignalHistInt / recoTotalHistInt * 100) + "%")

purityHistList = [trueSignalHist, trueCCHist, trueOutFiducialHist, truePiPlusHist, trueNoProtonHist, truePluralProtonHist, trueNoPhotonHist, trueOnePhotonHist, trueManyPhotonHist]

purityCanvas, purityStack, purityLegend, purityInt = \
    histStackFill("NC 2g + 1p Event Reconstruction Purity", purityHistList, "All Reconstructed Events: (", "vtxDistToTrue", "Events per " + targetPOTstring + " POT", ntuplePOTsum)

outFile = rt.TFile(args.outfile, "RECREATE")
purityCanvas.Write("PurityHist")