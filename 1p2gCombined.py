import sys, argparse
import numpy as np
import ROOT as rt

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, histStackFill, trueCCCutLoose, recoCCCutLoose, recoInvariantMassCalculations

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="13Nov1p2gCombinedEFOutput_Invariant.root", help="output root file name")
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

# Define histograms to be filled through loop
trueTotalHist = rt.TH1F("trueTotalHist", "All True Signal Events", 60, 0, 1500)
recoNoVertexHist = rt.TH1F("NoVertexFound", "Vertex not reconstructed", 60, 0, 1500)
recoOutFiducialHist = rt.TH1F("outsideFiducial", "Vertex reconstructed outside fiducial volume", 60, 0, 1500)
recoCCHist = rt.TH1F("recoCC", "Event reconstructed as charged-current", 60, 0, 1500)
recoPiPlusHist = rt.TH1F("recoPiPlus", "Charged pion reconstructed", 60, 0, 1500)
recoNoProtonHist = rt.TH1F("recoNoProton", "No proton reconstructed", 60, 0, 1500)
recoPluralProtonHist = rt.TH1F("recoPluralProton", "2+ protons reconstructed", 60, 0, 1500)
recoWrongProtonHist = rt.TH1F("recoWrongProton", "Reconstructed proton failed TID-matching", 60, 0, 1500)
recoNoPhotonHist = rt.TH1F("recoNoPhoton", "No photon reconstructed", 60, 0, 1500)
recoOnePhotonHist = rt.TH1F("recoOnePhoton", "One photon reconstructed", 60, 0, 1500)
recoManyPhotonHist = rt.TH1F("recoManyPhoton", "3+ photons reconstructed", 60, 0, 1500)
recoSignalHist = rt.TH1F("recoSignal", "1p2g successfully reconstructed", 60, 0, 1500)


purityxMax = 1900
purityBins = 100
recoTotalHist = rt.TH1F("recoTotalHist", "All Reco Signal Events", purityBins, 0, purityxMax)
trueOutFiducialHist = rt.TH1F("trueOutsideFiducial", "True vertex outside fiducial volume", purityBins, 0, purityxMax)
trueCCHist = rt.TH1F("trueCC", "True event charged-current", purityBins, 0, purityxMax)
truePiPlusHist = rt.TH1F("truePiPlus", "True charged pion", purityBins, 0, purityxMax)
trueNoProtonHist = rt.TH1F("trueNoProton", "No true proton", purityBins, 0, purityxMax)
truePluralProtonHist = rt.TH1F("truePluralProton", "2+ true protons", purityBins, 0, purityxMax)
trueWrongProtonHist = rt.TH1F("trueWrongProton", "Reconstructed proton failed TID-matching", purityBins, 0, purityxMax)
trueNoPhotonHist = rt.TH1F("trueNoPhoton", "No true photon", purityBins, 0, purityxMax)
trueOnePhotonHist = rt.TH1F("trueOnePhoton", "One true photon", purityBins, 0, purityxMax)
trueManyPhotonHist = rt.TH1F("trueManyPhoton", "3+ true photons", purityBins, 0, purityxMax)
trueSignalHist = rt.TH1F("trueSignal", "1p2g successfully truth-matched", purityBins, 0, purityxMax)

# Define min, max, fiducial and data function
fiducialWidth = 10
fiducialDict = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}

#---------------- Efficiency Loop ----------------#
#---------------- Efficiency Loop ----------------#
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
# True charged-current cut
    if trueCCCut(eventTree):
        continue
# True loose cc cut
#    if trueCCCutLoose(eventTree):
#        continue
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
#------------ Reco Matching ------------#

# start with whether vertex was found
    if eventTree.foundVertex == 0:
        recoNoVertexHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
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
# reco loose cc check
#    if recoCCCutLoose(eventTree):
#        recoCCHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue
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
    nRecoProtons, recoProtonTID, recoProtonIndex = recoProtonSelection(eventTree)
    if nRecoProtons != 1:
        continue
# Reco photon selection
    recoPhotonTIDList, recoLeadingPhotonEnergy, recoPhotonIndexList = recoPhotonSelection(eventTree, fiducialWidth) 
    if len(recoPhotonTIDList) != 2:
        continue
    
# Reco invariant mass calc
    recoProtonInv, recoPiZeroInv, recoDeltaInv = recoInvariantMassCalculations(eventTree, recoProtonIndex, recoPhotonIndexList)

    recoTotalHist.Fill(recoDeltaInv, eventTree.xsecWeight)

#------------ TruthMatching ------------#
# true cc check
    if trueCCCut(eventTree):
        trueCCHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue

# true loose cc check
#    if trueCCCutLoose(eventTree):
#        trueCCHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue
#true fiducial check
    if trueFiducialCut(eventTree, fiducialWidth):
        trueOutFiducialHist.Fill(recoDeltaInv, eventTree.xsecWeight)
        continue
# true pi+ check
    if truePiPlusCut(eventTree):
        truePiPlusHist.Fill(recoDeltaInv, eventTree.xsecWeight)
        continue
# true proton selection
    nTrueProtons, trueProtonTID, trueProtonIndex = trueProtonSelection(eventTree)
    if nTrueProtons == 0:
        trueNoProtonHist.Fill(recoDeltaInv, eventTree.xsecWeight)
        continue
    if nTrueProtons >= 2:
        truePluralProtonHist.Fill(recoDeltaInv, eventTree.xsecWeight)
        continue
# uncomment below for proton tid-matching
#    if trueProtonTID != trueProtonTID:
#        trueWrongProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue    
# true photon selection up top
    truePhotonTIDList, trueLeadingPhotonEnergy, truePhotonIndexList = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) == 0:
        trueNoPhotonHist.Fill(recoDeltaInv, eventTree.xsecWeight)
        continue
    if len(truePhotonTIDList) == 1:
        trueOnePhotonHist.Fill(recoDeltaInv, eventTree.xsecWeight)
        continue
    if len(truePhotonTIDList) >= 3:
        trueManyPhotonHist.Fill(recoDeltaInv, eventTree.xsecWeight)
        continue

    trueSignalHist.Fill(recoDeltaInv, eventTree.xsecWeight)
 
#------------------ End of Loops ------------------#

trueTotalHistInt = trueTotalHist.Integral(1, 60)
recoTotalHistInt = recoTotalHist.Integral(1, purityBins)

recoSignalHistInt = recoSignalHist.Integral(1, 60)
trueSignalHistInt = trueSignalHist.Integral(1, purityBins)

print("Efficiency: " + str(recoSignalHistInt / trueTotalHistInt * 100) + "%")
print("Purity: " + str(trueSignalHistInt / recoTotalHistInt * 100) + "%")

efficiencyHistList = [recoSignalHist, recoNoVertexHist, recoCCHist, recoNoPhotonHist, recoOnePhotonHist, recoManyPhotonHist, recoOutFiducialHist, recoPiPlusHist, recoNoProtonHist, recoPluralProtonHist]
purityHistList = [trueSignalHist, trueCCHist, trueOutFiducialHist, truePiPlusHist, trueNoProtonHist, truePluralProtonHist, trueNoPhotonHist, trueOnePhotonHist, trueManyPhotonHist,]

efficiencyCanvas, efficiencyStack, efficiencyLegend, efficiencyInt = \
    histStackFill("Reconstruction of True NC 1p2g Events", efficiencyHistList, "Total Truth Signal: (", "True Leading Photon Energy (MeV)", "Events per 6.67e+20 POT", ntuplePOTsum)

purityCanvas, purityStack, purityLegend, purityInt = \
    histStackFill("Truth-Matching of Reconstructed NC 1p2g Events", purityHistList, "All Reconstructed Events: (", "Reconstructed 1p+2g Invariant Mass (MeV/c^2)", "Events per 6.67e+20 POT", ntuplePOTsum)

# Write Efficiency and Purity Canvases to Outfile
outFile = rt.TFile(args.outfile, "RECREATE")
efficiencyCanvas.Write("EfficiencyHist")
purityCanvas.Write("PurityHist")