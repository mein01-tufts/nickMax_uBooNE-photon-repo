import sys, argparse
import numpy as np
import ROOT as rt

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelectionPiZero, recoPhotonSelection, histStackFill, trueCCCutLoose, recoCCCutLoose

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="NCPiZeroSignal.root", help="output root file name")
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


purityxMax = 1500
recoTotalHist = rt.TH1F("recoTotalHist", "All Reco Signal Events", 60, 0, purityxMax)
trueOutFiducialHist = rt.TH1F("trueOutsideFiducial", "True vertex outside fiducial volume", 60, 0, purityxMax)
trueCCHist = rt.TH1F("trueCC", "True event charged-current", 60, 0, purityxMax)
truePiPlusHist = rt.TH1F("truePiPlus", "True charged pion", 60, 0, purityxMax)
trueNoProtonHist = rt.TH1F("trueNoProton", "No true proton", 60, 0, purityxMax)
truePluralProtonHist = rt.TH1F("truePluralProton", "2+ true protons", 60, 0, purityxMax)
trueWrongProtonHist = rt.TH1F("trueWrongProton", "Reconstructed proton failed TID-matching", 60, 0, purityxMax)
trueNoPhotonHist = rt.TH1F("trueNoPhoton", "No true photon", 60, 0, purityxMax)
trueOnePhotonHist = rt.TH1F("trueOnePhoton", "One true photon", 60, 0, purityxMax)
trueManyPhotonHist = rt.TH1F("trueManyPhoton", "3+ true photons", 60, 0, purityxMax)
trueSignalHist = rt.TH1F("trueSignal", "1p2g successfully reconstructed", 60, 0, purityxMax)

# Define min, max, fiducial and data function
fiducialWidth = 10
fiducialDict = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}

trueSignal = 0
noProton, oneProton, manyProton = 0, 0, 0
onePhoton, twoPhoton, threePhoton = 0, 0, 0

final0p1g, final1p1g, final0p2g, final1p2g = 0, 0, 0, 0
#---------------- Efficiency Loop ----------------#
#---------------- Efficiency Loop ----------------#
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)   
    if trueCCCut(eventTree):
        continue
# true loose cc check
#    if trueCCCutLoose(eventTree):
#        trueCCHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue
#true fiducial check
    if trueFiducialCut(eventTree, fiducialWidth):
        continue
# true photon selection up top
    truePhotonTIDList, trueLeadingPhotonEnergy = truePhotonSelectionPiZero(eventTree, fiducialWidth)
    if len(truePhotonTIDList) == 0:
        continue

    nTrueProtons, trueProtonTID = trueProtonSelection(eventTree)
    trueSignal += 1
    if len(truePhotonTIDList) == 1 and nTrueProtons == 0:
        final0p1g += 1
    if len(truePhotonTIDList) == 1 and nTrueProtons == 1:
        final1p1g += 1
    if len(truePhotonTIDList) == 2 and nTrueProtons == 0:
        final0p2g += 1
    if len(truePhotonTIDList) == 2 and nTrueProtons == 1:
        final1p2g += 1    
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
    nRecoProtons, recoProtonTID = recoProtonSelection(eventTree)
    if nRecoProtons != 1:
        continue
# Reco photon selection
    recoPhotonTIDList, recoLeadingPhotonEnergy = recoPhotonSelection(eventTree, fiducialWidth) 
    if len(recoPhotonTIDList) != 2:
        continue
    
#    deltaX = abs(eventTree.trueVtxX - eventTree.vtxX)
#    deltaY = abs(eventTree.trueVtxY - eventTree.vtxY)
#    deltaZ = abs(eventTree.trueVtxZ - eventTree.vtxZ)

#    vtxDistance = np.sqrt(np.square(deltaX) + np.square(deltaY) + np.square(deltaZ))

    recoTotalHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)

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
        trueOutFiducialHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
# true photon selection up top
    truePhotonTIDList, trueLeadingPhotonEnergy = truePhotonSelectionPiZero(eventTree, fiducialWidth)
    if len(truePhotonTIDList) == 0:
        trueNoPhotonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
    trueSignalHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
 
#------------------ End of Loops ------------------#

purityHistList = [trueSignalHist, trueCCHist, trueOutFiducialHist, trueNoPhotonHist]

purityCanvas, purityStack, purityLegend, purityInt = \
    histStackFill("Truth-Matching of Reconstructed NC 1p2g Events", purityHistList, "Total Reconstruction Signal: (", "Reconstructed Leading Photon Energy (MeV)", "Events per 6.67e+20 POT", ntuplePOTsum)

# Write Efficiency and Purity Canvases to Outfile
outFile = rt.TFile(args.outfile, "RECREATE")
purityCanvas.Write("PurityHist")

print(str(trueSignal * targetPOT/ntuplePOTsum) + " pi0 per 6.67e+20")
print(str(final0p1g/trueSignal*100) + " percent 0p1g")
print(str(final0p2g/trueSignal*100) + " percent 0p2g")
print(str(final1p1g/trueSignal*100) + " percent 1p1g")
print(str(final1p2g/trueSignal*100) + " percent 1p2g")
