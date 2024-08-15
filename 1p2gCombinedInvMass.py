import sys, argparse
import numpy as np
import ROOT as rt

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, histStackFill

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gCombinedOutputInvariant3.root", help="output root file name")
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
efficiencyXmax = 2500
trueTotalHist = rt.TH1F("trueTotalHist", "All True Signal Events", 60, 0, efficiencyXmax)
recoNoVertexHist = rt.TH1F("NoVertexFound", "Vertex not reconstructed", 60, 0, efficiencyXmax)
recoOutFiducialHist = rt.TH1F("outsideFiducial", "Vertex reconstructed outside fiducial volume", 60, 0, efficiencyXmax)
recoCCHist = rt.TH1F("recoCC", "Event reconstructed as charged-current", 60, 0, efficiencyXmax)
recoPiPlusHist = rt.TH1F("recoPiPlus", "Charged pion reconstructed", 60, 0, efficiencyXmax)
recoNoProtonHist = rt.TH1F("recoNoProton", "No proton reconstructed", 60, 0, efficiencyXmax)
recoPluralProtonHist = rt.TH1F("recoPluralProton", "2+ protons reconstructed", 60, 0, efficiencyXmax)
recoWrongProtonHist = rt.TH1F("recoWrongProton", "Reconstructed proton failed TID-matching", 60, 0, efficiencyXmax)
recoNoPhotonHist = rt.TH1F("recoNoPhoton", "No photon reconstructed", 60, 0, efficiencyXmax)
recoOnePhotonHist = rt.TH1F("recoOnePhoton", "One photon reconstructed", 60, 0, efficiencyXmax)
recoManyPhotonHist = rt.TH1F("recoManyPhoton", "3+ photons reconstructed", 60, 0, efficiencyXmax)
recoSignalHist = rt.TH1F("recoSignal", "1p2g successfully reconstructed", 60, 0, efficiencyXmax)

recoTotalHist = rt.TH1F("recoTotalHist", "All Reco Signal Events", 60, 0, 1500)
trueOutFiducialHist = rt.TH1F("trueOutsideFiducial", "True vertex outside fiducial volume", 60, 0, 1500)
trueCCHist = rt.TH1F("trueCC", "True event charged-current", 60, 0, 1500)
truePiPlusHist = rt.TH1F("truePiPlus", "True charged pion", 60, 0, 1500)
trueNoProtonHist = rt.TH1F("trueNoProton", "No true proton", 60, 0, 1500)
truePluralProtonHist = rt.TH1F("truePluralProton", "2+ true protons", 60, 0, 1500)
trueWrongProtonHist = rt.TH1F("trueWrongProton", "Reconstructed proton failed TID-matching", 60, 0, 1500)
trueNoPhotonHist = rt.TH1F("trueNoPhoton", "No true photon", 60, 0, 1500)
trueOnePhotonHist = rt.TH1F("trueOnePhoton", "One true photon", 60, 0, 1500)
trueManyPhotonHist = rt.TH1F("trueManyPhoton", "3+ true photons", 60, 0, 1500)
trueSignalHist = rt.TH1F("trueSignal", "1p2g successfully reconstructed", 60, 0, 1500)


eventSpace = 0
sameMID = 0
diffMID = 0
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
    nTrueProtons, trueProtonTID = trueProtonSelection(eventTree)
    if nTrueProtons != 1:
        continue
# True photon selection
    truePhotonTIDList, trueLeadingPhotonEnergy = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) != 2:
        continue

    scaledEnergy = []
    truePhotonIndices = []
    invariantMass = 0
    for i in range(eventTree.nTrueSimParts):
        if truePhotonTIDList[0] == eventTree.trueSimPartTID[i]:
            truePhotonIndices.append(i)
        if truePhotonTIDList[1] == eventTree.trueSimPartTID[i]:
            truePhotonIndices.append(i)
    for x in truePhotonIndices:
        scaledEnergy.append(eventTree.trueSimPartE[x])

    eventSpace += 1
    MID1 = eventTree.trueSimPartMID[truePhotonIndices[0]]
    MID2 = eventTree.trueSimPartMID[truePhotonIndices[1]]
    
    if MID1 == MID2:
        sameMID += 1
    if MID1 != MID2:
        diffMID += 1

    a = truePhotonIndices[0]
    b = truePhotonIndices[1]
    lengthA = np.sqrt(np.square(eventTree.trueSimPartPx[a]) + np.square(eventTree.trueSimPartPy[a]) + np.square(eventTree.trueSimPartPz[a]))
    lengthB = np.sqrt(eventTree.trueSimPartPx[b]**2 + eventTree.trueSimPartPy[b]**2 + eventTree.trueSimPartPz[b]**2)
    aDotB = eventTree.trueSimPartPx[a]*eventTree.trueSimPartPx[b] + eventTree.trueSimPartPy[a]*eventTree.trueSimPartPy[b] + eventTree.trueSimPartPz[a]*eventTree.trueSimPartPz[b]
    invariantMass = np.sqrt((scaledEnergy[0]*scaledEnergy[1]) - (scaledEnergy[0]*scaledEnergy[1])*((aDotB)/(lengthA*lengthB)))
    invariantMass = np.sqrt((2*lengthA*lengthB) - (2*aDotB))

    trueTotalHist.Fill(invariantMass, eventTree.xsecWeight)
#------------ Reco Matching ------------#

# start with whether vertex was found
    if eventTree.foundVertex == 0:
        recoNoVertexHist.Fill(invariantMass, eventTree.xsecWeight)
        continue
# reco cc check
    if recoCCCut(eventTree):
        recoCCHist.Fill(invariantMass, eventTree.xsecWeight)
        continue
# reco fiducial check
    if recoFiducialCut(eventTree, fiducialWidth):
        recoOutFiducialHist.Fill(invariantMass, eventTree.xsecWeight)
        continue
# reco pi+ check
    if recoPiPlusCut(eventTree):
        recoPiPlusHist.Fill(invariantMass, eventTree.xsecWeight)
        continue
# reco proton selection and filling
    nRecoProtons, recoProtonTID = recoProtonSelection(eventTree)
    if nRecoProtons == 0:
        recoNoProtonHist.Fill(invariantMass, eventTree.xsecWeight)
        continue
    if nRecoProtons >= 2:
        recoPluralProtonHist.Fill(invariantMass, eventTree.xsecWeight)
        continue

# uncomment below for tid-matching
#    if recoProtonTID != trueProtonTID:
#        recoWrongProtonHist.Fill(trueLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue

    recophotonTIDList, recoLeadingPhotonEnergy = recoPhotonSelection(eventTree, fiducialWidth)
    if len(recophotonTIDList) == 0:
        recoNoPhotonHist.Fill(invariantMass, eventTree.xsecWeight)
        continue
    if len(recophotonTIDList) == 1:
        recoOnePhotonHist.Fill(invariantMass, eventTree.xsecWeight)
        continue
    if len(recophotonTIDList) >= 3:
        recoManyPhotonHist.Fill(invariantMass, eventTree.xsecWeight)
        continue

    recoSignalHist.Fill(invariantMass, eventTree.xsecWeight)

#------------------ Purity Loop ------------------#
#------------------ Purity Loop ------------------#

for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
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
    nRecoProtons, recoProtonTID = recoProtonSelection(eventTree)
    if nRecoProtons != 1:
        continue
# Reco photon selection
    recoPhotonTIDList, recoLeadingPhotonEnergy = recoPhotonSelection(eventTree, fiducialWidth) 
    if len(recoPhotonTIDList) != 2:
        continue

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
    nTrueProtons, trueProtonTID = trueProtonSelection(eventTree)
    if nTrueProtons == 0:
        trueNoProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if nTrueProtons >= 2:
        truePluralProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
        continue

# uncomment below for tid-matching
#    if trueProtonTID != trueProtonTID:
#        trueWrongProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue

# true photon selection
    truePhotonTIDList, trueLeadingPhotonEnergy = truePhotonSelection(eventTree, fiducialWidth)
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
 
#------------------ End of Loops ------------------#

trueTotalHistInt = trueTotalHist.Integral(1, 60)
recoTotalHistInt = recoTotalHist.Integral(1, 60)

recoSignalHistInt = recoSignalHist.Integral(1, 60)
trueSignalHistInt = trueSignalHist.Integral(1, 60)

print("Efficiency: " + str(recoSignalHistInt / trueTotalHistInt * 100) + "%")
print("Purity: " + str(trueSignalHistInt / recoTotalHistInt * 100) + "%")

print(str(eventSpace) + " total true 2g events")
print(str(sameMID) + " events with shared photon MID")
print(str(diffMID) + " events with different photon MID")

print(str(ntuplePOTsum))

efficiencyHistList = [recoSignalHist, recoNoVertexHist, recoCCHist, \
                      recoOutFiducialHist, recoPiPlusHist, recoNoProtonHist, recoPluralProtonHist, \
                        recoNoPhotonHist, recoOnePhotonHist, recoManyPhotonHist]
purityHistList = [trueSignalHist, trueCCHist, trueOutFiducialHist, truePiPlusHist, \
                  trueNoProtonHist, truePluralProtonHist, trueNoPhotonHist, trueOnePhotonHist, \
                    trueManyPhotonHist]

efficiencyCanvas, efficiencyStack, efficiencyLegend, efficiencyInt = \
    histStackFill("Reconstruction of True NC 1p2g Events", efficiencyHistList, "Total Truth Signal: (", "True 2-Photon Invariant Mass (MeV)", "Events per 6.67e+20 POT", ntuplePOTsum)

purityCanvas, purityStack, purityLegend, purityInt = \
    histStackFill("Truth-Matching of Reconstructed NC 1p2g Events", purityHistList, "Total Reconstruction Signal: (", "Leading Reconstructed Photon Energy (MeV)", "Events per 6.67e+20 POT", ntuplePOTsum)

# Write Efficiency and Purity Canvases to Outfile
outFile = rt.TFile(args.outfile, "RECREATE")
efficiencyCanvas.Write("EfficiencyHist")
purityCanvas.Write("PurityHist")