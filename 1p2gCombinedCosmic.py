import sys, argparse
import numpy as np
import ROOT as rt

#import selection commands from cuts.py
from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, histStackFill, trueCCCutLoose, recoCCCutLoose, recoInvariantMassCalculations

#takes imput arguments re: file locations
parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-c", "--cosmicFile", type=str, required=True, help="input cosmic ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="27Oct1p2gCombinedCosmic.root", help="output root file name")
args = parser.parse_args()

# Get relevant information from input ntuple file
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

# Get relevant info from cosmic ntuple file
cosmic_ntuple_file = rt.TFile(args.cosmicFile)
cosmicTree = cosmic_ntuple_file.Get("EventTree")
cosmicPotTree = cosmic_ntuple_file.Get("potTree")

# Scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

cosmicPOTsum = 3.2974607516739725e+20

# Define histograms to be filled
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

#------------------ End of Loop ------------------#

#---------------- cosmicData Loop ----------------#
# Perform reco cuts on cosmic data to determine number of cosmics to be tagged as signal
cosmicEventsTally = 0
for i in range(cosmicTree.GetEntries()):
    cosmicTree.GetEntry(i)

# Keep a tally of all cosmic events    
    cosmicEventsTally += 1
# Cosmic CC Cut:
    if recoCCCut(cosmicTree):
        continue

# Cosmic Fiducial Cut:
    if recoFiducialCut(cosmicTree):
        continue
    