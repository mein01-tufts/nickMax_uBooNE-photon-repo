import sys, argparse
import numpy as np
import ROOT as rt
import math

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, recoPhotonSelectionInvMass

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gCombinedOutputInvariant4.root", help="output root file name")
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

recoTotalHist = rt.TH1F("recoTotalHist", "All Reco Signal Events", 60, 0, 1000)
trueOutFiducialHist = rt.TH1F("trueOutsideFiducial", "True vertex outside fiducial volume", 60, 0, 1000)
trueCCHist = rt.TH1F("trueCC", "True event charged-current", 60, 0, 1000)
truePiPlusHist = rt.TH1F("truePiPlus", "True charged pion", 60, 0, 1000)
trueNoProtonHist = rt.TH1F("trueNoProton", "No true proton", 60, 0, 1000)
truePluralProtonHist = rt.TH1F("truePluralProton", "2+ true protons", 60, 0, 1000)
trueWrongProtonHist = rt.TH1F("trueWrongProton", "Reconstructed proton failed TID-matching", 60, 0, 1000)
trueNoPhotonHist = rt.TH1F("trueNoPhoton", "No true photon", 60, 0, 1000)
trueOnePhotonHist = rt.TH1F("trueOnePhoton", "One true photon", 60, 0, 1000)
trueManyPhotonHist = rt.TH1F("trueManyPhoton", "3+ true photons", 60, 0, 1000)
trueSignalHist = rt.TH1F("trueSignal", "1p2g successfully reconstructed", 60, 0, 1000)


eventSpace = 0
sameMID = 0
diffMID = 0
# Define min, max, fiducial and data function
fiducialWidth = 10
fiducialDict = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}

#---------------- Efficiency Loop ----------------#
#---------------- Efficiency Loop ----------------#

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
    recoPhotonTIDList, recoInvMass = recoPhotonSelectionInvMass(eventTree, fiducialWidth) 
    if len(recoPhotonTIDList) != 2:
        continue

    recoTotalHist.Fill(recoInvMass, eventTree.xsecWeight)


#------------ TruthMatching ------------#

# true cc check
    if trueCCCut(eventTree):
        trueCCHist.Fill(recoInvMass, eventTree.xsecWeight)
        continue
#true fiducial check
    if trueFiducialCut(eventTree, fiducialWidth):
        trueOutFiducialHist.Fill(recoInvMass, eventTree.xsecWeight)
        continue
# true pi+ check
    if truePiPlusCut(eventTree):
        truePiPlusHist.Fill(recoInvMass, eventTree.xsecWeight)
        continue
# true proton selection
    nTrueProtons, trueProtonTID = trueProtonSelection(eventTree)
    if nTrueProtons == 0:
        trueNoProtonHist.Fill(recoInvMass, eventTree.xsecWeight)
        continue
    if nTrueProtons >= 2:
        truePluralProtonHist.Fill(recoInvMass, eventTree.xsecWeight)
        continue

# uncomment below for tid-matching
#    if trueProtonTID != trueProtonTID:
#        trueWrongProtonHist.Fill(recoLeadingPhotonEnergy, eventTree.xsecWeight)
#        continue

# true photon selection
    truePhotonTIDList, trueLeadingPhotonEnergy = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) == 0:
        trueNoPhotonHist.Fill(recoInvMass, eventTree.xsecWeight)
        continue
    if len(truePhotonTIDList) == 1:
        trueOnePhotonHist.Fill(recoInvMass, eventTree.xsecWeight)
        continue
    if len(truePhotonTIDList) >= 3:
        trueManyPhotonHist.Fill(recoInvMass, eventTree.xsecWeight)
        continue
    
    trueSignalHist.Fill(recoInvMass, eventTree.xsecWeight)
 
#------------------ End of Loops ------------------#

#efficiencyHistList = [recoSignalHist, recoNoVertexHist, recoCCHist, \
#                      recoOutFiducialHist, recoPiPlusHist, recoNoProtonHist, recoPluralProtonHist, \
#                        recoNoPhotonHist, recoOnePhotonHist, recoManyPhotonHist]
purityHistList = [trueSignalHist, trueCCHist, trueOutFiducialHist, truePiPlusHist, \
                  trueNoProtonHist, truePluralProtonHist, trueNoPhotonHist, trueOnePhotonHist, \
                    trueManyPhotonHist]

def histStackFill(title, histList, legendTitle, xTitle, yTitle, ntuplePOTSum):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.35, 0.5, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt. kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kBlack, rt.kOrange]
  targetPOT = 3.2974607516739725e+20
  colors = [rt.kGreen+2, rt.kRed, rt. kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kBlack, rt.kYellow, rt.kGreen]
  targetPOT = 6.67e+20
  integralSum = 0
  sum = 0
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    hist.Scale(targetPOT/ntuplePOTSum)
    hist.SetFillColor(colors[x%7])
    hist.SetMarkerStyle(21)
    hist.SetMarkerColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    integralSum += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round((histInt), 2))+" events per 6.67e+20POT", "f")
    stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas(str(title)) 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle(str(xTitle))
  stack.GetYaxis().SetTitle(str(yTitle))
  legendHeaderString = str(str(legendTitle) + str(round((integralSum),1)) + " Events per 6.67e+20 POT)") 
  legend.SetHeader(str(legendHeaderString), "C")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

#efficiencyCanvas, efficiencyStack, efficiencyLegend, efficiencyInt = \
#    histStackFill("Reconstruction of True NC 1p2g Events", efficiencyHistList, "Total Truth Signal: (", "True 2-Photon Invariant Mass (MeV)", "Events per 6.67e+20 POT", ntuplePOTsum)

purityCanvas, purityStack, purityLegend, purityInt = \
    histStackFill("Truth-Matching of Reconstructed NC 1p2g Events", purityHistList, "Total Reconstruction Signal: (", "Reconstructed 2-Photon Invariant Mass (MeV)", "Events per 6.67e+20 POT", ntuplePOTsum)

# Write Efficiency and Purity Canvases to Outfile
outFile = rt.TFile(args.outfile, "RECREATE")
#efficiencyCanvas.Write("EfficiencyHist")
purityCanvas.Write("PurityHist")
trueTotalHist.Write("TrueTotal")