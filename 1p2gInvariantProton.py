import sys, argparse
import numpy as np
import ROOT as rt
import math

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, trueInvariantMassCalculations, recoInvariantMassCalculations

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="9OctInvMassCutsNewDeltaTrueReco.root", help="output root file name")
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

# Define Histograms to be filled
histMin = 750
histMax = 1500

recoProtonMassHist = rt.TH1F("recoProtonInv", "reco signal of true nc 2p1g", 100, histMin, histMax)
trueProtonMassHist = rt.TH1F("trueProtonInv", "true nc 1p2g signal", 100, histMin, histMax)

truePhotonMassHist = rt.TH1F("truePhotonInv", "true nc 1p2g signal", 75, 0, 700)
recoPhotonMassHist = rt.TH1F("reco2photonInv", "reco signal of true nc 2p1g", 75, 0, 700)

trueDeltaMassHist = rt.TH1F("trueDeltaInv", "true nc 2p1g signal", 110, 900, 2000)
recoDeltaMassHist = rt.TH1F("recoDeltaInv", "reco signal of true nc 2p1g", 100, 900, 1900)

# Define min, max, fiducial and data function
fiducialWidth = 10
fiducialDict = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}

# Goals: This program should output reconstructed invariant mass of protons in true 1p2g events.
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
    truePhotonTIDList, trueLeadingPhotonEnergy, truePhotonIndexList = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) != 2:
        continue

    trueProtonInvMass, truePiZeroInvMass, trueDeltaInvMass = trueInvariantMassCalculations(eventTree, trueProtonIndex, truePhotonIndexList)

    trueProtonMassHist.Fill(trueProtonInvMass, eventTree.xsecWeight)
    truePhotonMassHist.Fill(truePiZeroInvMass, eventTree.xsecWeight)
    trueDeltaMassHist.Fill(trueDeltaInvMass, eventTree.xsecWeight)

#---------------- Purity Loop ----------------#
#---------------- Purity Loop ----------------#
eventNumber = 0
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
    eventNumber +=1 
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

#------------ TruthMatching ------------#

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
# TID-Matching of proton
    if recoProtonTID != trueProtonTID:
        continue
# True photon selection
    truePhotonTIDList, trueLeadingPhotonEnergy, truePhotonIndexList = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) != 2:
        continue

# Reconstructed proton invariant mass:
    recoProtonInvariantMass, recoPiZeroInvariantMass, recoDeltaInvariantMass =  recoInvariantMassCalculations(eventTree, recoProtonIndex, recoPhotonIndexList)
#    protonRestMass = 938.272089
#    protonEtot = eventTree.trackRecoE[recoProtonIndex] + protonRestMass
#    nx, ny, nz = eventTree.trackStartDirX[recoProtonIndex], eventTree.trackStartDirY[recoProtonIndex], eventTree.trackStartDirZ[recoProtonIndex]
#    momentumNorm = np.sqrt(np.square(protonEtot) - np.square(protonRestMass))
#    Px, Py, Pz = momentumNorm * nx, momentumNorm * ny, momentumNorm * nz

#    recoProtonInvariantMass = np.sqrt(np.square(protonEtot) - (np.square(Px) + np.square(Py) + np.square(Pz)))
    recoProtonMassHist.Fill(recoProtonInvariantMass, eventTree.xsecWeight)

# Reconstructed π0 invariant mass:
#    photonIndex1, photonIndex2 = recoPhotonIndexList[0], recoPhotonIndexList[1]
#    energy1, energy2 = eventTree.showerRecoE[photonIndex1], eventTree.showerRecoE[photonIndex2]
#    x1, y1, z1 = energy1*eventTree.showerStartDirX[photonIndex1], energy1*eventTree.showerStartDirY[photonIndex1], energy1*eventTree.showerStartDirZ[photonIndex1]
#    x2, y2, z2 = energy2*eventTree.showerStartDirX[photonIndex2], energy2*eventTree.showerStartDirY[photonIndex2], energy2*eventTree.showerStartDirZ[photonIndex2]
#    x, y, z = (x1+x2), (y1+y2), (z1+z2)
#    recoPhotonInvariantMass = np.sqrt(np.square(energy1 + energy2 ) - (((np.square(x))+(np.square(y))+np.square(z))))
    recoPhotonMassHist.Fill(recoPiZeroInvariantMass, eventTree.xsecWeight)

# Reconstructed delta invariant mass:
#    recoDeltaMass = np.sqrt(np.square(energy1 + energy2 + protonEtot) - (np.square(x+Px) + np.square(y+Py) + np.square(z+Pz)))
    recoDeltaMassHist.Fill(recoDeltaInvariantMass, eventTree.xsecWeight)
#------------------ End of Loops ------------------#

def histStackFill(title, histList, legendTitle, xTitle, yTitle, ntuplePOTSum):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.35, 0.5, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt. kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kBlack, rt.kOrange]
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

tProtonList, rProtonList = [trueProtonMassHist], [recoProtonMassHist]
tPhotonList, rPhotonList = [truePhotonMassHist], [recoPhotonMassHist]
tDeltaList, rDeltaList = [trueDeltaMassHist], [recoDeltaMassHist]

trueCanvas, trueStack, trueLegend, trueInt = histStackFill("Proton Invariant Mass in True 1p2g Events", tProtonList, "Truth Signal (", "Invariant Mass (MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

recoCanvas, recoStack, recoLegend, recoInt = histStackFill("Proton Invariant Mass in Reco. Truth-Matched 1p2g Events", rProtonList, "Reco Signal: (", "Invariant Mass(MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

tPhotonCanvas, tPhotonStack, tPhotonLegend, tPhotonInt = histStackFill("2-Photon Invariant Mass in True 1p2g Events", tPhotonList, "Truth Signal: (", "Invariant Mass(MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)
rPhotonCanvas, rPhotonStack, rPhotonLegend, rPhotonInt = histStackFill("2-Photon Invariant Mass in Reco. Truth-Matched 1p2g Events", rPhotonList, "Reco Signal: (", "Invariant Mass(MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

tDeltaCanvas, tDeltaStack, tDeltaLegend, tDeltaInt = histStackFill("Delta Invariant Mass in True 1p2g Events", tDeltaList, "Truth Signal: (", "Invariant Mass(MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

rDeltaCanvas, rDeltaStack, rDeltaLegend, rDeltaInt = histStackFill("Delta Invariant Mass in Reco. Truth-Matched 1p2g Events", rDeltaList, "Reco Signal: (", "Invariant Mass(MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

outFile = rt.TFile(args.outfile, "RECREATE")
trueCanvas.Write("True Proton")
recoCanvas.Write("Reco Proton")
tPhotonCanvas.Write("True Pi0")
rPhotonCanvas.Write("Reco Pi0")
tDeltaCanvas.Write("True Delta+")
rDeltaCanvas.Write("Reco Delta+")