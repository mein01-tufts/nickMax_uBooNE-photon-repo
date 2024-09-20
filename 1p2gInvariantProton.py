import sys, argparse
import numpy as np
import ROOT as rt
import math

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="20SepInvMassProtonPhotonDeltaTrueReco.root", help="output root file name")
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
histMax = 1250

recoProtonMassHist = rt.TH1F("recoProtonInv", "reco signal of true nc 2p1g", 100, histMin, histMax)
trueProtonMassHist = rt.TH1F("trueProtonInv", "true nc 1p2g signal", 100, histMin, histMax)
recoPhotonMassHist = rt.TH1F("reco2photonInv", "reco signal of true nc 2p1g", 75, 0, 700)
trueDeltaMassHist = rt.TH1F("trueDeltaInv", "true nc 2p1g signal", 100, 900, 1900)
recoDeltaMassHist = rt.TH1F("recoDeltaInv", "reco signal of true nc 2p1g", 100, 900, 1900)

# Define min, max, fiducial and data function
fiducialWidth = 10
fiducialDict = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}

def recoPhotonSelection(eventTree, fiducialWidth):
    reco = 0
    recoPhotonTIDList = []
    recoPhotonIndexList = []
    xMin, xMax = 0, 256
    yMin, yMax = -116.5, 116.5
    zMin, zMax = 0, 1036
    for i in range(eventTree.nShowers):
        if eventTree.showerPID[i] == 22:
            if eventTree.showerStartPosX[i] <= (xMin + fiducialWidth) or eventTree.showerStartPosX[i] >= (xMax - fiducialWidth) or \
                eventTree.showerStartPosY[i] <= (yMin + fiducialWidth) or eventTree.showerStartPosY[i] >= (yMax - fiducialWidth) or \
                    eventTree.showerStartPosZ[i] <= (zMin + fiducialWidth) or eventTree.showerStartPosZ[i] >= (zMax - fiducialWidth):
                reco += 1
            else:
                recoPhotonTIDList.append(eventTree.showerTrueTID[i])
                recoPhotonIndexList.append(i)

    return recoPhotonTIDList, recoPhotonIndexList


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
    truePhotonTIDList, trueLeadingPhotonEnergy = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) != 2:
        continue
# True proton invariant mass:
    for i in range(eventTree.nTrueSimParts):
        if eventTree.trueSimPartTID[i] == trueProtonTID:
            pNum = i
    
    protonE = eventTree.trueSimPartE[pNum]
    protonPX, protonPY, protonPZ = eventTree.trueSimPartPx[pNum], eventTree.trueSimPartPy[pNum], eventTree.trueSimPartPz[pNum]
    protonTrueInv = np.sqrt(np.square(protonE) - (np.square(protonPX) + np.square(protonPY) + np.square(protonPZ)))

    trueProtonMassHist.Fill(protonTrueInv, eventTree.xsecWeight)
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
    recoPhotonTIDList, recoPhotonIndexList = recoPhotonSelection(eventTree, fiducialWidth) 
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
    truePhotonTIDList, trueLeadingPhotonEnergy = truePhotonSelection(eventTree, fiducialWidth)
    if len(truePhotonTIDList) != 2:
        continue

# True proton invariant mass:
    protonE = eventTree.trueSimPartE[trueProtonIndex]
    protonPX, protonPY, protonPZ = eventTree.trueSimPartPx[trueProtonIndex], eventTree.trueSimPartPy[trueProtonIndex], eventTree.trueSimPartPz[trueProtonIndex]
    protonTrueInv = np.sqrt(np.square(protonE) - (np.square(protonPX) + np.square(protonPY) + np.square(protonPZ)))

# Reconstructed proton invariant mass:
# I think the last one was throwing errors like this: index cant be found -> sentinel value for E -> negative square root
    if eventTree.trackRecoE[recoProtonIndex] < 0:
        print("event " + str(eventNumber) + " " + str(eventTree.trackRecoE[recoProtonIndex]))
    protonRestMass = 938.272089
    protonEtot = eventTree.trackRecoE[recoProtonIndex] + protonRestMass
    nx, ny, nz = eventTree.trackStartDirX[recoProtonIndex], eventTree.trackStartDirY[recoProtonIndex], eventTree.trackStartDirZ[recoProtonIndex]
    if nx == -9 or ny == -9 or nz == -9:
        print("event " + str(eventNumber) + ": sentinel value for proton track")
    momentumNorm = np.sqrt(np.square(protonEtot) - np.square(protonRestMass))
    Px, Py, Pz = momentumNorm * nx, momentumNorm * ny, momentumNorm * nz

    recoProtonInvariantMass = np.sqrt(np.square(protonEtot) - (np.square(Px) + np.square(Py) + np.square(Pz)))

    recoProtonMassHist.Fill(recoProtonInvariantMass, eventTree.xsecWeight)

# Reconstructed 2-photon invariant mass:
    photonIndex1, photonIndex2 = recoPhotonIndexList[0], recoPhotonIndexList[1]
    energy1, energy2 = eventTree.showerRecoE[photonIndex1], eventTree.showerRecoE[photonIndex2]
    x1, y1, z1 = energy1*eventTree.showerStartDirX[photonIndex1], energy1*eventTree.showerStartDirY[photonIndex1], energy1*eventTree.showerStartDirZ[photonIndex1]
    x2, y2, z2 = energy2*eventTree.showerStartDirX[photonIndex2], energy2*eventTree.showerStartDirY[photonIndex2], energy2*eventTree.showerStartDirZ[photonIndex2]
    x, y, z = (x1+x2), (y1+y2), (z1+z2)
    recoPhotonInvariantMass = np.sqrt(np.square(energy1 + energy2 ) - (((np.square(x))+(np.square(y))+np.square(z))))
    recoPhotonMassHist.Fill(recoPhotonInvariantMass, eventTree.xsecWeight)

    recoDeltaMass = np.sqrt(np.square(energy1 + energy2 + protonEtot) - (np.square(x+Px) + np.square(y+Py) + np.square(z+Pz)))
    recoDeltaMassHist.Fill(recoDeltaMass, eventTree.xsecWeight)
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


truthList = [trueProtonMassHist]
recoList = [recoProtonMassHist]
photonList = [recoPhotonMassHist]
deltaList = [recoDeltaMassHist]

trueCanvas, trueStack, trueLegend, trueInt = histStackFill("Proton Invariant Mass in True 1p2g Events", truthList, "Truth Signal (", "Invariant Mass (MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

recoCanvas, recoStack, recoLegend, recoInt = histStackFill("Proton Invariant Mass in Reco. Truth-Matched 1p2g Events", recoList, "Reco Signal: (", "Invariant Mass(MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

photonCanvas, photonStack, photonLegend, photonInt = histStackFill("2-Photon Invariant Mass in Reco. Truth-Matched 1p2g Events", photonList, "Reco Signal: (", "Invariant Mass(MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

deltaCanvas, deltaStack, deltaLegend, deltaInt = histStackFill("Delta Invariant Mass in Reco. Truth-Matched 1p2g Events", deltaList, "Reco Signal: (", "Invariant Mass(MeV/c^2)", "Events per 6.67e+20POT", ntuplePOTsum)

outFile = rt.TFile(args.outfile, "RECREATE")
trueCanvas.Write("True Proton")
recoCanvas.Write("Reco Proton")
photonCanvas.Write("Reco Pi0")
deltaCanvas.Write("Reco Delta+")