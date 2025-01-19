import sys, argparse
import numpy as np
import ROOT as rt

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, histStackFill, trueCCCutLoose, recoCCCutLoose, recoInvariantMassCalculations

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="17JanDataPlots.root", help="output root file name")
args = parser.parse_args()

# Grab ntuple file, eventTree, and potTree for reference
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

dataPhotonHist = rt.TH1F("DataEventsGamma", "Data Events with 2g + 1p", 60, 0, 1200)

dataDeltaHist = rt.TH1F("DataEventsDelta", "Data Events with 2g + 1p, delta", 600, 800, 1800)

dataPiZeroHist = rt.TH1F("DataEventsPiZero", "Data Events with 2g + 1p. pi zero", 60, 0, 800)

# Define min, max, fiducial and data function
fiducialWidth = 10
fiducialDict = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}

ccCut = 0
fiducialCut = 0
piPlusCut = 0
protonCut = 0
photonCut = 0
dataTally = 0
signalTally = 0

# Start loop for cuts
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
    dataTally += 1
# Charged-current cut
    if recoCCCut(eventTree):
        ccCut += 1
        continue
# Fiducial volume cut
    if recoFiducialCut(eventTree, fiducialWidth):
        fiducialCut += 1
        continue
# "Cosmic cut"
#    if eventTree.vtxFracHitsOnCosmic >= 1.:
#        continue
# Pi+ cut
    if recoPiPlusCut(eventTree):
        piPlusCut += 1
        continue
# Proton selection
    nRecoProtons = 0
    recoProtonIndex = 0
    for i in range(eventTree.nTracks):
        if eventTree.trackPID[i] == 2212 and eventTree.trackRecoE[i] >= 60:
            nRecoProtons += 1
            recoProtonIndex = i
    if nRecoProtons != 1:
        protonCut += 1
        continue
# Photon selection
    reco = 0
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
                recoPhotonIndexList.append(i)

    recoPhotonEnergyList = [0]
    recoLeadingPhotonEnergy = 0.
  
    for i in range(len(recoPhotonIndexList)):
        recoPhotonEnergyMeV = eventTree.showerRecoE[recoPhotonIndexList[i]] 
        recoPhotonEnergyList.append(recoPhotonEnergyMeV)
        recoLeadingPhotonEnergy = max(recoPhotonEnergyList)
 
    if len(recoPhotonIndexList) != 2:
        photonCut += 1
        continue

    signalTally += 1
    
    dataPhotonHist.Fill(recoLeadingPhotonEnergy, 1)

    recoProtonInv, recoPiZeroInv, recoDeltaInv = recoInvariantMassCalculations(eventTree, recoProtonIndex, recoPhotonIndexList)
    
    dataDeltaHist.Fill(recoDeltaInv, 1)
    dataPiZeroHist.Fill(recoPiZeroInv, 1)

def histStackFill(title, histList, legendTitle, xTitle, yTitle):
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
    hist.SetFillColor(colors[x%7])
    hist.SetMarkerStyle(21)
    hist.SetMarkerColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    integralSum += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round((histInt), 2))+" events per 4.4e+19POT", "f")
    stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas(str(title)) 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle(str(xTitle))
  stack.GetYaxis().SetTitle(str(yTitle))
  legendHeaderString = str(str(legendTitle) + str(round((integralSum),1)) + " Events per 4.4e+19 POT)") 
  legend.SetHeader(str(legendHeaderString), "C")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

dataList = [dataPhotonHist]
deltaList = [dataDeltaHist]
piZeroList = [dataPiZeroHist]

dataCanvas, dataStack, dataLegend, dataInt = histStackFill("Data Events with 2g + 1p", dataList, "Total: ", "Reconstructed Leading Photon Energy (MeV)", "events per 4.4e+19 POT", )

deltaCanvas, deltaStack, deltaLegend, deltaInt = histStackFill("Data Events with 2g + 1p, delta mass", deltaList, "Total: ", "Reconstructed 2g + 1p Invariant Mass (MeV)", "events per 4.4e+19 POT", )
piZeroCanvas, piZeroStack, piZeroLegend, piZeroInt = histStackFill("Data Events with 2g + 1p, pi0 mass", piZeroList, "Total: ", "Reconstructed 2g Invariant Mass (MeV)", "events per 4.4e+19 POT", )

outFile = rt.TFile(args.outfile, "RECREATE")
dataCanvas.Write("LeadingPhoton")
deltaCanvas.Write("DeltaMass")
piZeroCanvas.Write("PiZeroMass")

print("Of " + str(dataTally) + " total data events: ")
print(str(ccCut) + " were dropped by the charged-current cut")
print(str(fiducialCut) + " were dropped by the fiducial cut")
print(str(piPlusCut) + " were dropped by the pi+ cut")
print(str(protonCut) + " were dropped by the proton cut")
print(str(photonCut) + " were dropped by the photon cut")
print(str(signalTally) + " are reconstructed as signal events")