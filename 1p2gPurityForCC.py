import sys, argparse
import numpy as np
import ROOT as rt

from cuts import histStackFill, kineticEnergyCalculator, sStackFillNS, particleDistancesCalculator

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gCC725.root", help="output root file name")
args = parser.parse_args()

# open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

# scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles

# calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

# define histograms to be filled
muonEnergyHist = rt.TH1F("muonEnergyHist", "Primary Muons in True CC Reco NC Events",60,0,2500)
allPrimaryMuonHist = rt.TH1F("allPrimaryMuonHist", "Primary Muons in True CC Events",60,0,2500)
muonTrackLengthHist = rt.TH1F("muonTrackHist", "Primary Muon Track Length in True CC, Reco NC Events",60,0,1100)
allMuonTrackHist = rt.TH1F("allMuonTrackHist", "Primary Muon Track Length in True CC Events",60,0,1100)
muonDetectorDistanceHist = rt.TH1F("muonDistHist", "Primary Muon End Distance to Detector Wall in True CC Events",60,0,120)
allMuonDetectorDistanceHist = rt.TH1F("allMuonDistHist", "Primary Muon Track End Distance to Detector Wall in True CC Events",60,0,120)

xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

chargeCurrentEvents = 0
muonsObservedEnergy = 0
muonsObservedLength = 0
electronsObserved = 0
recoProtonIsProton = 0
recoProtonIsMuon = 0
recoProtonIsOther = 0
nonClassifiedTracks = 0
nonClassifiedTrackIsMuon = 0
recoMuons = 0
# begin loop over events in ntuple file:
# start with culling for reco-defined 1p2g signal, then reiterate using true to
# determine what reco actually saw
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

# find and fill typical primary muon energy of events
#    if eventTree.trueNuCCNC != 1:
#        for i in range(eventTree.nTruePrimParts):
#           if abs(eventTree.truePrimPartPDG[i]) == 13:
#                muonsObservedEnergy += 1
#                muonKE_MeV = (kineticEnergyCalculator(eventTree, i))*1000
#                allPrimaryMuonHist.Fill(muonKE_MeV, 1)

# find and fill typical primary muon track length of cc events
    if eventTree.trueNuCCNC != 1:
        for i in range(eventTree.nTrueSimParts):
            if eventTree.trueSimPartProcess[i] == 0 and eventTree.trueSimPartPDG[i] == 13:
                muonsObservedLength += 1
                muonTrackLength, muonDetectorDistance = particleDistancesCalculator(eventTree, i)
                allMuonTrackHist.Fill(muonTrackLength, eventTree.xsecWeight)
                allMuonDetectorDistanceHist.Fill(muonDetectorDistance, eventTree.xsecWeight)

# event vertex reconstructed cut
    if eventTree.foundVertex != 1:
        continue

# fiducial volume cut
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
        eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
        eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
        continue

# reco "cosmic cut"
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        continue

# reco charged-current cut:        
# iterate through all tracks/showers in event,
# look for non-secondary tracks identified as muons or electrons
# and for non-secondary showers identified as electrons
    recoPrimaryMuonFound = False
    recoPrimaryElectronTrackFound = False
    recoPrimaryElectronShowerFound = False
    for i in range(eventTree.nTracks):
        if eventTree.trackIsSecondary[i] == 0:
            if abs(eventTree.trackPID[i]) == 13:
                recoPrimaryMuonFound = True
            if abs(eventTree.trackPID[i]) == 11:
                recoPrimaryElectronTrackFound = True
    for i in range(eventTree.nShowers):
        if eventTree.showerIsSecondary[i] == 0:
            if abs(eventTree.showerPID[i]) == 11:
                recoPrimaryElectronShowerFound == True
# cut all events with primary electrons/muons
    if recoPrimaryMuonFound or recoPrimaryElectronTrackFound or recoPrimaryElectronShowerFound:    
        continue

# reco pi+ cut: cut all events w/ pi+ of KE >= 30MeV
    recoPiPlusPresent = False
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 211:
                if eventTree.trackRecoE[i] >= 30:
                    recoPiPlusPresent = True
                    break
    if recoPiPlusPresent:
        continue

# reco protons cut: all protons of KE >= 60MeV 
    nRecoProtons = 0
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 2212:
            if eventTree.trackRecoE[i] >= 60.0:
                nRecoProtons += 1
                recoProtonPID = eventTree.trackTruePID[i]
                recoProtonTID = eventTree.trackTrueTID[i]
    if nRecoProtons != 1:
        continue

# reco photon cut: find all photons, then append truth-matched TID to a list
    recoPhotonTIDList = []
    recoPhotonIndexList = []
    for i in range(eventTree.nShowers):
        if eventTree.showerPID[i] == 22:
            recoPhotonTIDList.append(eventTree.showerTrueTID[i])
            recoPhotonIndexList.append(i)
# cut events w/o two reco photons    
    if len(recoPhotonTIDList) != 2:
        continue
    #-------------- end reco selection --------------#
    #------------- begin truth matching -------------#

# true charged-current selection: Find true cc events, then fill primary muon's KE to histogram
    if eventTree.trueNuCCNC == 1:
        continue
    chargeCurrentEvents += 1

# find track length of all primary muons
    for i in range(eventTree.nTrueSimParts):
        if eventTree.trueSimPartProcess[i] == 0 and eventTree.trueSimPartPDG[i] == 13:
            muonTrackLength, muonDetectorDistance = particleDistancesCalculator(eventTree, i)
            muonTrackLengthHist.Fill(muonTrackLength, eventTree.xsecWeight)
            muonDetectorDistanceHist.Fill(muonDetectorDistance, eventTree.xsecWeight)

#    for i in range(eventTree.nTruePrimParts):
#        if abs(eventTree.truePrimPartPDG[i]) == 13:
#            muonsObserved += 1
#            muonKE_MeV = (kineticEnergyCalculator(eventTree, i))*1000
#            muonEnergyHist.Fill(muonKE_MeV, 1)
#        if abs(eventTree.truePrimPartPDG[i]) == 11:
#            electronsObserved += 1

#----- end of event loop ---------------------------------------------#

def sStackFillS(title, hist, kColor, canvasTitle, xAxis):
#Forms a filled stacked histogram based on only one. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.45, 0.8, 0.9, 0.9)
  targetPOT = 6.67e+20
  ntuplePOTSum = 4.675690535431973e+20
  bins = hist.GetNbinsX()
  hist.Scale(targetPOT/ntuplePOTSum)
  hist.SetFillColor(kColor)
  hist.SetMarkerStyle(21)
  hist.SetMarkerColor(kColor)
  histInt = hist.Integral(1, int(bins))
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1))+" muons per 6.67e+20 POT", "f")
  stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas(str(canvasTitle)) 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle(str(xAxis))
  stack.GetYaxis().SetTitle("Muons per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

def sStackFillS2(title, hist, kColor, canvasTitle, xAxis):
#Forms a filled stacked histogram based on only one. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.45, 0.8, 0.9, 0.9)
  targetPOT = 6.67e+20
  ntuplePOTSum = 4.675690535431973e+20
  bins = hist.GetNbinsX()
  hist.Scale(targetPOT/ntuplePOTSum)
  hist.SetFillColor(kColor)
  hist.SetMarkerStyle(21)
  hist.SetMarkerColor(kColor)
  histInt = hist.Integral(1, int(bins))
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1))+" muons per 6.67e+20 POT", "f")
  stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas(str(canvasTitle)) 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle(str(xAxis))
  stack.GetYaxis().SetTitle("Muons per 6.67e+20 POT")
  legend.Draw()
  rt.gPad.SetLogy(1)
  rt.gPad.Update()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt


#muonEnergyCanvas, muonEnergyStack, muonEnergyLegend, muonEnergyInt = sStackFillS("Muon Energy of Reco NC, True CC Events", muonEnergyHist, rt.kBlue, "muon energy")

#allMuonCanvas, allMuonStack, allMuonLegend, allMuonInt = sStackFillS("Muon Energy in True CC Events", allPrimaryMuonHist, rt.kBlue, "all muon energy")

#outFile = rt.TFile(args.outfile, "RECREATE")
#muonEnergyCanvas.Write()
#allMuonCanvas.Write()

muonTrackCanvas, muonTrackStack, muonTrackLegend, muonTrackInt = sStackFillS("Muon Track Length of True CC, Reco NC Events", muonTrackLengthHist, rt.kBlue, "muon track length", "Muon Track Length (cm)")

allMuonTrackCanvas, allMuonTrackStack, allMuonTrackLegend, allMuonTrackInt = sStackFillS("Muon Track Length in True CC Events", allMuonTrackHist, rt.kBlue, "all muon track length", "Muon Track Length (cm)")

muonDistCanvas, muonDistStack, muonDistLegend, muonDistkInt = sStackFillS2("Muon-Detector Distance of True CC, Reco NC Events", muonDetectorDistanceHist, rt.kBlue, "muon distance", "Muon-Detector Distance (cm)")

allMuonDistCanvas, allMuonDistStack, allMuonDistLegend, allMuonDistkInt = sStackFillS2("Muon-Detector Distance of True CC Events", allMuonDetectorDistanceHist, rt.kBlue, "all muon distance", "Muon-Detector Distance (cm)")


outFile = rt.TFile(args.outfile, "RECREATE")
muonTrackCanvas.Write()
muonDistCanvas.Write()
allMuonTrackCanvas.Write()
allMuonDistCanvas.Write()


print(str(muonsObservedEnergy) + "muons counted from trueprimparts")
print(str(muonsObservedLength) + "muons counted from truesimparts")
print(str(chargeCurrentEvents) + " total cc mis-ID events")

