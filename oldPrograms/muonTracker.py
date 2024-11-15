import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)

from cuts import trueCutNC, trueCutFiducials, trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
args = parser.parse_args()

#needed for proper scaling of error bars:
#rt.TH1.SetDefaultSumw2(rt.kTRUE)

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20"

#calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

#HISTOGRAMS DEFINED AND PREPARED HERE
effHistE1 = rt.TH1F("EnergyEfficiency1", "Vertex < 1 cm from true (energy)",30,0,6)
effHistE2 = rt.TH1F("EnergyEfficiency2", "Vertex 1-3 cm from true (energy)",30,0,6)
effHistE3 = rt.TH1F("EnergyEfficiency3", "Vertex 3-5 cm from true (energy)",30,0,6)
effHistE4 = rt.TH1F("EnergyEfficiency4", "Vertex 5+ cm from true (energy)",30,0,6)

effListE = [effHistE1, effHistE2, effHistE3, effHistE4]

effTotalE = rt.TH1F("TotalEfficiencyE", "Total Efficiency",30,0,6)

foundPhotonE1 = rt.TH1F("FoundE1", "Found Muon",30,0,6)
lostPhotonE1 = rt.TH1F("LostE1", "Track not found",30,0,6)
weirdE1 = rt.TH1F("WeirdE1", "Misclassified track",30,0,6)
unclassifiedE1 = rt.TH1F("UnclassifiedE1", "Unclassified track",30,0,6)

foundPhotonE2 = rt.TH1F("FoundE2", "Found Muon",30,0,6)
lostPhotonE2 = rt.TH1F("LostE2", "Track not found",30,0,6)
weirdE2 = rt.TH1F("WeirdE2", "Misclassified track",30,0,6)
unclassifiedE2 = rt.TH1F("UnclassifiedE2", "Unclassified track",30,0,6)

foundPhotonE3 = rt.TH1F("FoundE3", "Found Muon",30,0,6)
lostPhotonE3 = rt.TH1F("LostE3", "Track not found",30,0,6)
weirdE3 = rt.TH1F("WeirdE3", "Misclassified track",30,0,6)
unclassifiedE3 = rt.TH1F("UnclassifiedE3", "Unclassified track",30,0,6)

foundPhotonE4 = rt.TH1F("FoundE4", "Found Photons",30,0,6)
lostPhotonE4 = rt.TH1F("LostE4", "Track not found",30,0,6)
weirdE4 = rt.TH1F("WeirdE4", "Misclassified track",30,0,6)
unclassifiedE4 = rt.TH1F("UnclassifiedE4", "Unclassified track",30,0,6)

energyList1 = [lostPhotonE1, weirdE1, unclassifiedE1, foundPhotonE1]
energyList2 = [lostPhotonE2, weirdE2, unclassifiedE2, foundPhotonE2]
energyList3 = [lostPhotonE3, weirdE3, unclassifiedE3, foundPhotonE3]
energyList4 = [lostPhotonE4, weirdE4, unclassifiedE4, foundPhotonE4]

energyListList = [energyList1, energyList2, energyList3, energyList4]

foundListE = [foundPhotonE1, foundPhotonE2, foundPhotonE3, foundPhotonE4]
lostListE = [lostPhotonE1, lostPhotonE2, lostPhotonE3, lostPhotonE4]
weirdListE = [weirdE1, weirdE2, weirdE3, weirdE4]
unclassifiedListE = [unclassifiedE1, unclassifiedE2, unclassifiedE3, unclassifiedE4]

#Functions for making histograms
def addHist(ntuple, energy, HistList1, weight):
  if ntuple.vtxDistToTrue <= 1:
    HistList1[0].Fill(energy, weight)
  elif ntuple.vtxDistToTrue <= 3:
    HistList1[1].Fill(energy, weight)
  elif ntuple.vtxDistToTrue <= 5:
    HistList1[2].Fill(energy, weight)
  else:
    HistList1[3].Fill(energy, weight)

def addEffHist(effHist, histList):
  for x in range(1, 31):
    if histList[3].GetBinContent(x)+histList[0].GetBinContent(x)+histList[1].GetBinContent(x)+histList[2].GetBinContent(x) > 0:
      efficiency = histList[3].GetBinContent(x)/(histList[3].GetBinContent(x)+histList[0].GetBinContent(x)+histList[1].GetBinContent(x)+histList[2].GetBinContent(x))
    else:
      efficiency = 0
    effHist.SetBinContent(x, efficiency)


#Variables for program review
muonCount = 0
recoCount = 0
matchCount = 0
recoRight = 0
noMatch = 0

#Variables for program function
recoPhotonIDs = []
truePhotonIDs = []
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}

#Event loop begins
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #RECO CUTS
  if recoNoVertex(eventTree) == False:
    continue

  #Make sure the event is within the fiducial volume
  if recoFiducials(eventTree, fiducialData) == False:
    continue

  #Cut events with suitably energetic protons or charged pions
  if recoPionProton(eventTree) == False:
    continue

  #Since we're investigating CC events, we want to check for photons to ensure we're getting an otherwise accurate representation of our dataset
  recoPhotons = recoPhotonList(eventTree)
  if len(recoPhotons) == 0:
    continue
  
  #FIND MUON EVENTS
  #Get a list of true photons, scaled energy, and invariant mass
  muonList = []
  for x in range(eventTree.nTrueSimParts):
    if eventTree.trueSimPartPDG[x] == 13 and eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
      muonList.append(x)

  #Remove events unless they have a primary muon
  if len(muonList) == 0:
    continue
  muonCount += 1
  
  #List photons the reco can find
  recoMuons = []
  for x in range(eventTree.nTracks):
    if eventTree.trackPID[x] == 13 and eventTree.trackRecoE[x] > 20:
      recoMuons.append(x)

  if len(recoMuons) > 0:
    recoCount += 1
  #See how the algorithm treats each reco photon
  matches = 0
  distVtx = eventTree.vtxDistToTrue
  for x in muonList:
    recoMatch = False
    #Quick setup
    for y in range(eventTree.nTracks):
      #Match the true event to a reco event, if we can
      if eventTree.trackTrueTID[y] == eventTree.trueSimPartTID[x]:
        recoMatch = True
        matchCount += 1
        #If the track was correctly identified by the reco, we store it here
        if y in recoMuons:
          addHist(eventTree, eventTree.trueSimPartE[x]/1000, foundListE, eventTree.xsecWeight)
          recoRight += 1
        #If the track wasn't classified, we store it here
        elif eventTree.trackClassified[y] == 0:
          addHist(eventTree, eventTree.trueSimPartE[x]/1000, unclassifiedListE, eventTree.xsecWeight)
        #If it somehow got identified, but not as an electron or a photon, we should mark that down too
        else:
          addHist(eventTree, eventTree.trueSimPartE[x]/1000, weirdListE, eventTree.xsecWeight)

    if recoMatch == True:
      matches += 1
    #If we didn't find any track at all, we store it here
    else:
      addHist(eventTree, eventTree.trueSimPartE[x]/1000, lostListE, eventTree.xsecWeight)
      noMatch += 1

#Fill efficiency histograms
for x in range(len(effListE)):
  addEffHist(effListE[x], energyListList[x])

#Fill histograms for total efficiency
for x in range(1, 31):
  totalFound = foundPhotonE1.GetBinContent(x) + foundPhotonE2.GetBinContent(x) + foundPhotonE3.GetBinContent(x) + foundPhotonE4.GetBinContent(x)
  totalLost = lostPhotonE1.GetBinContent(x) + lostPhotonE2.GetBinContent(x) + lostPhotonE3.GetBinContent(x) + lostPhotonE4.GetBinContent(x)
  totalOther = weirdE1.GetBinContent(x) + weirdE2.GetBinContent(x) + weirdE3.GetBinContent(x) + weirdE4.GetBinContent(x) + unclassifiedE1.GetBinContent(x) + unclassifiedE2.GetBinContent(x) + unclassifiedE3.GetBinContent(x) +unclassifiedE4.GetBinContent(x)
  if totalFound + totalLost + totalOther > 0:
    totalE = totalFound/(totalFound + totalLost + totalOther)
  else:
    totalE = 0
  effTotalE.SetBinContent(x, totalE)
  
#effCanvas, effStack, effLegend, effInt = histStack("Efficiencies of photon showers over distance", effList)
#effCanvasE, effStackE, effLegendE, effIntE = histStack("Efficiencies of photon showers over energy", effListE)

effTotalE.SetLineColor(rt.kRed)

#Create canvases for efficiency histograms
effCanvasE1 = rt.TCanvas("effCanvasE1","effCanvasE1")
effHistE1.SetLineColor(rt.kGreen+2)
effHistE1.Draw()
effTotalE.Draw("Same")

effCanvasE2 = rt.TCanvas("effCanvasE2","effCanvasE2")
effHistE2.SetLineColor(rt.kGreen+2)
effHistE2.Draw()
effTotalE.Draw("Same")

effCanvasE3 = rt.TCanvas("effCanvasE3","effCanvasE3")
effHistE3.SetLineColor(rt.kGreen+2)
effHistE3.Draw()
effTotalE.Draw("Same")

effCanvasE4 = rt.TCanvas("effCanvasE4","effCanvasE4")
effHistE4.SetLineColor(rt.kGreen+2)
effHistE4.Draw()
effTotalE.Draw("Same")


#Create canvases for energy histograms
histCanvasE1, stackE1, legendE1, histIntE1 = histStack("Energy - Vertex within 1 cm of true vertex", energyList1, ntuplePOTsum)
histCanvasE2, stackE2, legendE2, histIntE2 = histStack("Energy - Vertex 1-3 cm from true vertex", energyList2, ntuplePOTsum)
histCanvasE3, stackE3, legendE3, histIntE3 = histStack("Energy - Vertex 3-5 cm from true vertex", energyList3, ntuplePOTsum)
histCanvasE4, stackE4, legendE4, histIntE4 = histStack("Energy - Vertex more than 5 cm from true vertex", energyList4, ntuplePOTsum)

#Make canvas list so we can export efficiently
canvasList = [effCanvasE1, effCanvasE2, effCanvasE3, effCanvasE4, histCanvasE1, histCanvasE2, histCanvasE3, histCanvasE4]

#Change axes
for hist in effListE:
  hist.GetXaxis().SetTitle("True Energy (GeV)")
  hist.GetYaxis().SetTitle("Efficiency")

stackListE = [stackE1, stackE2, stackE3, stackE4]
for stack in stackListE:
  stack.GetXaxis().SetTitle("True Muon Energy (GeV)")
  stack.GetYaxis().SetTitle("Number of Muons (Scaled to 6.67e+20)")

#Save to file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in canvasList:
  canvas.Write()


print("Total with true muons:", muonCount)
print("Total with reconstructed muons:", recoCount)
print("Total where true muon matched with any track", matchCount)
print("Total where the reco correctly reconstructed the muon", recoRight)
print("Total where the reco missed the track completely", noMatch)