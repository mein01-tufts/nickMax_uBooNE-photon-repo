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
effHist1 = rt.TH1F("DistanceEfficiency1", "Vertex < 1 cm from true (distance)",30,0,200)
effHist2 = rt.TH1F("DistanceEfficiency2", "Vertex 1-3 cm from true (distance)",30,0,200)
effHist3 = rt.TH1F("DistanceEfficiency3", "Vertex 3-5 cm from true (distance)",30,0,200)
effHist4 = rt.TH1F("DistanceEfficiency4", "Vertex 5+ cm from true (distance)",30,0,200)

effHistE1 = rt.TH1F("EnergyEfficiency1", "Vertex < 1 cm from true (energy)",30,0,1.6)
effHistE2 = rt.TH1F("EnergyEfficiency2", "Vertex 1-3 cm from true (energy)",30,0,1.6)
effHistE3 = rt.TH1F("EnergyEfficiency3", "Vertex 3-5 cm from true (energy)",30,0,1.6)
effHistE4 = rt.TH1F("EnergyEfficiency4", "Vertex 5+ cm from true (energy)",30,0,1.6)


effList = [effHist1, effHist2, effHist3, effHist4]
effListE = [effHistE1, effHistE2, effHistE3, effHistE4]

effTotal = rt.TH1F("TotalEfficiency", "Total Efficiency",30,0,200)
effTotalE = rt.TH1F("TotalEfficiencyE", "Total Efficiency",30,0,1.6)

foundPhotonHist1 = rt.TH1F("Found1", "Found Photons",30,0,200)
lostPhotonHist1 = rt.TH1F("Lost1", "Shower not found",30,0,200)
weirdHist1 = rt.TH1F("Weird1", "Misclassified shower",30,0,200)
unclassifiedHist1 = rt.TH1F("Unclassified1", "Unclassified shower",30,0,200)

foundPhotonHist2 = rt.TH1F("Found2", "Found Photons",30,0,200)
lostPhotonHist2 = rt.TH1F("Lost2", "Shower not found",30,0,200)
weirdHist2 = rt.TH1F("Weird2", "Misclassified shower",30,0,200)
unclassifiedHist2 = rt.TH1F("Unclassified2", "Unclassified shower",30,0,200)

foundPhotonHist3 = rt.TH1F("Found3", "Found Photons",30,0,200)
lostPhotonHist3 = rt.TH1F("Lost3", "Shower not found",30,0,200)
weirdHist3 = rt.TH1F("Weird3", "Misclassified shower",30,0,200)
unclassifiedHist3 = rt.TH1F("Unclassified3", "Unclassified shower",30,0,200)

foundPhotonHist4 = rt.TH1F("Found4", "Found Photons",30,0,200)
lostPhotonHist4 = rt.TH1F("Lost4", "Shower not found",30,0,200)
weirdHist4 = rt.TH1F("Weird4", "Misclassified shower",30,0,200)
unclassifiedHist4 = rt.TH1F("Unclassified4", "Unclassified shower",30,0,200)

distanceList1 = [lostPhotonHist1, weirdHist1, unclassifiedHist1, foundPhotonHist1]
distanceList2 = [lostPhotonHist2, weirdHist2, unclassifiedHist2, foundPhotonHist2]
distanceList3 = [lostPhotonHist3, weirdHist3, unclassifiedHist3, foundPhotonHist3]
distanceList4 = [lostPhotonHist4, weirdHist4, unclassifiedHist4, foundPhotonHist4]

distanceListList = [distanceList1, distanceList2, distanceList3, distanceList4]

foundList = [foundPhotonHist1, foundPhotonHist2, foundPhotonHist3, foundPhotonHist4]
lostList = [lostPhotonHist1, lostPhotonHist2, lostPhotonHist3, lostPhotonHist4]
weirdList = [weirdHist1, weirdHist2, weirdHist3, weirdHist4]
unclassifiedList = [unclassifiedHist1, unclassifiedHist2, unclassifiedHist3, unclassifiedHist4]

foundPhotonE1 = rt.TH1F("FoundE1", "Found Photons",30,0,1.6)
lostPhotonE1 = rt.TH1F("LostE1", "Shower not found",30,0,1.6)
weirdE1 = rt.TH1F("WeirdE1", "Misclassified shower",30,0,1.6)
unclassifiedE1 = rt.TH1F("UnclassifiedE1", "Unclassified shower",30,0,1.6)

foundPhotonE2 = rt.TH1F("FoundE2", "Found Photons",30,0,1.6)
lostPhotonE2 = rt.TH1F("LostE2", "Shower not found",30,0,1.6)
weirdE2 = rt.TH1F("WeirdE2", "Misclassified shower",30,0,1.6)
unclassifiedE2 = rt.TH1F("UnclassifiedE2", "Unclassified shower",30,0,1.6)

foundPhotonE3 = rt.TH1F("FoundE3", "Found Photons",30,0,1.6)
lostPhotonE3 = rt.TH1F("LostE3", "Shower not found",30,0,1.6)
weirdE3 = rt.TH1F("WeirdE3", "Misclassified shower",30,0,1.6)
unclassifiedE3 = rt.TH1F("UnclassifiedE3", "Unclassified shower",30,0,1.6)

foundPhotonE4 = rt.TH1F("FoundE4", "Found Photons",30,0,1.6)
lostPhotonE4 = rt.TH1F("LostE4", "Shower not found",30,0,1.6)
weirdE4 = rt.TH1F("WeirdE4", "Misclassified shower",30,0,1.6)
unclassifiedE4 = rt.TH1F("UnclassifiedE4", "Unclassified shower",30,0,1.6)

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
def addHist(distVtx, distance, energy, HistList1, HistList2, weight):
  if eventTree.vtxDistToTrue <= 1:
    HistList1[0].Fill(distance, weight)
    HistList2[0].Fill(energy, weight)
  elif eventTree.vtxDistToTrue <= 3:
    HistList1[1].Fill(distance, weight)
    HistList2[1].Fill(energy, weight)
  elif eventTree.vtxDistToTrue <= 5:
    HistList1[2].Fill(distance, weight)
    HistList2[2].Fill(energy, weight)
  elif eventTree.vtxDistToTrue > 5:
    HistList1[3].Fill(distance, weight)
    HistList2[3].Fill(energy, weight)

def addEffHist(effHist, histList):
  for x in range(1, 31):
    if histList[3].GetBinContent(x)+histList[0].GetBinContent(x)+histList[1].GetBinContent(x)+histList[2].GetBinContent(x) > 0:
      efficiency = histList[3].GetBinContent(x)/(histList[3].GetBinContent(x)+histList[0].GetBinContent(x)+histList[1].GetBinContent(x)+histList[2].GetBinContent(x))
    else:
      efficiency = 0
    effHist.SetBinContent(x, efficiency)


#Variables for program review


#Variables for program function
recoPhotonIDs = []
truePhotonIDs = []
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}
firstTime = True
firstTimeE = True

#Event loop begins
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #RECO CUTS
  if recoNoVertex(eventTree) == False:
    continue

  #See if the event is neutral current
  if recoNeutralCurrent(eventTree) == False:
    continue

  #Make sure the event is within the fiducial volume
  if recoFiducials(eventTree, fiducialData) == False:
    continue

  #Cut events with suitably energetic protons or charged pions
  if recoPionProton(eventTree) == False:
    continue
  
  #FIND TWO-PHOTON EVENTS
  #Get a list of true photons, scaled energy, and invariant mass
  truePhotonIDs = truePhotonList(eventTree, truePhotonIDs, fiducialData)
  #Remove events unless they have exactly two true photons
  if len(truePhotonIDs) != 2:
    continue

  #Calculate leading photon
  leadingPhoton, invariantMass = scaleTrueEnergy(eventTree, truePhotonIDs)

  #List photons the reco can find
  recoList = recoPhotonList(eventTree)
  #See how the algorithm treats each reco photon
  matches = 0
  distVtx = eventTree.vtxDistToTrue
  for x in truePhotonIDs:
    #Quick setup
    trueNo = x
    distancexy = np.sqrt((eventTree.trueSimPartEDepX[trueNo] - eventTree.trueVtxX)**2 + (eventTree.trueSimPartEDepY[trueNo] - eventTree.trueVtxY)**2)
    distancexyz = np.sqrt((distancexy)**2 + (eventTree.trueSimPartEDepZ[trueNo] - eventTree.trueVtxZ)**2)
    recoMatch = False
    for y in range(eventTree.nShowers):
      #Match the true event to a reco event, if we can
      if eventTree.showerTrueTID[y] == eventTree.trueSimPartTID[x]:
        recoMatch = True
        #If the track was correctly identified by the reco, we store it here
        if y in recoList:
          addHist(distVtx, distancexyz, eventTree.trueSimPartE[x]/1000, foundList, foundListE, eventTree.xsecWeight)
        #If the track wasn't classified, we store it here
        elif eventTree.showerClassified[y] == 0:
          addHist(distVtx, distancexyz, eventTree.trueSimPartE[x]/1000, unclassifiedList, unclassifiedListE, eventTree.xsecWeight)
        #If it somehow got identified, but not as an electron or a photon, we should mark that down too
        else:
          addHist(distVtx, distancexyz, eventTree.trueSimPartE[x]/1000, weirdList, weirdListE, eventTree.xsecWeight)

    if recoMatch == True:
      matches += 1
    #If we didn't find any track at all, we store it here
    else:
      addHist(distVtx, distancexyz, eventTree.trueSimPartE[x]/1000, lostList, lostListE, eventTree.xsecWeight)


#Fill efficiency histograms
for x in range(len(effList)):
  addEffHist(effList[x], distanceListList[x])
  addEffHist(effListE[x], energyListList[x])

for x in range(effHist1.GetNbinsX()):
  total = (effHist1.GetBinContent(x) + effHist2.GetBinContent(x) + effHist3.GetBinContent(x) + effHist4.GetBinContent(x))/4
  effTotal.SetBinContent(x, total)

  totalE = (effHistE1.GetBinContent(x) + effHistE2.GetBinContent(x) + effHistE3.GetBinContent(x) + effHistE4.GetBinContent(x))/4
  effTotalE.SetBinContent(x, totalE)
  
#effCanvas, effStack, effLegend, effInt = histStack("Efficiencies of photon showers over distance", effList)
#effCanvasE, effStackE, effLegendE, effIntE = histStack("Efficiencies of photon showers over energy", effListE)

effTotal.SetLineColor(rt.kRed)
effTotalE.SetLineColor(rt.kRed)

#Create canvases for efficiency histograms
effCanvas1 = rt.TCanvas("effCanvas1","effCanvas1")
effHist1.SetLineColor(rt.kGreen+2)
effHist1.Draw()
effTotal.Draw("Same")


effCanvas2 = rt.TCanvas("effCanvas2","effCanvas2")
effHist2.SetLineColor(rt.kGreen+2)
effHist2.Draw()
effTotal.Draw("Same")

effCanvas3 = rt.TCanvas("effCanvas3","effCanvas3")
effHist3.SetLineColor(rt.kGreen+2)
effHist3.Draw()
effTotal.Draw("Same")

effCanvas4 = rt.TCanvas("effCanvas4","effCanvas4")
effHist4.SetLineColor(rt.kGreen+2)
effHist4.Draw()
effTotal.Draw("Same")

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


#Create canvases for distance histograms
histCanvas1, stack1, legend1, histInt1 = histStack("Distances - Vertex within 1 cm of true vertex", distanceList1)
histCanvas2, stack2, legend2, histInt2 = histStack("Distances - Vertex 1-3 cm from true vertex", distanceList2)
histCanvas3, stack3, legend3, histInt3 = histStack("Distances - Vertex 3-5 cm from true vertex", distanceList3)
histCanvas4, stack4, legend4, histInt4 = histStack("Distances - Vertex more than 5 cm from true vertex", distanceList4)

#Create canvases for energy histograms
histCanvasE1, stackE1, legendE1, histIntE1 = histStack("Energy - Vertex within 1 cm of true vertex", energyList1)
histCanvasE2, stackE2, legendE2, histIntE2 = histStack("Energy - Vertex 1-3 cm from true vertex", energyList2)
histCanvasE3, stackE3, legendE3, histIntE3 = histStack("Energy - Vertex 3-5 cm from true vertex", energyList3)
histCanvasE4, stackE4, legendE4, histIntE4 = histStack("Energy - Vertex more than 5 cm from true vertex", energyList4)

#Make canvas list so we can export efficiently
canvasList = [effCanvas1, effCanvas2, effCanvas3, effCanvas4, effCanvasE1, effCanvasE2, effCanvasE3, effCanvasE4, histCanvas1, histCanvas2, histCanvas3, histCanvas4, histCanvasE1, histCanvasE2, histCanvasE3, histCanvasE4]

#Change axes
effHist1.GetXaxis().SetTitle("True Shower Distance to Vertex (cm)")
effHistE1.GetXaxis().SetTitle("True Energy (GeV)")
effHist1.GetYaxis().SetTitle("Efficiency")
effHistE1.GetYaxis().SetTitle("Efficiency")

stackList = [stack1, stack2, stack3, stack4]
for stack in stackList:
  stack.GetXaxis().SetTitle("True Distance to Vertex (cm)")
  stack.GetYaxis().SetTitle("Number of Photons")

stackListE = [stackE1, stackE2, stackE3, stackE4]
for stack in stackListE:
  stack.GetXaxis().SetTitle("True Photon Energy GeV")
  stack.GetYaxis().SetTitle("Number of Photons")

#Save to file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in canvasList:
  canvas.Write()
