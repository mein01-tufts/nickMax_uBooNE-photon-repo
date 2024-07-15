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
noneFound = rt.TH1F("noneFound", "None Found",60,0,1.6)
oneFound = rt.TH1F("oneFound", "One Found",60,0,1.6)
twoFound = rt.TH1F("twoFound", "Two Found",60,0,1.6)

totalList2 = [noneFound, oneFound, twoFound]

effHist1 = rt.TH1F("DistanceEfficiency", "Reconstructed Vertex < 1 cm from true",60,0,200)
effHist2 = rt.TH1F("DistanceEfficiency", "Reconstructed Vertex < 3 cm from true",60,0,200)
effHist3 = rt.TH1F("DistanceEfficiency", "Reconstructed Vertex < 5 cm from true",60,0,200)
effHist4 = rt.TH1F("DistanceEfficiency", "Reconstructed Vertex > 5 cm from true",60,0,200)

effHistE1 = rt.TH1F("DistanceEfficiency", "Reconstructed Vertex < 1 cm from true",60,0,200)
effHistE2 = rt.TH1F("DistanceEfficiency", "Reconstructed Vertex < 3 cm from true",60,0,200)
effHistE3 = rt.TH1F("DistanceEfficiency", "Reconstructed Vertex < 5 cm from true",60,0,200)
effHistE4 = rt.TH1F("DistanceEfficiency", "Reconstructed Vertex > 5 cm from true",60,0,200)

effList = [effHist1, effHist2, effHist3, effHist4]
effListE = [effHistE1, effHistE2, effHistE3, effHistE4]

foundPhotonHist1 = rt.TH1F("Found2", "Found Photons",60,0,200)
lostPhotonHist1 = rt.TH1F("Lost2", "Shower not found",60,0,200)
weirdHist1 = rt.TH1F("Weird2", "Misclassified shower",60,0,200)
unclassifiedHist1 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,200)

foundPhotonHist2 = rt.TH1F("Found2", "Found Photons",60,0,200)
lostPhotonHist2 = rt.TH1F("Lost2", "Shower not found",60,0,200)
weirdHist2 = rt.TH1F("Weird2", "Misclassified shower",60,0,200)
unclassifiedHist2 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,200)

foundPhotonHist3 = rt.TH1F("Found2", "Found Photons",60,0,200)
lostPhotonHist3 = rt.TH1F("Lost2", "Shower not found",60,0,200)
weirdHist3 = rt.TH1F("Weird2", "Misclassified shower",60,0,200)
unclassifiedHist3 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,200)

foundPhotonHist4 = rt.TH1F("Found2", "Found Photons",60,0,200)
lostPhotonHist4 = rt.TH1F("Lost2", "Shower not found",60,0,200)
weirdHist4 = rt.TH1F("Weird2", "Misclassified shower",60,0,200)
unclassifiedHist4 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,200)

distanceList1 = [lostPhotonHist1, weirdHist1, unclassifiedHist1, foundPhotonHist1]
distanceList2 = [lostPhotonHist2, weirdHist2, unclassifiedHist2, foundPhotonHist2]
distanceList3 = [lostPhotonHist3, weirdHist3, unclassifiedHist3, foundPhotonHist3]
distanceList4 = [lostPhotonHist4, weirdHist4, unclassifiedHist4, foundPhotonHist4]

distanceListList = [distanceList1, distanceList2, distanceList3, distanceList4]

foundList = [foundPhotonHist1, foundPhotonHist2, foundPhotonHist3, foundPhotonHist4]
lostList = [lostPhotonHist1, lostPhotonHist2, lostPhotonHist3, lostPhotonHist4]
weirdList = [weirdHist2, weirdHist3, weirdHist4, weirdHist1]
unclassifiedList = [unclassifiedHist1, unclassifiedHist2, unclassifiedHist3, unclassifiedHist4]

foundPhotonE1 = rt.TH1F("Found2", "Found Photons",60,0,1.6)
lostPhotonE1 = rt.TH1F("Lost2", "Shower not found",60,0,1.6)
weirdE1 = rt.TH1F("Weird2", "Misclassified shower",60,0,1.6)
unclassifiedE1 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,1.6)

foundPhotonE2 = rt.TH1F("Found2", "Found Photons",60,0,1.6)
lostPhotonE2 = rt.TH1F("Lost2", "Shower not found",60,0,1.6)
weirdE2 = rt.TH1F("Weird2", "Misclassified shower",60,0,1.6)
unclassifiedE2 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,1.6)

foundPhotonE3 = rt.TH1F("Found2", "Found Photons",60,0,1.6)
lostPhotonE3 = rt.TH1F("Lost2", "Shower not found",60,0,1.6)
weirdE3 = rt.TH1F("Weird2", "Misclassified shower",60,0,1.6)
unclassifiedE3 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,1.6)

foundPhotonE4 = rt.TH1F("Found2", "Found Photons",60,0,1.6)
lostPhotonE4 = rt.TH1F("Lost2", "Shower not found",60,0,1.6)
weirdE4 = rt.TH1F("Weird2", "Misclassified shower",60,0,1.6)
unclassifiedE4 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,1.6)

energyList1 = [lostPhotonE1, weirdE1, unclassifiedE1, foundPhotonE1]
energyList2 = [lostPhotonE2, weirdE2, unclassifiedE2, foundPhotonE2]
energyList3 = [lostPhotonE3, weirdE3, unclassifiedE3, foundPhotonE3]
energyList4 = [lostPhotonE4, weirdE4, unclassifiedE4, foundPhotonE4]

energyListList = [energyList1, energyList2, energyList3, energyList4]

foundListE = [foundPhotonE1, foundPhotonE2, foundPhotonE3, foundPhotonE4]
lostListE = [lostPhotonE1, lostPhotonE2, lostPhotonE3, lostPhotonE4]
weirdListE = [weirdE2, weirdE3, weirdE4, weirdE1]
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
  else:
    HistList1[3].Fill(distance, weight)
    HistList2[3].Fill(energy, weight)

def addEffHist(effHist, histList):
  for x in range(1, 61):
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
    #Now we store the event in one of the total histograms depending on how many of its photons we found
#  if matches == 0:
#    noneFound.Fill(invariantMass, eventTree.xsecWeight)
#  elif matches == 1:
#    oneFound.Fill(invariantMass, eventTree.xsecWeight)
#  elif matches == 2:
#    twoFound.Fill(invariantMass, eventTree.xsecWeight)
#  else:
#    print("Something wrong with the number of matches here in event", i)

#Fill efficiency histograms
for x in range(len(effList)):
  addEffHist(effList[x], distanceListList[x])
  addEffHist(effListE[x], energyListList[x])

effCanvas, effStack, effLegend, effInt = histStack("Efficiencies of photon showers over distance", effList)
effCanvasE, effStackE, effLegendE, effIntE = histStack("Efficiencies of photon showers over energy", effListE)

histCanvas1, stack1, legend1, histInt1 = histStack("Distances of Showers in Two Gamma + 0 True Events - reconstructed vertex within 1 cm of true vertex", distanceList1)
histCanvas2, stack2, legend2, histInt2 = histStack("Distances of Showers in Two Gamma + 0 True Events - reconstructed vertex within 3 cm of true vertex", distanceList2)
histCanvas3, stack3, legend3, histInt3 = histStack("Distances of Showers in Two Gamma + 0 True Events - reconstructed vertex qithin 5 cm of true vertex", distanceList3)
histCanvas4, stack4, legend4, histInt4 = histStack("Distances of Showers in Two Gamma + 0 True Events - reconstructed vertex more than 5 cm from true vertex", distanceList4)

histCanvasE1, stackE1, legendE1, histIntE1 = histStack("Energy of Showers in Two Gamma + 0 True Events - reconstructed vertex within 1 cm of true vertex", energyList1)
histCanvasE2, stackE2, legendE2, histIntE2 = histStack("Energy of Showers in Two Gamma + 0 True Events - reconstructed vertex within 3 cm of true vertex", energyList2)
histCanvasE3, stackE3, legendE3, histIntE3 = histStack("Energy of Showers in Two Gamma + 0 True Events - reconstructed vertex within 5 cm of true vertex", energyList3)
histCanvasE4, stackE4, legendE4, histIntE4 = histStack("Energy of Showers in Two Gamma + 0 True Events - reconstructed vertex more than 5 cm from true vertex", energyList4)

canvasList = [effCanvas, effCanvasE, histCanvas1, histCanvas2, histCanvas3, histCanvas4, histCanvasE1, histCanvasE2, histCanvasE3, histCanvasE4]

effStack.GetXaxis().SetTitle("True Shower Distance to Vertex (cm)")
effStackE.GetXaxis().SetTitle("True Energy (GeV)")
effStack.GetYaxis().SetTitle("Efficiency")
effStackE.GetYaxis().SetTitle("Efficiency")

stackList = [stack1, stack2, stack3, stack4]
for stack in stackList:
  stack.GetXaxis().SetTitle("True Distance to Vertex (cm)")
  stack.GetYaxis().SetTitle("Number of Photons")

stackListE = [stackE1, stackE2, stackE3, stackE4]
for stack in stackListE:
  stack.GetXaxis().SetTitle("True Photon Energy GeV")
  stack.GetYaxis().SetTitle("Number of Photons")
  
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in canvasList:
  canvas.Write()
