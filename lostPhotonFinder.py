import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

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

foundPhotonDist2 = rt.TH1F("Found2", "Found Photons",60,0,200)
lostPhotonDist2 = rt.TH1F("Lost2", "Shower not found",60,0,200)
electronHist2 = rt.TH1F("Electron2", "Misclassified as electron",60,0,200)
weirdHist2 = rt.TH1F("Weird2", "Misclassified, not electron",60,0,200)
unclassifiedHist2 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,200)

distanceList2 = [lostPhotonDist2, electronHist2, weirdHist2, unclassifiedHist2, foundPhotonDist2]

foundPhotonE2 = rt.TH1F("Found2", "Found Photons",60,0,1.6)
lostPhotonE2 = rt.TH1F("Lost2", "Shower not found",60,0,1.6)
electronE2 = rt.TH1F("Electron2", "Misclassified as electron",60,0,1.6)
weirdE2 = rt.TH1F("Weird2", "Misclassified, not electron",60,0,1.6)
unclassifiedE2 = rt.TH1F("Unclassified2", "Unclassified shower",60,0,1.6)

energyList2 = [lostPhotonE2, electronE2, weirdE2, unclassifiedE2, foundPhotonE2]

#Functions for efficiency
def addHist(distance, energy, Hist1, Hist2, weight):
  Hist1.Fill(distance, weight)
  Hist2.Fill(energy, weight)

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
          addHist(distancexyz, eventTree.trueSimPartE[x]/1000, foundPhotonDist2, foundPhotonE2, eventTree.xsecWeight)

        #If the track wasn't classified, we store it here
        elif eventTree.showerClassified[y] == 0:
          addHist(distancexyz, eventTree.trueSimPartE[x]/1000, unclassifiedHist2, unclassifiedE2, eventTree.xsecWeight)

        #If the track was mischaracterized as an electron, we store it here
        elif eventTree.showerPID[y] == 11:
          addHist(distancexyz, eventTree.trueSimPartE[x]/1000, electronHist2, electronE2, eventTree.xsecWeight)

        #If it somehow got identified, but not as an electron or a photon, we should mark that down too
        else:
          addHist(distancexyz, eventTree.trueSimPartE[x]/1000, weirdHist2, weirdE2, eventTree.xsecWeight)

        #Whatever the case, we found the corresponding track, so we can end the loop now
        break
    if recoMatch == True:
      matches += 1
    #If we didn't find any track at all, we store it here
    else:
      addHist(distancexyz, eventTree.trueSimPartE[x]/1000, lostPhotonDist2, lostPhotonE2, eventTree.xsecWeight)
    
    #Now we store the event in one of the total histograms depending on how many of its photons we found
    if matches == 0:
      noneFound.Fill(invariantMass, eventTree.xsecWeight)
    elif matches == 1:
      oneFound.Fill(invariantMass, eventTree.xsecWeight)
    elif matches == 2:
      twoFound.Fill(invariantMass, eventTree.xsecWeight)
    else:
      print("Something wrong with the number of matches here in event", i)
        


histCanvas1, stack1, legend1, histInt1 = histStack("Distances of all Photons in Two-Photon True Events", distanceList2)
histCanvas2, stack2, legend2, histInt2 = histStack("Energy of all Photons in Two-Photon True Events", energyList2)
histCanvas3, stack3, legend3, histInt3 = histStack("Invariant Mass of All Events", totalList2)


stack1.GetXaxis().SetTitle("True Distance to Vertex (cm)")
stack3.GetXaxis().SetTitle("True Invariant Mass")
outFile = rt.TFile(args.outfile, "RECREATE")
histCanvas3.Write()
histCanvas1.Write()
histCanvas2.Write()
