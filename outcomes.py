import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials,trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoPion, recoProton, recoCutElectronScore, recoCutShowerFromChargeScore, recoCutLongTracks, recoPhotonListFiducial, recoCutPrimary, recoCutShortTracks, recoPhotonListTracks, recoCutFarShowers

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")

args = parser.parse_args()

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20"
ntuplePOTsum = 0

for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

#Hists created and organized here
foundHist = rt.TH1F("foundHist", "Found Photon",60,0,600)
unclassifiedHist = rt.TH1F("unclassifiedHist", "Unclassified Shower",60,0,600)
foundTrackHist = rt.TH1F("foundTrackHist", "Reconstructed as Track",60,0,600)
unclassifiedTrackHist = rt.TH1F("unclassifiedTrackHist", "Reconstructed as Unclassified Track",60,0,600)
lostHist = rt.TH1F("lostHist", "Lost Photon",60,0,600)

histList1 = [foundHist, unclassifiedHist, foundTrackHist, unclassifiedTrackHist, lostHist]

recoEff = rt.TH1F("recoEff", "Efficiency of Reconstruction",60,0,600)

#Built-in functions here
def addHist(eventTree, photonList, photonList2, histList, variable, weight):
  if len(photonList) + len(photonList2) == 1:
    histList[0].Fill(variable, weight)
  elif len(photonList) + len(photonList2) == 2:
    histList[1].Fill(variable, weight)
  else:
    histList[2].Fill(variable, weight)


#Variables for program review
showerMatch = 0
trackMatch = 0
rightCount = 0

#We put this into the addHist function for truth-based graphs
emptyList = []

#Variables for program function
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":30}


#BEGINNING EVENT LOOP FOR EFFICIENCY
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #Selecting events using truth
  if trueCutNC(eventTree) == False:
    continue

  if trueCutFiducials(eventTree, fiducialData) == False:
    continue
    
  if trueCutCosmic(eventTree) == False:
    continue
    
  if trueCutPionProton(eventTree) == False:
    continue
  
  truePhotonIDs = truePhotonList(eventTree, fiducialData)

  if len(truePhotonIDs) == 0:
    continue
  
  #GRAPHING BASED ON RECO

  #Iterate over photons to see which ones match to showers
  for x in truePhotonIDs:
    #Calculate the sum of the pixels; this is our new variables, which will hopefully let us determine what showers are actually reconstructable
    energySum = eventTree.trueSimPartPixelSumYplane[x]*0.0126
    matched = False
    for y in range(eventTree.nShowers):
      if eventTree.showerTrueTID[y] == eventTree.trueSimPartTID[x]:
        showerMatch += 1
        matched = True
        showerID = y
        break
    if matched == True:
      #Graph based on whether or not it's classified
      if eventTree.showerClassified[showerID] == 1:
        foundHist.Fill(energySum, eventTree.xsecWeight)
      else:
        unclassifiedHist.Fill(energySum, eventTree.xsecWeight)
    #If we can't get a match, we continue
    else:
      #Check to see if it reconstructed as a track for some reason
      for y in range(eventTree.nTracks):
        trackMatched = False
        if eventTree.trackTrueTID[y] == eventTree.trueSimPartTID[x]:
          trackMatch += 1
          trackMatched = True
          if eventTree.trackTruePID[y] == 22:
            rightCount += 1
          else:
            print("Track TID:", eventTree.trackTrueTID[y], "Particle PDG:", eventTree.trueSimPartPDG[x])
          trackID = y
          break
        if trackMatched == True:
          if eventTree.trackClassified[trackID] == 1:
            foundTrackHist.Fill(energySum, eventTree.xsecWeight)
          else:
            unclassifiedTrackHist.Fill(energySum, eventTree.xsecWeight)
            print("Event Number", i, "went to the unclassified track list")
        else:
          lostHist.Fill(energySum, eventTree.xsecWeight)

for x in range(1, 61):
  denominator = histList1[0].GetBinContent(x) + histList1[1].GetBinContent(x) + histList1[2].GetBinContent(x) + histList1[3].GetBinContent(x) + histList1[4].GetBinContent(x)
  if denominator > 0:
    efficiency = histList1[0].GetBinContent(x)/denominator
  else:
    efficiency = 0
  recoEff.SetBinContent(x, efficiency)

#Stacking histograms
Canvas1, Stack1, Legend1, Int1 = histStack("Photon Reconstruction", histList1, ntuplePOTsum)

Stack1.GetXaxis().SetTitle("Total pixels")
recoEff.GetXaxis().SetTitle("Total pixels")
recoEff.GetYaxis().SetTitle("Efficiency")

#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
Canvas1.Write()
recoEff.Write()

print(showerMatch, "particles matched with showers.")
print(trackMatch, "particles matched with tracks. It was correct in", rightCount, "cases.")