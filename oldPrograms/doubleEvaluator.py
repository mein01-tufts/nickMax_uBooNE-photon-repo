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

cosmicPOTsum = 3.2974607516739725e+20

#Histograms
EDepHist = rt.TH1F("EDep", "Successful deposit in Fiducial", 60, 0, 1036)
EDepFailHist = rt.TH1F("EDepFail", "Failed to deposit in Fiducial", 60, 0, 1036)
EDepList = [EDepHist, EDepFailHist]

EDepZHist = rt.TH2F("EDepZHist", "No. of Pixels (y) vs EDep z-coordinate", 60, 0, 10000, 60, 0, 1036)

#Variables for program review

#Variables for program function
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":30}

#BEGINNING EVENT LOOP FOR DEFAULT PURITY
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #Selecting events with reco
  #See if the event has a vertex
  if recoNoVertex(eventTree) == False:
    continue

  #See if the event is neutral current
  if recoNeutralCurrent(eventTree) == False:
    continue

  #Use Matt's Cosmic Cut
  if trueCutCosmic(eventTree) == False:
    continue 

  #Make sure the event is within the fiducial volume
  if recoFiducials(eventTree, fiducialData) == False:
    continue

  #Cut events with suitably energetic protons or charged pions 
  if recoPion(eventTree) == False:
    continue

  for x in range(eventTree.nTrueSimParts):
    if eventTree.trueSimPartPDG[x] == 22:
      EDepZHist.Fill(eventTree.trueSimPartPixelSumYplane[x], eventTree.trueSimPartEDepZ[x], eventTree.xsecWeight)

#canvas1, stack1, legend1, histInt1 = histStack("Energy Deposition based on Z-axis position", EDepList, ntuplePOTsum)

#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
EDepZHist.Write()