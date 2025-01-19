import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials,trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoPion, recoProton, recoCutShowerFromChargeScore, recoCutLongTracks, recoPhotonListFiducial, recoCutPrimary, recoCutShortTracks, recoPhotonListTracks, recoCutFarShowers, trueCutMuons, trueCutElectrons, recoCutMuons, recoCutElectrons, recoCutManyTracks, recoCutCompleteness, recoCutMuonCompleteness, histStackTwoSignal, histStackData

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
#targetPOT = 6.67e+20
#targetPOTstring = "6.67e+20"
ntuplePOTsum = 4.4e19

#Hists created and organized here
#PURITY HISTOGRAMS

#1 GAMMA + 0 HISTS
noProtonSignal = rt.TH1F("1 Gamma + 0 Events", "1 Gamma + 0 Events",60,0,2)

#1 GAMMA + 1P HISTS
protonSignal = rt.TH1F("1 Gamma + 1P Events", "1 Gamma + 1P Events",60,0,2)

#HISTLIST
dataHists = [noProtonSignal, protonSignal]

#Built-in functions here
def addHist(ntuple, protonNumber, histList, variable, weight):
  if protonNumber == 0:
    histList[0].Fill(variable, weight)
  elif protonNumber == 1:
    histList[1].Fill(variable, weight)


def histScale(hist, POTSum):
  POTTarget = 6.67e+20
  hist.Scale(POTTarget/POTSum)
  return hist

#Variables for program review

#Variables for program function
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":30}
classificationThreshold = 0

#BEGINNING EVENT LOOP FOR DEFAULT PURITY
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #Selecting events with reco
  #See if the event has a vertex
  if recoNoVertex(eventTree) == False:
    continue

  #See if the event is neutral current
  if recoCutMuons(eventTree, classificationThreshold) == False:
    continue

  if recoCutElectrons(eventTree, classificationThreshold) == False:
    continue

  #Use Matt's Cosmic Cut
  if trueCutCosmic(eventTree) == False:
    continue

  #Make sure the event is within the fiducial volume
  if recoFiducials(eventTree, fiducialData) == False:
    continue

  #Cut events with suitably energetic charged pions 
  if recoPion(eventTree, classificationThreshold) == False:
    continue

  #Cut events with far-travelling protons
  recoProtonCount = recoProton(eventTree, classificationThreshold)
  if recoProtonCount > 1:
    continue

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, eventTree, classificationThreshold)

  #List all photons classified as tracks
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree, classificationThreshold)

  if len(recoList) + len(recoTrackList) != 1:
    continue

  #Try cutting based on data for Shower from Charged
  if recoCutShowerFromChargeScore(eventTree, recoList, recoTrackList) == False:
    continue
  
  #Cut based on primary score
  if recoCutPrimary(eventTree, recoList, recoTrackList) == False:
    continue

  #Cut based on the presence of tracks over 20 cm
  if recoCutLongTracks(eventTree, fiducialData) == False:
    continue

  #Cut based on the completeness of known Muons
  if recoCutMuonCompleteness(eventTree) == False:
    continue

  #We made it! Time to graph
  leadingPhoton = scaleRecoEnergy(eventTree, recoList, recoTrackList)

  addHist(eventTree, recoProtonCount, dataHists, leadingPhoton, 1)


#Stacking histograms
dataCanvas1, dataStack1, dataLegend1, dataInt1 = histStackTwoSignal("All events in the 1 Gamma + 0 and 1 Gamma + 1P Analysis", dataHists, ntuplePOTsum)

writeList = [dataCanvas1]

for stack in [dataStack1]:
  stack.GetXaxis().SetTitle("Reconstructed Leading Photon Energy (GeV)")

for stack in [dataStack1]:
  stack.GetYaxis().SetTitle("Events per 4.4e19 POT")

#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in writeList:
  canvas.Write()