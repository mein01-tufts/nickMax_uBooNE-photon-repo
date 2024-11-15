import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials, trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, recoCutLowEnergy

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
parser.add_argument("-c", "--cosmicFile", required = True, help = "input cosmic background file")

args = parser.parse_args()

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


def addHist(recoList, leadingPhoton, Hist1, Hist2, Hist3, weight):
  #Adds to one of three histograms depending on reco photon
  if len(recoList) == 1:
    Hist1.Fill(leadingPhoton, weight)
  elif len(recoList) == 2:
    Hist2.Fill(invariantMass, weight)
  else:
    Hist3.Fill(leadingPhoton, weight)

  return Hist1, Hist2, Hist3

  
#HISTOGRAMS DEFINED AND PREPARED HERE
totalHist1 = rt.TH1F("total1", "One photon",60,0,2)
totalHist2 = rt.TH1F("total2", "Two photons",60,0,2)
totalHist3 = rt.TH1F("total3", "Three photons",60,0,2)

totalList = [totalHist1, totalHist2, totalHist3]

muonHist = rt.TH1F("Muon1", "Charged Current - Missed Muon",60,0,1.2)
electronHist = rt.TH1F("Electron1", "Electron flagged as photon",60,0,1.2)
chargeCurrentHist = rt.TH1F("CC1", "Other Charged Current",60,0,1.2)
fiducialHist = rt.TH1F("fiducial1", "Out of Fiducial",60,0,1.2)
cosmicHist = rt.TH1F("cosmic1", "Cosmic Overlap",60,0,1.2)
pionProtonHist = rt.TH1F("pionProton1", "Charged Pion or Proton",60,0,1.2)
noPhotonHist = rt.TH1F("noPhotons1", "No Photons",60,0,1.2)
twoPhotonHist = rt.TH1F("twoPhotons1", "Two Photons",60,0,1.2)
manyPhotonHist = rt.TH1F("manyPhotons1", "3+ Photons",60,0,1.2)
successHist = rt.TH1F("success1", "Signal",60,0,1.2)

histList = [successHist, electronHist, muonHist, chargeCurrentHist, fiducialHist, cosmicHist, pionProtonHist, noPhotonHist, twoPhotonHist, manyPhotonHist]

muonHist2 = rt.TH1F("Muon2", "Charged Current - Missed Muon",60,0,0.8)
electronHist2 = rt.TH1F("Electron2", "Electron flagged as photon",60,0,0.8)
chargeCurrentHist2 = rt.TH1F("CC2", "Other Charged current",60,0,0.8)
fiducialHist2 = rt.TH1F("fiducial2", "Out of fiducial",60,0,0.8)
cosmicHist2 = rt.TH1F("cosmic2", "Cosmic overlap",60,0,0.8)
pionProtonHist2 = rt.TH1F("pionProton2", "Charged Pion or Proton",60,0,0.8)
noPhotonHist2 = rt.TH1F("noPhotons2", "No Photons",60,0,0.8)
onePhotonHist2 = rt.TH1F("twoPhotons2", "One Photon",60,0,0.8)
manyPhotonHist2 = rt.TH1F("manyPhotons2", "3+ Photons",60,0,0.8)
successHist2 = rt.TH1F("success2", "Signal",60,0,0.8)

histList2 = [successHist2, electronHist2, muonHist2, chargeCurrentHist2, fiducialHist2, cosmicHist2, pionProtonHist2, noPhotonHist2, onePhotonHist2, manyPhotonHist2]

muonHist3 = rt.TH1F("Muon3", "Charged Current - Missed Muon",60,0,2)
electronHist3 = rt.TH1F("Electron3", "Electron flagged as photon",60,0,2)
chargeCurrentHist3 = rt.TH1F("CC3", "Other Charged current",60,0,2)
fiducialHist3 = rt.TH1F("fiducial3", "Out of fiducial",60,0,2)
cosmicHist3 = rt.TH1F("cosmic3", "Cosmic overlap",60,0,2)
pionProtonHist3 = rt.TH1F("pionProton3", "Charged Pion or Proton",60,0,2)
noPhotonHist3 = rt.TH1F("noPhotons3", "No Photons",60,0,2)
twoPhotonHist3 = rt.TH1F("twoPhotons3", "Two Photons",60,0,2)
onePhotonHist3 = rt.TH1F("onePhotons", "One Photon",60,0,2)
successHist3 = rt.TH1F("success3", "Signal",60,0,2)

histList3 = [successHist3, electronHist3, muonHist3, chargeCurrentHist3, fiducialHist3, cosmicHist3, pionProtonHist3, noPhotonHist3, twoPhotonHist3, onePhotonHist3]

#Variables for program review
recoCount = 0

#Variables for program function
recoPhotonIDs = []
truePhotonIDs = []
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}

#Event loop begins
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #CUT WITH RECO
  #See if the event has a vertex
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

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonList(eventTree)

  #This thing reduces cosmic background by ~85%
#  if recoCutLowEnergy(recoList, eventTree) == False:
#    continue
  
  if len(recoList) == 0:
    continue
  else:
    recoCount += 1

  #SORT USING TRUTH
  #Scale energy, find leading photon
  leadingPhoton, invariantMass = scaleRecoEnergy(eventTree, recoList)
  
  
  #Fill total histogram
  addHist(recoList, leadingPhoton, totalHist1, totalHist2, totalHist3, eventTree.xsecWeight)
  
  #Fill for charge current
  if trueCutNC(eventTree) == False:
    added = False
    for x in recoList:
      if eventTree.showerTruePID[x] == 11:
        addHist(recoList, leadingPhoton, electronHist, electronHist2, electronHist3, eventTree.xsecWeight)
        added = True
      elif 13 in eventTree.truePrimPartPDG:
        addHist(recoList, leadingPhoton, muonHist, muonHist2, muonHist3, eventTree.xsecWeight)
        added = True
    if added == False:
      addHist(recoList, leadingPhoton, chargeCurrentHist, chargeCurrentHist2, chargeCurrentHist3, eventTree.xsecWeight)
    continue

  #Fill for fiducial failures
  if trueCutFiducials(eventTree, fiducialData) == False:
    addHist(recoList, leadingPhoton, fiducialHist, fiducialHist2, fiducialHist3, eventTree.xsecWeight)
    continue
  
  #Fill for cosmic failures
  if trueCutCosmic(eventTree) == False:
    addHist(recoList, leadingPhoton, cosmicHist, cosmicHist2, cosmicHist3, eventTree.xsecWeight)
    continue

  #Fill for pions and protons slipping in
  if trueCutPionProton(eventTree) == False:
    addHist(recoList, leadingPhoton, pionProtonHist, pionProtonHist2, pionProtonHist3, eventTree.xsecWeight)
    continue

  #Find out the actual number of photons, and sort according to whether we got it right
  truePhotonIDs = truePhotonList(eventTree, fiducialData)
  if len(truePhotonIDs) == 0:
    addHist(recoList, leadingPhoton, noPhotonHist, noPhotonHist2, noPhotonHist3, eventTree.xsecWeight)

  elif len(truePhotonIDs) == 1:
    addHist(recoList, leadingPhoton, successHist, onePhotonHist2, onePhotonHist3, eventTree.xsecWeight)

  elif len(truePhotonIDs) == 2:
    addHist(recoList, leadingPhoton, twoPhotonHist, successHist2, twoPhotonHist3, eventTree.xsecWeight)
  else: 
    addHist(recoList, leadingPhoton, manyPhotonHist, manyPhotonHist2, successHist3, eventTree.xsecWeight)
  
#END OF LOOP
#Stack Histograms
totalHistCanvas, stack0, legend0, histInt0 = histStack("Chart of all reconstructed events", totalList, ntuplePOTsum)
histCanvas, stack, legend, histInt = histStack("Outcome of single-photon reco events", histList, ntuplePOTsum)
histCanvas2, stack2, legend2, histInt2 = histStack("Outcome of two-photon reco events", histList2, ntuplePOTsum)
histCanvas3, stack3, legend3, histInt3 = histStack("Outcome of three-photon reco events", histList3, ntuplePOTsum)

stack.GetXaxis().SetTitle("Reconstructed Leading Photon Energy (GeV)")
stack2.GetXaxis().SetTitle("Reconstructed Invariant Mass (GeV)")

#


#Export canvases to file
outFile = rt.TFile(args.outfile, "RECREATE")
totalHistCanvas.Write()
histCanvas.Write()
histCanvas2.Write()
histCanvas3.Write()
