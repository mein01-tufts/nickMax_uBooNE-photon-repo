import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials,trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoPion, recoProton, recoCutShowerFromChargeScore, recoCutLongTracks, recoPhotonListFiducial, recoCutPrimary, recoCutShortTracks, recoPhotonListTracks, recoCutFarShowers, trueCutMuons, trueCutElectrons, recoCutMuons, recoCutElectrons, recoCutManyTracks

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
#parser.add_argument("-c", "--cosmicFile", type=str, required=True, help="cosmic input file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")

args = parser.parse_args()

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

#cosmic_file = rt.TFile(args.cosmicFile)
#cosmicTree = cosmic_file.Get("EventTree")
#cosmicPotTree = ntuple_file.Get("cosmicPotTree")

#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20"
ntuplePOTsum = 0

for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

cosmicPOTsum = 3.2974607516739725e+20

#Hists created and organized here
#TRACK HISTOGRAMS
#Muons:
histMuonPhScore = rt.TH1F("MuonPhScore", "Track Photon Score",30,0,10)
histMuonElScore = rt.TH1F("MuonElScore", "Track Electron Score",30,0,10)
histMuonMuScore = rt.TH1F("MuonMuScore", "Track Muon Score",30,0,10)
histMuonPiScore = rt.TH1F("MuonPiScore", "Track Pion Score",30,0,10)
histMuonPrScore = rt.TH1F("MuonPrScore", "Track Proton Score",30,0,10)
histMuonCompPurity = rt.TH2F("MuonCompPurity", "Correctly Identified Muon Comp vs. Purity",30,0,1,30,0,1)

histDeadlyMuonPhScore = rt.TH1F("DeadlyMuonPhScore", "Track Photon Score",30,0,10)
histDeadlyMuonElScore = rt.TH1F("DeadlyMuonElScore", "Track Electron Score",30,0,10)
histDeadlyMuonMuScore = rt.TH1F("DeadlyMuonMuScore", "Track Muon Score",30,0,10)
histDeadlyMuonPiScore = rt.TH1F("DeadlyMuonPiScore", "Track Pion Score",30,0,10)
histDeadlyMuonPrScore = rt.TH1F("DeadlyMuonPrScore", "Track Proton Score",30,0,10)
histDeadlyMuonCompPurity = rt.TH2F("DeadlyMuonCompPurity", "Misidentified, Over-threshold Muon Comp vs. Purity",30,0,1,30,0,1)

histUnderMuonPhScore = rt.TH1F("UnderMuonPhScore", "Track Photon Score",30,0,10)
histUnderMuonElScore = rt.TH1F("UnderMuonElScore", "Track Electron Score",30,0,10)
histUnderMuonMuScore = rt.TH1F("UnderMuonMuScore", "Track Muon Score",30,0,10)
histUnderMuonPiScore = rt.TH1F("UnderMuonPiScore", "Track Pion Score",30,0,10)
histUnderMuonPrScore = rt.TH1F("UnderMuonPrScore", "Track Proton Score",30,0,10)
histUnderMuonCompPurity = rt.TH2F("UnderMuonCompPurity", "Correctly Identified, Over-threshold Muon Comp vs. Purity",30,0,1,30,0,1)

histBenignMuonPhScore = rt.TH1F("NonMuonPhScore", "Track Photon Score",30,0,10)
histBenignMuonElScore = rt.TH1F("NonMuonElScore", "Track Electron Score",30,0,10)
histBenignMuonMuScore = rt.TH1F("NonMuonMuScore", "Track Muon Score",30,0,10)
histBenignMuonPiScore = rt.TH1F("NonMuonPiScore", "Track Pion Score",30,0,10)
histBenignMuonPrScore = rt.TH1F("NonMuonPrScore", "Track Proton Score",30,0,10)
histBenignMuonCompPurity = rt.TH2F("MuonNonCompPurity", "Incorrectly Identified, Under-Threshold Muon Comp vs. Purity",30,0,1,30,0,1)
#Pions:
histPionPhScore = rt.TH1F("PionPhScore", "Track Photon Score",30,0,10)
histPionElScore = rt.TH1F("PionElScore", "Track Electron Score",30,0,10)
histPionMuScore = rt.TH1F("PionMuScore", "Track Muon Score",30,0,10)
histPionPiScore = rt.TH1F("PionPiScore", "Track Pion Score",30,0,10)
histPionPrScore = rt.TH1F("PionPrScore", "Track Proton Score",30,0,10)
histPionCompPurity = rt.TH2F("PionCompPurity", "Correctly Identified Pion Comp vs. Purity",30,0,1,30,0,1)

histDeadlyPionPhScore = rt.TH1F("DeadlyPionPhScore", "Track Photon Score",30,0,10)
histDeadlyPionElScore = rt.TH1F("DeadlyPionElScore", "Track Electron Score",30,0,10)
histDeadlyPionMuScore = rt.TH1F("DeadlyPionMuScore", "Track Muon Score",30,0,10)
histDeadlyPionPiScore = rt.TH1F("DeadlyPionPiScore", "Track Pion Score",30,0,10)
histDeadlyPionPrScore = rt.TH1F("DeadlyPionPrScore", "Track Proton Score",30,0,10)
histDeadlyPionCompPurity = rt.TH2F("DeadlyPionCompPurity", "Incorrectly Identified, Over-Threshold Pion Comp vs. Purity",30,0,1,30,0,1)

histUnderPionPhScore = rt.TH1F("UnderPionPhScore", "Track Photon Score",30,0,10)
histUnderPionElScore = rt.TH1F("UnderPionElScore", "Track Electron Score",30,0,10)
histUnderPionMuScore = rt.TH1F("UnderPionMuScore", "Track Muon Score",30,0,10)
histUnderPionPiScore = rt.TH1F("UnderPionPiScore", "Track Pion Score",30,0,10)
histUnderPionPrScore = rt.TH1F("UnderPionPrScore", "Track Proton Score",30,0,10)
histUnderPionCompPurity = rt.TH2F("UnderPionCompPurity", "Correctly Identified, Over-Threshold Pion Comp vs. Purity",30,0,1,30,0,1)

histBenignPionPhScore = rt.TH1F("NonPionPhScore", "Track Photon Score",30,0,10)
histBenignPionElScore = rt.TH1F("NonPionElScore", "Track Electron Score",30,0,10)
histBenignPionMuScore = rt.TH1F("NonPionMuScore", "Track Muon Score",30,0,10)
histBenignPionPiScore = rt.TH1F("NonPionPiScore", "Track Pion Score",30,0,10)
histBenignPionPrScore = rt.TH1F("NonPionPrScore", "Track Proton Score",30,0,10)
histBenignPionCompPurity = rt.TH2F("PionNonCompPurity", "Incorrectly Identified, Low-Energy Pion Comp vs. Purity",30,0,1,30,0,1)
#Protons:
histProtonPhScore = rt.TH1F("ProtonPhScore", "Track Photon Score",30,0,10)
histProtonElScore = rt.TH1F("ProtonElScore", "Track Electron Score",30,0,10)
histProtonMuScore = rt.TH1F("ProtonMuScore", "Track Muon Score",30,0,10)
histProtonPiScore = rt.TH1F("ProtonPiScore", "Track Pion Score",30,0,10)
histProtonPrScore = rt.TH1F("ProtonPrScore", "Track Proton Score",30,0,10)
histProtonCompPurity = rt.TH2F("ProtonCompPurity", "Correctly Identified, Under Threshold Proton Comp vs. Purity",30,0,1,30,0,1)

histUnderProtonPhScore = rt.TH1F("UnderProtonPhScore", "Track Photon Score",30,0,10)
histUnderProtonElScore = rt.TH1F("UnderProtonElScore", "Track Electron Score",30,0,10)
histUnderProtonMuScore = rt.TH1F("UnderProtonMuScore", "Track Muon Score",30,0,10)
histUnderProtonPiScore = rt.TH1F("UnderProtonPiScore", "Track Pion Score",30,0,10)
histUnderProtonPrScore = rt.TH1F("UnderProtonPrScore", "Track Proton Score",30,0,10)
histUnderProtonCompPurity = rt.TH2F("UnderProtonCompPurity", "Correctly Identified, Over-Threshold Proton Comp vs. Purity",30,0,1,30,0,1)

histDeadlyProtonPhScore = rt.TH1F("DeadlyProtonPhScore", "Track Photon Score",30,0,10)
histDeadlyProtonElScore = rt.TH1F("DeadlyProtonElScore", "Track Electron Score",30,0,10)
histDeadlyProtonMuScore = rt.TH1F("DeadlyProtonMuScore", "Track Muon Score",30,0,10)
histDeadlyProtonPiScore = rt.TH1F("DeadlyProtonPiScore", "Track Pion Score",30,0,10)
histDeadlyProtonPrScore = rt.TH1F("DeadlyProtonPrScore", "Track Proton Score",30,0,10)
histDeadlyProtonCompPurity = rt.TH2F("DeadlyProtonCompPurity", "Incorrectly Identified, Over-Threshold Proton Comp vs. Purity",30,0,1,30,0,1)

histBenignProtonPhScore = rt.TH1F("BenignProtonPhScore", "Track Photon Score",30,0,10)
histBenignProtonElScore = rt.TH1F("BenignProtonElScore", "Track Electron Score",30,0,10)
histBenignProtonMuScore = rt.TH1F("BenignProtonMuScore", "Track Muon Score",30,0,10)
histBenignProtonPiScore = rt.TH1F("BenignProtonPiScore", "Track Pion Score",30,0,10)
histBenignProtonPrScore = rt.TH1F("BenignProtonPrScore", "Track Proton Score",30,0,10)
histBenignProtonCompPurity = rt.TH2F("BenignProtonCompPurity", "Incorrectly Identified, Under-Threshold Proton Comp vs. Purity",30,0,1,30,0,1)

histNotTrackPhScore = rt.TH1F("nothingTrackPhScore", "Track Photon Score",30,0,10)
histNotTrackElScore = rt.TH1F("nothingTrackElScore", "Track Electron Score",30,0,10)
histNotTrackMuScore = rt.TH1F("nothingTrackMuScore", "Track Muon Score",30,0,10)
histNotTrackPiScore = rt.TH1F("nothingTrackPiScore", "Track Pion Score",30,0,10)
histNotTrackPrScore = rt.TH1F("nothingTruckPrScore", "Track Proton Score",30,0,10)
histNotTrackCompPurity = rt.TH2F("NoTrackCompPurity", "Fake Track Comp vs. Purity",30,0,1,30,0,1)



#TRACK HISTLISTS
#Muons:
MuonHistList = [histMuonPhScore, histMuonElScore, histMuonMuScore, histMuonPiScore, histMuonPrScore]
UnderMuonHistList = [histUnderMuonPhScore, histUnderMuonElScore, histUnderMuonMuScore, histUnderMuonPiScore, histUnderMuonPrScore]
DeadlyMuonHistList = [histDeadlyMuonPhScore, histDeadlyMuonElScore, histDeadlyMuonMuScore, histDeadlyMuonPiScore, histDeadlyMuonPrScore]
BenignMuonHistList = [histBenignMuonPhScore, histBenignMuonElScore, histBenignMuonMuScore, histBenignMuonPiScore, histBenignMuonPrScore]
#Pions:
PionHistList = [histPionPhScore, histPionElScore, histPionMuScore, histPionPiScore, histPionPrScore]
UnderPionHistList = [histUnderPionPhScore, histUnderPionElScore, histUnderPionMuScore, histUnderPionPiScore, histUnderPionPrScore]
DeadlyPionHistList = [histDeadlyPionPhScore, histDeadlyPionElScore, histDeadlyPionMuScore, histDeadlyPionPiScore, histDeadlyPionPrScore]
BenignPionHistList = [histBenignPionPhScore, histBenignPionElScore, histBenignPionMuScore, histBenignPionPiScore, histBenignPionPrScore]
#Protons:
ProtonHistList = [histProtonPhScore, histProtonElScore, histProtonMuScore, histProtonPiScore, histProtonPrScore]
UnderProtonHistList = [histUnderProtonPhScore, histUnderProtonElScore, histUnderProtonMuScore, histUnderProtonPiScore, histUnderProtonPrScore]
DeadlyProtonHistList = [histDeadlyProtonPhScore, histDeadlyProtonElScore, histDeadlyProtonMuScore, histDeadlyProtonPiScore, histDeadlyProtonPrScore]
BenignProtonHistList = [histBenignProtonPhScore, histBenignProtonElScore, histBenignProtonMuScore, histBenignProtonPiScore, histBenignProtonPrScore]
#Figments:
NotTrackHistList = [histNotTrackPhScore, histNotTrackElScore, histNotTrackMuScore, histNotTrackPiScore, histNotTrackPrScore]

#SHOWER HISTOGRAMS
#Photons:
histPhotonPhScore = rt.TH1F("PhotonPhScore", "Shower Photon Score",30,0,10)
histPhotonElScore = rt.TH1F("PhotonPurity", "Shower Electron Score",30,0,10)
histPhotonMuScore = rt.TH1F("PhotonMuScore", "Shower Muon Score",30,0,10)
histPhotonPiScore = rt.TH1F("PhotonPiScore", "Shower Pion Score",30,0,10)
histPhotonPrScore = rt.TH1F("PhotonPrScore", "Shower Proton Score",30,0,10)
histPhotonCompPurity = rt.TH2F("PhotonCompPurity", "Correctly Identified, Under-Threshold Photon Comp vs. Purity",30,0,1,30,0,1)

histUnderPhotonPhScore = rt.TH1F("UnderPhotonPhScore", "Shower Photon Score",30,0,10)
histUnderPhotonElScore = rt.TH1F("UnderPhotonPurity", "Shower Electron Score",30,0,10)
histUnderPhotonMuScore = rt.TH1F("UnderPhotonMuScore", "Shower Muon Score",30,0,10)
histUnderPhotonPiScore = rt.TH1F("UnderPhotonPiScore", "Shower Pion Score",30,0,10)
histUnderPhotonPrScore = rt.TH1F("UnderPhotonPrScore", "Shower Proton Score",30,0,10)
histUnderPhotonCompPurity = rt.TH2F("UnderPhotonCompPurity", "Correctly Identified, Over-Threshold Photon Comp vs. Purity",30,0,1,30,0,1)

histDeadlyPhotonPhScore = rt.TH1F("DeadlyPhotonPhScore", "Shower Photon Score",30,0,10)
histDeadlyPhotonElScore = rt.TH1F("DeadlyPhotonPurity", "Shower Electron Score",30,0,10)
histDeadlyPhotonMuScore = rt.TH1F("DeadlyPhotonMuScore", "Shower Muon Score",30,0,10)
histDeadlyPhotonPiScore = rt.TH1F("DeadlyPhotonPiScore", "Shower Pion Score",30,0,10)
histDeadlyPhotonPrScore = rt.TH1F("DeadlyPhotonPrScore", "Shower Proton Score",30,0,10)
histDeadlyPhotonCompPurity = rt.TH2F("DeadlyPhotonCompPurity", "Incorrectly Identified, Over-Threshold Photon Comp vs. Purity",30,0,1,30,0,1)

histBenignPhotonPhScore = rt.TH1F("BenignPhotonPhScore", "Shower Photon Score",30,0,10)
histBenignPhotonElScore = rt.TH1F("BenignPhotonPurity", "Shower Electron Score",30,0,10)
histBenignPhotonMuScore = rt.TH1F("BenignPhotonMuScore", "Shower Muon Score",30,0,10)
histBenignPhotonPiScore = rt.TH1F("BenignPhotonPiScore", "Shower Pion Score",30,0,10)
histBenignPhotonPrScore = rt.TH1F("BenignPhotonPrScore", "Shower Proton Score",30,0,10)
histBenignPhotonCompPurity = rt.TH2F("BenignPhotonCompPurity", "Incorrectly Identified, Under-Threshold Photon Comp vs. Purity",30,0,1,30,0,1)
#Electrons:
histElectronPhScore = rt.TH1F("ElectronPhScore", "Shower Photon Score",30,0,10)
histElectronElScore = rt.TH1F("ElectronElScore", "Shower Electron Score",30,0,10)
histElectronMuScore = rt.TH1F("ElectronMuScore", "Shower Muon Score",30,0,10)
histElectronPiScore = rt.TH1F("ElectronPiScore", "Shower Pion Score",30,0,10)
histElectronPrScore = rt.TH1F("ElectronPrScore", "Shower Proton Score",30,0,10)
histElectronCompPurity = rt.TH2F("ElectronCompPurity", "Correctly Identified, Under Threshold Electron Comp vs. Purity",30,0,1,30,0,1)

histUnderElectronPhScore = rt.TH1F("UnderElectronPhScore", "Shower Photon Score",30,0,10)
histUnderElectronElScore = rt.TH1F("UnderElectronElScore", "Shower Electron Score",30,0,10)
histUnderElectronMuScore = rt.TH1F("UnderElectronMuScore", "Shower Muon Score",30,0,10)
histUnderElectronPiScore = rt.TH1F("UnderElectronPiScore", "Shower Pion Score",30,0,10)
histUnderElectronPrScore = rt.TH1F("UnderElectronPrScore", "Shower Proton Score",30,0,10)
histUnderElectronCompPurity = rt.TH2F("UnderElectronCompPurity", "Correctly Identified, Over-Threshold Electron Comp vs. Purity",30,0,1,30,0,1)

histDeadlyElectronPhScore = rt.TH1F("DeadlyElectronPhScore", "Shower Photon Score",30,0,10)
histDeadlyElectronElScore = rt.TH1F("DeadlyElectronElScore", "Shower Electron Score",30,0,10)
histDeadlyElectronMuScore = rt.TH1F("DeadlyElectronMuScore", "Shower Muon Score",30,0,10)
histDeadlyElectronPiScore = rt.TH1F("DeadlyElectronPiScore", "Shower Pion Score",30,0,10)
histDeadlyElectronPrScore = rt.TH1F("DeadlyElectronPrScore", "Shower Proton Score",30,0,10)
histDeadlyElectronCompPurity = rt.TH2F("DeadlyElectronCompPurity", "Incorrectly Identified, Over-Threshold Electron Comp vs. Purity",30,0,1,30,0,1)

histBenignElectronPhScore = rt.TH1F("ElectronBenignPhScore", "Shower Photon Score",30,0,10)
histBenignElectronElScore = rt.TH1F("ElectronBenignElScore", "Shower Electron Score",30,0,10)
histBenignElectronMuScore = rt.TH1F("ElectronBenignMuScore", "Shower Muon Score",30,0,10)
histBenignElectronPiScore = rt.TH1F("ElectronBenignPiScore", "Shower Pion Score",30,0,10)
histBenignElectronPrScore = rt.TH1F("ElectronBenignPrScore", "Shower Proton Score",30,0,10)
histBenignElectronCompPurity = rt.TH2F("BenignElectronCompPurity", "Incorrectly Identified, Under-Threshold Electron Comp vs. Purity",30,0,1,30,0,1)
#Figments:
histNotShowerPhScore = rt.TH1F("NotShowerPhScore", "Shower Photon Score",30,0,10)
histNotShowerElScore = rt.TH1F("NotShowerElScore", "Shower Electron Score",30,0,10)
histNotShowerMuScore = rt.TH1F("NotShowerMuScore", "Shower Muon Score",30,0,10)
histNotShowerPiScore = rt.TH1F("NotShowerPiScore", "Shower Pion Score",30,0,10)
histNotShowerPrScore = rt.TH1F("NotShowerPrScore", "Shower Proton Score",30,0,10)
histNotShowerCompPurity = rt.TH2F("NotShowerCompPurity", "Fake Shower Comp vs. Purity",30,0,1,30,0,1)

#SHOWER HISTLISTS
PhotonHistList = [histPhotonPhScore, histPhotonElScore, histPhotonMuScore, histPhotonPiScore, histPhotonPrScore]
UnderPhotonHistList = [histUnderPhotonPhScore, histUnderPhotonElScore, histUnderPhotonMuScore, histUnderPhotonPiScore, histUnderPhotonPrScore]
DeadlyPhotonHistList = [histDeadlyPhotonPhScore, histDeadlyPhotonElScore, histDeadlyPhotonMuScore, histDeadlyPhotonPiScore, histDeadlyPhotonPrScore]
BenignPhotonHistList = [histBenignPhotonPhScore, histBenignPhotonElScore, histBenignPhotonMuScore, histBenignPhotonPiScore, histBenignPhotonPrScore]

ElectronHistList = [histElectronPhScore, histElectronElScore, histElectronMuScore, histElectronPiScore, histElectronPrScore]
UnderElectronHistList = [histUnderElectronPhScore, histUnderElectronElScore, histUnderElectronMuScore, histUnderElectronPiScore, histUnderElectronPrScore]
DeadlyElectronHistList = [histDeadlyElectronPhScore, histDeadlyElectronElScore, histDeadlyElectronMuScore, histDeadlyElectronPiScore, histDeadlyElectronPrScore]
BenignElectronHistList = [histBenignElectronPhScore, histBenignElectronElScore, histBenignElectronMuScore, histBenignElectronPiScore, histBenignElectronPrScore]

NotShowerHistList = [histNotShowerPhScore, histNotShowerElScore, histNotShowerMuScore, histNotShowerPiScore, histNotShowerPrScore]

#Built-in functions here
def addHistTracks(ntuple, trackNo, histList, compPurHist, weight):
  #Adds the track's data to all relevant hists
  histList[0].Fill(abs(ntuple.trackPhScore[trackNo]), weight)
  histList[1].Fill(abs(ntuple.trackElScore[trackNo]), weight)
  histList[2].Fill(abs(ntuple.trackMuScore[trackNo]), weight)
  histList[3].Fill(abs(ntuple.trackPiScore[trackNo]), weight)
  histList[4].Fill(abs(ntuple.trackPrScore[trackNo]), weight)
  compPurHist.Fill(ntuple.trackComp[trackNo], ntuple.trackPurity[trackNo], weight)

def addHistShowers(ntuple, showerNo, histList, compPurHist, weight):
  #Adds the shower's data to all relevant hists
  histList[0].Fill(ntuple.showerComp[showerNo], weight)
  histList[1].Fill(ntuple.showerPurity[showerNo], weight)
  histList[2].Fill(ntuple.showerPhScore[showerNo], weight)
  histList[3].Fill(ntuple.showerElScore[showerNo], weight)
  compPurHist.Fill(ntuple.showerComp[showerNo], ntuple.showerPurity[showerNo], weight)

def histScale(hist, POTSum):
  POTTarget = 6.67e+20
  hist.Scale(POTTarget/POTSum)
  return hist

#Variables for program review
weirdEvent = 0

#We put this into the addHist function for truth-based graphs
emptyList = []

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
  if recoCutMuons(eventTree) == False:
    continue

  if recoCutElectrons(eventTree) == False:
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

  #Cut events with far-travelling protons
  recoProtonCount = recoProton(eventTree)
  if recoProtonCount > 1:
    continue

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, eventTree)

  #List all photons classified as tracks
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree)

  if len(recoList) + len(recoTrackList) == 0:
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

  #Cut based on number of tracks
  #if recoCutManyTracks(eventTree) == False:
  #  continue

  #Cut unclassified tracks too short to be trusted
  #if recoCutShortTracks(eventTree) == False:
  #  continue

  #Now we try to match reconstructed tracks to real ones
  for x in range(eventTree.nTracks):
    particleEnergy = 0
    for y in range(eventTree.nTrueSimParts):
      #Calculate particle energy using the truth-matched particle
      if eventTree.trackTrueTID[x] == eventTree.trueSimPartTID[y]:
        particleEnergy = eventTree.trueSimPartE[y] - np.sqrt(abs(eventTree.trueSimPartE[y]**2 - (eventTree.trueSimPartPx[y]**2+eventTree.trueSimPartPy[y]**2+eventTree.trueSimPartPz[y]**2)))
      #distxyz = np.sqrt((eventTree.trackStartPosX[x] - eventTree.trackEndPosX[x])**2 + (eventTree.trackStartPosY[x] - eventTree.trackEndPosY[x])**2 + (eventTree.trackStartPosZ[x] - eventTree.trackEndPosZ[x])**2)

      #See if the track is a muon
      if eventTree.trackTruePID[x] == 13:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(len(eventTree.truePrimPartPDG)):
            if eventTree.truePrimPartPDG[y] == 13:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        #See if Reco was able to identify it correctly, and if so, if it got the energy right (to the extent that it matters)
        if eventTree.trackPID[x] == 13:
          if particleEnergy > 100: 
            addHistTracks(eventTree, x, UnderMuonHistList, histUnderMuonCompPurity, eventTree.xsecWeight)
          else:
            addHistTracks(eventTree, x, MuonHistList, histMuonCompPurity, eventTree.xsecWeight)

        elif eventTree.trackPID[x] != 0:
          if particleEnergy > 100:
            addHistTracks(eventTree, x, DeadlyMuonHistList, histDeadlyMuonCompPurity, eventTree.xsecWeight)
          else:
            addHistTracks(eventTree, x, BenignMuonHistList, histBenignMuonCompPurity, eventTree.xsecWeight)


      #See if the track is a pion
      elif abs(eventTree.trackTruePID[x]) == 211:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(len(eventTree.truePrimPartPDG)):
            if abs(eventTree.truePrimPartPDG[y]) == 211:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        #See if Reco was able to identify it correctly, and if so, if it got the energy right (to the extent that it matters)
        if eventTree.trackPID[x] == 211:
          if particleEnergy >= 50: 
            addHistTracks(eventTree, x, UnderPionHistList, histUnderPionCompPurity, eventTree.xsecWeight)
          else:
            addHistTracks(eventTree, x, PionHistList, histPionCompPurity, eventTree.xsecWeight)

        elif eventTree.trackPID[x] != 0:
          if particleEnergy >= 50:
            addHistTracks(eventTree, x, DeadlyPionHistList, histDeadlyPionCompPurity, eventTree.xsecWeight)
          else:
            addHistTracks(eventTree, x, BenignPionHistList, histBenignPionCompPurity, eventTree.xsecWeight)

      #See if the track is a proton
      elif eventTree.trackTruePID[x] == 2212:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(eventTree.nTruePrimParts):
            if eventTree.truePrimPartPDG[y] == 2212:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        #See if Reco was able to identify it correctly, and if so, if it got the energy right (to the extent that it matters)
        if eventTree.trackPID[x] == 2212:
          if particleEnergy >= 100: 
            addHistTracks(eventTree, x, UnderProtonHistList, histUnderProtonCompPurity, eventTree.xsecWeight)
          else:
            addHistTracks(eventTree, x, ProtonHistList, histProtonCompPurity, eventTree.xsecWeight)

        elif eventTree.trackPID[x] != 0:
          if particleEnergy >= 100:
            addHistTracks(eventTree, x, DeadlyProtonHistList, histDeadlyProtonCompPurity, eventTree.xsecWeight)
          else:
            addHistTracks(eventTree, x, BenignProtonHistList, histBenignProtonCompPurity, eventTree.xsecWeight)
      
      #See if the track is just a figment of the reco's imagination
      elif eventTree.trackTruePID[x] == 0: 
        addHistTracks(eventTree, x, NotTrackHistList, histNotTrackCompPurity, eventTree.xsecWeight)

      #Otherwise, it's something weird; we're just going to track the number of these, for now
      else:
        weirdEvent += 1


  #NOW WE ITERATE OVER SHOWERS
  for x in range(eventTree.nShowers):
    particleEnergy = 0
    for y in range(eventTree.nTrueSimParts):
      #Calculate particle energy using the truth-matched particle
      if eventTree.showerTrueTID[x] == eventTree.trueSimPartTID[y]:
        particleEnergy = eventTree.trueSimPartE[y]
      #distxyz = np.sqrt((eventTree.trackStartPosX[x] - eventTree.trackEndPosX[x])**2 + (eventTree.trackStartPosY[x] - eventTree.trackEndPosY[x])**2 + (eventTree.trackStartPosZ[x] - eventTree.trackEndPosZ[x])**2)

      #See if the track is a photon
      if eventTree.showerTruePID[x] == 22:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(len(eventTree.truePrimPartPDG)):
            if eventTree.truePrimPartPDG[y] == 22:
              particleEnergy = eventTree.truePrimPartE[y]

        #See if Reco was able to identify it correctly, and if so, if it got the energy right (to the extent that it matters)
        if eventTree.showerPID[x] == 22:
          if particleEnergy >= 20: 
            addHistShowers(eventTree, x, UnderPhotonHistList, histUnderPhotonCompPurity, eventTree.xsecWeight)
          else:
            addHistShowers(eventTree, x, PhotonHistList, histPhotonCompPurity, eventTree.xsecWeight)

        elif eventTree.showerPID[x] != 0:
          if particleEnergy >= 20:
            addHistShowers(eventTree, x, DeadlyPhotonHistList, histDeadlyPhotonCompPurity, eventTree.xsecWeight)
          else:
            addHistShowers(eventTree, x, BenignPhotonHistList, histBenignPhotonCompPurity, eventTree.xsecWeight)


      #See if the track is a electron
      elif abs(eventTree.showerTruePID[x]) == 11:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(len(eventTree.truePrimPartPDG)):
            if abs(eventTree.truePrimPartPDG[y]) == 211:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        
        #See if Reco was able to identify it correctly, and if so, if it got the energy right (to the extent that it matters)
        if eventTree.showerPID[x] == 11:
          if particleEnergy >= 10: 
            addHistShowers(eventTree, x, UnderElectronHistList, histUnderElectronCompPurity, eventTree.xsecWeight)
          else:
            addHistShowers(eventTree, x, ElectronHistList, histElectronCompPurity, eventTree.xsecWeight)

        elif eventTree.showerPID[x] != 0:
          if particleEnergy >= 10:
            addHistShowers(eventTree, x, DeadlyElectronHistList, histDeadlyElectronCompPurity, eventTree.xsecWeight)
          else:
            addHistShowers(eventTree, x, BenignElectronHistList, histBenignElectronCompPurity, eventTree.xsecWeight)


      #See if the shower is just a figment of the reco's imagination
      elif eventTree.showerTruePID[x] == 0: 
        addHistShowers(eventTree, x, NotShowerHistList, histNotShowerCompPurity, eventTree.xsecWeight)

      #Otherwise, it's something weird; we're just going to track the number of these, for now
      else:
        weirdEvent += 1



#LOOPS OVER - HISTOGRAM ORGANIZING TIME
print("There were", weirdEvent, "events that defied current categorization.")

#Stacking histograms
muonCanvas, muonStack, muonLegend, muonInt = histStack("Correctly Identified, Under-Threshold Muons", MuonHistList, ntuplePOTsum)
muonUnderCanvas, muonUnderStack, muonUnderLegend, muonUnderInt = histStack("Correctly Identified, Over-Threshold Muons", UnderMuonHistList, ntuplePOTsum)
muonDeadlyCanvas, muonDeadlyStack, muonDeadlyLegend, muonDeadlyInt = histStack("Incorrectly Identified, Over-Threshold Muons", DeadlyMuonHistList, ntuplePOTsum)
muonBenignCanvas, muonBenignStack, muonBenignLegend, muonBenignInt = histStack("Incorrectly Identified, Under-Threshold Muons", BenignMuonHistList, ntuplePOTsum)

pionCanvas, pionStack, pionLegend, pionInt = histStack("Correctly Identified, Under-Threshold Pions", PionHistList, ntuplePOTsum)
pionUnderCanvas, pionUnderStack, pionUnderLegend, pionUnderInt = histStack("Correctly Identified, Over-Threshold Pions", UnderPionHistList, ntuplePOTsum)
pionDeadlyCanvas, pionDeadlyStack, pionDeadlyLegend, pionDeadlyInt = histStack("Incorrectly Identified, Over-Threshold Pions", DeadlyPionHistList, ntuplePOTsum)
pionBenignCanvas, pionBenignStack, pionBenignLegend, pionBenignInt = histStack("Incorrectly Identified, Under-Threshold Pions", BenignPionHistList, ntuplePOTsum)

protonCanvas, protonStack, protonLegend, protonInt = histStack("Correctly Identified, Under-Threshold Protons", ProtonHistList, ntuplePOTsum)
protonUnderCanvas, protonUnderStack, protonUnderLegend, protonUnderInt = histStack("Correctly Identified, Over-Threshold Protons", UnderProtonHistList, ntuplePOTsum)
protonDeadlyCanvas, protonDeadlyStack, protonDeadlyLegend, protonDeadlyInt = histStack("Incorrectly Identified, Over-Threshold Protons", DeadlyProtonHistList, ntuplePOTsum)
protonBenignCanvas, protonBenignStack, protonBenignLegend, protonBenignInt = histStack("Incorrectly Identified, Under-Threshold Protons", BenignProtonHistList, ntuplePOTsum)


notTrackCanvas, notTrackStack, notTrackLegend, notTrackInt = histStack("False Tracks (no particle identified)", NotTrackHistList, ntuplePOTsum)

photonCanvas, photonStack, photonLegend, photonInt = histStack("Correctly Identified, Under-Threshold Photons", PhotonHistList, ntuplePOTsum)
photonUnderCanvas, photonUnderStack, photonUnderLegend, photonUnderInt = histStack("Correctly Identified, Over-Threshold Photons", UnderPhotonHistList, ntuplePOTsum)
photonDeadlyCanvas, photonDeadlyStack, photonDeadlyLegend, photonDeadlyInt = histStack("Incorrectly Identified, Over-Threshold Photons", DeadlyPhotonHistList, ntuplePOTsum)
photonBenignCanvas, photonBenignStack, photonBenignLegend, photonBenignInt = histStack("Incorrectly Identified, Under-Threshold Photons", BenignPhotonHistList, ntuplePOTsum)

electronCanvas, electronStack, electronLegend, electronInt = histStack("Correctly Identified, Under-Threshold Electrons", ElectronHistList, ntuplePOTsum)
electronUnderCanvas, electronUnderStack, electronUnderLegend, electronUnderInt = histStack("Correctly Identified, Over-Threshold Electrons", UnderElectronHistList, ntuplePOTsum)
electronDeadlyCanvas, electronDeadlyStack, electronDeadlyLegend, electronDeadlyInt = histStack("Incorrectly Identified, Over-Threshold Electrons", DeadlyElectronHistList, ntuplePOTsum)
electronBenignCanvas, electronBenignStack, electronBenignLegend, electronBenignInt = histStack("Incorrectly Identified, Under-Threshold Electrons", BenignElectronHistList, ntuplePOTsum)

notShowerCanvas, notShowerStack, notShowerLegend, notShowerInt = histStack("False Showers (no particle identified)", NotShowerHistList, ntuplePOTsum)

canvasList = [muonCanvas, muonUnderCanvas, muonDeadlyCanvas, muonBenignCanvas, pionCanvas, pionUnderCanvas, pionDeadlyCanvas, pionBenignCanvas, protonCanvas, protonUnderCanvas, protonDeadlyCanvas, protonBenignCanvas, photonCanvas, photonUnderCanvas, photonDeadlyCanvas, photonBenignCanvas, electronCanvas, electronUnderCanvas, electronDeadlyCanvas, electronBenignCanvas, notShowerCanvas, notTrackCanvas]
stackList = [muonStack, muonUnderStack, muonDeadlyStack, muonBenignStack, pionStack, pionUnderStack, pionDeadlyStack, pionBenignStack, protonStack, protonUnderStack, protonDeadlyStack, protonBenignStack, photonStack, photonUnderStack, photonDeadlyStack, photonBenignStack, electronStack, electronUnderStack, electronDeadlyStack, electronBenignStack, notShowerStack, notTrackStack]
Hist2dList = [histMuonCompPurity, histUnderMuonCompPurity, histDeadlyMuonCompPurity, histBenignMuonCompPurity, histPionCompPurity, histUnderPionCompPurity, histDeadlyPionCompPurity, histBenignPionCompPurity, histProtonCompPurity, histUnderProtonCompPurity, histDeadlyProtonCompPurity, histBenignProtonCompPurity, histPhotonCompPurity, histUnderPhotonCompPurity, histDeadlyPhotonCompPurity, histBenignPhotonCompPurity, histElectronCompPurity, histUnderElectronCompPurity, histDeadlyElectronCompPurity, histBenignElectronCompPurity, histNotTrackCompPurity, histNotShowerCompPurity]

#Label axes correctly
for stack in stackList:
  stack.GetXaxis().SetTitle("Particle Scores")

for hist in Hist2dList:
  hist.GetXaxis().SetTitle("Completeness")
  hist.GetYaxis().SetTitle("Purity")

#Now all that's left to do is write the canvases and histograms to the file
outFile = rt.TFile(args.outfile, "RECREATE")
for hist in Hist2dList:
  hist.Write()

for canvas in canvasList:
  canvas.Write()