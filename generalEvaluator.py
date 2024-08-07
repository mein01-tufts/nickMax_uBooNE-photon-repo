import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials,trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoPion, recoProton, recoCutElectronScore, recoCutShowerFromChargeScore, recoCutLongTracks, recoPhotonListFiducial, recoCutPrimary, recoCutShortTracks, recoPhotonListTracks, recoCutFarShowers

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-c", "--cosmicFile", type=str, required=True, help="cosmic input file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")

args = parser.parse_args()

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

cosmic_file = rt.TFile(args.cosmicFile)
cosmicTree = cosmic_file.Get("EventTree")
cosmicPotTree = ntuple_file.Get("cosmicPotTree")

#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20"
ntuplePOTsum = 0

for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

cosmicPOTsum = 3.2974607516739725e+20

#Hists created and organized here
#PURITY HISTOGRAMS
purityTotal1 = rt.TH1F("PTotal1", "One Photon",60,0,2)
puritySignal1 =  rt.TH1F("PSignal1", "Signal",60,0,2)
purityCC1 = rt.TH1F("PCC1", "Actually Charged Current",60,0,2)
purityFiducials1 = rt.TH1F("PFiducial1", "Out of Fiducial",60,0,2)
purityCosmic1 = rt.TH1F("PCosmic1", "Failed Cosmic",60,0,2)
purityPionProton1 = rt.TH1F("PPionProton1", "Charged Pion or Proton",60,0,2)
purityNoPhotons1 = rt.TH1F("PNoPhoton1", "No Real Photons",60,0,2)
purityTwoPhotons1 = rt.TH1F("PTwoPhoton1", "2 Real Photons",60,0,2)
purityManyPhotons1 = rt.TH1F("PMorePhoton1", "3+ Real Photons",60,0,2)

twoPhotonsLost1 = rt.TH1F("TwoPhotonLost", "2 Real Photons (retrievable)",60,0,2)
twoPhotonsUnclassified1 = rt.TH1F("PTwoPhotonUnclassified", "2 Real Photons (irretrievable)",60,0,2)

purityTotal2 = rt.TH1F("PTotal2", "One Photon",60,0,2)
puritySignal2 =	 rt.TH1F("PSignal2", "Signal",60,0,2)
purityCC2 = rt.TH1F("PCC2", "Actually Charged Current",60,0,2)
purityFiducials2 = rt.TH1F("PFiducial2", "Out of Fiducial",60,0,2)
purityCosmic2 = rt.TH1F("PCosmic2", "Failed Cosmic",60,0,2)
purityPionProton2 = rt.TH1F("PPionProton2", "Charged Pion or Proton",60,0,2)
purityNoPhotons2 = rt.TH1F("PNoPhoton2", "No Real Photons",60,0,2)
purityOnePhoton2 = rt.TH1F("POnePhoton2", "1 Real Photon",60,0,2)
purityManyPhotons2 = rt.TH1F("PManyPhoton2", "3+ Real Photons",60,0,2)

purityTotal3 = rt.TH1F("PTotal3", "One Photon",60,0,2)
puritySignal3 =	 rt.TH1F("PSignal3", "Signal",60,0,2)
purityCC3 = rt.TH1F("PCC3", "Actually Charged Current",60,0,2)
purityFiducials3 = rt.TH1F("PFiducial3", "Out of Fiducial",60,0,2)
purityCosmic3 = rt.TH1F("PCosmic3", "Failed Cosmic",60,0,2)
purityPionProton3 = rt.TH1F("PPionProton3", "Charged Pion or Proton",60,0,2)
purityNoPhotons3 = rt.TH1F("PNoPhoton3", "No Real Photons",60,0,2)
purityOnePhoton3 = rt.TH1F("POnePhoton3", "2 Real Photons",60,0,2)
purityTwoPhotons3 = rt.TH1F("PTwoPhotons33", "3+ Real Photons",60,0,2)

#PURITY HISTLISTS
totalPHists = [purityTotal1, purityTotal2, purityTotal3]
signalPHists = [puritySignal1, puritySignal2, puritySignal3]
CCPHists = [purityCC1, purityCC2, purityCC3]
fiducialPHists = [purityFiducials1, purityFiducials2, purityFiducials3]
cosmicPHists = [purityCosmic1, purityCosmic2, purityCosmic3]
pionProtonPHists = [purityPionProton1, purityPionProton2, purityPionProton3]
noPhotonPHists = [purityNoPhotons1, purityNoPhotons2, purityNoPhotons3]
onePhotonPHists = [puritySignal1, purityOnePhoton2, purityOnePhoton3]
twoPhotonPHists = [purityTwoPhotons1, puritySignal2, purityTwoPhotons3]
manyPhotonPHists = [purityManyPhotons1, purityManyPhotons2, puritySignal3]

#Cosmics go here, so we can put them on purity
cosmicOnePhoton = rt.TH1F("cBackground1", "One Photon Cosmic",60,0,2)
cosmicTwoPhotons = rt.TH1F("cBackground2", "Two Photon Cosmic",60,0,2)
cosmicThreePhotons = rt.TH1F("cBackground3", "3+ Photon Cosmic",60,0,2)
cosmicList = [cosmicOnePhoton, cosmicTwoPhotons, cosmicThreePhotons]
cosmicThreshold = rt.TH1F("cosmicThreshold", "Necessary Signal to Beat Cosmics (Single Photon)",60,0,2)

#Big Lists, for Big Plots
pList1 = [puritySignal1, purityCC1, purityFiducials1, purityCosmic1, purityPionProton1, purityNoPhotons1, purityTwoPhotons1, purityManyPhotons1]
pList2 = [puritySignal2, purityCC2, purityFiducials2, purityCosmic2, purityPionProton2, purityNoPhotons2, purityOnePhoton2, purityManyPhotons2] 
pList3 = [puritySignal3, purityCC3, purityFiducials3, purityCosmic3, purityPionProton3,	purityNoPhotons3, purityOnePhoton3, purityTwoPhotons3]


#EFFICIENCY HISTOGRAMS
effTotal1 = rt.TH1F("effTotal1", "One Photon",60,0,2)
effNoVertex1 = rt.TH1F("effNoVertex1", "No Vertex Found",60,0,2)
effCC1 = rt.TH1F("effCC1", "CC False Positive",60,0,2)
effFiducial1 = rt.TH1F("effFiducial1", "Placed out of Fiducial",60,0,2)
effPion1 = rt.TH1F("effPion1", "Pion False Positive",60,0,2)
effProton1 =  rt.TH1F("effProton1", "Proton False Positive",60,0,2)
effNoPhotons1 = rt.TH1F("effNoPhotons1", "No Photons Found",60,0,2)
effSignal1 = rt.TH1F("effSignal 1", "Signal",60,0,2)
effTwoPhotons1 = rt.TH1F("effTwoPhotons1", "Two Photons Found",60,0,2)
effManyPhotons1 = rt.TH1F("effManyPhotons1", "Many Photons Found",60,0,2)
effShowerCharge1 = rt.TH1F("effShowerCharge1", "Shower from Charged Cut",60,0,2)
effPrimary1 = rt.TH1F("effPrimary1", "Primary Score Cut",60,0,2)
effLongTracks1 = rt.TH1F("effLongTracks1", "Tracks with length > 20 cm",60,0,2)
effShortTrack1 = rt.TH1F("effShortTrack1", "Unclassified tracks with length < 10 cm",60,0,2)
#effLongShowers1 = rt.TH1F("effLongShowers1", "Non-Photon showers > 80 cm from vertex",60,0,2)

effTotal2 = rt.TH1F("effTotal2", "Two Photons",60,0,2)
effNoVertex2 = rt.TH1F("effNoVertex2", "No Vertex Found",60,0,2)
effCC2 = rt.TH1F("effCC2", "CC False Positive",60,0,2)
effFiducial2 = rt.TH1F("effFiducial2", "Placed out of Fiducial",60,0,2)
effPion2 = rt.TH1F("effPion2", "Pion False Positive",60,0,2)
effProton2 =  rt.TH1F("effProton2", "Proton False Positive",60,0,2)
effNoPhotons2 =	rt.TH1F("effNoPhotons2", "No Photons Found",60,0,2)
effSignal2 = rt.TH1F("effSignal2", "Signal",60,0,2)
effOnePhoton2 = rt.TH1F("effOnePhoton2", "One Photon Found",60,0,2)
effManyPhotons2 = rt.TH1F("effManyPhotons2", "Many Photons Found",60,0,2)
effShowerCharge2 = rt.TH1F("effShowerCharge2", "Shower from Charged Cut",60,0,2)
effPrimary2 = rt.TH1F("effPrimary2", "Primary Score Cut",60,0,2)
effLongTracks2 = rt.TH1F("effLongTracks2", "Tracks with length > 20 cm",60,0,2)
effShortTrack2 = rt.TH1F("effShortTrack2", "Unclassified tracks with length < 10 cm",60,0,2)
#effLongShowers2 = rt.TH1F("effLongShowers2", "Non-Photon showers > 80 cm from vertex",60,0,2)

effTotal3 = rt.TH1F("effTotal3", "3+ Photons",60,0,2)
effNoVertex3 = rt.TH1F("effNoVertex3", "No Vertex Found",60,0,2)
effCC3 = rt.TH1F("effCC3", "CC False Positive",60,0,2)
effFiducial3 = rt.TH1F("effFiducial3", "Placed out of Fiducial",60,0,2)
effPion3 = rt.TH1F("effPion3", "Pion False Positive",60,0,2)
effProton3 =  rt.TH1F("effProton3", "Proton False Positive",60,0,2)
effNoPhotons3 =	rt.TH1F("effNoPhotons3", "No Photons Found",60,0,2)
effSignal3 = rt.TH1F("effSignal 3", "Signal",60,0,2)
effOnePhoton3 = rt.TH1F("effManyPhotons3", "Many Photons Found",60,0,2)
effTwoPhotons3 = rt.TH1F("effTwoPhotons3", "Two Photons Found",60,0,2)
effShowerCharge3 = rt.TH1F("effShowerCharge3", "Shower from Charged Cut",60,0,2)
effPrimary3 = rt.TH1F("effPrimary3", "Primary Score Cut",60,0,2)
effLongTracks3 = rt.TH1F("effLongTracks3", "Tracks with length > 20 cm",60,0,2)
effShortTrack3 = rt.TH1F("effShortTrack3", "Unclassified tracks with length < 10 cm",60,0,2)
#effLongShowers3 = rt.TH1F("effLongShowers3", "Non-Photon showers > 80 cm from vertex",60,0,2)

#Histogram Lists!
effTotalList = [effTotal1, effTotal2, effTotal3]
effNoVertexHists = [effNoVertex1, effNoVertex2, effNoVertex3]
effCCHists = [effCC1, effCC2, effCC3]
effFiducialHists = [effFiducial1, effFiducial2, effFiducial3]
effPionHists = [effPion1, effPion2, effPion3]
effProtonHists = [effProton1, effProton2, effProton3]
effNoPhotonHists = [effNoPhotons1, effNoPhotons2, effNoPhotons3]
effOnePhotonHists = [effSignal1, effOnePhoton2, effOnePhoton3] 
effTwoPhotonHists = [effTwoPhotons1, effSignal2, effTwoPhotons3]
effManyPhotonHists = [effManyPhotons1, effManyPhotons2, effSignal3]
effShowerChargeHists = [effShowerCharge1, effShowerCharge2, effShowerCharge3]
effLongTrackHists = [effLongTracks1, effLongTracks2, effLongTracks3]
effPrimaryHists = [effPrimary1, effPrimary2, effPrimary3]
effShortTrackHists = [effShortTrack1, effShortTrack2, effShortTrack3]
#effFarShowers = [effLongShowers1, effLongShowers2, effLongShowers3]

effList1 = [effSignal1, effNoVertex1, effCC1, effPion1, effProton1, effShowerCharge1, effNoPhotons1, effTwoPhotons1, effManyPhotons1, effPrimary1, effLongTracks1, effShortTrack1]
effList2 = [effSignal2, effNoVertex2, effCC2, effPion2, effProton2, effShowerCharge1, effNoPhotons2, effOnePhoton2, effManyPhotons2, effPrimary2, effLongTracks2, effShortTrack1]
effList3 = [effSignal3, effNoVertex3, effCC3, effPion3, effProton3, effShowerCharge1, effNoPhotons3, effOnePhoton3, effTwoPhotons3, effPrimary3, effLongTracks3, effShortTrack1]

#Built-in functions here
def addHist(eventTree, photonList, photonList2, histList, variable, weight):
  if len(photonList) + len(photonList2) == 1:
    histList[0].Fill(variable, weight)
  elif len(photonList) + len(photonList2) == 2:
    histList[1].Fill(variable, weight)
  else:
    histList[2].Fill(variable, weight)

def purityStack(title, purityList, cosmicHist, POTSum, cosmicSum):
  #Create Component Variables
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt.kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kGreen+4, rt.kOrange+1]
  POTTarget = 6.67e+20
  histIntTotal = 0
  #Organize other histograms
  for x in range(len(purityList)):
    hist = purityList[x]
    bins = hist.GetNbinsX()
    hist.Scale(POTTarget/POTSum)
    hist.SetLineColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    histIntTotal += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
    stack.Add(hist)

  #Format and add the cosmic histogram
  hist = cosmicHist
  bins = hist.GetNbinsX()
  hist.Scale(POTTarget/cosmicSum)
  hist.SetLineColor(rt.kBlack)
  histInt = hist.Integral(1, int(bins))
  histIntTotal += histInt
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  legendHeaderString = "Total: " + str(round((histIntTotal),1)) 
  legend.SetHeader(str(legendHeaderString), "C")
  stack.Add(hist)
  #Finish working on the Canvas and return necessary components
  histCanvas = rt.TCanvas() 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle("Photon Energy (GeV)")
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()
  return histCanvas, stack, legend, histInt

def histScale(hist, POTSum):
  POTTarget = 6.67e+20
  hist.Scale(POTTarget/POTSum)
  return hist

#Variables for program review
initialCount = 0
vertexCount = 0
NCCount = 0
fiducialCount = 0
noPionCount = 0
recoCount = 0
onePhoton = 0
twoPhotons = 0
threePhotons = 0
count = 0

totalCosmics = 0
vertexCosmics = 0 
NCCosmics = 0
MattProofCosmics = 0
fiducialCosmics = 0
pionlessCosmics = 0
photonCosmics = 0
uncutCosmics = 0
untrackedCosmics = 0

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

  #Cut unclassified tracks too short to be trusted
  #if recoCutShortTracks(eventTree) == False:
  #  continue


  #PURITY - GRAPHING BASED ON TRUTH
  
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(eventTree, recoList, recoTrackList)

  #Fill totals
  addHist(eventTree, recoList, recoTrackList, totalPHists, leadingPhoton, eventTree.xsecWeight)
  initialCount += 1
  #Neutral current!
  if trueCutNC(eventTree) == False:
    addHist(eventTree, recoList, recoTrackList, CCPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  NCCount += 1

  #I suppose we can pretend that this is doing something
  if trueCutCosmic(eventTree) == False:
    addHist(eventTree, recoList, recoTrackList, cosmicPHists, leadingPhoton, eventTree.xsecWeight)
    continue

  if trueCutFiducials(eventTree, fiducialData) == False:
    addHist(eventTree, recoList, recoTrackList, fiducialPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  fiducialCount += 1

  #pions and protons!
  pionCount, protonCount = trueCutPionProton(eventTree)
  if pionCount > 0 or protonCount > 1:
    addHist(eventTree, recoList, recoTrackList, pionProtonPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  noPionCount += 1

  #Now we make a list of the actual photons!
  truePhotonIDs = truePhotonList(eventTree, fiducialData)

  #Are there actually any photons?
  if len(truePhotonIDs) == 0:
    addHist(eventTree, recoList, recoTrackList, noPhotonPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  recoCount += 1
  #Is there one Photon?
  if len(truePhotonIDs) == 1:
    addHist(eventTree, recoList, recoTrackList, onePhotonPHists, leadingPhoton, eventTree.xsecWeight)
    onePhoton += 1
  #Are there two?
  elif len(truePhotonIDs) == 2:
    addHist(eventTree, recoList, recoTrackList, twoPhotonPHists, leadingPhoton, eventTree.xsecWeight)
    twoPhotons += 1
  
  #In that case, there should be at least three
  else:
    addHist(eventTree, recoList, recoTrackList, manyPhotonPHists, leadingPhoton, eventTree.xsecWeight)
    threePhotons += 1

#BEGINNING EVENT LOOP FOR COSMICS
for i in range(cosmicTree.GetEntries()):
  cosmicTree.GetEntry(i)

  #COSMICS - SELECTING EVENTS BASED ON RECO
  #See if the event has a vertex
  totalCosmics += 1

  if recoNoVertex(cosmicTree) == False:
    continue
  vertexCosmics += 1
  #See if the event is neutral current                                                                                                  
  if recoNeutralCurrent(cosmicTree) == False:
    continue
  NCCosmics += 1
  #Use Matt's Cosmic Cut
  if trueCutCosmic(cosmicTree) == False:
    continue 
  MattProofCosmics += 1
  #Make sure the event is within the fiducial volume
  if recoFiducials(cosmicTree, fiducialData) == False:
    continue
  fiducialCosmics += 1
  #Cut events with suitably energetic protons or charged pions
  if recoPion(cosmicTree) == False:
    continue
  pionlessCosmics += 1
  #Cut events with far-travelling protons
  recoProtonCount = recoProton(cosmicTree)
  if recoProtonCount > 1:
    continue

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, cosmicTree)

  #List all photons classified as tracks
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree)
  
  if len(recoList) + len(recoTrackList) == 0:
    continue
  photonCosmics += 1
  
  #Try cutting based on data for Shower from Charged
  if recoCutShowerFromChargeScore(cosmicTree, recoList, recoTrackList) == False:
    continue

  #Cut based on primary score
  if recoCutPrimary(cosmicTree, recoList, recoTrackList) == False:
    continue
  uncutCosmics += 1

  #Longtracks cut
  if recoCutLongTracks(cosmicTree, fiducialData) == False:
    continue

  #Cut unclassified tracks too short to be trusted
  #if recoCutShortTracks(eventTree) == False:
  #  continue

  untrackedCosmics += 1

  #graphing based on photon count
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(cosmicTree, recoList, recoTrackList)
  addHist(cosmicTree, recoList, recoTrackList, cosmicList, leadingPhoton, 1)

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
    
  pionCount, protonCount = trueCutPionProton(eventTree)
  if pionCount > 0 or protonCount > 1:
    continue
  
  truePhotonIDs = truePhotonList(eventTree, fiducialData)

  if len(truePhotonIDs) == 0:
    continue

  
  #EFFICIENCY - GRAPHING BASED ON RECO

  leadingPhoton = scaleTrueEnergy(eventTree, truePhotonIDs)

  if recoNoVertex(eventTree) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effNoVertexHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #See if the event is neutral current                                                                                                  
  if recoNeutralCurrent(eventTree) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effCCHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut events with vertexes outside the fiducial
  if recoFiducials(eventTree, fiducialData) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effFiducialHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut events with too many protons
  recoProtonCount = recoProton(eventTree)
  if recoProtonCount > 1:
    addHist(eventTree, truePhotonIDs, emptyList, effProtonHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut events with pions present
  if recoPion(eventTree) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effPionHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, eventTree)
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree)

  #Try cutting for Shower from Charged Score 
  if recoCutShowerFromChargeScore(eventTree, recoList, recoTrackList) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effShowerChargeHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut based on Electron Score
  if recoCutPrimary(eventTree, recoList, recoTrackList) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effPrimaryHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Try cutting based on data for Track Lengths
  if recoCutLongTracks(eventTree, fiducialData) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effLongTrackHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut unclassified tracks too short to be trusted
  #if recoCutShortTracks(eventTree) == False:
  #  addHist(eventTree, truePhotonIDs, emptyList, effShortTrackHists, leadingPhoton, eventTree.xsecWeight)
  #  continue

  #Now we're pretty sure the event is legitimate, so we go ahead and graph based on the number of photons
  if len(recoList) + len(recoTrackList) == 0:
    addHist(eventTree, truePhotonIDs, emptyList, effNoPhotonHists, leadingPhoton, eventTree.xsecWeight)
  elif len(recoList) + len(recoTrackList) == 1:
    addHist(eventTree, truePhotonIDs, emptyList, effOnePhotonHists, leadingPhoton, eventTree.xsecWeight)
  elif len(recoList) + len(recoTrackList) == 2:
    addHist(eventTree, truePhotonIDs, emptyList, effTwoPhotonHists, leadingPhoton, eventTree.xsecWeight)
  else:
    addHist(eventTree, truePhotonIDs, emptyList, effManyPhotonHists, leadingPhoton, eventTree.xsecWeight)

#LOOPS OVER - HISTOGRAM ORGANIZING TIME

#Making our threshold histogram for the cosmics
for x in range(1, 61):
  value = np.sqrt(cosmicOnePhoton.GetBinContent(x))
  cosmicThreshold.SetBinContent(x, value)
cosmicThreshold.GetXaxis().SetTitle("Leading Photon Energy (GeV)")
cosmicThreshold.GetYaxis().SetTitle("Square Root of N")

#Stacking histograms
purityCanvas1, purityStack1, purityLegend1, purityInt1 = purityStack("Single-Photon Purity", pList1, cosmicList[0], ntuplePOTsum, cosmicPOTsum)
purityCanvas2, purityStack2, purityLegend2, purityInt2 = purityStack("Two-Photon Purity", pList2, cosmicList[1], ntuplePOTsum, cosmicPOTsum)
purityCanvas3, purityStack3, purityLegend3, purityInt3 = purityStack("3+ Photon Purity", pList3, cosmicList[2], ntuplePOTsum, cosmicPOTsum)
pTotalCanvas, pTotalStack, pTotalLegend, pTotalINt = histStack("Total purity", totalPHists, ntuplePOTsum)
effCanvas1, effStack1, effLegend1, effInt1 = histStack("Single-Photon Efficiency", effList1, ntuplePOTsum)
effCanvas2, effStack2, effLegend2, effInt2 = histStack("Two-Photon Efficiency", effList2, ntuplePOTsum)
effCanvas3, effStack3, effLegend3, effInt3 = histStack("3+ Photon Efficiency", effList3, ntuplePOTsum)

cosmicCanvas, cosmicStack, cosmicLegend, cosmicInt = histStack("Cosmic Background", cosmicList, cosmicPOTsum)

writeList = [purityCanvas1, purityCanvas2, purityCanvas3, effCanvas1, effCanvas2, effCanvas3, cosmicCanvas, cosmicThreshold]
legendList = [purityLegend1, purityLegend2, purityLegend3, pTotalLegend, effLegend1, effLegend2, effLegend3, cosmicLegend]


#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in writeList:
  canvas.Write()

print("Vertex reconstructed:", initialCount)
print("Neutral Current:", NCCount)
print("In Fiducial:", fiducialCount)
print("No pions, protons:", noPionCount)
print("Fully reconstructed:", recoCount)
print(onePhoton, "events had one photon")
print(twoPhotons, "events had two photons")
print(threePhotons, "events had 3+ photons")

print("POTsum:", ntuplePOTsum)