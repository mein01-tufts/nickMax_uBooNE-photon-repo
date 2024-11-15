import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials,trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoPion, recoProton, recoCutShowerFromChargeScore, recoCutLongTracks, recoPhotonListFiducial, recoCutPrimary, recoCutShortTracks, recoPhotonListTracks, recoCutFarShowers, trueCutMuons, trueCutElectrons, recoCutMuons, recoCutElectrons, recoCutManyTracks, recoCutCompleteness, recoCutMuonCompleteness, histStackTwoSignal

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
#1 GAMMA + 0 HISTS
puritySignal = rt.TH1F("PSignal1", "Signal",60,0,2)
purityOffSignal = rt.TH1F("POffSignal1", "Signal, One Undetected Proton",60,0,2)
purityMuon = rt.TH1F("PMuon1", "Over-threshold Muon",60,0,2)
purityElectron = rt.TH1F("PElectron1", "Over-threshold Electron",60,0,2)
purityFiducials = rt.TH1F("PFiducial1", "Out of Fiducial",60,0,2)
purityPion = rt.TH1F("PPion1", "Charged Pion",60,0,2)
purityProton = rt.TH1F("PProton1", "2+ Protons",60,0,2)
purityNoPhotons = rt.TH1F("PNoPhoton1", "No Real Photons",60,0,2)
purityTwoPhotons = rt.TH1F("PTwoPhoton1", "2 Real Photons",60,0,2)
purityManyPhotons = rt.TH1F("PMorePhoton1", "3+ Real Photons",60,0,2)

#1 GAMMA + 1P HISTS
protonPuritySignal = rt.TH1F("ProtonPSignal1", "Signal",60,0,2)
protonPurityOffSignal = rt.TH1F("ProtonPOffSignal1", "Signal, No Real Protons",60,0,2)
protonPurityMuon = rt.TH1F("ProtonPMuon1", "Over-threshold Muon",60,0,2)
protonPurityElectron = rt.TH1F("ProtonPElectron1", "Over-threshold Electron",60,0,2)
protonPurityFiducials = rt.TH1F("ProtonPFiducial1", "Out of Fiducial",60,0,2)
protonPurityPion = rt.TH1F("ProtonPPion1", "Charged Pion",60,0,2)
protonPurityProton = rt.TH1F("ProtonPProton1", "2+ Protons",60,0,2)
protonPurityNoPhotons = rt.TH1F("ProtonPNoPhoton1", "No Real Photons",60,0,2)
protonPurityTwoPhotons = rt.TH1F("ProtonPTwoPhoton1", "2 Real Photons",60,0,2)
protonPurityManyPhotons = rt.TH1F("ProtonPMorePhoton1", "3+ Real Photons",60,0,2)


#PURITY HISTLISTS
signalPHistsNoProton = [puritySignal, purityOffSignal]
signalPHistsProton = [protonPurityOffSignal, protonPuritySignal]
muonPHists = [purityMuon, protonPurityMuon]
electronPHists = [purityElectron, protonPurityElectron]
fiducialPHists = [purityFiducials, protonPurityFiducials]
pionPHists = [purityPion, protonPurityPion]
protonPHists = [purityProton, protonPurityProton]
noPhotonPHists = [purityNoPhotons, protonPurityNoPhotons]
twoPhotonPHists = [purityTwoPhotons, protonPurityTwoPhotons]
manyPhotonPHists = [purityManyPhotons, protonPurityManyPhotons]

#Cosmics go here, so we can put them on purity
cosmicOnePhoton = rt.TH1F("cBackground1", "Cosmic Background",60,0,2)
protonCosmicOnePhoton = rt.TH1F("cProtonBackground1", "Cosmic Background",60,0,2)
cosmicList = [cosmicOnePhoton, protonCosmicOnePhoton]

#Big Lists, for Big Plots
pList1 = [puritySignal, purityOffSignal, purityMuon, purityElectron, purityFiducials, purityPion, purityProton, purityNoPhotons, purityTwoPhotons, purityManyPhotons, cosmicOnePhoton]

pProtonList1 = [protonPuritySignal, protonPurityOffSignal, protonPurityMuon, protonPurityElectron, protonPurityFiducials, protonPurityPion, protonPurityProton, protonPurityNoPhotons, protonPurityTwoPhotons, protonPurityManyPhotons, protonCosmicOnePhoton]

#EFFICIENCY HISTOGRAMS
effNoVertex = rt.TH1F("effNoVertex1", "No Vertex Found",60,0,2)
effMuon = rt.TH1F("effMuon1", "Muon False Positive",60,0,2)
effElectron = rt.TH1F("effElectron1", "Electron False Positive",60,0,2)
effFiducial = rt.TH1F("effFiducial1", "Placed out of Fiducial",60,0,2)
effPion = rt.TH1F("effPion1", "Pion False Positive",60,0,2)
effProton =  rt.TH1F("effProton1", "Proton False Positive",60,0,2)
effNoPhotons = rt.TH1F("effNoPhotons1", "No Photons Found",60,0,2)
effSignal = rt.TH1F("effSignal 1", "Signal",60,0,2)
effTwoPhotons = rt.TH1F("effTwoPhotons1", "Two Photons Found",60,0,2)
effManyPhotons = rt.TH1F("effManyPhotons1", "Many Photons Found",60,0,2)
effShowerCharge = rt.TH1F("effShowerCharge1", "Shower from Charged Cut",60,0,2)
effPrimary = rt.TH1F("effPrimary1", "Primary Score Cut",60,0,2)
effLongTracks = rt.TH1F("effLongTracks1", "Tracks with length > 20 cm",60,0,2)
effMuonComp = rt.TH1F("effMuonComp1", "Muons with too-low Efficiency",60,0,2)

effProtonNoVertex = rt.TH1F("effProtonNoVertex1", "No Vertex Found",60,0,2)
effProtonMuon = rt.TH1F("effProtonMuon1", "Muon False Positive",60,0,2)
effProtonElectron = rt.TH1F("effProtonElectron1", "Electron False Positive",60,0,2)
effProtonFiducial = rt.TH1F("effProtonFiducial1", "Placed out of Fiducial",60,0,2)
effProtonPion = rt.TH1F("effProtonPion1", "Pion False Positive",60,0,2)
effProtonProton =  rt.TH1F("effProtonProton1", "Proton False Positive",60,0,2)
effProtonNoPhotons = rt.TH1F("effProtonNoPhotons1", "No Photons Found",60,0,2)
effProtonSignal = rt.TH1F("effProtonSignal 1", "Signal",60,0,2)
effProtonTwoPhotons = rt.TH1F("effProtonTwoPhotons1", "Two Photons Found",60,0,2)
effProtonManyPhotons = rt.TH1F("effProtonManyPhotons1", "Many Photons Found",60,0,2)
effProtonShowerCharge = rt.TH1F("effProtonShowerCharge1", "Shower from Charged Cut",60,0,2)
effProtonPrimary = rt.TH1F("effProtonPrimary1", "Primary Score Cut",60,0,2)
effProtonLongTracks = rt.TH1F("effProtonLongTracks1", "Tracks with length > 20 cm",60,0,2)
effProtonMuonComp = rt.TH1F("effProtonMuonComp1", "Muons with too-low Efficiency",60,0,2)


#Histogram Lists!
effNoVertexHists = [effNoVertex, effProtonNoVertex]
effMuonHists = [effMuon, effProtonMuon]
effElectronHists = [effElectron, effProtonElectron]
effFiducialHists = [effFiducial, effProtonFiducial]
effPionHists = [effPion, effProtonPion]
effProtonHists = [effProton, effProtonProton]
effNoPhotonHists = [effNoPhotons, effProtonNoPhotons]
effOnePhotonHists = [effSignal, effProtonSignal] 
effTwoPhotonHists = [effTwoPhotons, effProtonTwoPhotons]
effManyPhotonHists = [effManyPhotons, effProtonManyPhotons]
effShowerChargeHists = [effShowerCharge, effProtonShowerCharge]
effPrimaryHists = [effPrimary, effProtonPrimary]
effLongTrackHists = [effLongTracks, effProtonLongTracks]
effMuonCompHists = [effMuonComp, effProtonMuonComp]

#Lists of histograms for stacking. Place with signal first, then put the others in backwards in order of application (so the last one to apply would immediately follow the signal, then the second to last, all the way down to the first)
effList1 = [effSignal, effMuonComp, effLongTracks, effPrimary, effManyPhotons, effTwoPhotons, effNoPhotons, effShowerCharge, effProton, effPion, effElectron, effMuon, effNoVertex]
effProtonList1 = [effProtonSignal, effProtonMuonComp, effProtonLongTracks, effProtonPrimary, effProtonManyPhotons, effProtonTwoPhotons, effProtonNoPhotons, effProtonShowerCharge, effProtonProton, effProtonPion, effProtonElectron, effProtonMuon, effProtonNoVertex]


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

passTruth = 0
hasVertex = 0
hasNC = 0
inFiducial = 0
pionProtonFine = 0
survivesCuts = 0
noEffPhotons = 0
oneEffPhoton = 0
twoEffPhotons = 0
manyEffPhotons = 0

#We put this into the addHist function for truth-based graphs
emptyList = []

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

  #Cut based on the completeness of known Muons
  if recoCutMuonCompleteness(eventTree) == False:
    continue

  #PURITY - GRAPHING BASED ON TRUTH
  
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(eventTree, recoList, recoTrackList)
  
  #Cut muons and electrons
  if trueCutMuons(eventTree) == False:
    addHist(eventTree, recoProtonCount, muonPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  if trueCutElectrons(eventTree) == False:
    addHist(eventTree, recoProtonCount, electronPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  NCCount += 1

  if trueCutFiducials(eventTree, fiducialData) == False:
    addHist(eventTree, recoProtonCount, fiducialPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  fiducialCount += 1

  #pions and protons!
  pionCount, protonCount = trueCutPionProton(eventTree)
  if pionCount > 0:
    addHist(eventTree, recoProtonCount, pionPHists, leadingPhoton, eventTree.xsecWeight)
    continue
    
  elif protonCount > 1:
    addHist(eventTree, recoProtonCount, protonPHists, leadingPhoton, eventTree.xsecWeight)
    continue

  noPionCount += 1

  #Now we make a list of the actual photons!
  truePhotonIDs = truePhotonList(eventTree, fiducialData)

  #Are there actually any photons?
  if len(truePhotonIDs) == 0:
    addHist(eventTree, recoProtonCount, noPhotonPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  recoCount += 1
  #Is there one Photon?
  if len(truePhotonIDs) == 1:
    if protonCount == 0:
      addHist(eventTree, recoProtonCount, signalPHistsNoProton, leadingPhoton, eventTree.xsecWeight)
    elif protonCount == 1:
      addHist(eventTree, recoProtonCount, signalPHistsProton, leadingPhoton, eventTree.xsecWeight)
    onePhoton += 1
    continue
  #Are there two?
  elif len(truePhotonIDs) == 2:
    addHist(eventTree, recoProtonCount, twoPhotonPHists, leadingPhoton, eventTree.xsecWeight)
    twoPhotons += 1
    continue
  
  #In that case, there should be at least three
  else:
    addHist(eventTree, recoProtonCount, manyPhotonPHists, leadingPhoton, eventTree.xsecWeight)
    threePhotons += 1
    continue

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
  if recoCutMuons(cosmicTree) == False:
    continue

  if recoCutElectrons(cosmicTree) == False:
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
  #Cut events with too many protons
  recoProtonCount = recoProton(cosmicTree)
  if recoProtonCount > 1:
    continue

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, cosmicTree)

  #List all photons classified as tracks
  recoTrackList = recoPhotonListTracks(fiducialData, cosmicTree)
  
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

  #Cut based on the completeness of known Muons
  if recoCutMuonCompleteness(eventTree) == False:
    continue

  untrackedCosmics += 1

  #graphing based on photon count
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(cosmicTree, recoList, recoTrackList)
  addHist(cosmicTree, recoProtonCount, cosmicList, leadingPhoton, 1)

#BEGINNING EVENT LOOP FOR EFFICIENCY
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #Selecting events using truth
  if trueCutMuons(eventTree) == False:
    continue

  if trueCutElectrons(eventTree) == False:
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
  passTruth += 1

  if recoNoVertex(eventTree) == False:
    addHist(eventTree, protonCount, effNoVertexHists, leadingPhoton, eventTree.xsecWeight)
    continue
  hasVertex += 1
  #See if the event is neutral current                                                                                                  
  #if recoNeutralCurrent(eventTree) == False:
  #  addHist(eventTree, truePhotonIDs, emptyList, effCCHists, leadingPhoton, eventTree.xsecWeight)
  #  continue

  #Check for above-threshold muons and electrons 
  if recoCutMuons(eventTree, classificationThreshold) == False:
    addHist(eventTree, protonCount, effMuonHists, leadingPhoton, eventTree.xsecWeight)
    continue

  if recoCutElectrons(eventTree, classificationThreshold) == False:
    addHist(eventTree, protonCount, effElectronHists, leadingPhoton, eventTree.xsecWeight)
    continue

  hasNC += 1
  #Cut events with vertexes outside the fiducial
  if recoFiducials(eventTree, fiducialData) == False:
    addHist(eventTree, protonCount, effFiducialHists, leadingPhoton, eventTree.xsecWeight)
    continue
  inFiducial += 1
  #Cut events with too many protons
  recoProtonCount = recoProton(eventTree, classificationThreshold)
  if recoProtonCount > 1:
    addHist(eventTree, protonCount, effProtonHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut events with pions present
  if recoPion(eventTree, classificationThreshold) == False:
    addHist(eventTree, protonCount, effPionHists, leadingPhoton, eventTree.xsecWeight)
    continue
  pionProtonFine += 1

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, eventTree, classificationThreshold)
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree, classificationThreshold)

  #Try cutting for Shower from Charged Score 
  if recoCutShowerFromChargeScore(eventTree, recoList, recoTrackList) == False:
    addHist(eventTree, protonCount, effShowerChargeHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut based on Primary Score
  if recoCutPrimary(eventTree, recoList, recoTrackList) == False:
    addHist(eventTree, protonCount, effPrimaryHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Try cutting based on data for Track Lengths
  if recoCutLongTracks(eventTree, fiducialData) == False:
    addHist(eventTree, protonCount, effLongTrackHists, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut based on the completeness of known Muons
  if recoCutMuonCompleteness(eventTree) == False:
    addHist(eventTree, protonCount, effMuonCompHists, leadingPhoton, eventTree.xsecWeight)
    continue

  survivesCuts += 1

  #Now we're pretty sure the event is legitimate, so we go ahead and graph based on the number of photons
  if len(recoList) + len(recoTrackList) == 0:
    addHist(eventTree, protonCount, effNoPhotonHists, leadingPhoton, eventTree.xsecWeight)
    noEffPhotons += 1
  elif len(recoList) + len(recoTrackList) == 1:
    addHist(eventTree, protonCount, effOnePhotonHists, leadingPhoton, eventTree.xsecWeight)
    oneEffPhoton += 1
  elif len(recoList) + len(recoTrackList) == 2:
    addHist(eventTree, protonCount, effTwoPhotonHists, leadingPhoton, eventTree.xsecWeight)
    twoEffPhotons += 1
  else:
    addHist(eventTree, protonCount, effManyPhotonHists, leadingPhoton, eventTree.xsecWeight)
    manyEffPhotons += 1

#LOOPS OVER - HISTOGRAM ORGANIZING TIME

#Scaling the Cosmic Histograms
for hist in cosmicList:
  hist.Scale(ntuplePOTsum/cosmicPOTsum)


#Stacking histograms
purityCanvas1, purityStack1, purityLegend1, purityInt1 = histStackTwoSignal("1 Gamma + 0 Sample", pList1, ntuplePOTsum)
purityProtonCanvas1, purityProtonStack1, purityProtonLegend1, purityProtonInt1 = histStackTwoSignal("1 Gamma + 1P Sample", pProtonList1, ntuplePOTsum)

effCanvas1, effStack1, effLegend1, effInt1 = histStack("True 1 Gamma + 0  Outcomes", effList1, ntuplePOTsum)
effProtonCanvas1, effProtonStack1, effProtonLegend1, effProtonInt1 = histStack("True 1 Gamma + 1P  Outcomes", effProtonList1, ntuplePOTsum)

writeList = [purityCanvas1, purityProtonCanvas1, effCanvas1, effProtonCanvas1]

for stack in [purityStack1, purityProtonStack1]:
  stack.GetXaxis().SetTitle("Reconstructed Leading Photon Energy (GeV)")

for stack in [effStack1, effProtonStack1]:
  stack.GetXaxis().SetTitle("True LeadingPhoton Energy (GeV)")

#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in writeList:
  canvas.Write()

print("PURITY STATS:")
print("Vertex reconstructed:", initialCount)
print("Neutral Current:", NCCount)
print("In Fiducial:", fiducialCount)
print("No pions, protons:", noPionCount)
print("Fully reconstructed:", recoCount)
print(onePhoton, "events had one photon")
print(twoPhotons, "events had two photons")
print(threePhotons, "events had 3+ photons")

print("EFFICIENCY STATS:")
print("Total signal space:", passTruth)
print("Total with vertex:", hasVertex)
print("Total NC:", hasNC)
print("Total reconstructed in Fiducial:", inFiducial)
print("Total with acceptable pion/protons:", pionProtonFine)
print("Total that survive our cuts:", survivesCuts)
print("Total with no photons:", noEffPhotons)
print("Total with one photon:", oneEffPhoton)
print("Total with two photons:", twoEffPhotons)
print("Total with three photons:", manyEffPhotons)


print("POTsum:", ntuplePOTsum)