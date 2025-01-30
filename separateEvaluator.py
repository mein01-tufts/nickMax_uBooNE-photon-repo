import sys, argparse
import numpy as np
import ROOT as rt
from math import sqrt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials,trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoPion, recoProton, recoCutShowerFromChargeScore, recoCutLongTracks, recoPhotonListFiducial, recoCutPrimary, recoCutShortTracks, recoPhotonListTracks, recoCutFarShowers, trueCutMuons, trueCutElectrons, recoCutMuons, recoCutElectrons, recoCutManyTracks, recoCutCompleteness, recoCutMuonCompleteness, histStackTwoSignal, recoCutMaxInTime, trueTwoPhotonOpeningAngle

from helpers.larflowreco_ana_funcs import getCosThetaGravVector
from selection_1g1p import run_1g1p_reco_selection_cuts
from truthdef import truthdef_1gamma_cuts

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-c", "--cosmicFile", type=str, required=True, help="cosmic input file")
parser.add_argument("-b", "--beamFile",type=str,default="__none__",help="[optional] beam data input file")
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
#cosmicPotTree = ntuple_file.Get("cosmicPotTree")

beamTree = None
nBeamEntries = 0
if args.beamFile != "__none__":
  beam_file = rt.TFile(args.beamFile)
  beamTree = beam_file.Get("EventTree")
  nBeamEntries = beamTree.GetEntries()

#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 4.4e+19
targetPOTstring = "4.4e+19"
ntuplePOTsum = 0

for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

#For the old cosmic file
#cosmicPOTsum = 3.2974607516739725e+20
#cosmicWeight = 1

#For the new cosmic file
#cosmicPOTsum = 5.28e+19
#cosmicWeight = 0.4
cosmicPOTsum = 1.1e20
cosmicWeight = 1.0

# ROOT FILE TO SAVE Histograms
outrootfile = rt.TFile("outhist_separateEvaluator.root","recreate")

cosmic_event_list = open("cosmic_event_list.txt",'w')

#Hists created and organized here
#PURITY HISTOGRAMS
#1 GAMMA + 0 HISTS

histNbins = 20
histXmax = 2.0

# intimepixelsummax
#histNbins = 50
#histXmax = 1000.0

#histNbins = 25
#histXmax = 1.0

puritySignal = rt.TH1F("PSignal1", "Signal",histNbins,0,histXmax)
purityOffSignal = rt.TH1F("POffSignal1", "Signal, One Undetected Proton",histNbins,0,histXmax)
purityMuon = rt.TH1F("PMuon1", "Over-threshold Muon",histNbins,0,histXmax)
purityElectron = rt.TH1F("PElectron1", "Over-threshold Electron",histNbins,0,histXmax)
purityFiducials = rt.TH1F("PFiducial1", "Out of Fiducial (1g)",histNbins,0,histXmax)
purityFiducials2g = rt.TH1F("PFiducials2g1", "Out of Fiducial (2g)",histNbins,0,histXmax)
purityPion = rt.TH1F("PPion1", "Charged Pion",histNbins,0,histXmax)
purityProton = rt.TH1F("PProton1", "2+ Protons",histNbins,0,histXmax)
purityNoPhotons = rt.TH1F("PNoPhoton1", "No Real Photons",histNbins,0,histXmax)
purityTwoPhotons = rt.TH1F("PTwoPhoton1", "2 Real Photons",histNbins,0,histXmax)
purityManyPhotons = rt.TH1F("PMorePhoton1", "3+ Real Photons",histNbins,0,histXmax)

#1 GAMMA + 1P HISTS
protonPuritySignal = rt.TH1F("ProtonPSignal1", "Signal",histNbins,0,histXmax)
protonPurityOffSignal = rt.TH1F("ProtonPOffSignal1", "Signal, No Real Protons",histNbins,0,histXmax)
protonPurityMuon = rt.TH1F("ProtonPMuon1", "Over-threshold Muon",histNbins,0,histXmax)
protonPurityElectron = rt.TH1F("ProtonPElectron1", "Over-threshold Electron",histNbins,0,histXmax)
protonPurityFiducials = rt.TH1F("ProtonPFiducial1", "Out of Fiducial (1g)",histNbins,0,histXmax)
protonPurityFiducials2g = rt.TH1F("ProtonPFiducial2g1", "Out of Fiducial (2g)",histNbins,0,histXmax)
protonPurityPion = rt.TH1F("ProtonPPion1", "Charged Pion",histNbins,0,histXmax)
protonPurityProton = rt.TH1F("ProtonPProton1", "2+ Protons",histNbins,0,histXmax)
protonPurityNoPhotons = rt.TH1F("ProtonPNoPhoton1", "No Real Photons",histNbins,0,histXmax)
protonPurityTwoPhotons = rt.TH1F("ProtonPTwoPhoton1", "2 Real Photons",histNbins,0,histXmax)
protonPurityManyPhotons = rt.TH1F("ProtonPMorePhoton1", "3+ Real Photons",histNbins,0,histXmax)


#PURITY HISTLISTS
signalPHistsNoProton = [puritySignal, purityOffSignal]
signalPHistsProton = [protonPurityOffSignal, protonPuritySignal]
muonPHists = [purityMuon, protonPurityMuon]
electronPHists = [purityElectron, protonPurityElectron]
fiducialPHists   = [purityFiducials,protonPurityFiducials]
fiducialPHists2g = [purityFiducials2g,protonPurityFiducials2g]
pionPHists = [purityPion, protonPurityPion]
protonPHists = [purityProton, protonPurityProton]
noPhotonPHists = [purityNoPhotons, protonPurityNoPhotons]
twoPhotonPHists = [purityTwoPhotons, protonPurityTwoPhotons]
manyPhotonPHists = [purityManyPhotons, protonPurityManyPhotons]

#Cosmics go here, so we can put them on purity
cosmicOnePhoton = rt.TH1F("cBackground1", "Cosmic Background",histNbins,0,histXmax)
protonCosmicOnePhoton = rt.TH1F("cProtonBackground1", "Cosmic Background",histNbins,0,histXmax)
cosmicList = [cosmicOnePhoton, protonCosmicOnePhoton]

#Beam Data
if beamTree is not None:
  beamOnePhoton = rt.TH1F("beam1g0X1", "Cosmic Background",histNbins,0,histXmax)
  protonBeamOnePhoton = rt.TH1F("beam1g1p1", "Cosmic Background",histNbins,0,histXmax)
else:
  beamOnePhoton = None
  protonBeamOnePhoton = None
beamList = [beamOnePhoton,protonBeamOnePhoton]

#Big Lists, for Big Plots
pList1 = [puritySignal,
          purityOffSignal,
          purityFiducials,
          purityFiducials2g,
          purityMuon,
          purityElectron,
          purityPion,
          purityProton,
          purityNoPhotons,
          purityTwoPhotons,
          purityManyPhotons,
          cosmicOnePhoton]

pProtonList1 = [protonPuritySignal,
                protonPurityOffSignal,
                protonPurityFiducials,
                protonPurityFiducials2g,                
                protonPurityMuon,
                protonPurityElectron,
                protonPurityPion,
                protonPurityProton,
                protonPurityNoPhotons,
                protonPurityTwoPhotons,
                protonPurityManyPhotons,
                protonCosmicOnePhoton]

#EFFICIENCY HISTOGRAMS
effNoVertex = rt.TH1F("effNoVertex1", "No Vertex Found",histNbins,0,histXmax)
effMuon = rt.TH1F("effMuon1", "Muon False Positive",histNbins,0,histXmax)
effElectron = rt.TH1F("effElectron1", "Electron False Positive",histNbins,0,histXmax)
effCosmicPixel = rt.TH1F("effCosmicPixel1", "No cosmic-tagged pixels",histNbins,0,histXmax)
effFiducial = rt.TH1F("effFiducial1", "Placed out of Fiducial",histNbins,0,histXmax)
effPion = rt.TH1F("effPion1", "Pion False Positive",histNbins,0,histXmax)
effProton =  rt.TH1F("effProton1", "Proton False Positive",histNbins,0,histXmax)
effNoPhotons = rt.TH1F("effNoPhotons1", "No Photons Found",histNbins,0,histXmax)
effSignal = rt.TH1F("effSignal 1", "Signal",histNbins,0,histXmax)
effTwoPhotons = rt.TH1F("effTwoPhotons1", "Two Photons Found",histNbins,0,histXmax)
effManyPhotons = rt.TH1F("effManyPhotons1", "Many Photons Found",histNbins,0,histXmax)
effShowerCharge = rt.TH1F("effShowerCharge1", "Shower from Charged Cut",histNbins,0,histXmax)
effPrimary = rt.TH1F("effPrimary1", "Primary Score Cut",histNbins,0,histXmax)
effLongTracks = rt.TH1F("effLongTracks1", "Tracks with length > 20 cm",histNbins,0,histXmax)
effMuonComp = rt.TH1F("effMuonComp1", "Muons with too-low Efficiency",histNbins,0,histXmax)
effMaxInTime = rt.TH1F("effMaxInTime1", "Likely Incorrect Vertex",histNbins,0,histXmax)
effShowerComp = rt.TH1F("effShowerComp1", "Shower Completeness",histNbins,0,histXmax)

effProtonNoVertex = rt.TH1F("effProtonNoVertex1", "No Vertex Found",histNbins,0,histXmax)
effProtonMuon = rt.TH1F("effProtonMuon1", "Muon False Positive",histNbins,0,histXmax)
effProtonElectron = rt.TH1F("effProtonElectron1", "Electron False Positive",histNbins,0,histXmax)
effProtonCosmicPixel = rt.TH1F("effProtonCosmicPixel1", "No cosmic-tagged pixels",histNbins,0,histXmax)
effProtonFiducial = rt.TH1F("effProtonFiducial1", "Placed out of Fiducial",histNbins,0,histXmax)
effProtonPion = rt.TH1F("effProtonPion1", "Pion False Positive",histNbins,0,histXmax)
effProtonProton =  rt.TH1F("effProtonProton1", "Proton False Positive",histNbins,0,histXmax)
effProtonNoPhotons = rt.TH1F("effProtonNoPhotons1", "No Photons Found",histNbins,0,histXmax)
effProtonSignal = rt.TH1F("effProtonSignal 1", "Signal",histNbins,0,histXmax)
effProtonTwoPhotons = rt.TH1F("effProtonTwoPhotons1", "Two Photons Found",histNbins,0,histXmax)
effProtonManyPhotons = rt.TH1F("effProtonManyPhotons1", "Many Photons Found",histNbins,0,histXmax)
effProtonShowerCharge = rt.TH1F("effProtonShowerCharge1", "Shower from Charged Cut",histNbins,0,histXmax)
effProtonPrimary = rt.TH1F("effProtonPrimary1", "Primary Score Cut",histNbins,0,histXmax)
effProtonLongTracks = rt.TH1F("effProtonLongTracks1", "Tracks with length > 20 cm",histNbins,0,histXmax)
effProtonMuonComp = rt.TH1F("effProtonMuonComp1", "Muons with too-low Efficiency",histNbins,0,histXmax)
effProtonMaxInTime = rt.TH1F("effProtonMaxInTime1", "Likely Incorrect Vertex",histNbins,0,histXmax)
effProtonShowerComp = rt.TH1F("effProtonShowerComp1", "Shower Completeness",histNbins,0,histXmax)

#Histogram Lists!
effNoVertexHists = [effNoVertex, effProtonNoVertex]
effMuonHists = [effMuon, effProtonMuon]
effElectronHists = [effElectron, effProtonElectron]
effCosmicPixelHists = [effCosmicPixel, effProtonCosmicPixel]
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
effMaxInTimeHists = [effMaxInTime, effProtonMaxInTime]
effShowerCompHists = [effShowerComp,effProtonShowerComp]

#Lists of histograms for stacking. Place with signal first, then put the others in backwards in order of application (so the last one to apply would immediately follow the signal, then the second to last, all the way down to the first)
effList1 = [effSignal,
            effShowerComp,
            effMaxInTime,
            effMuonComp,
            effLongTracks,
            effPrimary,
            effShowerCharge,            
            effManyPhotons,
            effTwoPhotons,
            effNoPhotons,
            effProton,
            effPion,
            effFiducial,
            effCosmicPixel,
            effElectron,
            effMuon,
            effNoVertex]
effProtonList1 = [effProtonSignal,
                  effProtonShowerComp,
                  effProtonMaxInTime,
                  effProtonMuonComp,
                  effProtonLongTracks,
                  effProtonPrimary,
                  effProtonShowerCharge,                  
                  effProtonManyPhotons,
                  effProtonTwoPhotons,
                  effProtonNoPhotons,
                  effProtonProton,
                  effProtonPion,
                  effProtonFiducial,
                  effProtonCosmicPixel,
                  effProtonElectron,
                  effProtonMuon,
                  effProtonNoVertex]


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
noVertexCosmics = 0
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

distList = []

#We put this into the addHist function for truth-based graphs
emptyList = []

#Variables for program function
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":15, "photonWidth":3}
classificationThreshold = 0
showerRecoThreshold = 10.0
photonEDepThreshold = 10.0

#BEGINNING EVENT LOOP FOR DEFAULT PURITY
for i in range(eventTree.GetEntries()):

  if False:
    break
  
  if i>0 and i%10000==0:
    print("eventTree Entry[",i,"]")
  
  eventTree.GetEntry(i)

  # #Selecting events with reco
  # #See if the event has a vertex
  # if recoNoVertex(eventTree) == False:
  #   continue

  # #See if the event is neutral current
  # if recoCutMuons(eventTree, classificationThreshold) == False:
  #   continue

  # if recoCutElectrons(eventTree, classificationThreshold) == False:
  #   continue

  # #Use Matt's Cosmic Cut
  # if trueCutCosmic(eventTree) == False:
  #   continue

  # #Make sure the event is within the fiducial volume
  # if recoFiducials(eventTree, fiducialData) == False:
  #   continue

  # #Cut events with suitably energetic charged pions 
  # if recoPion(eventTree, classificationThreshold) == False:
  #   continue

  # #Cut events with far-travelling protons
  # recoProtonCount = recoProton(eventTree, classificationThreshold)
  # if recoProtonCount > 1:
  #   continue

  # #See if there are any photons in the event - if so, list them
  # recoList = recoPhotonListFiducial(fiducialData, eventTree, classificationThreshold)

  # #List all photons classified as tracks
  # #recoTrackList = recoPhotonListTracks(fiducialData, eventTree, classificationThreshold)
  # recoTrackList = []

  # if len(recoList) + len(recoTrackList) != 1:
  #   continue

  # #Try cutting based on data for Shower from Charged
  # if recoCutShowerFromChargeScore(eventTree, recoList, recoTrackList) == False:
  #   continue
  
  # #Cut based on primary score
  # if recoCutPrimary(eventTree, recoList, recoTrackList) == False:
  #   continue

  # #Cut based on the presence of tracks over 20 cm
  # if recoCutLongTracks(eventTree, fiducialData) == False:
  #   continue

  # #Cut based on the completeness of known Muons
  # if recoCutMuonCompleteness(eventTree) == False:
  #   continue

  # if recoCutMaxInTime(eventTree, recoProtonCount) == False:    
  #   continue

  eventPassAllCuts, cuts_passed, recoList, recoTrackList, recoProtonCount = run_1g1p_reco_selection_cuts( eventTree, classificationThreshold, fiducialData )
  if eventPassAllCuts==False:
    #print("entry[",i,"] passes original cuts, but does not pass selection function")
    #for cutname in cuts_passed:
    #  print("  [",cutname,"] ",cuts_passed[cutname])
    continue

  #PURITY - GRAPHING BASED ON TRUTH
  
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(eventTree, recoList, recoTrackList)
  #leadingPhoton = eventTree.vtxMaxIntimePixelSum*0.0126 # hack
  #for recoIDX in recoList+recoTrackList:
  #  leadingPhoton = eventTree.showerComp[recoIDX]

  vertexDistX = eventTree.vtxX - eventTree.trueVtxX
  vertexDistY = eventTree.vtxY - eventTree.trueVtxY
  vertexDistZ = eventTree.vtxZ - eventTree.trueVtxZ
  vertexDist = np.sqrt(vertexDistX**2 + vertexDistY**2 + vertexDistZ**2)
  #if vertexDist < 3:
  #  continue

  #Now we make a list of the actual photons!
  truePhotonIDs = truePhotonList(eventTree, fiducialData, threshold=photonEDepThreshold)
  nEdepPhotons = len(truePhotonIDs)

  missingPhotonDist2vtx = -1.0
  missingPhotonEDep = -1.0
  missingPhotonE = -1.0
  for truePhotonIDX in truePhotonIDs:
    truePhotonTID = eventTree.trueSimPartTID[truePhotonIDX]
    idfound = False
    for recoIDX in recoList:
      if eventTree.showerTrueTID[recoIDX]==truePhotonTID:
        idfound = True
    for recoIDX in recoTrackList:
      if eventTree.trackTrueTID[recoIDX]==truePhotonTID:
        idfound = True
    if not idfound:
      # we missed this true photon
      dist2vtx = 0.0
      d2vx=eventTree.trueSimPartEDepX[truePhotonIDX]-eventTree.trueVtxX
      d2vy=eventTree.trueSimPartEDepY[truePhotonIDX]-eventTree.trueVtxY
      d2vz=eventTree.trueSimPartEDepZ[truePhotonIDX]-eventTree.trueVtxZ
      dist2vtx = np.sqrt( d2vx*d2vx + d2vy*d2vy + d2vz*d2vz )
      if missingPhotonDist2vtx<0.0 or dist2vtx < missingPhotonDist2vtx:
        missingPhotonDist2vtx = dist2vtx

      photonedep = eventTree.trueSimPartPixelSumYplane[truePhotonIDX]*0.0126
      if missingPhotonEDep<0.0 or photonedep<missingPhotonEDep:
        missingPhotonEDep = photonedep
      
  #leadingPhoton = missingPhotonEDep

  trueOpeningAngle = 0.0
  if nEdepPhotons==2:    
    id1 = truePhotonIDs[0]
    id2 = truePhotonIDs[1]
    trueOpeningAngle = trueTwoPhotonOpeningAngle( eventTree, truePhotonIDs[0], truePhotonIDs[1] )
    if trueOpeningAngle<20.0:
      nEdepPhotons=1
  #leadingPhoton = trueOpeningAngle

  # out of fiducial vertex
  if trueCutFiducials(eventTree, fiducialData) == False:
    if nEdepPhotons==1:
      addHist(eventTree, recoProtonCount, fiducialPHists, leadingPhoton, eventTree.xsecWeight)
    else:
      addHist(eventTree, recoProtonCount, fiducialPHists2g, leadingPhoton, eventTree.xsecWeight)
    continue
  fiducialCount += 1
  
  #Cut muons and electrons
  if trueCutMuons(eventTree) == False:
    addHist(eventTree, recoProtonCount, muonPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  if trueCutElectrons(eventTree) == False:
    addHist(eventTree, recoProtonCount, electronPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  NCCount += 1

  #pions and protons!
  pionCount, protonCount = trueCutPionProton(eventTree)
  if pionCount > 0:
    addHist(eventTree, recoProtonCount, pionPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  if protonCount > 1:
    addHist(eventTree, recoProtonCount, protonPHists, leadingPhoton, eventTree.xsecWeight)
    continue

  noPionCount += 1


  #Are there actually any photons?
  if nEdepPhotons==0:
    addHist(eventTree, recoProtonCount, noPhotonPHists, leadingPhoton, eventTree.xsecWeight)
    continue
  
  recoCount += 1
  #Is there one Photon?
  if nEdepPhotons==1:
    if protonCount == 0:
      addHist(eventTree, recoProtonCount, signalPHistsNoProton, leadingPhoton, eventTree.xsecWeight)
    elif protonCount == 1:
      addHist(eventTree, recoProtonCount, signalPHistsProton, leadingPhoton, eventTree.xsecWeight)
    onePhoton += 1
    continue
  #Are there two?
  elif nEdepPhotons == 2:
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
  if i>0 and i%10000==0:
    print("cosmicTree Entry[",i,"]")
  cosmicTree.GetEntry(i)

  cosmicEventPassAllCuts, cosmic_cuts_passed, recoList, recoTrackList, recoProtonCount = run_1g1p_reco_selection_cuts( cosmicTree,
                                                                                                                       classificationThreshold,
                                                                                                                       fiducialData )
  untrackedCosmics += 1

  if not cosmicEventPassAllCuts:
    continue

  #graphing based on photon count
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(cosmicTree, recoList, recoTrackList)
  #leadingPhoton = cosmicTree.vtxMaxIntimePixelSum*0.0126 # hack
  #for recoIDX in recoList+recoTrackList:
  #  leadingPhoton = cosmicTree.showerComp[recoIDX]  
  addHist(cosmicTree, recoProtonCount, cosmicList, leadingPhoton, 1.0)

  print(cosmicTree.run," ",cosmicTree.subrun," ",cosmicTree.event," ",cosmicTree.fileid, file=cosmic_event_list)


#BEGINNING EVENT LOOP FOR COSMICS
for i in range(nBeamEntries):
  if i>0 and i%10000==0:
    print("beamTree Entry[",i,"]")
  beamTree.GetEntry(i)

  beamPass, beam_cutspassed, recoList, recoTrackList, recoProtonCount = run_1g1p_reco_selection_cuts( beamTree, classificationThreshold, fiducialData )
  if not beamPass:
    continue

  #graphing based on photon count
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(beamTree, recoList, recoTrackList)
  #leadingPhoton = beamTree.vtxMaxIntimePixelSum*0.0126 # hack
  #for recoIDX in recoList+recoTrackList:
  #  leadingPhoton = beamTree.showerComp[recoIDX]  
  addHist(beamTree, recoProtonCount, beamList, leadingPhoton, 1.0)
  
#BEGINNING EVENT LOOP FOR EFFICIENCY
for i in range(eventTree.GetEntries()):

  if i>0 and i%10000==0:
    print("eventTree Entry[",i,"]")
  
  eventTree.GetEntry(i)

  # #Selecting events using truth
  # if trueCutFiducials(eventTree, fiducialData) == False:
  #   continue
  
  # if trueCutMuons(eventTree) == False:
  #   continue

  # if trueCutElectrons(eventTree) == False:
  #   continue

  # #if trueCutCosmic(eventTree) == False:
  # #  continue

  # pionCount, protonCount = trueCutPionProton(eventTree)
  # if pionCount > 0:
  #   continue
  # if protonCount > 1:
  #   continue

  # truePhotonIDs = truePhotonList(eventTree, fiducialData, threshold=photonEDepThreshold)
  # nEdepPhotons = len(truePhotonIDs)
  # trueOpeningAngle = 0.0
  # if nEdepPhotons==2:
  #   # WC inclusive opening angle acceptance
  #   id1 = truePhotonIDs[0]
  #   id2 = truePhotonIDs[1]
  #   trueOpeningAngle = trueTwoPhotonOpeningAngle( eventTree, truePhotonIDs[0], truePhotonIDs[1] )
  #   if trueOpeningAngle<20.0:
  #     nEdepPhotons=1

  # if nEdepPhotons != 1:
  #   continue
  passes_truthdef, truthcuts = truthdef_1gamma_cuts( eventTree, photonEDepThreshold, fiducialData, return_on_fail=True )
  if not passes_truthdef:
    continue

  nEdepPhotons = truthcuts["NEdepPhotons"]
  trueOpeningAngle = truthcuts["trueOpeningAngle"]
  truePhotonIDs = truthcuts["truePhotonIDs"]
  protonCount = truthcuts["protonCount"]
  
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

  #Use Matt's Cosmic Cut
  if trueCutCosmic(eventTree) == False:
    addHist(eventTree, protonCount, effCosmicPixelHists, leadingPhoton, eventTree.xsecWeight)
  
  #Cut events with vertexes outside the fiducial
  if recoFiducials(eventTree, fiducialData) == False:
    addHist(eventTree, protonCount, effFiducialHists, leadingPhoton, eventTree.xsecWeight)
    continue
  inFiducial += 1

  #Cut events with pions present
  if recoPion(eventTree, classificationThreshold) == False:
    addHist(eventTree, protonCount, effPionHists, leadingPhoton, eventTree.xsecWeight)
    continue
  pionProtonFine += 1
  
  #Cut events with too many protons
  recoProtonCount = recoProton(eventTree, classificationThreshold)
  if recoProtonCount > 1:
    addHist(eventTree, protonCount, effProtonHists, leadingPhoton, eventTree.xsecWeight)
    continue


  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, eventTree, showerRecoThreshold)
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree, classificationThreshold)
  #recoTrackList = []

  # only 1 photon reco cut happens here
  if len(recoList) + len(recoTrackList) != 1:
    if len(recoList) + len(recoTrackList) == 0:
      addHist(eventTree, protonCount, effNoPhotonHists, leadingPhoton, eventTree.xsecWeight)
      noEffPhotons += 1    
      continue
    elif len(recoList) + len(recoTrackList) == 2:
      addHist(eventTree, protonCount, effTwoPhotonHists, leadingPhoton, eventTree.xsecWeight)
      twoEffPhotons += 1
      continue
    else:
      addHist(eventTree, protonCount, effManyPhotonHists, leadingPhoton, eventTree.xsecWeight)
      manyEffPhotons += 1
      continue
    
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

  # meant to reduce vertices selected when there is a clear in-time muon
  if recoCutMaxInTime(eventTree, protonCount) == False:
    addHist(eventTree, protonCount, effMaxInTimeHists, leadingPhoton, eventTree.xsecWeight)
    continue

  if recoCutCompleteness(eventTree, recoList, recoTrackList)==False:    
    addHist(eventTree, protonCount, effShowerCompHists, leadingPhoton, eventTree.xsecWeight)
    continue

  survivesCuts += 1

  #Now we're pretty sure the event is legitimate, so we go ahead and graph based on the number of photons
  addHist(eventTree, protonCount, effOnePhotonHists, leadingPhoton, eventTree.xsecWeight)
  oneEffPhoton += 1

#LOOPS OVER - HISTOGRAM ORGANIZING TIME

#Scaling the Cosmic Histograms
for hist in cosmicList:
  hist.Scale((targetPOT/cosmicPOTsum)*(ntuplePOTsum/targetPOT)) # the second factor is there to get canceled in histstacktwosignal


#Stacking histograms
# Prediction with purity categories
purityCanvas1, purityStack1, purityLegend1, purityInt1 = histStackTwoSignal("1 Gamma + 0 Sample", pList1,
                                                                            ntuplePOTsum, targetPOT, beamOnePhoton)
purityProtonCanvas1, purityProtonStack1, purityProtonLegend1, purityProtonInt1 = histStackTwoSignal("1 Gamma + 1P Sample", pProtonList1,
                                                                                                    ntuplePOTsum, targetPOT, protonBeamOnePhoton)

# efficiency
effCanvas1, effStack1, effLegend1, effInt1 = histStack("TrueOutcomes1", "True 1 Gamma + 0  Outcomes", effList1, ntuplePOTsum)
effProtonCanvas1, effProtonStack1, effProtonLegend1, effProtonInt1 = histStack("TrueOutcomes1P", "True 1 Gamma + 1P  Outcomes", effProtonList1, ntuplePOTsum)

writeList = [purityCanvas1, purityProtonCanvas1, effCanvas1, effProtonCanvas1]

for stack in [purityStack1, purityProtonStack1]:
  stack.GetXaxis().SetTitle("Reconstructed Leading Photon Energy (GeV)")
  print(stack.GetName(),": ",stack.GetMaximum())

for stack in [effStack1, effProtonStack1]:
  stack.GetXaxis().SetTitle("True LeadingPhoton Energy (GeV)")

for stack in [effStack1, effProtonStack1, purityStack1, purityProtonStack1]:
  stack.GetYaxis().SetTitle("Events per 4.4e19 POT")

#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
for n,canvas in enumerate(writeList):
  canvas.cd()
  if n==0:
    beamOnePhoton.Draw("E1same")
  elif n==1:
    protonBeamOnePhoton.Draw("E1same")
  else:
    pass
  canvas.Update()
  canvas.Write()

outrootfile.Write()
  
print("PURITY STATS:")
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

print("There were", totalCosmics, "total cosmics")
print(vertexCosmics, "had a reconstructed vertex")
print(noVertexCosmics, "had no reconstructed vertex")


#print("Average vertex distance:", sum(distList)/len(distList))
