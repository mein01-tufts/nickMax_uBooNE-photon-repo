import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials, trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoCutLowEnergy, recoPion, recoProton, CCSeeker

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
purityTotal1 = rt.TH1F("PTotal1", "One Photon",60,0,10)
puritySignal1 =  rt.TH1F("PSignal1", "Signal",60,0,10)
purityCC1 = rt.TH1F("PCC1", "Actually Charged Current",60,0,10)
purityFiducials1 = rt.TH1F("PFiducial1", "Out of Fiducial",60,0,10)
purityCosmic1 = rt.TH1F("PCosmic1", "Failed Cosmic",60,0,10)
purityPionProton1 = rt.TH1F("PPionProton1", "Charged Pion or Proton",60,0,10)
purityNoPhotons1 = rt.TH1F("PNoPhoton1", "No Real Photons",60,0,10)
purityTwoPhotons1 = rt.TH1F("PTwoPhoton1", "2 Real Photons",60,0,10)
purityManyPhotons1 = rt.TH1F("PMorePhoton1", "3+ Real Photons",60,0,10)

purityTotal2 = rt.TH1F("PTotal2", "One Photon",60,0,10)
puritySignal2 =	 rt.TH1F("PSignal2", "Signal",60,0,10)
purityCC2 = rt.TH1F("PCC2", "Actually Charged Current",60,0,10)
purityFiducials2 = rt.TH1F("PFiducial2", "Out of Fiducial",60,0,10)
purityCosmic2 = rt.TH1F("PCosmic2", "Failed Cosmic",60,0,10)
purityPionProton2 = rt.TH1F("PPionProton2", "Charged Pion or Proton",60,0,10)
purityNoPhotons2 = rt.TH1F("PNoPhoton2", "No Real Photons",60,0,10)
purityOnePhoton2 = rt.TH1F("POnePhoton2", "1 Real Photon",60,0,10)
purityManyPhotons2 = rt.TH1F("PManyPhoton2", "3+ Real Photons",60,0,10)

purityTotal3 = rt.TH1F("PTotal3", "One Photon",60,0,10)
puritySignal3 =	 rt.TH1F("PSignal3", "Signal",60,0,10)
purityCC3 = rt.TH1F("PCC3", "Actually Charged Current",60,0,10)
purityFiducials3 = rt.TH1F("PFiducial3", "Out of Fiducial",60,0,10)
purityCosmic3 = rt.TH1F("PCosmic3", "Failed Cosmic",60,0,10)
purityPionProton3 = rt.TH1F("PPionProton3", "Charged Pion or Proton",60,0,10)
purityNoPhotons3 = rt.TH1F("PNoPhoton3", "No Real Photons",60,0,10)
purityOnePhoton3 = rt.TH1F("POnePhoton3", "2 Real Photons",60,0,10)
purityTwoPhotons3 = rt.TH1F("PTwoPhotons33", "3+ Real Photons",60,0,10)

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
cosmicOnePhoton = rt.TH1F("cBackground1", "One Photon Cosmic",60,0,10)
cosmicTwoPhotons = rt.TH1F("cBackground2", "Two Photon Cosmic",60,0,10)
cosmicThreePhotons = rt.TH1F("cBackground3", "3+ Photon Cosmic",60,0,10)
cosmicList = [cosmicOnePhoton, cosmicTwoPhotons, cosmicThreePhotons]

#Overall purity hitograms
overallPurityE1 = rt.TH1F("overallPurityE1", "One Photon",60,0,10)
overallPurityE2 = rt.TH1F("overallPurityE2", "One Photon",60,0,10)
overallPurityE3 = rt.TH1F("overallPurityE3", "One Photon",60,0,10)
overallPurityList = [overallPurityE1, overallPurityE2, overallPurityE3]

#Big Lists, for Big Plots
pList1 = [puritySignal1, purityCC1, purityFiducials1, purityCosmic1, purityPionProton1, purityNoPhotons1, purityTwoPhotons1, purityManyPhotons1]
pList2 = [puritySignal2, purityCC2, purityFiducials2, purityCosmic2, purityPionProton2, purityNoPhotons2, purityOnePhoton2, purityManyPhotons2] 
pList3 = [puritySignal3, purityCC3, purityFiducials3, purityCosmic3, purityPionProton3,	purityNoPhotons3, purityOnePhoton3, purityTwoPhotons3]

#Built-in functions here
def addHist(eventTree, photonList, histList, variable, weight):
  if len(photonList) == 1:
    histList[0].Fill(variable, weight)
  elif len(photonList) == 2:
    histList[1].Fill(variable, weight)
  else:
    histList[2].Fill(variable, weight)

def purityStack(title, purityList, cosmicHist, POTSum, cosmicSum):
  #Create Component Variables
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt.kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kGreen+4, rt.kOrange+1]
  POTTarget = 6.67e+20
  #Organize other histograms
  for x in range(len(purityList)):
    hist = purityList[x]
    bins = hist.GetNbinsX()
    hist.Scale(POTTarget/POTSum)
    hist.SetLineColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
    stack.Add(hist)

  #Format and add the cosmic histogram
  hist = cosmicHist
  bins = hist.GetNbinsX()
  hist.Scale(POTTarget/cosmicSum)
  hist.SetLineColor(rt.kBlack)
  histInt = hist.Integral(1, int(bins))
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
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

#Variables for program function
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10}


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

  #Make sure the event is within the fiducial volume
  if recoFiducials(eventTree, fiducialData) == False:
    continue

  #Cut events with suitably energetic protons or charged pions 
  if recoPionProton(eventTree) == False:
    continue
  
  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonList(eventTree)
  if len(recoList) == 0:
    continue
  
  #Cut events where the reco thinks the photon came from a charged parent
  if CCSeeker(eventTree, recoList) == False:
    continue
  
  #PURITY - GRAPHING BASED ON TRUTH
  
  #Calculating graphing values
  leadingPhoton = recoList[0]
  for x in range(len(recoList)):
    if eventTree.showerRecoE[x] > eventTree.showerRecoE[leadingPhoton]:
      leadingPhoton = x
  ElScore = abs(eventTree.showerFromNeutralScore[leadingPhoton])

  #if eventTree.nTracks > 0:
  #  leadingMuonTrack = eventTree.trackMuScore[0]
  #  for x in range(eventTree.nTracks):
  #    if eventTree.trackMuScore[x] > leadingMuonTrack:
  #      leadingMuonTrack = eventTree.trackMuScore[x]
  #else:
  #  leadingMuonTrack = 13
  #ElScore = abs(leadingMuonTrack)
  #addHist(cosmicTree, recoList, cosmicList, ElScore, 1)

  #Fill totals
  addHist(eventTree, recoList, totalPHists, ElScore, eventTree.xsecWeight)
 
  #Neutral current!
  if trueCutNC(eventTree) == False:
    addHist(eventTree, recoList, CCPHists, ElScore, eventTree.xsecWeight)
    continue
 
  #I suppose we can pretend that this is doing something
  if trueCutCosmic(eventTree) == False:
    addHist(eventTree, recoList, cosmicPHists, ElScore, eventTree.xsecWeight)
    continue

  if trueCutFiducials(eventTree, fiducialData) == False:
    addHist(eventTree, recoList, fiducialPHists, ElScore, eventTree.xsecWeight)
    continue

  #pions and protons!
  if trueCutPionProton(eventTree) == False:
    addHist(eventTree, recoList, pionProtonPHists, ElScore, eventTree.xsecWeight)
    continue
  #Now we make a list of the actual photons!
  truePhotonIDs = truePhotonList(eventTree, fiducialData)

  #Are there actually any photons?
  if len(truePhotonIDs) == 0:
    addHist(eventTree, recoList, noPhotonPHists, ElScore, eventTree.xsecWeight)
    continue
 
  #Is there one Photon?
  if len(truePhotonIDs) == 1:
    addHist(eventTree, recoList, onePhotonPHists, ElScore, eventTree.xsecWeight)

 #Are there two?
  elif len(truePhotonIDs) == 2:
    addHist(eventTree, recoList, twoPhotonPHists, ElScore, eventTree.xsecWeight)
 
  #In that case, there should be at least three
  else:
    addHist(eventTree, recoList, manyPhotonPHists, ElScore, eventTree.xsecWeight)

#BEGINNING EVENT LOOP FOR COSMICS
for i in range(cosmicTree.GetEntries()):
  cosmicTree.GetEntry(i)

  #COSMICS - SELECTING EVENTS BASED ON RECO
  #See if the event has a vertex
  if recoNoVertex(cosmicTree) == False:
    continue

  #See if the event is neutral current                                                                                                  
  if recoNeutralCurrent(cosmicTree) == False:
    continue

  #Make sure the event is within the fiducial volume
  if recoFiducials(cosmicTree, fiducialData) == False:
    continue

  #Cut events with suitably energetic protons or charged pions
  if recoPionProton(cosmicTree) == False:
    continue

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonList(cosmicTree)
  if len(recoList) == 0:
    continue

  #Try and put the process checker to work
  if CCSeeker(cosmicTree, recoList) == False:
    continue

  #graphing based on photon count
  #Calculating graphing values
  leadingPhoton = recoList[0]
  for x in range(len(recoList)):
    if cosmicTree.showerRecoE[x] > cosmicTree.showerRecoE[leadingPhoton]:
      leadingPhoton = x
  ElScore = abs(cosmicTree.showerFromNeutralScore[leadingPhoton])
  #if eventTree.nTracks > 0:
  #  leadingMuonTrack = eventTree.trackMuScore[0]
  #  for x in range(eventTree.nTracks):
  #    if eventTree.trackMuScore[x] > leadingMuonTrack:
  #      leadingMuonTrack = eventTree.trackMuScore[x]
  #else:
  #  leadingMuonTrack = 15
  #ElScore = abs(leadingMuonTrack)
  #addHist(cosmicTree, recoList, cosmicList, ElScore, 1)


#LOOPS OVER - HISTOGRAM ORGANIZING TIME
purityCanvas1, purityStack1, purityLegend1, purityInt1 = purityStack("Single-Photon Purity", pList1, cosmicList[0], ntuplePOTsum, cosmicPOTsum)
purityCanvas2, purityStack2, purityLegend2, purityInt2 = purityStack("Two-Photon Purity", pList2, cosmicList[1], ntuplePOTsum, cosmicPOTsum)
purityCanvas3, purityStack3, purityLegend3, purityInt3 = purityStack("3+ Photon Purity", pList3, cosmicList[2], ntuplePOTsum, cosmicPOTsum)
cosmicCanvas, cosmicStack, cosmicLegend, cosmicInt = histStack("Cosmic Background", cosmicList, cosmicPOTsum)

stackList = [purityStack1, purityStack2, purityStack3,cosmicStack]
for stack in stackList:
  stack.GetXaxis().SetTitle("Shower-From-Charged Score")

#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
purityCanvas1.Write()
purityCanvas2.Write()
purityCanvas3.Write()
cosmicCanvas.Write()
#for hist in overallPurityList:
#  hist.Write()