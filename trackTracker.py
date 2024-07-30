import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)

from cuts import trueCutNC, trueCutFiducials, trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
args = parser.parse_args()

#needed for proper scaling of error bars:
rt.TH1.SetDefaultSumw2(rt.kTRUE)

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
foundHistMuon = rt.TH1F("foundHistMuon", "Muon Correctly Reconstructed",60,-0.5,10)
unclassifiedHistMuon = rt.TH1F("unclassifiedHistMuon", "Muon Track Unclassified",60,-0.5,10)
weirdHistMuon = rt.TH1F("weirdHistMuon", "Muon Reconstructed as Something Else",60,-0.5,10)
lostHistMuon = rt.TH1F("lostHistMuon", "Muon Track not Reconstructed",60,-0.5,10)

muonHistList = [foundHistMuon, unclassifiedHistMuon, weirdHistMuon, lostHistMuon]

foundHistPion = rt.TH1F("foundHistPion", "Pion Correctly Reconstructed",60,0,3)
unclassifiedHistPion = rt.TH1F("unclassifiedHistPion", "Pion Track Unclassified",60,0,3)
weirdHistPion = rt.TH1F("weirdHistPion", "Pion Reconstructed as Something Else",60,0,3)
lostHistPion = rt.TH1F("lostHistPion", "Pion Track not Reconstructed",60,0,3)

pionHistList = [foundHistPion, unclassifiedHistPion, weirdHistPion, lostHistPion]

foundHistProton = rt.TH1F("foundHistProton", "Proton Correctly Reconstructed",60,0,3)
unclassifiedHistProton = rt.TH1F("unclassifiedHistProton", "Proton Track Unclassified",60,0,3)
weirdHistProton = rt.TH1F("weirdHistProton", "Proton Reconstructed as Something Else",60,0,3)
lostHistProton = rt.TH1F("lostHistProton", "Proton Track not Reconstructed",60,0,3)

protonHistList = [foundHistProton, unclassifiedHistProton, weirdHistProton, lostHistProton]

effMuons = rt.TH1F("effMuons", "Muon Reconstruction Efficiency",60,-0.5,10)
effPions = rt.TH1F("effPions", "Pion Reconstruction Efficiency",60,0,3)
effProtons = rt.TH1F("effProtons", "Proton Reconstruction Efficiency",60,0,3)

foundList = [foundHistMuon, foundHistPion, foundHistProton]
unclassifiedList = [unclassifiedHistMuon, unclassifiedHistPion, unclassifiedHistProton]
weirdList = [weirdHistMuon, weirdHistPion, weirdHistProton]
lostList = [lostHistMuon, lostHistPion, lostHistProton]
effList = [effMuons, effPions, effProtons]

#Functions for making histograms


#Variables for program review
tracked = 0
unclassified = 0
weird = 0
lost = 0

pi1 = 0
pi2 = 0

pr1 = 0
pr2 = 0

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

  #Make sure the event is within the fiducial volume
  if recoFiducials(eventTree, fiducialData) == False:
    continue
  
  #Check for photons to ensure we're getting an otherwise accurate representation of our dataset
  recoPhotons = recoPhotonList(eventTree)
  if len(recoPhotons) == 0:
    continue

  #For each event, use truth to find and list relevant particles
  protonList = []
  pionList = []
  muonList = []

  #Find and list primary track particles in the event
  for x in range(eventTree.nTrueSimParts):
    if eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
      particleEnergy = eventTree.trueSimPartE[x] - np.sqrt(abs(eventTree.trueSimPartE[x]**2 - (eventTree.trueSimPartPx[x]**2+eventTree.trueSimPartPy[x]**2+eventTree.trueSimPartPz[x]**2)))
      if eventTree.trueSimPartPDG[x] == 13: #and particleEnergy >= 20:
        muonList.append(x)
      elif abs(eventTree.trueSimPartPDG[x]) == 211: #and particleEnergy >= 30:
        pionList.append(x)
      elif eventTree.trueSimPartPDG[x] == 2212: #and particleEnergy >= 60:
        protonList.append(x)

  #Make sure there are actually real tracks in this event
  if len(muonList) == 0 and len(pionList) == 0 and len(protonList) == 0:
    continue

  #Now we try to match reconstructed tracks to real ones
  for y in muonList:
    particleEnergy = eventTree.trueSimPartE[y] - np.sqrt(abs(eventTree.trueSimPartE[y]**2 - (eventTree.trueSimPartPx[y]**2+eventTree.trueSimPartPy[y]**2+eventTree.trueSimPartPz[y]**2)))
    trackFound = False
    for x in range(eventTree.nTracks):
      if eventTree.trackTrueTID[x] == eventTree.trueSimPartTID[y]:
        trackFound = True
        #distxyz = np.sqrt(abs(eventTree.trackStartPosX[x] - eventTree.trackEndPosX[x])**2 + (eventTree.trackStartPosY[x] - eventTree.trackEndPosY[x])**2 + (eventTree.trackStartPosZ[x] - eventTree.trackEndPosZ[x])**2)
        #muonScore = abs(eventTree.trackMuScore[x])
        if eventTree.trackPID[x] == 13:
          foundHistMuon.Fill(abs(eventTree.trackMuScore[x]), eventTree.xsecWeight)
          tracked += 1
          break
        elif eventTree.trackClassified[x] == 0:
          unclassifiedHistMuon.Fill(abs(eventTree.trackMuScore[x]), eventTree.xsecWeight)
          unclassified += 1
          break
        else:
          weirdHistMuon.Fill(abs(eventTree.trackMuScore[x]), eventTree.xsecWeight)
          weird += 1
          break
    if trackFound == False:
      distxyz = -0.1
      lostHistMuon.Fill(distxyz, eventTree.xsecWeight)
      lost += 1

  for y in pionList:
    trackFound = False
    particleEnergy = eventTree.trueSimPartE[y] - np.sqrt(abs(eventTree.trueSimPartE[y]**2 - (eventTree.trueSimPartPx[y]**2+eventTree.trueSimPartPy[y]**2+eventTree.trueSimPartPz[y]**2)))
    for x in range(eventTree.nTracks):
      if eventTree.trackTrueTID[x] == eventTree.trueSimPartTID[y]:
        trackFound = True
        #pionScore = abs(eventTree.trackPiScore[x])
        #distxyz = np.sqrt(abs(eventTree.trackStartPosX[x] - eventTree.trackEndPosX[x])**2 + (eventTree.trackStartPosY[x] - eventTree.trackEndPosY[x])**2 + (eventTree.trackStartPosZ[x] - eventTree.trackEndPosZ[x])**2)
        if abs(eventTree.trackPID[x]) == 211:
          foundHistPion.Fill(particleEnergy/1000, eventTree.xsecWeight)
          pi1 += 1
          break
        elif eventTree.trackClassified[x] == 0:
          unclassifiedHistPion.Fill(particleEnergy/1000, eventTree.xsecWeight)
          pi2 += 1
          break
        else:
          weirdHistPion.Fill(particleEnergy/1000, eventTree.xsecWeight)
          pi2 += 1
          break
    if trackFound == False:
      distxyz = -0.1
      lostHistPion.Fill(particleEnergy/1000, eventTree.xsecWeight)
      pi2 += 1


  for y in protonList:
    trackFound = False
    particleEnergy = eventTree.trueSimPartE[y] - np.sqrt(abs(eventTree.trueSimPartE[y]**2 - (eventTree.trueSimPartPx[y]**2+eventTree.trueSimPartPy[y]**2+eventTree.trueSimPartPz[y]**2)))
    for x in range(eventTree.nTracks):
      if eventTree.trackTrueTID[x] == eventTree.trueSimPartTID[y]:
        trackFound = True
        #distxyz = np.sqrt(abs(eventTree.trackStartPosX[x] - eventTree.trackEndPosX[x])**2 + (eventTree.trackStartPosY[x] - eventTree.trackEndPosY[x])**2 + (eventTree.trackStartPosZ[x] - eventTree.trackEndPosZ[x])**2)
        #protonScore = abs(eventTree.trackPrScore[x])
        if eventTree.trackPID[x] == 2212:
          foundHistProton.Fill(particleEnergy/1000, eventTree.xsecWeight)
          pr1 += 1
          break
        elif eventTree.trackClassified[x] == 0:
          unclassifiedHistProton.Fill(particleEnergy/1000, eventTree.xsecWeight)
          pr2 += 1
          break
        else:
          weirdHistProton.Fill(particleEnergy/1000, eventTree.xsecWeight)
          pr2 += 1
          break
    if trackFound == False:
      #distxyz = -0.1
      lostHistProton.Fill(particleEnergy/1000, eventTree.xsecWeight)
      pr2 += 1


#END OF LOOP - STACKING HISTOGRAMS

for y in range(3):
  for x in range(1, 61):
    signal = foundList[y].GetBinContent(x)
    total = foundList[y].GetBinContent(x) + unclassifiedList[y].GetBinContent(x) + weirdList[y].GetBinContent(x) + lostList[y].GetBinContent(x)
    if total > 0:
      efficiency = signal/total
    else:
      efficiency = 0
    effList[y].SetBinContent(x, efficiency)

muonCanvas, muonStack, muonLegend, muonInt = histStack("Reconstructions of Muon Tracks of Various Lengths", muonHistList, ntuplePOTsum)
pionCanvas, pionStack, pionLegend, pionInt = histStack("Reconstructions of Pion Tracks of Various Lengths", pionHistList, ntuplePOTsum)
protonCanvas, protonStack, protonLegend, protonInt = histStack("Reconstructions of Proton Tracks of Various Lengths", protonHistList, ntuplePOTsum)

canvasList =[muonCanvas, pionCanvas, protonCanvas]

stackList = [muonStack, pionStack, protonStack]
for stack in stackList:
  stack.GetXaxis().SetTitle("Particle True Kinetic Energy (GeV)")
  stack.GetYaxis().SetTitle("No. of Particles per 6.67e+20 POT")

for hist in effList:
  hist.GetXaxis().SetTitle("Particle True Kinetic Energy (GeV)")
  hist.GetYaxis().SetTitle("Efficiency")

muonStack.GetXaxis().SetTitle("Muon Score")

#Save to file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in canvasList:
  canvas.Write()
for hist in effList:
  hist.Write()

#Calculate overall Efficiencies
muonEff = tracked/(tracked + unclassified + weird + lost)
pionEff = pi1/(pi1 + pi2)
protonEff = pr1/(pr1 + pr2)

#Report results
print("Muons found:", tracked)
print("Muons unclassified:", unclassified)
print("Muons misclassified:", weird)
print("Muons lost:", lost)
print("Muon Efficiency:", muonEff)
print("Pion Efficiency:", pionEff)
print("Proton Efficiency:", protonEff)