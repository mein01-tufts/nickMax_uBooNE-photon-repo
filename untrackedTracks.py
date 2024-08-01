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
muonHist = rt.TH1F("muonHist", "Under-Threshold Muons",60,0,80)
protonHist = rt.TH1F("protonHist", "Under-Threshold Protons",60,0,80)
pionHist = rt.TH1F("pionHist", "Under-Threshold Pions",60,0,80)
muonHistBackground = rt.TH1F("muonHistBackground", "Over-Threshold Muons",60,0,80)
protonHistBackground = rt.TH1F("protonHistBackground", "Over-Threshold Protons",60,0,80)
pionHistBackground = rt.TH1F("pionHistBackground", "Over-Threshold Pions",60,0,80)
nothingHist = rt.TH1F("nothingHist", "Not a Real Track",60,0,80)
otherHist = rt.TH1F("otherHist", "Something else",60,0,80)

muonHistList = [muonHist, muonHistBackground]
protonHistList = [protonHist, protonHistBackground]
pionHistList = [pionHist, pionHistBackground]
otherHistList = [nothingHist, otherHist]

#Functions for making histograms


#Variables for program review

#Variables for program function
recoPhotonIDs = []
truePhotonIDs = []
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":30}

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

  #Now we try to match reconstructed tracks to real ones
  for x in range(eventTree.nTracks):
    if eventTree.trackClassified[x] == 0:
      particleEnergy = 0
      for y in range(eventTree.nTrueSimParts):
        if eventTree.trackTrueTID[x] == eventTree.trueSimPartTID[y]:
          particleEnergy = eventTree.trueSimPartE[y] - np.sqrt(abs(eventTree.trueSimPartE[y]**2 - (eventTree.trueSimPartPx[y]**2+eventTree.trueSimPartPy[y]**2+eventTree.trueSimPartPz[y]**2)))
          break
      distxyz = np.sqrt((eventTree.trackStartPosX[x] - eventTree.trackEndPosX[x])**2 + (eventTree.trackStartPosY[x] - eventTree.trackEndPosY[x])**2 + (eventTree.trackStartPosZ[x] - eventTree.trackEndPosZ[x])**2)

      #See if the track is a muon
      if eventTree.trackTruePID[x] == 13:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(len(eventTree.truePrimPartPDG)):
            if eventTree.truePrimPartPDG[y] == 13:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        if particleEnergy > 100:
          muonHistBackground.Fill(distxyz, eventTree.xsecWeight)
        else:
          muonHist.Fill(distxyz, eventTree.xsecWeight)

      #See if the track is a pion
      elif abs(eventTree.trackTruePID[x]) == 211:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(len(eventTree.truePrimPartPDG)):
            if abs(eventTree.truePrimPartPDG[y]) == 211:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        if particleEnergy > 50:
          pionHistBackground.Fill(distxyz, eventTree.xsecWeight)
        else:
          pionHist.Fill(distxyz, eventTree.xsecWeight)

      #See if the track is a proton
      elif eventTree.trackTruePID[x] == 2212:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(eventTree.nTruePrimParts):
            if eventTree.truePrimPartPDG[y] == 2212:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        if particleEnergy > 100:
          protonHistBackground.Fill(distxyz, eventTree.xsecWeight)
        else:
          protonHist.Fill(distxyz, eventTree.xsecWeight)
      
      #See if the track is just a figment of the reco's imagination
      elif eventTree.trackTruePID[x] == 0: 
        nothingHist.Fill(distxyz, eventTree.xsecWeight)

      else:
        #Fill histogram for other particles
        otherHist.Fill(distxyz, eventTree.xsecWeight)

#END OF LOOP - STACKING HISTOGRAMS
muonList = [muonHist, muonHistBackground]

muonCanvas1, muonStack1, muonLegend1, muonInt1 = histStack("Unclassified muon tracks in events with photons", muonHistList, ntuplePOTsum)
protonCanvas1, protonStack1, protonLegend1, protonInt1 = histStack("Unclassified proton tracks in events with photons", protonHistList, ntuplePOTsum)
pionCanvas1, pionStack1, pionLegend1, pionInt1 = histStack("Unclassified pion tracks in events with photons", pionHistList, ntuplePOTsum)
otherCanvas1, otherStack1, otherLegend1, otherInt1 = histStack("Other unclassified tracks in events with photons", otherHistList, ntuplePOTsum)

canvasList =[muonCanvas1, protonCanvas1, pionCanvas1, otherCanvas1]

stackList = [muonStack1, protonStack1, pionStack1, otherStack1]

for stack in stackList:
  stack.GetXaxis().SetTitle("Track Length (cm)")
  stack.GetYaxis().SetTitle("No. of Particles per 6.67e+20 POT")

#Save to file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in canvasList:
  canvas.Write()