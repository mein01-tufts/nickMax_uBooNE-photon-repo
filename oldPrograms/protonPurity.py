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
#Histograms for proton tracks reco'd above 100 MeV
noTrackHighE = rt.TH1F("noTrackHigh", "No particle present",60,0,200)
muonHighE = rt.TH1F("muonHigh", "Actually a Muon",60,0,200)
notProtonHighE = rt.TH1F("notProtonHigh", "Some other Particle",60,0,200)
overEstimateHighE = rt.TH1F("overEstimate", "Actually < 100 MeV",60,0,200)
signalProtonHighE = rt.TH1F("signalProtonHigh", "Signal",60,0,200)

overTrackHists = [noTrackHighE, muonHighE, notProtonHighE, overEstimateHighE, signalProtonHighE]

#Histograms for proton tracks below 100 MeV
noTrackLowE = rt.TH1F("noTrackLow", "No particle present",60,0,15)
muonLowE = rt.TH1F("muonHighE", "Actually a Muon",60,0,15)
notProtonLowE = rt.TH1F("notProtonLow", "Some other Particle",60,0,15)
underEstimateLowE = rt.TH1F("underEstimate", "Actually > 100 MeV",60,0,15)
signalProtonLowE = rt.TH1F("signalProtonLow", "Correctly reconstructed",60,0,15)

underTrackHists = [noTrackLowE, muonLowE, notProtonLowE, underEstimateLowE, signalProtonLowE]

#Histograms for proton showers above 100 MeV
noShowerHighE = rt.TH1F("noShowerHigh", "No particle present",60,0,15)
muonShowerHighE = rt.TH1F("muonShowerHighE", "Actually a Muon",60,0,15)
notProtonShowerHighE = rt.TH1F("notProtonShowerHighE", "Some other Particle",60,0,15)
overEstimateShowerHighE = rt.TH1F("overEstimateShowerHighE", "Actually < 100 MeV",60,0,15)
signalProtonShowerHighE = rt.TH1F("signalProtonShowerHighE", "Signal",60,0,15)

overShowerHists = [noShowerHighE, muonShowerHighE, notProtonShowerHighE, overEstimateShowerHighE, signalProtonShowerHighE]

#Histograms for proton showers below 100 MeV
noShowerLowE = rt.TH1F("noShowerLowE", "No particle present",60,0,20)
muonShowerLowE = rt.TH1F("muonShowerLowE", "Actually a Muon",60,0,20)
notProtonShowerLowE = rt.TH1F("notProtonShowerLowE", "Some other Particle",60,0,20)
underEstimateShowerLowE = rt.TH1F("overEstimateShowerLowE", "Actually > 100 MeV",60,0,20)
signalProtonShowerLowE = rt.TH1F("signalProtonShowerLowE", "Signal",60,0,20)

underShowerHists = [noShowerLowE, muonShowerLowE, notProtonShowerLowE, underEstimateShowerLowE, signalProtonShowerLowE]

#Functions for making histograms


#Variables for program review
reconstructedPhotons = 0
correctPhotons = 0

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

  #Iterate over reconstructed tracks to see what we're identifying as a proton
  for x in range(eventTree.nTracks):
    if eventTree.trackPID[x] == 2212:
      #If possible, match to a proton in True Sim Parts
      match = 0
      for y in range(eventTree.nTrueSimParts):
        if eventTree.trackTrueTID[x] == eventTree.trueSimPartTID[y]:
          match = y
          break
      #Calculate reco track length, since that can be used to make cuts
      distxyz = np.sqrt((eventTree.trackStartPosX[x] - eventTree.trackEndPosX[x])**2 + (eventTree.trackStartPosY[x] - eventTree.trackEndPosY[x])**2 + (eventTree.trackStartPosZ[x] - eventTree.trackEndPosZ[x])**2)
      #See if the proton reconstructs over the threshold
      if eventTree.trackRecoE[x] > 100:
        #Make sure there's actually a particle present
        if eventTree.trackTruePID[x] == 0 or match == 0:
          noTrackHighE.Fill(distxyz, eventTree.xsecWeight)
        #Check to see if that particle is a muon (I suspect this might be a common issue)
        elif eventTree.trackTruePID[x] == 13:
          muonHighE.Fill(distxyz, eventTree.xsecWeight)
        #Check to see if it's anything else
        elif eventTree.trackTruePID[x] != 2212:
          notProtonHighE.Fill(distxyz, eventTree.xsecWeight)
        #Make sure the proton is actually over the energy threshold
        elif eventTree.trueSimPartE[match] - np.sqrt(eventTree.trueSimPartE[match]**2 - ((eventTree.trueSimPartPx[match]**2)+eventTree.trueSimPartPy[match]**2+eventTree.trueSimPartPz[match]**2)) < 100:
          overEstimateHighE.Fill(distxyz, eventTree.xsecWeight)
        #If it passes all of this, it's signal
        else:
          signalProtonHighE.Fill(distxyz, eventTree.xsecWeight)
      
      #Now we look at proton tracks that reconstructed below threshold
      else:
        #Make sure there's actually a particle present
        if eventTree.trackTruePID[x] == 0 or match == 0:
          noTrackLowE.Fill(distxyz, eventTree.xsecWeight)
        #Check to see if that particle is a muon (I suspect this might be a common issue)
        elif eventTree.trackTruePID[x] == 13:
          muonLowE.Fill(distxyz, eventTree.xsecWeight)
        #Check to see if it's anything else
        elif eventTree.trackTruePID[x] != 2212:
          notProtonLowE.Fill(distxyz, eventTree.xsecWeight)
        #Make sure the proton is actually over the energy threshold
        elif eventTree.trueSimPartE[match] - np.sqrt(eventTree.trueSimPartE[match]**2 - ((eventTree.trueSimPartPx[match]**2)+eventTree.trueSimPartPy[match]**2+eventTree.trueSimPartPz[match]**2)) > 100:
          underEstimateLowE.Fill(distxyz, eventTree.xsecWeight)
        #If it passes all of this, it's signal
        else:
          signalProtonLowE.Fill(distxyz, eventTree.xsecWeight)

  #Now we do the same thing - but for SHOWERS!
  for x in range(eventTree.nShowers):
    if eventTree.showerPID[x] == 2212:
      #If possible, match to a proton in True Sim Parts
      match = 0
      for y in range(eventTree.nTrueSimParts):
        if eventTree.showerTrueTID[x] == eventTree.trueSimPartTID[y]:
          match = y
          break
      #Calculate reco shower distance, since that can be used to make cuts
      distxyz = np.sqrt((eventTree.showerStartPosX[x] - eventTree.vtxX)**2 + (eventTree.showerStartPosY[x] - eventTree.vtxY)**2 + (eventTree.showerStartPosZ[x] - eventTree.vtxZ)**2)
      #See if the proton reconstructs over the threshold
      if eventTree.showerRecoE[x] > 100:
        #Make sure there's actually a particle present
        if eventTree.showerTruePID[x] == 0:
          noShowerHighE.Fill(distxyz, eventTree.xsecWeight)
        #Check to see if that particle is a muon (I suspect this might be a common issue)
        elif eventTree.showerTruePID[x] == 13:
          muonShowerHighE.Fill(distxyz, eventTree.xsecWeight)
        #Check to see if it's anything else
        elif eventTree.showerTruePID[x] != 2212:
          notProtonShowerHighE.Fill(distxyz, eventTree.xsecWeight)
        #Make sure the proton is actually over the energy threshold
        elif eventTree.trueSimPartE[match] - np.sqrt(eventTree.trueSimPartE[match]**2 - ((eventTree.trueSimPartPx[match]**2)+eventTree.trueSimPartPy[match]**2+eventTree.trueSimPartPz[match]**2)) < 100:
          overEstimateShowerHighE.Fill(distxyz, eventTree.xsecWeight)
        #If it passes all of this, it's signal
        else:
          signalProtonShowerHighE.Fill(distxyz, eventTree.xsecWeight)
      
      #Now we look at proton tracks that reconstructed below threshold
      else:
        #Make sure there's actually a particle present
        if eventTree.showerTruePID[x] == 0:
          noShowerLowE.Fill(distxyz, eventTree.xsecWeight)
        #Check to see if that particle is a muon (I suspect this might be a common issue)
        elif eventTree.showerTruePID[x] == 13:
          muonShowerLowE.Fill(distxyz, eventTree.xsecWeight)
        #Check to see if it's anything else
        elif eventTree.showerTruePID[x] != 2212:
          notProtonShowerLowE.Fill(distxyz, eventTree.xsecWeight)
        #Make sure the proton is actually over the energy threshold
        elif eventTree.trueSimPartE[match] - np.sqrt(eventTree.trueSimPartE[match]**2 - ((eventTree.trueSimPartPx[match]**2)+eventTree.trueSimPartPy[match]**2+eventTree.trueSimPartPz[match]**2)) > 100:
          underEstimateShowerLowE.Fill(distxyz, eventTree.xsecWeight)
        #If it passes all of this, it's signal
        else:
          signalProtonShowerLowE.Fill(distxyz, eventTree.xsecWeight)
            

#END OF LOOP - STACKING HISTOGRAMS
highTrackCanvas, highTrackStack, highTrackLegend, highTrackInt = histStack("Tracks reconstructed as protons Over 100 MeV", overTrackHists, ntuplePOTsum)
lowTrackCanvas, lowTrackStack, lowTrackLegend, lowTrackInt = histStack("Tracks reconstructed as protons Under 100 MeV", underTrackHists, ntuplePOTsum)
highShowerCanvas, highShowerStack, highShowerLegend, highShowerInt = histStack("Showers reconstructed as protons Over 100 MeV", overShowerHists, ntuplePOTsum)
lowShowerCanvas, lowShowerStack, lowShowerLegend, lowShowerInt = histStack("Showers reconstructed as protons Under 100 MeV", underShowerHists, ntuplePOTsum)

canvasList =[highTrackCanvas, lowTrackCanvas, highShowerCanvas, lowShowerCanvas]

trackStackList = [highTrackStack, lowTrackStack]
showerStackList = [highShowerStack, lowShowerStack]

for stack in trackStackList:
  stack.GetXaxis().SetTitle("Track Length (cm)")
  stack.GetYaxis().SetTitle("No. of Particles per 6.67e+20 POT")

for stack in showerStackList:
  stack.GetXaxis().SetTitle("Shower Distance (cm)")
  stack.GetYaxis().SetTitle("No. of Particles per 6.67e+20 POT")


#Save to file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in canvasList:
  canvas.Write()