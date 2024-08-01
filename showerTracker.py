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
photonUnderHist = rt.TH1F("photonUnder", "Under-Threshold Photons",60,0,200)
photonOverHist = rt.TH1F("photonOver", "Over-Threshold Photons",60,0,200)
electronUnderHist = rt.TH1F("electronUnder", "Under-Threshold Electrons",60,0,200)
electronOverHist = rt.TH1F("electronOver", "Over-Threshold Electrons",60,0,200)
nothingHist = rt.TH1F("nothingHist", "Not a Real Shower",60,0,200)
positronHist = rt.TH1F("positronHist", "Positron",60,0,200)
protonHist = rt.TH1F("protonHist", "Protons",60,0,200)
muonHist = rt.TH1F("muonHist", "Muons",60,0,200)
pionHist = rt.TH1F("pionHist", "Pions",60,0,200)
reallyWeirdHist = rt.TH1F("reallyWeirdHist", "Something Weird",60,0,200)

photonHistList = [photonUnderHist, photonOverHist]
electronHistList = [electronUnderHist, electronOverHist]
otherHistList = [nothingHist, positronHist, muonHist, protonHist, pionHist, reallyWeirdHist]

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

  #Now we try to match reconstructed showers to real ones
  for x in range(eventTree.nShowers):
    if eventTree.showerClassified[x] == 1:
      particleEnergy = 0
      for y in range(eventTree.nTrueSimParts):
        if eventTree.showerTrueTID[x] == eventTree.trueSimPartTID[y]:
          particleEnergy = eventTree.trueSimPartE[y] - np.sqrt(abs(eventTree.trueSimPartE[y]**2 - (eventTree.trueSimPartPx[y]**2+eventTree.trueSimPartPy[y]**2+eventTree.trueSimPartPz[y]**2)))
          break
      distxyz = np.sqrt((eventTree.showerStartPosX[x] - eventTree.vtxX)**2 + (eventTree.showerStartPosY[x] - eventTree.vtxY)**2 + (eventTree.showerStartPosZ[x] - eventTree.vtxZ)**2)

      #See if the shower is a photon
      if eventTree.showerTruePID[x] == 22:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(len(eventTree.truePrimPartPDG)):
            if eventTree.truePrimPartPDG[y] == 22:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        if particleEnergy > 20:
          photonOverHist.Fill(distxyz, eventTree.xsecWeight)
        else:
          photonUnderHist.Fill(distxyz, eventTree.xsecWeight)

      #See if the shower is an electron
      elif eventTree.showerTruePID[x] == 11:
        #Get around a bug where the truth-matched particle only appears in truePrimParts
        if particleEnergy == 0:
          for y in range(len(eventTree.truePrimPartPDG)):
            if eventTree.truePrimPartPDG[y] == 11:
              particleEnergy = eventTree.truePrimPartE[y] - np.sqrt(abs(eventTree.truePrimPartE[y]**2 - (eventTree.truePrimPartPx[y]**2+eventTree.truePrimPartPy[y]**2+eventTree.truePrimPartPz[y]**2)))
              particleEnergy = particleEnergy/1000
        if particleEnergy > 10:
          electronOverHist.Fill(distxyz, eventTree.xsecWeight)
        else:
          electronUnderHist.Fill(distxyz, eventTree.xsecWeight)

      #See if the shower is just a figment of the reco's imagination
      elif eventTree.showerTruePID[x] == 0:
        notReal += 1
        if eventTree.showerPID[x] == 22:
          notRealPhoton += 1
        nothingHist.Fill(distxyz, eventTree.xsecWeight)

      else:
        #Fill histogram for other particles
        print("True:", eventTree.showerTruePID[x], "Reco:", eventTree.showerPID[x])
        if eventTree.showerTruePID[x] == -11:
          positronHist.Fill(distxyz, eventTree.xsecWeight)
        elif abs(eventTree.showerTruePID[x]) == 211:
          pionHist.Fill(distxyz, eventTree.xsecWeight)
        elif eventTree.showerTruePID[x] == 2212:
          protonHist.Fill(distxyz, eventTree.xsecWeight)
        elif eventTree.showerTruePID[x] == 13:
          muonHist.Fill(distxyz, eventTree.xsecWeight)
        else:
          reallyWeirdHist.Fill(distxyz, eventTree.xsecWeight)
        
            

#END OF LOOP - STACKING HISTOGRAMS
photonCanvas1, photonStack1, photonLegend1, photonInt1 = histStack("Classified photon shower in events with reconstructed photons", photonHistList, ntuplePOTsum)
electronCanvas1, electronStack1, electronLegend1, electronInt1 = histStack("Classified electron showers in events with reconstructed photons", electronHistList, ntuplePOTsum)
otherCanvas1, otherStack1, otherLegend1, otherInt1 = histStack("Other Classified showers in events with reconstructed photons", otherHistList, ntuplePOTsum)

canvasList =[photonCanvas1, electronCanvas1, otherCanvas1]

stackList = [photonStack1, electronStack1, otherStack1]

for stack in stackList:
  stack.GetXaxis().SetTitle("Shower start distance from vertex (cm)")
  stack.GetYaxis().SetTitle("No. of Particles per 6.67e+20 POT")

#Save to file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in canvasList:
  canvas.Write()