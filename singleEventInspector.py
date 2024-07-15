#Import the usual suspects
import sys, argparse
import numpy as np
import ROOT as rt
from helpers.larflowreco_ana_funcs import getCosThetaGravVector


parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-n", "--eventNumber", type=int, required=True, help="the number of the event you want to search")
parser.add_argument("-listTrue", "--listTrueParticles", action="store_true", help="printthe lists of tracks and particles with corresponding information")
parser.add_argument("-listReco", "--listRecoParticles", action="store_true", help="printthe lists of tracks and particles with corresponding information")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
args = parser.parse_args()

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

#calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

#Set detector volumes
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

#Set variables for program review

#Open specified entry
eventTree.GetEntry(args.eventNumber)

#Initial announcment
print("NOW READING: EVENT NUMBER", args.eventNumber)
print("RUN NUMBER:", eventTree.run, "SUBRUN NUMBER:", eventTree.subrun, "EVENT:", eventTree.event)
#RECO EVALUATION
#Check for photon, charge current
photonFound = False
chargeCurrent = False
photonIDs = []
for x in range(eventTree.nShowers):
  if eventTree.showerClassified[x] == 1:
    if eventTree.showerPID[x] == 22:
      photonFound = True
      photonIDs.append(x)
    elif eventTree.showerPID[x] == 13:
      chargeCurrent = True

print("This collision had a neutrino energy of", eventTree.trueNuE, "GeV.")
      
if photonFound == True:
  print("This event appears to have", len(photonIDs), " photons present in the reco.")
else:
  print("This event does not appear to contain any photons based on the reco.")

#Check tracks to check for Neutral Current
for x in range(eventTree.nTracks):
  if eventTree.trackClassified[x] == 1:
    if eventTree.trackPID[x] == 13:
      chargeCurrent = True
      
if chargeCurrent == True:
  if eventTree.trueNuCCNC == 1:
    print("This event has a charged current, but the reco thinks it's neutral.")
  else:
    print("This event has a charge current, and the reco is able to detect it.")
  
else:
  if eventTree.trueNuCCNC == 0:
    print("This event has a neutral current, but the reco thinks it's charged.")
  else:
    print("This event has a charged current,  and the reco is able to detect it.")

  
#EVALUATION USING TRUTH
truePhotonIDs = []
for x in range(eventTree.nTrueSimParts):
  if eventTree.trueSimPartPDG[x] == 22:
    truePhotonIDs.append(x)

trueSecondaryIDs = []
for x in range(len(truePhotonIDs)):
  if eventTree.trueSimPartPDG[x] == 22:
      if abs(eventTree.trueSimPartX[x] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[x] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[x] - eventTree.trueVtxZ) <= 0.15:
        trueSecondaryIDs.append(x)

if len(trueSecondaryIDs) > 0:
  print("Based on truth, this event actually contains", len(truePhotonIDs), "photons, and based on vertex information", len(trueSecondaryIDs), "of them are secondary particles.")
elif len(truePhotonIDs) > 0:
  print("Based on truth, this event actually contains", len(truePhotonIDs), "photons. Based on vertex information, none of them are secondary particles.")
else:  
  print("Based on truth, this event contains no photons")
  
#PRINTING LISTS
if args.listRecoParticles:
  print("SHOWERS AND ENERGY BASED ON RECO")
  for x in range(eventTree.nShowers):
    if eventTree.showerClassified[x] == 1:
      print("Shower ID:", eventTree.showerPID[x], "Shower Energy (MeV):", eventTree.showerRecoE[x])

  print("TRACKS AND ENERGY BASED ON RECO")
  for x in range(eventTree.nTracks):
    if eventTree.trackClassified[x] == 1:
      print("Track ID:", eventTree.trackPID[x], "Track Energy (MeV):", eventTree.trackRecoE[x])

        
if args.listTrueParticles:
  print("PRIMARY PARTICLE PDG VALUES BASED ON TRUTH")
  for x in range(eventTree.nTruePrimParts):
    print(eventTree.truePrimPartPDG[x])

  print("TRUE SIM PARTICLE TID, PDG,  MOTHER TID VALUES, AND ENERGY:")
  for x in range(eventTree.nTrueSimParts):
    print("TID:", eventTree.trueSimPartTID[x], "PDG:", eventTree.trueSimPartPDG[x], "Mother TID:", eventTree.trueSimPartMID[x], "TOTAL ENERGY(MeV):", eventTree.trueSimPartE[x])
