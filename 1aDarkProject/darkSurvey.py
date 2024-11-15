import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials,trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoPion, recoProton, recoCutShowerFromChargeScore, recoCutLongTracks, recoPhotonListFiducial, recoCutPrimary, recoCutShortTracks, recoPhotonListTracks, recoCutFarShowers, trueCutMuons, trueCutElectrons, recoCutMuons, recoCutElectrons, recoCutManyTracks, recoCutCompleteness, recoCutMuonCompleteness

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
args = parser.parse_args()

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")


#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20"
ntuplePOTsum = 0

for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT


showerHits =  rt.TH1F("showerHits", "ShowerNHits Values",60,0,4000)


#Variables for program review

#Variables for program function
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":30}
classificationThreshold = 0

trueElectronCountList = [0,0,0,0]
truePhotonCountList = [0,0,0,0]
trueMuonCountList = [0,0,0,0]
trueProtonCountList = [0,0,0,0]
truePionCountList = [0,0,0,0]
trueOtherNonsenseList = [0,0,0,0]

recoElectronCountList = [0,0,0,0]
recoPhotonCountList = [0,0,0,0]
recoMuonCountList = [0,0,0,0]
recoProtonCountList = [0,0,0,0]
recoPionCountList = [0,0,0,0]
recoOtherNonsenseList = [0,0,0,0]

vertexFound = 0
showerClassified = 0

#BEGINNING EVENT LOOP FOR DEFAULT PURITY
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)
  if eventTree.foundVertex == 1:
    vertexFound += 1

  #Reset looped variables
  trueElectronCount = 0
  truePhotonCount = 0
  trueMuonCount = 0
  trueProtonCount = 0
  truePionCount = 0
  trueOtherNonsense = 0

  recoElectronCount = 0
  recoPhotonCount = 0
  recoMuonCount = 0
  recoProtonCount = 0
  recoPionCount = 0
  recoOtherNonsense = 0

  for x in range(eventTree.nTrueSimParts):
    if abs(eventTree.trueSimPartPDG[x]) == 11:
      trueElectronCount += 1

    elif abs(eventTree.trueSimPartPDG[x]) == 22:
      truePhotonCount += 1

    elif abs(eventTree.trueSimPartPDG[x]) == 13:
      trueMuonCount += 1

    elif abs(eventTree.trueSimPartPDG[x]) == 211:
      truePionCount += 1

    elif abs(eventTree.trueSimPartPDG[x]) == 2212:
      trueProtonCount += 1

    else:
      trueOtherNonsense += 1
    
  for x in range(eventTree.nShowers):
    #print("Shower classified:", eventTree.showerClassified[x], "Shower size:", eventTree.showerSize[x], "reconstructed energy:", eventTree.showerRecoE[x])
    #print("ShowerNHits:", eventTree.showerNHits[x])
    if eventTree.showerClassified[x] != 0:
      showerClassified += 1
      print(eventTree.showerPID[x])
    showerHits.Fill(eventTree.showerNHits[x], 1)
    if abs(eventTree.showerPID[x]) == 11:
      recoElectronCount += 1

    elif abs(eventTree.showerPID[x]) == 22:
      recoPhotonCount += 1

    elif abs(eventTree.showerPID[x]) == 13:
      recoMuonCount += 1

    elif abs(eventTree.showerPID[x]) == 211:
      recoPionCount += 1

    elif abs(eventTree.showerPID[x]) == 2212:
      recoProtonCount += 1

    else:
      recoOtherNonsense += 1

  for x in range(eventTree.nTracks):
    if abs(eventTree.trackPID[x]) == 11:
      recoElectronCount += 1

    if abs(eventTree.trackPID[x]) == 22:
      recoPhotonCount += 1

    if abs(eventTree.trackPID[x]) == 13:
      recoMuonCount += 1

    if abs(eventTree.trackPID[x]) == 211:
      recoPionCount += 1

    if abs(eventTree.trackPID[x]) == 2212:
      recoProtonCount += 1

    else:
      recoOtherNonsense += 1

  #Now we adjust the tallies to get data

  #FOR TRUE VALUES
  if truePionCount == 0:
    truePionCountList[0] += 1
  elif truePionCount == 1:
    truePionCountList[1] += 1
  elif truePionCount == 2:
    truePionCountList[2] += 1
  else:
    truePionCountList[3] += 1

  if trueProtonCount == 0:
    trueProtonCountList[0] += 1
  elif trueProtonCount == 1:
    trueProtonCountList[1] += 1
  elif trueProtonCount == 2:
    trueProtonCountList[2] += 1
  else:
    trueProtonCountList[3] += 1

  if trueMuonCount == 0:
    trueMuonCountList[0] += 1
  elif trueMuonCount == 1:
    trueMuonCountList[1] += 1
  elif trueMuonCount == 2:
    trueMuonCountList[2] += 1
  else:
    trueMuonCountList[3] += 1

  if truePhotonCount == 0:
    truePhotonCountList[0] += 1
  elif truePhotonCount == 1:
    truePhotonCountList[1] += 1
  elif truePhotonCount == 2:
    truePhotonCountList[2] += 1
  else:
    truePhotonCountList[3] += 1

  if trueElectronCount == 0:
    trueElectronCountList[0] += 1
  elif trueElectronCount == 1:
    trueElectronCountList[1] += 1
  elif trueElectronCount == 2:
    trueElectronCountList[2] += 1
  else:
    trueElectronCountList[3] += 1

  #FOR RECO VALUES
  if recoPionCount == 0:
    recoPionCountList[0] += 1
  elif recoPionCount == 1:
    recoPionCountList[1] += 1
  elif recoPionCount == 2:
    recoPionCountList[2] += 1
  else:
    recoPionCountList[3] += 1

  if recoProtonCount == 0:
    recoProtonCountList[0] += 1
  elif recoProtonCount == 1:
    recoProtonCountList[1] += 1
  elif recoProtonCount == 2:
    recoProtonCountList[2] += 1
  else:
    recoProtonCountList[3] += 1

  if recoMuonCount == 0:
    recoMuonCountList[0] += 1
  elif recoMuonCount == 1:
    recoMuonCountList[1] += 1
  elif recoMuonCount == 2:
    recoMuonCountList[2] += 1
  else:
    recoMuonCountList[3] += 1

  if recoPhotonCount == 0:
    recoPhotonCountList[0] += 1
  elif recoPhotonCount == 1:
    recoPhotonCountList[1] += 1
  elif recoPhotonCount == 2:
    recoPhotonCountList[2] += 1
  else:
    recoPhotonCountList[3] += 1

  if recoElectronCount == 0:
    recoElectronCountList[0] += 1
  elif recoElectronCount == 1:
    recoElectronCountList[1] += 1
  elif recoElectronCount == 2:
    recoElectronCountList[2] += 1
  else:
    recoElectronCountList[3] += 1


#END OF LOOP
print("There were", eventTree.GetEntries(), "events. We found vertices for", vertexFound, "of them")
print("Of these,", trueElectronCountList[0], "contained no electrons,", trueElectronCountList[1], "contained one electron,", trueElectronCountList[2], "contained two electrons, and", trueElectronCountList[3], "contained three or more electrons.")
print(truePhotonCountList[0], "contained no photons,", truePhotonCountList[1], "contained one photon,", truePhotonCountList[2], "contained two photons, and", truePhotonCountList[3], "contained three or more photons.")
print(trueMuonCountList[0], "contained no muons,", trueMuonCountList[1], "contained one muon,", trueMuonCountList[2], "contained two muons, and", trueMuonCountList[3], "contained three or more muons.")
print(truePionCountList[0], "contained no pions,", truePionCountList[1], "contained one pion,", truePionCountList[2], "contained two pions, and", truePionCountList[3], "contained three or more pions.")
print(trueProtonCountList[0], "contained no protons,", trueProtonCountList[1], "contained one proton,", trueProtonCountList[2], "contained two protons, and", trueProtonCountList[3], "contained three or more protons.")

print("In terms of Reco:")
print(recoElectronCountList[0], "were reconstructed with no electrons,", recoElectronCountList[1], "were reconstructed with one electron,", recoElectronCountList[2], "were reconstructed with two electrons, and", recoElectronCountList[3], "were reconstructed with three or more electrons.")
print(recoPhotonCountList[0], "were reconstructed with no photons,", recoPhotonCountList[1], "were reconstructed with one photon,", recoPhotonCountList[2], "were reconstructed with two photons, and", recoPhotonCountList[3], "were reconstructed with three or more photons.")
print(recoMuonCountList[0], "were reconstructed with no muons,", recoMuonCountList[1], "were reconstructed with one muon,", recoMuonCountList[2], "were reconstructed with two muons, and", recoMuonCountList[3], "were reconstructed with three or more muons.")
print(recoPionCountList[0], "were reconstructed with no pions,", recoPionCountList[1], "were reconstructed with one pion,", recoPionCountList[2], "were reconstructed with two pions, and", recoPionCountList[3], "were reconstructed with three or more pions.")
print(recoProtonCountList[0], "were reconstructed with no protons,", recoProtonCountList[1], "were reconstructed with one proton,", recoProtonCountList[2], "were reconstructed with two protons, and", recoProtonCountList[3], "were reconstructed with three or more protons.")
print("There were", showerClassified, "Classified showers")

outFile = rt.TFile(args.outfile, "RECREATE")
showerHits.Write()