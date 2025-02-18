import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials, trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, trickyPionProtonCuts

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

lowECount = 0

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
args = parser.parse_args()

#needed for proper scaling of error bars:
#rt.TH1.SetDefaultSumw2(rt.kTRUE)

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
totalHist1 = rt.TH1F("total1", "One photon",60,0,2)
totalHist2 = rt.TH1F("total2", "Two photons",60,0,2)
totalHist3 = rt.TH1F("total3", "Three photons",60,0,2)
totalList = [totalHist1, totalHist2, totalHist3]

noVertexHist = rt.TH1F("noVertex1", "No vertex",60,0,2)
NCHist = rt.TH1F("NC1", "Charged current detected",60,0,2)
pionProtonHist = rt.TH1F("pionProton1", "Pion/Proton detected",60,0,2)
fiducialHist = rt.TH1F("fiducial1", "Out of fiducial",60,0,2)
noPhotonHist = rt.TH1F("noPhotons1", "No Photons detected",60,0,2)
tooManyHist = rt.TH1F("tooMany1", "Excess Photons detected",60,0,2)
tooFewHist = rt.TH1F("tooFew1", "Too few photons detected",60,0,2)
successHist = rt.TH1F("success1", "Signal",60,0,2)

histList = [tooManyHist, pionProtonHist, fiducialHist, NCHist, noVertexHist, noPhotonHist, tooFewHist, successHist]

noVertexHist2 = rt.TH1F("noVertex2", "No vertex",60,0,2)
NCHist2 = rt.TH1F("NC1", "Charged current detected",60,0,2)
pionProtonHist2 = rt.TH1F("pionProton2", "Pion/Proton detected",60,0,2)
fiducialHist2 = rt.TH1F("fiducial2", "Out of fiducial",60,0,2)
noPhotonHist2 = rt.TH1F("noPhotons2", "No Photons detected",60,0,2)
tooManyHist2 = rt.TH1F("tooMany2", "Excess Photons detected",60,0,2)
tooFewHist2 = rt.TH1F("tooFew2", "Too few photons detected",60,0,2)
successHist2 = rt.TH1F("success2", "Signal",60,0,2)

histList2 = [tooManyHist2, pionProtonHist2, fiducialHist2, NCHist2, noVertexHist2, noPhotonHist2, tooFewHist2, successHist2]

noVertexHist3 = rt.TH1F("noVertex3", "No vertex",60,0,2)
NCHist3 = rt.TH1F("NC1", "Charged current detected",60,0,2)
pionProtonHist3 = rt.TH1F("pionProton3", "Pion/Proton detected",60,0,2)
fiducialHist3 = rt.TH1F("fiducial3", "Out of fiducial",60,0,2)
noPhotonHist3 = rt.TH1F("noPhotons3", "No Photons detected",60,0,2)
tooManyHist3 = rt.TH1F("tooMany3", "Excess Photons detected",60,0,2)
tooFewHist3 = rt.TH1F("tooFew3", "Too few photons detected",60,0,2)
successHist3 = rt.TH1F("success3", "Signal",60,0,2)

histList3 = [tooManyHist3, pionProtonHist3, fiducialHist3, NCHist3, noVertexHist3, noPhotonHist3, tooFewHist3, successHist3]


#Variables for program review
trueSuccess = 0
NCSuccess = 0
fiducialSuccess = 0
cosmicSuccess = 0
pionProtonSuccess = 0

recoSuccess = 0


#Variables for program function
truePhotonIDs = []
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":10} 

#Event loop begins
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #Iterate with truth variables
  if trueCutNC(eventTree) == False:
    continue
  else:
    NCSuccess += 1

  if trueCutFiducials(eventTree, fiducialData) == False:
    continue
  else:
    fiducialSuccess += 1
    
  if trueCutCosmic(eventTree) == False:
    continue
  else:
    cosmicSuccess += 1
    
  if trueCutPionProton(eventTree) == False:
    continue
  else:
    pionProtonSuccess += 1
  
  truePhotonIDs = truePhotonList(eventTree, truePhotonIDs, fiducialData)
  if len(truePhotonIDs) > 0:
    trueSuccess += 1
  else:
    continue

  if trickyPionProtonCuts(eventTree) == False:
    continue
  else:
    lowECount += 1
    
  #SORT USING RECO VARIABLES AND FILL HISTOGRAMS
  #Scale and determine leading photon
  scaledEnergy = []
  for x in range(len(truePhotonIDs)):
    energyGeV = eventTree.trueSimPartE[truePhotonIDs[x]]/1000
    scaledEnergy.append(energyGeV)
    
  leadingPhoton = scaledEnergy[0]
  for x in range(len(scaledEnergy)):
    if scaledEnergy[x] > leadingPhoton:
      leadingPhoton = scaledEnergy[x]

  #Fill total graph
  if len(truePhotonIDs) == 1:
    totalHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

  elif len(truePhotonIDs) == 2:
    totalHist2.Fill(leadingPhoton, eventTree.xsecWeight)

  else:
    totalHist3.Fill(leadingPhoton, eventTree.xsecWeight)    
  
  #Checking the number of photons the reco finds
  recoList = recoPhotonList(eventTree)

  #Too many photons
  if len(recoList) > len(truePhotonIDs):
    if len(truePhotonIDs) == 1:
      tooManyHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonIDs) == 2:
      tooManyHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      tooManyHist3.Fill(leadingPhoton, eventTree.xsecWeight)
    continue

  #Fill for events with suitably energetic pions or protons                                                                                    
  if recoPionProton(eventTree) == False:
    if len(truePhotonIDs) == 1:
      pionProtonHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonIDs) == 2:
      pionProtonHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      pionProtonHist3.Fill(leadingPhoton, eventTree.xsecWeight)
    continue

  #Fill for events outside the fiducial                                                                                                        
  if recoFiducials(eventTree, fiducialData) == False:
    if len(truePhotonIDs) == 1:
      fiducialHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonIDs) == 2:
      fiducialHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      fiducialHist3.Fill(leadingPhoton, eventTree.xsecWeight)
    continue
  
  #Fill for events with charged current                                                                                                       
  if recoNeutralCurrent(eventTree) == False:
    if len(truePhotonIDs) == 1:
      NCHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonIDs) == 2:
      NCHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      NCHist3.Fill(leadingPhoton, eventTree.xsecWeight)

    continue
  
  #Fill graph for events with no vertex                                                                                                        
  if recoNoVertex(eventTree) == False:
    if len(truePhotonIDs) == 1:
      noVertexHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonIDs) == 2:
      noVertexHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      noVertexHist3.Fill(leadingPhoton, eventTree.xsecWeight)

    continue

  #Filling for no photons                                                                                                                       
  if len(recoList) == 0:
    if len(truePhotonIDs) == 1:
      noPhotonHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonIDs) == 2:
      noPhotonHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      noPhotonHist3.Fill(leadingPhoton, eventTree.xsecWeight)

    continue

  #Too few photons                                                                                                                              
  if len(recoList) < len(truePhotonIDs):
    if len(truePhotonIDs) == 1:
      tooFewHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonIDs) == 2:
      tooFewHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      tooFewHist3.Fill(leadingPhoton, eventTree.xsecWeight)
    continue

  #If an event has made it this far, it's signal
  recoSuccess += 1
  if len(truePhotonIDs) == 1:
    successHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

  elif len(truePhotonIDs) == 2:
    successHist2.Fill(leadingPhoton, eventTree.xsecWeight)

  else:
    successHist3.Fill(leadingPhoton, eventTree.xsecWeight)
  

histCanvas, stack, legend, histInt = histStack("Outcome of single-photon true events (reco reversed)", histList)
histCanvas2, stack2, legend2, histInt2 = histStack("Outcome of two-photon true events (reco reversed", histList2)
histCanvas3, stack3, legend3, histInt3 = histStack("Outcome of three-photon true events (reco reversed)", histList3)
totalHistCanvas, stack4, legend4, histInt4 = histStack("Chart of all true events", totalList)



outFile = rt.TFile(args.outfile, "RECREATE")
totalHistCanvas.Write()
histCanvas.Write()
histCanvas2.Write()
histCanvas3.Write()


print(NCSuccess, "events passed the NC Check")
print(fiducialSuccess, "events passed the fiducial Check")
print(cosmicSuccess, "events passed the Cosmic Check")
print(pionProtonSuccess, "events passed the charged pion and proton Check")
print(trueSuccess, "events passed the entire truth check")
print(lowECount, "events passed the low energy pion/proton check")
print(recoSuccess, "events passed the reco")
