import sys, argparse
import numpy as np
import ROOT as rt

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

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

#define histograms to fill
#we will write histograms to output file for:
successHist1 = rt.TH1F("SuccessHist1", "Only photons" ,60,0,2)
NCFailHist1 = rt.TH1F("NCFailHist1", "Photons and protons",60,0,2)
tooManyPhotonHist1 = rt.TH1F("photonFailHist1", "Photons and charged pions",60,0,2)
tooFewPhotonHist1 = rt.TH1F("photonFailHist1", "Photons and charged pions",60,0,2)
fiducialFailHist1 = rt.TH1F("fiducialFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
pionProtonFailHist1 = rt.TH1F("pionProtonFailHist1", "Photons, protons, and charged pions", 60, 0, 2)

successHist2 = rt.TH1F("SuccessHist1", "Only photons" ,60,0,2)
NCFailHist2 = rt.TH1F("NCFailHist1", "Photons and protons",60,0,2)
tooManyPhotonHist2 = rt.TH1F("photonFailHist2", "Photons and charged pions",60,0,2)
tooFewPhotonHist2 = rt.TH1F("photonFailHist2", "Photons and charged pions",60,0,2)
fiducialFailHist2 = rt.TH1F("fiducialFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
pionProtonFailHist2 = rt.TH1F("pionProtonFailHist1", "Photons, protons, and charged pions", 60, 0, 2)

successHist3 = rt.TH1F("SuccessHist1", "Only photons" ,60,0,2)
NCFailHist3 = rt.TH1F("NCFailHist1", "Photons and protons",60,0,2)
tooManyPhotonHist3 = rt.TH1F("photonFailHist3", "Photons and charged pions",60,0,2)
tooFewPhotonHist3 = rt.TH1F("photonFailHist3", "Photons and charged pions",60,0,2)
fiducialFailHist3 = rt.TH1F("fiducialFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
pionProtonFailHist3 = rt.TH1F("pionProtonFailHist1", "Photons, protons, and charged pions", 60, 0, 2)

#set histogram axis titles and increase line width
def configureHist(h):
  h.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h.GetXaxis().SetTitle("energy (GeV)")
  h.SetLineWidth(2)
  return h

#Scale the histograms
successHist1 = configureHist(successHist1)
NCFailHist1  = configureHist(NCFailHist1)
tooManyPhotonHist1 = configureHist(tooManyPhotonHist1)
tooFewPhotonHist1 = configureHist(tooFewPhotonHist1)
fiducialFailHist1 = configureHist(fiducialFailHist1)
pionProtonFailHist1 = configureHist(pionProtonFailHist1)

successHist2 = configureHist(successHist3)
NCFailHist2  = configureHist(NCFailHist3)
tooManyPhotonHist2 = configureHist(tooManyPhotonHist2)
tooFewPhotonHist2 = configureHist(tooFewPhotonHist2)
fiducialFailHist2 = configureHist(fiducialFailHist3)
pionProtonFailHist2 = configureHist(pionProtonFailHist3)

successHist3 = configureHist(successHist3)
NCFailHist3 = configureHist(NCFailHist3)
tooManyPhotonHist3 = configureHist(tooManyPhotonHist3)
tooFewPhotonHist3 = configureHist(tooFewPhotonHist3)
fiducialFailHist3 = configureHist(fiducialFailHist3)
pionProtonFailHist3 = configureHist(pionProtonFailHist3)

#Set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

#Establish variables for program evaluation
signalEvents = 0
NCCuts = 0
CosmicCuts = 0
FiducialCuts = 0
DepositCuts = 0
PionProtonCuts = 0
CountCuts = 0

#Establish firmly held beliefs
SantaReal = True


#begin loop over events in ntuple file

for i in range(eventTree.GetEntries()):

  eventTree.GetEntry(i)
  
  #Defining detector dimensions for fiducial
  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  #Defining fiducial width
  fiducialWidth = 10

  #Defining/resetting our variables
  badEvent = False
  primList = []
  photonIDList = []
  protonPresent = False
  PionPresent = False

  #checking for neutral current
  if eventTree.trueNuCCNC != 1:
    NCCuts += 1
    continue
    
  #fiducial cut removes events within fiducial
  if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
    FiducialCuts += 1
    continue
    
  #skip events where all hits overlap with tagged cosmic rays
  if eventTree.vtxFracHitsOnCosmic >= 1.:
    continue


  #Determining presence of suitably energetic protons and charged pions - if they're present, the event is beyond the scope of our investigation
  if badEvent == False:
    for x in range(len(eventTree.truePrimPartPDG)):
      if eventTree.truePrimPartPDG[x] == 211:
        pionEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
        if pionEnergy >= 0.03:
          pionPresent = True
          badEvent = True
          PionProtonCuts += 1
          break
      
      elif eventTree.truePrimPartPDG[x] == 2212:
        protonEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
        if protonEnergy >= 0.06:
          protonPresent = True
          badEvent = True
          PionProtonCuts += 1
          break

  #Establish a list of TIDs for primary particles (we'll use this in a moment)
  for x in range(len(eventTree.trueSimPartTID)):
    if eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
      primList.append(eventTree.trueSimPartTID[x])
  
  #Iterate through to find photons
  for x in range(len(eventTree.trueSimPartPDG)):
    if eventTree.trueSimPartPDG[x] == 22:
    #Check for parent particle in the primary list - if it has one, the photon is secondary. This mainly finds edge cases
      if eventTree.trueSimPartMID[x] in primList:
          photonIDList.append(x)

          #Failing that, check if the photon has coordinates within 1.5 mm of those of the event vertex
      elif abs(eventTree.trueSimPartX[x] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[x] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[x] -eventTree.trueVtxZ) <= 0.15:
        photonIDList.append(x)
    #Flag the event unless the number of secondary photons matches the number of photons found in reco
  if len(photonIDList) == 0:
    CountCuts += 1
    badEvent = True
    
  #Discard the photon unless it begins to deposit energy within the fiducial volume
  if badEvent == False:
    truePhotonList = []
    for x in photonIDList:
      if eventTree.trueSimPartEDepX[x] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[x] >= (xMax - fiducialWidth) or eventTree.trueSimPartEDepY[x] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[x] >= (yMax - fiducialWidth) or eventTree.trueSimPartEDepZ[x] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[x] >= (zMax - fiducialWidth):
        SantaReal = False
      else:
        truePhotonList.append(x)
    if len(truePhotonList) == 0:
      badEvent = True

  if badEvent == True:
    continue
  
  #HERE IS WHERE WE WILL COMPARE THE RESULTS TO RECONSTRUCTED DATA
  signalEvents += 1
  
  #Defining our reco variables
  photonFound = False
  protonFound = False
  chargedPionFound = False
  chargeCurrent = False
  fiducialFail = False
  recoPhotonIDs = []
  scaledEnergy = []


  for x in range(len(truePhotonList)):
    energyGeV = eventTree.trueSimPartE[truePhotonList[x]]/1000
    scaledEnergy.append(energyGeV)

  for x in range(eventTree.nTracks):
    if eventTree.trackClassified[x] == 1:
      if eventTree.trackPID[x] == 13:
        chargeCurrent = True
        break
      elif eventTree.trackPID[x] == 2212:
        if eventTree.trackRecoE[x] >= 60:
          protonFound = True
          break
      elif eventTree.trackPID[x] == 211:
        if eventTree.trackRecoE[x] >= 30:
          chargedPionFound = True
          break

  for x in range(eventTree.nShowers):
    if eventTree.showerClassified[x] == 1:
      if eventTree.showerPID[x] == 22:
        recoPhotonIDs.append(x)
        photonFound = True
      elif eventTree.showerPID[x] == 13:
        chargeCurrent = True

  if eventTree.foundVertex == 1:
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
      fiducialFail = True
  
  #NOW WE FILL THE HISTOGRAMS

  #Filling graphs for the erroneous detection of charged current
  if chargeCurrent == True:
    if len(truePhotonList) == 1:
      NCFailHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonList) == 2:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      NCFailHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      NCFailHist3.Fill(leadingPhoton, eventTree.xsecWeight)

  #Filling graphs for the erroneous detection of protons or pions among the primary particles
  if chargedPionFound == True or protonFound == True:
    if len(truePhotonList) == 1:
      pionProtonFailHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonList) == 2:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      pionProtonFailHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      pionProtonFailHist3.Fill(leadingPhoton, eventTree.xsecWeight)

  #Filling graphs for erroneously low photon counts
  if len(truePhotonList) > len(recoPhotonIDs):
    if len(truePhotonList) == 1:
      tooFewPhotonHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonList) == 2:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      tooFewPhotonHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      tooFewPhotonHist3.Fill(leadingPhoton, eventTree.xsecWeight)

  #Filling graphs for erroneously high photon counts
  if len(truePhotonList) < len(recoPhotonIDs):
    if len(truePhotonList) == 1:
      tooManyPhotonHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonList) == 2:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      tooManyPhotonHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      tooManyPhotonHist3.Fill(leadingPhoton, eventTree.xsecWeight)

  #Filling graphs for vertices erroneously placed outside of fiducial
  if fiducialFail == True:
    if len(truePhotonList) == 1:
      fiducialFailHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonList) == 2:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      fiducialFailHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      fiducialFailHist3.Fill(leadingPhoton, eventTree.xsecWeight)


  #Filling the success graph:
  if fiducialFail == False and len(truePhotonList) == len(recoPhotonIDs) and chargedPionFound == False and protonFound == False and chargeCurrent == False:
    if len(truePhotonList) == 1:
      successHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonList) == 2:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      successHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      successHist3.Fill(leadingPhoton, eventTree.xsecWeight)

#----- end of event loop ---------------------------------------------#

#scale histograms to target POT
successHist1.Scale(targetPOT/ntuplePOTsum)
NCFailHist1.Scale(targetPOT/ntuplePOTsum)
tooManyPhotonHist1.Scale(targetPOT/ntuplePOTsum)
tooFewPhotonHist1.Scale(targetPOT/ntuplePOTsum)
fiducialFailHist1.Scale(targetPOT/ntuplePOTsum)
pionProtonFailHist1.Scale(targetPOT/ntuplePOTsum)

successHist2.Scale(targetPOT/ntuplePOTsum)
NCFailHist2.Scale(targetPOT/ntuplePOTsum)
tooManyPhotonHist2.Scale(targetPOT/ntuplePOTsum)
tooFewPhotonHist2.Scale(targetPOT/ntuplePOTsum)
fiducialFailHist2.Scale(targetPOT/ntuplePOTsum)
pionProtonFailHist2.Scale(targetPOT/ntuplePOTsum)

successHist3.Scale(targetPOT/ntuplePOTsum)
NCFailHist3.Scale(targetPOT/ntuplePOTsum)
tooManyPhotonHist3.Scale(targetPOT/ntuplePOTsum)
tooFewPhotonHist3.Scale(targetPOT/ntuplePOTsum)
fiducialFailHist3.Scale(targetPOT/ntuplePOTsum)
pionProtonFailHist3.Scale(targetPOT/ntuplePOTsum)

#Create stack histograms and add the filled graphs to them
#Too many photons                                                                                              
successStack = rt.THStack("successStack", "Events in which the reco matched the signal")

successHist1.SetLineColor(rt.kRed)
successStack.Add(successHist1)

successHist2.SetLineColor(rt.kBlue)
successStack.Add(successHist2)

successHist3.SetLineColor(rt.kOrange)
successStack.Add(successHist3)

successLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates                    

successLegend.AddEntry(successHist1, "One Photon", "l")
successLegend.AddEntry(successHist2, "Two Photons", "l")
successLegend.AddEntry(successHist3, "Three Photons", "l")

successCanvas = rt.TCanvas()
successStack.Draw("HIST")
successLegend.Draw()
rt.gPad.Update()

#Too many photons
tooManyPhotonStack = rt.THStack("tooManyPhotonStack", "Events in which the reco found too many photons")

tooManyPhotonHist1.SetLineColor(rt.kRed)
tooManyPhotonStack.Add(tooManyPhotonHist1)

tooManyPhotonHist2.SetLineColor(rt.kBlue)
tooManyPhotonStack.Add(tooManyPhotonHist2)

tooManyPhotonHist3.SetLineColor(rt.kOrange)
tooManyPhotonStack.Add(tooManyPhotonHist3)

tooManyPhotonLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates

tooManyPhotonLegend.AddEntry(tooManyPhotonHist1, "One Photon", "l")
tooManyPhotonLegend.AddEntry(tooManyPhotonHist2, "Two Photons", "l")
tooManyPhotonLegend.AddEntry(tooManyPhotonHist3, "Three Photons", "l")

tooManyPhotonCanvas = rt.TCanvas()
tooManyPhotonStack.Draw("HIST")
tooManyPhotonLegend.Draw()
rt.gPad.Update()

#Too few photons
tooFewPhotonStack = rt.THStack("tooFewPhotonStack", "Events in which the reco did not find all photons")

tooFewPhotonHist1.SetLineColor(rt.kRed)
tooFewPhotonStack.Add(tooFewPhotonHist1)

tooFewPhotonHist2.SetLineColor(rt.kBlue)
tooFewPhotonStack.Add(tooFewPhotonHist2)

tooFewPhotonHist3.SetLineColor(rt.kOrange)
tooFewPhotonStack.Add(tooFewPhotonHist3)

tooFewPhotonLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates                    

tooFewPhotonLegend.AddEntry(tooFewPhotonHist1, "One Photon", "l")
tooFewPhotonLegend.AddEntry(tooFewPhotonHist2, "Two Photons", "l")
tooFewPhotonLegend.AddEntry(tooFewPhotonHist3, "Three Photons", "l")

tooFewPhotonCanvas = rt.TCanvas()
tooFewPhotonStack.Draw("HIST")
tooFewPhotonLegend.Draw()
rt.gPad.Update()

#Protons or Pions detected                                                                                              
pionProtonStack = rt.THStack("pionProtonStack", "Events in which the reco found pions and protons, but they did not in fact exist")

pionProtonFailHist1.SetLineColor(rt.kRed)
pionProtonStack.Add(tooManyPhotonHist1)

pionProtonFailHist2.SetLineColor(rt.kBlue)
pionProtonStack.Add(tooManyPhotonHist2)

pionProtonFailHist3.SetLineColor(rt.kOrange)
pionProtonStack.Add(tooManyPhotonHist3)

pionProtonFailLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates

pionProtonFailLegend.AddEntry(pionProtonFailHist1, "One Photon", "l")
pionProtonFailLegend.AddEntry(pionProtonFailHist2, "Two Photons", "l")
pionProtonFailLegend.AddEntry(pionProtonFailHist3, "Three Photons", "l")

pionProtonCanvas = rt.TCanvas()
pionProtonStack.Draw("HIST")
pionProtonFailLegend.Draw()
rt.gPad.Update()

#Outside of fiducial detected
                                                                                                              
fiducialStack = rt.THStack("fiducialStack", "Events in which the reco thought the event was out of the fiducial, but it was not")

fiducialFailHist1.SetLineColor(rt.kRed)
fiducialStack.Add(fiducialFailHist1)

fiducialFailHist2.SetLineColor(rt.kBlue)
fiducialStack.Add(fiducialFailHist2)

fiducialFailHist3.SetLineColor(rt.kOrange)
fiducialStack.Add(fiducialFailHist3)

fiducialFailLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates                   

fiducialFailLegend.AddEntry(fiducialFailHist1, "One Photon", "l")
fiducialFailLegend.AddEntry(fiducialFailHist2, "Two Photons", "l")
fiducialFailLegend.AddEntry(fiducialFailHist3, "Three Photons", "l")

fiducialFailCanvas = rt.TCanvas()
fiducialStack.Draw("HIST")
fiducialFailLegend.Draw()
rt.gPad.Update()


#CC detected

NCStack = rt.THStack("fiducialStack", "Events in which the reco thought the event was charge current")

NCFailHist1.SetLineColor(rt.kRed)
NCStack.Add(NCFailHist1)

NCFailHist2.SetLineColor(rt.kBlue)
NCStack.Add(NCFailHist2)

NCFailHist3.SetLineColor(rt.kOrange)
NCStack.Add(NCFailHist3)

NCFailLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates                    

NCFailLegend.AddEntry(NCFailHist1, "One Photon", "l")
NCFailLegend.AddEntry(NCFailHist2, "Two Photons", "l")
NCFailLegend.AddEntry(NCFailHist3, "Three Photons", "l")

NCFailCanvas = rt.TCanvas()
NCStack.Draw("HIST")
NCFailLegend.Draw()
rt.gPad.Update()


#create output root file and write histograms to file
outFile = rt.TFile(args.outfile, "RECREATE")
NCFailCanvas.Write()
fiducialFailCanvas.Write()
tooManyPhotonCanvas.Write()
tooFewPhotonCanvas.Write()
successCanvas.Write()
pionProtonCanvas.Write()


print("There were", signalEvents, "total signal events")
