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
photonFailHist1 = rt.TH1F("photonFailHist1", "Photons and charged pions",60,0,2)
fiducialFailHist1 = rt.TH1F("fiducialFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
pionProtonFailHist1 = rt.TH1F("pionProtonFailHist1", "Photons, protons, and charged pions", 60, 0, 2)

successHist2 = rt.TH1F("SuccessHist1", "Only photons" ,60,0,2)
NCFailHist2 = rt.TH1F("NCFailHist1", "Photons and protons",60,0,2)
photonFailHist2 = rt.TH1F("photonFailHist1", "Photons and charged pions",60,0,2)
fiducialFailHist2 = rt.TH1F("fiducialFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
pionProtonFailHist2 = rt.TH1F("pionProtonFailHist1", "Photons, protons, and charged pions", 60, 0, 2)

successHist3 = rt.TH1F("SuccessHist1", "Only photons" ,60,0,2)
NCFailHist3 = rt.TH1F("NCFailHist1", "Photons and protons",60,0,2)
photonFailHist3 = rt.TH1F("photonFailHist1", "Photons and charged pions",60,0,2)
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
photonFailHist1 = configureHist(photonFailHist1)
fiducialFailHist1 = configureHist(fiducialFailHist1)
pionProtonFailHist1 = configureHist(pionProtonFailHist1)

successHist2 = configureHist(successHist3)
NCFailHist2  = configureHist(NCFailHist3)
photonFailHist2 = configureHist(photonFailHist3)
fiducialFailHist2 = configureHist(fiducialFailHist3)
pionProtonFailHist2 = configureHist(pionProtonFailHist3)

successHist3 = configureHist(successHist3)
NCFailHist3 = configureHist(NCFailHist3)
photonFailHist3 = configureHist(photonFailHist3)
fiducialFailHist3 = configureHist(fiducialFailHist3)
pionProtonFailHist3 = configureHist(pionProtonFailHist3)

#Set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

#Establish variables for program evaluation


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
    continue
    
  #fiducial cut removes events within fiducial
  if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
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
          break
      
      elif eventTree.truePrimPartPDG[x] == 2212:
        protonEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
        if protonEnergy >= 0.06:
          protonPresent = True
          badEvent = True
          break

  #Establish a list of TIDs for primary particles (we'll use this in a moment)
  if badEvent == False:
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
    if len(photonIDList) != len(recoPhotonIDs):
      motherVertexBackground += 1

  #Discard the photon unless it begins to deposit energy within the fiducial volume
  if badEvent == False:
    truePhotonList = []
    for x in range(len(photonIDList)):
      if eventTree.trueSimPartEDepX[x] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[x] >= (xMax - fiducialWidth) or eventTree.trueSimPartEDepY[x] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[x] >= (yMax - fiducialWidth) or eventTree.trueSimPartEDepZ[x] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[x] >= (zMax - fiducialWidth):
        SantaReal = False
      else:
        truePhotonList.append(x)
    if len(truePhotonList) == 0:
      badEvent = True
      break

  if badEvent == True:
    continue
  
  #HERE IS WHERE WE WILL COMPARE THE RESULTS TO RECONSTRUCTED DATA

  #Defining our reco variables
  photonFound = False
  protonFound = False
  chargedPionFound = False
  chargeCurret = False
  recoPhotonIDs = []
  scaledEnergy = []


  for x in range(len(truePhotonIDs)):
    energyGeV = eventTree.trueSimPartE[truePhotonIDs[x]]/1000
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

  #Filling graphs for erroneous photon counts
  if len(truePhotonList) != len(recoPhotonIDs):
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

  
#----- end of event loop ---------------------------------------------#

#scale histograms to target POT
successHist1(targetPOT/ntuplePOTsum)
NCFailHist1(targetPOT/ntuplePOTsum)
photonFailHist1(targetPOT/ntuplePOTsum)
fiducialFailHist1(targetPOT/ntuplePOTsum)
pionProtonFailHist1(targetPOT/ntuplePOTsum)

successHist2(targetPOT/ntuplePOTsum)
NCFailHist2(targetPOT/ntuplePOTsum)
photonFailHist2(targetPOT/ntuplePOTsum)
fiducialFailHist2(targetPOT/ntuplePOTsum)
pionProtonFailHist2(targetPOT/ntuplePOTsum)

successHist3(targetPOT/ntuplePOTsum)
NCFailHist3(targetPOT/ntuplePOTsum)
photonFailHist3(targetPOT/ntuplePOTsum)
fiducialFailHist3(targetPOT/ntuplePOTsum)
pionProtonFailHist3(targetPOT/ntuplePOTsum)

#Create stack histogram, add others to it
histStack = rt.THStack("histStack", "NC Histograms with Secondary Photons")

soleGammaHist.SetLineColor(rt.kGreen)
histStack.Add(soleGammaHist)

protonGammaHist.SetLineColor(rt.kRed)
histStack.Add(protonGammaHist)

pionGammaHist.SetLineColor(rt.kMagenta)
histStack.Add(pionGammaHist)

protonPionGammaHist.SetLineColor(rt.kCyan)
histStack.Add(protonPionGammaHist)

legend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates

legend.AddEntry(soleGammaHist, "Only photons", "l")
legend.AddEntry(protonGammaHist, "Photons and protons", "l")
legend.AddEntry(pionGammaHist, "Photons and charged pions", "l")
legend.AddEntry(protonPionGammaHist, "Photons, charged pions, and protons", "l")

histCanvas = rt.TCanvas()
histStack.Draw("HIST")
legend.Draw()
rt.gPad.Update()

#create output root file and write histograms to file
outFile = rt.TFile(args.outfile, "RECREATE")
histCanvas.Write()

