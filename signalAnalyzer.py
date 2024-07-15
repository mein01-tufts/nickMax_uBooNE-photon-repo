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
noVertexHist1 = rt.TH1F("NoVertexHist1", "Only photons" ,60,0,2)
successHist1 = rt.TH1F("SuccessHist1", "Only photons" ,60,0,2)
NCFailHist1 = rt.TH1F("NCFailHist1", "Photons and protons",60,0,2)
tooManyPhotonHist1 = rt.TH1F("photonFailHist1", "Photons and charged pions",60,0,2)
tooFewPhotonHist1 = rt.TH1F("photonFailHist1", "Photons and charged pions",60,0,2)
fiducialFailHist1 = rt.TH1F("fiducialFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
pionProtonFailHist1 = rt.TH1F("pionProtonFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
noPhotonHist1 =	rt.TH1F("noPhotonHist1", "Photons and charged pions",60,0,2)
totalHist1 = rt.TH1F("totalHist1", "Photons and charged pions",60,0,2)

noVertexHist2 =	rt.TH1F("NoVertexHist2", "Only photons" ,60,0,2)
successHist2 = rt.TH1F("SuccessHist1", "Only photons" ,60,0,2)
NCFailHist2 = rt.TH1F("NCFailHist1", "Photons and protons",60,0,2)
tooManyPhotonHist2 = rt.TH1F("photonFailHist2", "Photons and charged pions",60,0,2)
tooFewPhotonHist2 = rt.TH1F("photonFailHist2", "Photons and charged pions",60,0,2)
fiducialFailHist2 = rt.TH1F("fiducialFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
pionProtonFailHist2 = rt.TH1F("pionProtonFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
noPhotonHist2 =	rt.TH1F("noPhotonHist3", "Photons and charged pions",60,0,2)
totalHist2 = rt.TH1F("totalHist2", "Photons and charged pions",60,0,2)

noVertexHist3 =	rt.TH1F("NoVertexHist3", "Only photons" ,60,0,2)
successHist3 = rt.TH1F("SuccessHist1", "Only photons" ,60,0,2)
NCFailHist3 = rt.TH1F("NCFailHist1", "Photons and protons",60,0,2)
tooManyPhotonHist3 = rt.TH1F("photonFailHist3", "Photons and charged pions",60,0,2)
tooFewPhotonHist3 = rt.TH1F("photonFailHist3", "Photons and charged pions",60,0,2)
noPhotonHist3 = rt.TH1F("noPhotonHist3", "Photons and charged pions",60,0,2)
fiducialFailHist3 = rt.TH1F("fiducialFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
pionProtonFailHist3 = rt.TH1F("pionProtonFailHist1", "Photons, protons, and charged pions", 60, 0, 2)
totalHist3 = rt.TH1F("totalHist3", "Photons and charged pions",60,0,2)

#set histogram axis titles and increase line width
def configureHist(h):
  h.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h.GetXaxis().SetTitle("energy (GeV)")
  h.SetLineWidth(2)
  return h

#Scale the histograms
noVertexHist1 = configureHist(noVertexHist1)
totalHist1 = configureHist(totalHist1)
successHist1 = configureHist(successHist1)
NCFailHist1  = configureHist(NCFailHist1)
tooManyPhotonHist1 = configureHist(tooManyPhotonHist1)
tooFewPhotonHist1 = configureHist(tooFewPhotonHist1)
fiducialFailHist1 = configureHist(fiducialFailHist1)
pionProtonFailHist1 = configureHist(pionProtonFailHist1)
noPhotonHist1 = configureHist(noPhotonHist1)

noVertexHist2 =	configureHist(noVertexHist2)
totalHist2 = configureHist(totalHist2)
successHist2 = configureHist(successHist2)
NCFailHist2  = configureHist(NCFailHist2)
tooManyPhotonHist2 = configureHist(tooManyPhotonHist2)
tooFewPhotonHist2 = configureHist(tooFewPhotonHist2)
fiducialFailHist2 = configureHist(fiducialFailHist2)
pionProtonFailHist2 = configureHist(pionProtonFailHist2)
noPhotonHist2 =	configureHist(noPhotonHist2)

noVertexHist3 =	configureHist(noVertexHist3)
totalHist3 = configureHist(totalHist3)
successHist3 = configureHist(successHist3)
NCFailHist3 = configureHist(NCFailHist3)
tooManyPhotonHist3 = configureHist(tooManyPhotonHist3)
tooFewPhotonHist3 = configureHist(tooFewPhotonHist3)
fiducialFailHist3 = configureHist(fiducialFailHist3)
pionProtonFailHist3 = configureHist(pionProtonFailHist3)
noPhotonHist3 =	configureHist(noPhotonHist3)

#Set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

#Establish variables for program evaluation
signalEvents = 0
successCount = 0
NCCuts = 0
CosmicCuts = 0
FiducialCuts = 0
motherVertexCuts = 0
DepositCuts = 0
PionProtonCuts = 0
CountCuts = 0
oddSuccessCount = 0
makesSense = 0
noPhotons = 0
noVertex = 0
noVertexNoPhotons = 0

#begin loop over events in ntuple file

for i in range(eventTree.GetEntries()):

  eventTree.GetEntry(i)
  
  #Defining/resetting our variables
  badEvent = False
  primList = []
  photonIDList = []
  protonPresent = False
  PionPresent = False
  notAllPass = False

  #checking for neutral current
  if eventTree.trueNuCCNC != 1:
    continue
  else:
    NCCuts += 1
    
  #fiducial cut removes events within fiducial
  if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
    continue
  else:
    FiducialCuts += 1
  
  #skip events where all hits overlap with tagged cosmic rays
  if eventTree.vtxFracHitsOnCosmic >= 1.:
    continue
  else:
    CosmicCuts += 1

  #Determining presence of suitably energetic protons and charged pions - if they're present, the event is beyond the scope of our investigation
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

  if badEvent == True:
    continue
  else:
    PionProtonCuts += 1
    
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
  else:
    motherVertexCuts += 1
  #Discard the photon unless it begins to deposit energy within the fiducial volume
  if badEvent == False:
    truePhotonList = []
    for x in photonIDList:
      if eventTree.trueSimPartEDepX[x] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[x] >= (xMax - fiducialWidth) or eventTree.trueSimPartEDepY[x] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[x] >= (yMax - fiducialWidth) or eventTree.trueSimPartEDepZ[x] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[x] >= (zMax - fiducialWidth):
        notAllPass = True
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
  #Fillin the histogram displaying total events
  if len(truePhotonList) == 1:
    totalHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

  elif len(truePhotonList) == 2:
    leadingPhoton = scaledEnergy[0]
    for x in range(len(scaledEnergy)):
      if scaledEnergy[x] > leadingPhoton:
        leadingPhoton = scaledEnergy[x]
    totalHist2.Fill(leadingPhoton, eventTree.xsecWeight)

  else:
    leadingPhoton = scaledEnergy[0]
    for x in range(len(scaledEnergy)):
      if scaledEnergy[x] > leadingPhoton:
        leadingPhoton = scaledEnergy[x]
    totalHist3.Fill(leadingPhoton, eventTree.xsecWeight)

  #Events with no vertex!                                                                                      
  if eventTree.foundVertex != 1:
    if len(truePhotonList) == 1:
      noVertex += 1
      if len(recoPhotonIDs) == 0:
        makesSense += 1
    if len(truePhotonList) == 1:
      noVertexHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)
    elif len(truePhotonList) == 2:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      noVertexHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      noVertexHist3.Fill(leadingPhoton, eventTree.xsecWeight)


  #Filling in the graph for no photons
  if len(recoPhotonIDs) == 0:
    if len(truePhotonList) == 1:
      noPhotons += 1
      if eventTree.foundVertex != 1:
        noVertexNoPhotons += 1
    if len(truePhotonList) == 1:
      noPhotonHist1.Fill(scaledEnergy[0], eventTree.xsecWeight)

    elif len(truePhotonList) == 2:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      noPhotonHist2.Fill(leadingPhoton, eventTree.xsecWeight)

    else:
      leadingPhoton = scaledEnergy[0]
      for x in range(len(scaledEnergy)):
        if scaledEnergy[x] > leadingPhoton:
          leadingPhoton = scaledEnergy[x]
      noPhotonHist3.Fill(leadingPhoton, eventTree.xsecWeight)
    
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

  #Filling graphs for erroneously low nonzero photon counts
  if len(truePhotonList) > len(recoPhotonIDs) and len(recoPhotonIDs) != 0:
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
    successCount += 1
    if notAllPass == True:
      oddSuccessCount += 1
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
totalHist1.Scale(targetPOT/ntuplePOTsum)
noVertexHist1.Scale(targetPOT/ntuplePOTsum)

successHist2.Scale(targetPOT/ntuplePOTsum)
NCFailHist2.Scale(targetPOT/ntuplePOTsum)
tooManyPhotonHist2.Scale(targetPOT/ntuplePOTsum)
tooFewPhotonHist2.Scale(targetPOT/ntuplePOTsum)
fiducialFailHist2.Scale(targetPOT/ntuplePOTsum)
pionProtonFailHist2.Scale(targetPOT/ntuplePOTsum)
totalHist2.Scale(targetPOT/ntuplePOTsum)
noVertexHist2.Scale(targetPOT/ntuplePOTsum)

successHist3.Scale(targetPOT/ntuplePOTsum)
NCFailHist3.Scale(targetPOT/ntuplePOTsum)
tooManyPhotonHist3.Scale(targetPOT/ntuplePOTsum)
tooFewPhotonHist3.Scale(targetPOT/ntuplePOTsum)
fiducialFailHist3.Scale(targetPOT/ntuplePOTsum)
pionProtonFailHist3.Scale(targetPOT/ntuplePOTsum)
totalHist3.Scale(targetPOT/ntuplePOTsum)
noVertexHist3.Scale(targetPOT/ntuplePOTsum)

#Create stack histograms and add the filled graphs to them
#Total Histogram
totalStack = rt.THStack("totalStack", "All events that passed true selection")
totalHist1.SetLineColor(rt.kRed)
totalStack.Add(totalHist1)

totalHist2.SetLineColor(rt.kBlue)
totalStack.Add(totalHist2)

totalHist3.SetLineColor(rt.kOrange)
totalStack.Add(totalHist3)

totalInt1 = totalHist1.Integral(0,60)
totalInt2 = totalHist2.Integral(0,60)
totalInt3 = totalHist3.Integral(0,60)


totalLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates                          

totalLegend.AddEntry(totalHist1, "One Photon "+str(round(totalInt1, 1))+" Events", "l")
totalLegend.AddEntry(totalHist2, "Two Photons "+str(round(totalInt2, 1))+" Events", "l")
totalLegend.AddEntry(totalHist3, "Three or More Photons "+str(round(totalInt3, 1))+" Events", "l")

totalCanvas = rt.TCanvas()
totalStack.Draw("HIST")
totalLegend.Draw()
rt.gPad.Update()

#One, two, and three photon stacks
onePhotonStack = rt.THStack("onePhotonStack", "All one-photon events that passed true selection")
twoPhotonStack = rt.THStack("twoPhotonStack", "All two-photon events that passed true selection")
threePhotonStack = rt.THStack("threePhotonStack", "All three-photon events that passed true selection")

onePhotonLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)
twoPhotonLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)
threePhotonLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)

successHist1.SetLineColor(rt.kGreen)
onePhotonStack.Add(successHist1)

successHist2.SetLineColor(rt.kGreen)
twoPhotonStack.Add(successHist2)

successHist3.SetLineColor(rt.kGreen)
threePhotonStack.Add(successHist3)

successInt1 = successHist1.Integral(0,60)
successInt2 = successHist2.Integral(0,60)
successInt3 = successHist3.Integral(0,60)

onePhotonLegend.AddEntry(successHist1, "Successful reconstruction: "+str(round(successInt1, 1))+" Events", "l")
twoPhotonLegend.AddEntry(successHist2, "Successful reconstruction: "+str(round(successInt2, 1))+" Events", "l")
threePhotonLegend.AddEntry(successHist3, "Successful reconstruction: "+str(round(successInt3, 1))+" Events", "l")

#Too many photons
tooManyPhotonHist1.SetLineColor(rt.kCyan)
onePhotonStack.Add(tooManyPhotonHist1)

tooManyPhotonHist2.SetLineColor(rt.kCyan)
twoPhotonStack.Add(tooManyPhotonHist2)

tooManyPhotonHist3.SetLineColor(rt.kCyan)
threePhotonStack.Add(tooManyPhotonHist3)

tooManyPhotonInt1 = tooManyPhotonHist1.Integral(0,60)
tooManyPhotonInt2 = tooManyPhotonHist2.Integral(0,60)
tooManyPhotonInt3 = tooManyPhotonHist3.Integral(0,60)

tooManyPhotonLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates

onePhotonLegend.AddEntry(tooManyPhotonHist1, "Too many photons: "+str(round(tooManyPhotonInt1, 1))+" Events", "l")
twoPhotonLegend.AddEntry(tooManyPhotonHist2, "Too many photons: "+str(round(tooManyPhotonInt2, 1))+" Events", "l")
threePhotonLegend.AddEntry(tooManyPhotonHist3, "Too many photons: "+str(round(tooManyPhotonInt3, 1))+" Events", "l")

#Too few photons
tooFewPhotonHist1.SetLineColor(rt.kOrange)
onePhotonStack.Add(tooFewPhotonHist1)

tooFewPhotonHist2.SetLineColor(rt.kOrange)
twoPhotonStack.Add(tooFewPhotonHist2)

tooFewPhotonHist3.SetLineColor(rt.kOrange)
threePhotonStack.Add(tooFewPhotonHist3)

tooFewPhotonInt1 = tooFewPhotonHist1.Integral(0,60)
tooFewPhotonInt2 = tooFewPhotonHist2.Integral(0,60)
tooFewPhotonInt3 = tooFewPhotonHist3.Integral(0,60)

tooFewPhotonLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates                    

onePhotonLegend.AddEntry(tooFewPhotonHist1, "Too few photons: "+str(round(tooFewPhotonInt1, 1))+" Events", "l")
twoPhotonLegend.AddEntry(tooFewPhotonHist2, "Too few photons: "+str(round(tooFewPhotonInt2, 1))+" Events", "l")
threePhotonLegend.AddEntry(tooFewPhotonHist3, "Too few photons: "+str(round(tooFewPhotonInt3, 1))+" Events", "l")

#Protons or Pions detected                                                                                              
pionProtonFailHist1.SetLineColor(rt.kBlue)
onePhotonStack.Add(pionProtonFailHist1)

pionProtonFailHist2.SetLineColor(rt.kBlue)
twoPhotonStack.Add(pionProtonFailHist2)

pionProtonFailHist3.SetLineColor(rt.kBlue)
threePhotonStack.Add(pionProtonFailHist3)

pionProtonFailInt1 = pionProtonFailHist1.Integral(0,60)
pionProtonFailInt2 = pionProtonFailHist2.Integral(0,60)
pionProtonFailInt3 = pionProtonFailHist3.Integral(0,60)

onePhotonLegend.AddEntry(pionProtonFailHist1, "Charged Pion/Proton false positive: "+str(round(pionProtonFailInt1, 1))+" Events", "l")
twoPhotonLegend.AddEntry(pionProtonFailHist2, "Charged Pion/Proton false positive: "+str(round(pionProtonFailInt2, 1))+" Events", "l")
threePhotonLegend.AddEntry(pionProtonFailHist3, "Charged Pion/Proton false positive: "+str(round(pionProtonFailInt3, 1))+" Events", "l")

#Outside of fiducial detected
                                                                                                              
fiducialFailHist1.SetLineColor(rt.kYellow+2)
onePhotonStack.Add(fiducialFailHist1)

fiducialFailHist2.SetLineColor(rt.kYellow+2)
twoPhotonStack.Add(fiducialFailHist2)

fiducialFailHist3.SetLineColor(rt.kYellow+2)
threePhotonStack.Add(fiducialFailHist3)

fiducialFailInt1 = fiducialFailHist1.Integral(0,60)
fiducialFailInt2 = fiducialFailHist2.Integral(0,60)
fiducialFailInt3 = fiducialFailHist3.Integral(0,60)

fiducialFailLegend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates                   

onePhotonLegend.AddEntry(fiducialFailHist1, "Misplaced out of fiducial: "+str(round(fiducialFailInt1, 1))+" Events", "l")
twoPhotonLegend.AddEntry(fiducialFailHist2, "Misplaced out of fiducial: "+str(round(fiducialFailInt2, 1))+" Events", "l")
threePhotonLegend.AddEntry(fiducialFailHist3, "Misplaced out of fiducial: "+str(round(fiducialFailInt3, 1))+" Events", "l")

#CC detected

NCFailHist1.SetLineColor(rt.kRed)
onePhotonStack.Add(NCFailHist1)

NCFailHist2.SetLineColor(rt.kRed)
twoPhotonStack.Add(NCFailHist2)

NCFailHist3.SetLineColor(rt.kRed)
threePhotonStack.Add(NCFailHist3)


NCFailInt1 = NCFailHist1.Integral(0,60)
NCFailInt2 = NCFailHist2.Integral(0,60)
NCFailInt3 = NCFailHist3.Integral(0,60)


onePhotonLegend.AddEntry(NCFailHist1, "Charged current false positive: "+str(round(NCFailInt1, 1))+" Events", "l")
twoPhotonLegend.AddEntry(NCFailHist2, "Charged current false positive: "+str(round(NCFailInt2, 1))+" Events", "l")
threePhotonLegend.AddEntry(NCFailHist3, "Charged current false positive: "+str(round(NCFailInt3, 1))+" Events" "l")


#No Vertex Found

noVertexHist1.SetLineColor(rt.kMagenta)
onePhotonStack.Add(noVertexHist1)

noVertexHist2.SetLineColor(rt.kMagenta)
twoPhotonStack.Add(noVertexHist2)

noVertexHist3.SetLineColor(rt.kMagenta)
threePhotonStack.Add(noVertexHist3)


noVertexInt1 = noVertexHist1.Integral(0,60)
noVertexInt2 = noVertexHist2.Integral(0,60)
noVertexInt3 = noVertexHist3.Integral(0,60)


onePhotonLegend.AddEntry(noVertexHist1, "No vertex found: "+str(round(noVertexInt1, 1))+" Events", "l")
twoPhotonLegend.AddEntry(noVertexHist2, "No vertex found: "+str(round(noVertexInt2, 1))+" Events", "l")
threePhotonLegend.AddEntry(noVertexHist3, "No vertex found: "+str(round(noVertexInt3, 1))+" Events" "l")

#No photons found

noPhotonHist1.SetLineColor(rt.kBlack)
onePhotonStack.Add(noPhotonHist1)

noPhotonHist2.SetLineColor(rt.kBlack)
twoPhotonStack.Add(noPhotonHist2)

noPhotonHist3.SetLineColor(rt.kBlack)
threePhotonStack.Add(noPhotonHist3)


noPhotonInt1 = noPhotonHist1.Integral(0,60)
noPhotonInt2 = noPhotonHist2.Integral(0,60)
noPhotonInt3 = noPhotonHist3.Integral(0,60)


onePhotonLegend.AddEntry(noPhotonHist1, "No photons found: "+str(round(noPhotonInt1, 1))+" Events", "l")
twoPhotonLegend.AddEntry(noPhotonHist2, "No photons found: "+str(round(noPhotonInt2, 1))+" Events", "l")
threePhotonLegend.AddEntry(noPhotonHist3, "No photons found: "+str(round(noPhotonInt3, 1))+" Events" "l")

onePhotonCanvas = rt.TCanvas()
onePhotonStack.Draw("HIST")
onePhotonLegend.Draw()
rt.gPad.Update()

twoPhotonCanvas = rt.TCanvas()
twoPhotonStack.Draw("HIST")
twoPhotonLegend.Draw()
rt.gPad.Update()

threePhotonCanvas = rt.TCanvas()
threePhotonStack.Draw("HIST")
threePhotonLegend.Draw()
rt.gPad.Update()

#create output root file and write histograms to file
outFile = rt.TFile(args.outfile, "RECREATE")
totalCanvas.Write()
onePhotonCanvas.Write()
twoPhotonCanvas.Write()
threePhotonCanvas.Write()

print(NCCuts, "events passed the neutral current test")
print(FiducialCuts, "events passed the fiducial test")
print(CosmicCuts, "events passed the cosmic test")
print(PionProtonCuts, "events passed the pion/proton check")
print(motherVertexCuts, "events passed the mother particle/vertex check")
print(signalEvents, "events passed the energy deposit check and were declared signal")
