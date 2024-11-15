#Import the usual suspects
import sys, argparse
import numpy as np
import ROOT as rt
from helpers.larflowreco_ana_funcs import getCosThetaGravVector

#Import some functions of our own
from cuts import trueSignalFinder

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
signalHist = rt.TH2F("onePhotonSignalHist", "Signal events with one photon", 60, -1.2, 1.2, 60, 0, 1.5)
twoPhotonHist = rt.TH2F("oneRecoTwoTrueHist", "Signal events with one extra photon", 60, -1.2, 1.2, 60, 0, 1.5)
morePhotonHist = rt.TH2F("oneRecoManyTrueHist", "Signal events with two or more extra photons", 60, -1.2, 1.2, 60, 0, 1.5)
backgroundHist = rt.TH2F("onePhotonBackgroundHist", "Background events with one photon", 60, -1.2, 1.2, 60, 0, 1.5)

onePhotonHist2 = rt.TH2F("twoRecoOneTrueHist", "Signal events missing one photon", 60, -10, 10, 60, 0, 1.5)
signalHist2 = rt.TH2F("twoPhotonSignalHist", "Signal events with two photons", 60, -10, 10, 60, 0, 1.5)
morePhotonHist2 = rt.TH2F("twoRecoManyTrueHist", "Signal events with one or more extra photons", 60, -10, 10, 60, 0, 1.5)
backgroundHist2 = rt.TH2F("twoPhotonBackgroundHist", "Background events with two photons", 60, -10, 10, 60, 0, 1.5)

onePhotonHist3 = rt.TH2F("threeRecoOneTrueHist", "Signal events missing two photons", 60, -10, 10, 60, 0, 1.5)
twoPhotonHist3 = rt.TH2F("threeRecoTwoTrueHist", "Signal events missing one photon", 60, -10, 10, 60, 0, 1.5)
signalHist3 = rt.TH2F("threePhotonSignalHist", "Signal events with three or more photons", 60, -10, 10, 60, 0, 1.5)
backgroundHist3 = rt.TH2F("threePhotonBackgroundHist", "Background events with three or more photons", 60, -10, 10, 60, 0, 1.5)

#set histogram axis titles and increase line width
def configureHist(h):
  h.GetYaxis().SetTitle("")
  h.GetXaxis().SetTitle("Cosine of angle to beam")
  h.SetLineWidth(2)
  return h

#Scale the histograms
signalHist = configureHist(signalHist)
twoPhotonHist = configureHist(twoPhotonHist)
morePhotonHist = configureHist(morePhotonHist)
backgroundHist = configureHist(backgroundHist)


signalHist2 = configureHist(signalHist2)
onePhotonHist2 = configureHist(onePhotonHist2)
morePhotonHist2 = configureHist(morePhotonHist2)
backgroundHist2 = configureHist(backgroundHist2)


signalHist3 = configureHist(signalHist3)
onePhotonHist3 = configureHist(onePhotonHist3)
twoPhotonHist3 = configureHist(twoPhotonHist3)
backgroundHist3 = configureHist(backgroundHist3)

#Set detector volumes
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

#Set variables for program review
NCBackground = 0
CosmicBackground = 0
motherVertexBackground = 0
ProtonPionBackground = 0
FiducialBackground = 0
DepositBackground = 0
canItBeDone = 0

#Establish firmly held beliefs
SantaReal = True


#Now we loop through the events to form the histogram!
for i in range(eventTree.GetEntries()):

  eventTree.GetEntry(i)
  #PART 1 - FINDING RECO EVENTS

  #Iterate through the showers to see if we can find a photon (in which case it seems like signal) or a muon (in which case it's almost certainly charge current, and we need to throw it out)
  protonFound = False
  photonFound = False
  chargeCurrent = False
  chargedPionFound = False
  recoPhotonIDs = []
  for x in range(eventTree.nShowers):
    if eventTree.showerClassified[x] == 1:
      if eventTree.showerPID[x] == 22:
#        if eventTree.showerStartPosX[x] > (xMin + fiducialWidth) and eventTree.showerStartPosX[x] < (xMax - fiducialWidth) and eventTree.showerStartPosY[x] > (yMin + fiducialWidth) and eventTree.showerStartPosY[x] < (yMax - fiducialWidth) and eventTree.showerStartPosZ[x] > (zMin + fiducialWidth) and eventTree.vtxZ < (zMax - fiducialWidth):
        recoPhotonIDs.append(x)
        photonFound = True
      elif eventTree.showerPID[x] == 13:
        chargeCurrent = True

#Check tracks to check for Neutral Current, for Protons, and for Charged Pions
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

  if eventTree.foundVertex == 1:
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
      continue
    
  if photonFound == False:
    continue

  if chargeCurrent == True:
    continue

  if protonFound == True or chargedPionFound == True:
    continue

  
  
  #PART 2 - USING TRUTH TO SEPARATE SIGNAL FROM BACKGROUND
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
    badEvent = True
    NCBackground += 1
    
  #fiducial cut removes events within fiducial
  if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
    FiducialBackground += 1
    badEvent = True
    
  #skip events where all hits overlap with tagged cosmic rays
  if eventTree.vtxFracHitsOnCosmic >= 1.:
    CosmicBackground += 1
    badEvent = True


  #Determining presence of suitably energetic protons and charged pions - if they're present, the event is beyond the scope of our investigation
  if badEvent == False:
    for x in range(len(eventTree.truePrimPartPDG)):
      if eventTree.truePrimPartPDG[x] == 211:
        pionEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
        if pionEnergy >= 0.03:
          pionPresent = True
          badEvent = True
          ProtonPionBackground += 1
          break
      
      elif eventTree.truePrimPartPDG[x] == 2212:
        protonEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
        if protonEnergy >= 0.06:
          protonPresent = True
          badEvent = True
          ProtonPionBackground += 1
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
  truePhotonList = []
  for x in range(len(photonIDList)):
    if eventTree.trueSimPartEDepX[x] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[x] >= (xMax - fiducialWidth) or eventTree.trueSimPartEDepY[x] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[x] >= (yMax - fiducialWidth) or eventTree.trueSimPartEDepZ[x] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[x] >= (zMax - fiducialWidth):
      SantaReal = False
    else:
      truePhotonList.append(x)
  if len(truePhotonList) == 0:
    badEvent = True
    DepositBackground += 1
    
  #HERE IS WHERE WE WILL DIVIDE THE EVENTS INTO BINS
  scaledEnergy = []
  for x in range(len(recoPhotonIDs)):
    energyGeV = eventTree.showerRecoE[recoPhotonIDs[x]]/1000
    scaledEnergy.append(energyGeV)

  if len(recoPhotonIDs) == 1:
    if badEvent == True:
      backgroundHist.Fill(eventTree.showerCosTheta[recoPhotonIDs[0]], scaledEnergy[0], eventTree.xsecWeight)
    elif len(truePhotonList) == 1:
      signalHist.Fill(eventTree.showerCosTheta[recoPhotonIDs[0]], scaledEnergy[0], eventTree.xsecWeight)
    elif len(truePhotonList) == 2:
      twoPhotonHist.Fill(eventTree.showerCosTheta[recoPhotonIDs[0]], scaledEnergy[0], eventTree.xsecWeight)
    else:
      morePhotonHist.Fill(eventTree.showerCosTheta[recoPhotonIDs[0]], scaledEnergy[0], eventTree.xsecWeight)

#  elif len(recoPhotonIDs) == 2:
#    if badEvent == True:
#      highEnergy = scaledEnergy[0]
#      for x in len(scaledEnergy):
#        if scaledEnergy[x] > highEnergy:
#          highEnergy = scaledEnergy[x]
#      backgroundHist2.Fill(eventTree.trackCosTheta[x], highEnergy, eventTree.xsecWeight)
#    elif len(truePhotonList) ==1:
#      highEnergy = scaledEnergy[0]
#      for x in len(scaledEnergy):
#        if scaledEnergy[x] > highEnergy:
#          highEnergy = scaledEnergy[x]
#      onePhotonHist2.Fill(eventTree.trackCosTheta[x], scaledEnergy[0], eventTree.xsecWeight)
#    elif len(truePhotonList) == 2:
#      invariantMass = np.sqrt((scaledEnergy[0]*scaledEnergy[1]) - (scaledEnergy[0]*scaledEnergy[1]*eventTree.showerStartDirX[0])*(eventTree.showerStartDirX[1] + eventTree.showerStartDirY[0]*eventTree.showerStartDirY[1] + eventTree.showerStartDirZ[0]*eventTree.showerStartDirZ[1]))
#      signalHist2.Fill(eventTree.trackCosTheta[x], invariantMass, eventTree.xsecWeight)
#    else:
#      highEnergy = scaledEnergy[0]
#      for x in len(scaledEnergy):
#        if scaledEnergy[x] > highEnergy:
#          highEnergy = scaledEnergy[x]
#      morePhotonHist2.Fill(eventTree.trackCosTheta[x], scaledEnergy[0], eventTree.xsecWeight)

#  else:
#    if badEvent == True:
#      highEnergy = scaledEnergy[0]
#      for x in len(scaledEnergy):
#        if scaledEnergy[x] > highEnergy:
#          highEnergy = scaledEnergy[x]
#      backgroundHist3.Fill(eventTree.trackCosTheta[x], highEnergy, eventTree.xsecWeight)
      
#    elif len(truePhotonList) == 1:
#      highEnergy = scaledEnergy[0]
#      for x in len(scaledEnergy):
#        if scaledEnergy[x] > highEnergy:
#          highEnergy = scaledEnergy[x]
#      onePhotonHist3.Fill(eventTree.trackCosTheta[x], highEnergy, eventTree.xsecWeight)

#    elif len(truePhotonList) == 2:
#      highEnergy = scaledEnergy[0]
#      for x in len(scaledEnergy):
#        if scaledEnergy[x] > highEnergy:
#          highEnergy = scaledEnergy[x]
#      twoPhotonHist3.Fill(eventTree.trackCosTheta[x], highEnergy, eventTree.xsecWeight)

#    else:
#      highEnergy = scaledEnergy[0]
#      for x in len(scaledEnergy):
#        if scaledEnergy[x] > highEnergy:
#          highEnergy = scaledEnergy[x]
#      signalHist3.Fill(eventTree.trackCosTheta[x], highEnergy, eventTree.xsecWeight)

#----- end of event loop ---------------------------------------------#

#scale histograms to target POT
signalHist.Scale(targetPOT/ntuplePOTsum)
backgroundHist.Scale(targetPOT/ntuplePOTsum)
twoPhotonHist.Scale(targetPOT/ntuplePOTsum)
morePhotonHist.Scale(targetPOT/ntuplePOTsum)

#onePhotonHist2.Scale(targetPOT/ntuplePOTsum)
#morePhotonHist2.Scale(targetPOT/ntuplePOTsum)   
#signalHist2.Scale(targetPOT/ntuplePOTsum)
#backgroundHist2.Scale(targetPOT/ntuplePOTsum)

#onePhotonHist3.Scale(targetPOT/ntuplePOTsum)
#twoPhotonHist3.Scale(targetPOT/ntuplePOTsum)
#signalHist3.Scale(targetPOT/ntuplePOTsum)
#backgroundHist3.Scale(targetPOT/ntuplePOTsum)

#Integrate hists for data
#twoPhotonInt = twoPhotonHist.Integral(0,60)
#morePhotonInt = morePhotonHist.Integral(0,60)
#signalInt = signalHist.Integral(0,60)
#backgroundInt = backgroundHist.Integral(0,60)

#signalInt2 = signalHist2.Integral(0,60)
#backgroundInt2 = backgroundHist2.Integral(0,60)
#onePhotonInt2 = onePhotonHist2.Integral(0,60)
#morePhotonInt2 = morePhotonHist2.Integral(0,60)

#onePhotonInt3 = onePhotonHist3.Integral(0,60)
#twoPhotonInt3 = twoPhotonHist3.Integral(0,60)
#signalInt3 = signalHist3.Integral(0,60)
#backgroundInt3 = backgroundHist3.Integral(0,60)

#create output root file and write histograms to file
outFile = rt.TFile(args.outfile, "RECREATE")
signalHist.Write()
twoPhotonHist.Write()
morePhotonHist.Write()
backgroundHist.Write()

#signalHist2.Write()
#onePhotonHist2.Write()
#morePhotonHist2.Write()
#backgroundHist2.Write()

#signalHist3.Write()
#twoPhotonHist3.Write()
#onePhotonHist3.Write()
#backgroundHis3t.Write()



print("Neutral current background:", NCBackground)
print("Cosmic background:", CosmicBackground)
print("Mother particle or vertex background:", motherVertexBackground)
print("Proton or charged pion background:", ProtonPionBackground)
print("Fiducial background (event vertex-based):", FiducialBackground)
print("Fiducial background (energy deposit-based:", DepositBackground)
if SantaReal == True:
  print("We were right! Santa is real!")
else:
  print("This data does not suggest that Santa is real.")
