
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
onePhotonHist = rt.TH1F("onePhotonHist", "Events with just one photon" ,60,0,6)
twoPhotonHist = rt.TH1F("twoPhotonHist", "Events with two photons",60,0,6)
threePhotonHist = rt.TH1F("threePhotonHist", "Events with three photons",60,0,6)
manyPhotonHist = rt.TH1F("morePhotonHist", "Events with lots of photons", 60, 0, 6)

#set histogram axis titles and increase line width
def configureHist(h):
  h.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h.GetXaxis().SetTitle("energy (GeV)")
  h.SetLineWidth(2)
  return h

#Scale the histograms
onePhotonHist = configureHist(onePhotonHist)
twoPhotonHist = configureHist(twoPhotonHist)
threePhotonHist = configureHist(threePhotonHist)
manyPhotonHist = configureHist(manyPhotonHist)

#Set detector volumes
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

#Set variables for program review
problemCount = 0


#begin loop over events in ntuple file

for i in range(eventTree.GetEntries()):

  eventTree.GetEntry(i)

  #filter for neutral current
  if eventTree.trueNuCCNC != 1:
    continue

  #Remove photon and charged pion events
  pionPresent = False
  protonPresent = False                                    
  for x in range(len(eventTree.truePrimPartPDG)):
    if eventTree.truePrimPartPDG[x] == 211 and eventTree.truePrimPartE[x] >= 0.03:
      pionPresent = True
    elif eventTree.truePrimPartPDG[x] == 2212 and eventTree.truePrimPartE[x] >= 0.06:
      protonPresent = True

  if protonPresent == True or pionPresent == True:
    continue
  
  #filter for fiducial width
  if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or \
    eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or \
    eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
    continue
        
  #Filter out cosmic cuts, because that seems like the right thing to do
  #skip events where all hits overlap with tagged cosmic rays
  if eventTree.vtxFracHitsOnCosmic >= 1.:
    continue
  
  #Create sorting variables
  photonInSecondary = False
  primList = []
  photonList = []
  
  #Check if Neutral Pion and Kaon in primaries
  if 111 in eventTree.truePrimPartPDG or 311 in eventTree.truePrimPartPDG:
    #Check if there is actually a detectable photon
    for x in range(eventTree.nTrueSimParts):
      if eventTree.trueSimPartPDG[x] == 22:
        photonList.append(x)
        photonInSecondary = True
  else:
  #Create a list of prime particle Track IDs
    for x in range(len(eventTree.trueSimPartTID)):
      if eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
        primList.append(eventTree.trueSimPartTID[x])
    #Iterate through to find photons
    for x in range(len(eventTree.trueSimPartPDG)):
      if eventTree.trueSimPartPDG[x] == 22:
        photonList.append(x)
        #Check for parent particle in the primary list
        if eventTree.trueSimPartMID[x] in primList:
          photonInSecondary = True      
        #Failing that, check if the photon has coordinates within 15 mm of those of the event vertex
        elif abs(eventTree.trueSimPartX[x] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[x] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[x] -eventTree.trueVtxZ) <= 0.15:
          photonInSecondary = True
  #Discard event unless a secondary photon is found
  if photonInSecondary == False:
    continue

  #Discard event unless the photon begins to deposit energy within the fiducial volume
  for x in photonList:
    if eventTree.trueSimPartEDepX[x] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[x] >= (xMax - fiducialWidth) or eventTree.trueSimPartEDepY[x] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[x] >= (yMax - fiducialWidth) or eventTree.trueSimPartEDepZ[x] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[x] >= (zMax - fiducialWidth):
      continue
  
  #HERE IS WHERE WE WILL DIVIDE THE EVENTS INTO BINS  
  if len(photonList) == 0:
    problemCount += 1
  elif len(photonList) == 1:
    onePhotonHist.Fill(eventTree.trueNuE, eventTree.xsecWeight)
  elif len(photonList) == 2:
    twoPhotonHist.Fill(eventTree.trueNuE, eventTree.xsecWeight)
  elif len(photonList) == 3:
    threePhotonHist.Fill(eventTree.trueNuE, eventTree.xsecWeight)
  else:
    manyPhotonHist.Fill(eventTree.trueNuE, eventTree.xsecWeight)
#----- end of event loop ---------------------------------------------#

#scale histograms to target POT
onePhotonHist.Scale(targetPOT/ntuplePOTsum)
twoPhotonHist.Scale(targetPOT/ntuplePOTsum)
threePhotonHist.Scale(targetPOT/ntuplePOTsum)
manyPhotonHist.Scale(targetPOT/ntuplePOTsum)


#Create stack histogram, add others to it
histStack = rt.THStack("histStack", "NC Histograms with Secondary Photons")

onePhotonHist.SetLineColor(rt.kGreen)
histStack.Add(onePhotonHist)

twoPhotonHist.SetLineColor(rt.kRed)
histStack.Add(twoPhotonHist)

threePhotonHist.SetLineColor(rt.kMagenta)
histStack.Add(threePhotonHist)

manyPhotonHist.SetLineColor(rt.kCyan)
histStack.Add(manyPhotonHist)

legend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates

legend.AddEntry(onePhotonHist, "One photon", "l")
legend.AddEntry(twoPhotonHist, "Two photons", "l")
legend.AddEntry(threePhotonHist, "Three photons", "l")
legend.AddEntry(manyPhotonHist, "More than three photons", "l")

histCanvas = rt.TCanvas()
histStack.Draw("HIST")
legend.Draw()
rt.gPad.Update()

#create output root file and write histograms to file
outFile = rt.TFile(args.outfile, "RECREATE")
histCanvas.Write()


