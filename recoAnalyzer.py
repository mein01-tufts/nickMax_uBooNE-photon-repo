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
signalHist = rt.TH1F("threePhotonHist", "Events with three photons",60,0,2)
backgroundHist = rt.TH1F("morePhotonHist", "Events with lots of photons", 60, 0, 2)

#set histogram axis titles and increase line width
def configureHist(h):
  h.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h.GetXaxis().SetTitle("energy (GeV")
  h.SetLineWidth(2)
  return h

#Scale the histograms
signalHist = configureHist(signalHist)
backgroundHist = configureHist(backgroundHist)

#Set detector volumes
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

#Set variables for program review
NCBackground = 0
CosmicBackground = 0
VertexBackground = 0
trueEvents = 0
NCCaught = 0

#Now we loop through the events to form the histogram!
for i in range(eventTree.GetEntries()):

  eventTree.GetEntry(i)
  #PART 1 - FINDING RECO EVENTS

  #Iterate through the showers to see if we can find a photon (in which case it seems like signal) or a muon (in which case it's almost certainly charge current, and we need to throw it out)
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

#Check tracks to check for Neutral Current
  for x in range(eventTree.nTracks):
    if eventTree.trackClassified[x] == 1:
      if eventTree.trackPID[x] == 13:
        chargeCurrent = True
      
  if photonFound == False:
    continue

  if chargeCurrent == True:
    NCCaught += 1
    continue

  

  #PART 2 - USING TRUTH TO SEPARATE SIGNAL FROM BACKGROUND
  truthConsistent = True
  if trueSignalFinder(eventTree) == False:
    truthConsistent = False
    
#HERE IS WHERE WE WILL DIVIDE THE EVENTS INTO BINS
  scaledEnergy = []
  for x in range(len(photonIDs)):
    energyGeV = eventTree.showerRecoE[photonIDs[x]]/1000
    scaledEnergy.append(energyGeV)

  if truthConsistent == True:
    signalHist.Fill(scaledEnergy[0], eventTree.xsecWeight)

  else:
    backgroundHist.Fill(scaledEnergy[0], eventTree.xsecWeight)


#----- end of event loop ---------------------------------------------#

#scale histograms to target POT
signalHist.Scale(targetPOT/ntuplePOTsum)
backgroundHist.Scale(targetPOT/ntuplePOTsum)

#Integrate hists for data
signalInt = signalHist.Integral(0,60)
backgroundInt = backgroundHist.Integral(0,60)

#Create stack histogram, add others to it
histStack = rt.THStack("histStack", "NC Histograms with Secondary Photons")

signalHist.SetLineColor(rt.kGreen)
histStack.Add(signalHist)

backgroundHist.SetLineColor(rt.kRed)
histStack.Add(backgroundHist)

totalEvents = NCBackground + CosmicBackground + VertexBackground

legend = rt.TLegend(0.7, 0.7, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates

legend.AddEntry(signalHist, "Signal: "+str(signalInt)+" Events", "l")
legend.AddEntry(backgroundHist, "Background: "+str(backgroundInt)+" Events", "l")

histCanvas = rt.TCanvas()
histStack.Draw("HIST")
histStack.GetXaxis().SetTitle("Photon Energy (GeV)")
histStack.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
legend.Draw()
rt.gPad.Update()

#create output root file and write histograms to file
outFile = rt.TFile(args.outfile, "RECREATE")
histCanvas.Write()

print("Neutral current events detected using Reco:", NCCaught)
