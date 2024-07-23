import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials, trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoCutLowEnergy, recoPion, recoProton

from helpers.larflowreco_ana_funcs import getCosThetaGravVector

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-c", "--cosmicFile", type=str, required=True, help="cosmic input file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")

args = parser.parse_args()

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

cosmic_file = rt.TFile(args.cosmicFile)
cosmicTree = cosmic_file.Get("EventTree")
cosmicPotTree = ntuple_file.Get("cosmicPotTree")

#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20"
ntuplePOTsum = 0
for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

cosmicPOTsum = 3.2974607516739725e+20

#Hists created and organized here

#Signal Histograms (purity)
signal1 = rt.TH1F("signal1", "One-Photon Signal",50,0,50)
signal2 = rt.TH1F("signal2", "Two-Photon Signal",50,0,50)
signal3 = rt.TH1F("signal3", "3+ Photon Signal",50,0,50)
signalHistList = [signal1, signal2, signal3]

#Non-Cosmic Background Histograms (purity)
background1 = rt.TH1F("background1", "One-Photon Background",50,0,50)
background2 = rt.TH1F("background2", "Two-Photon Backround",50,0,50)
background3 = rt.TH1F("background3", "3+ Photon Background",50,0,50)
backgroundHistList = [background1, background2, background3]

#Cosmics go here, so we can put them on purity
cosmic1 = rt.TH1F("cBackground1", "One Photon Cosmic",50,0,50)
cosmic2 = rt.TH1F("cBackground2", "Two Photon Cosmic",50,0,50)
cosmic3 = rt.TH1F("cBackground3", "3+ Photon Cosmic",50,0,50)
cosmicHistList = [cosmic1, cosmic2, cosmic3]

histList1 = [signal1, background1]
histList2 = [signal2, background2]
histList3 = [signal3, background3]

purityHist1 = rt.TH1F("purity1", "Single-Photon Purity",50,0,50)

#Signal Histograms (efficiency)
effSignal1 = rt.TH1F("effsignal1", "One-Photon Signal",50,0,50)
effSignal2 = rt.TH1F("effsignal2", "Two-Photon Signal",50,0,50)
effSignal3 = rt.TH1F("effsignal3", "3+ Photon Signal",50,0,50)
effSignalHistList = [effSignal1, effSignal2, effSignal3]

#Non-Cosmic Background Histograms
effBackground1 = rt.TH1F("effbackground1", "One-Photon Background",50,0,50)
effBackground2 = rt.TH1F("effbackground2", "Two-Photon Backround",50,0,50)
effBackground3 = rt.TH1F("effbackground3", "3+ Photon Background",50,0,50)
effBackgroundHistList = [effBackground1, effBackground2, effBackground3]

effList1 = [effSignal1, effBackground1]
effList2 = [effSignal2, effBackground2]
effList3 = [effSignal3, effBackground3]

efficiencyHist1 = rt.TH1F("overallEfficiency1", "Single-Photon Efficiency",50,0,50)

#Built-in functions here
def addHist(eventTree, photonList, histList, variable, weight):
  if len(photonList) == 1:
    histList[0].Fill(variable, weight)
  elif len(photonList) == 2:
    histList[1].Fill(variable, weight)
  else:
    histList[2].Fill(variable, weight)

def purityStack(title, purityList, cosmicHist, POTSum, cosmicSum):
  #Create Component Variables
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt.kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kGreen+4, rt.kOrange+1]
  POTTarget = 6.67e+20
  #Organize other histograms
  for x in range(len(purityList)):
    hist = purityList[x]
    bins = hist.GetNbinsX()
    hist.Scale(POTTarget/POTSum)
    hist.SetLineColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
    stack.Add(hist)

  #Format and add the cosmic histogram
  hist = cosmicHist
  bins = hist.GetNbinsX()
  hist.Scale(POTTarget/cosmicSum)
  hist.SetLineColor(rt.kBlack)
  histInt = hist.Integral(1, int(bins))
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  stack.Add(hist)
  #Finish working on the Canvas and return necessary components
  histCanvas = rt.TCanvas() 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle("Photon Energy (GeV)")
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()
  return histCanvas, stack, legend, histInt

def histScale(hist, POTSum):
  POTTarget = 6.67e+20
  hist.Scale(POTTarget/POTSum)
  return hist

def countIncrease(photonList, countList):
  if len(photonList) == 1:
    countList[0] += 1
  elif len(photonList) == 2:
    countList[1] += 1
  else:
    countList[2] += 1

#Variables for program review


#Variables for program function


#We run the entire program once for each fiducial measurement from 1-50. The Login Node is going to love it
for z in range(1, 51):
  fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":z}
  
  #Defining/resetting the graphing variables
  signalList = [0, 0, 0]
  backgroundList = [0, 0, 0]
  cosmicList = [0, 0, 0]

  #BEGINNING EVENT LOOP FOR DEFAULT PURITY
  for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

    #Selecting events with reco
    #See if the event has a vertex
    if recoNoVertex(eventTree) == False:
      continue

    #See if the event is neutral current
    if recoNeutralCurrent(eventTree) == False:
      continue

    #Make sure the event is within the fiducial volume
    if recoFiducials(eventTree, fiducialData) == False:
      continue

    #Cut events with suitably energetic protons or charged pions 
    if recoPionProton(eventTree) == False:
      continue

    #See if there are any photons in the event - if so, list them
    recoList = recoPhotonList(eventTree)
    if len(recoList) == 0:
      continue
  
    #PURITY - GRAPHING BASED ON TRUTH
    #Neutral current!
    if trueCutNC(eventTree) == False:
      countIncrease(recoList, backgroundList)
      continue

    #Fiducials!
    if trueCutFiducials(eventTree, fiducialData) == False:
      countIncrease(recoList, backgroundList)
      continue

    #pions and protons!
    if trueCutPionProton(eventTree) == False:
      countIncrease(recoList, backgroundList)
      continue

    #Now we make a list of the actual photons!
    truePhotonIDs = truePhotonList(eventTree, fiducialData)

    #Are there actually any photons?
    if len(truePhotonIDs) == 0:
      countIncrease(recoList, backgroundList)
      continue
  
    #Is the number of photons correct?
    if len(recoList) == len(truePhotonIDs):
      countIncrease(recoList, signalList)
    else:
      countIncrease(recoList, backgroundList)

  #BEGINNING EVENT LOOP FOR COSMICS
  for i in range(cosmicTree.GetEntries()):
    cosmicTree.GetEntry(i)

    #COSMICS - SELECTING EVENTS BASED ON RECO
    #See if the event has a vertex
    if recoNoVertex(cosmicTree) == False:
      continue

    #See if the event is neutral current                                                                                                  
    if recoNeutralCurrent(cosmicTree) == False:
      continue

    #Make sure the event is within the fiducial volume
    if recoFiducials(cosmicTree, fiducialData) == False:
      continue

    #Cut events with suitably energetic protons or charged pions
    if recoPionProton(cosmicTree) == False:
      continue

    #See if there are any photons in the event - if so, list them
    recoList = recoPhotonList(cosmicTree)
    if len(recoList) == 0:
      continue

    #graphing based on photon count
    #Calculating graphing values
    countIncrease(recoList, cosmicList)

  #Fill histograms for purity
  for n in range(3):
    signalHistList[n].SetBinContent(z, signalList[n])
    backgroundHistList[n].SetBinContent(z, backgroundList[n])
    cosmicHistList[n].SetBinContent(z, cosmicList[n])

  #BEGINNING EVENT LOOP FOR EFFICIENCY
  #Reset variables 
  signalList = [0, 0, 0]
  backgroundList = [0, 0, 0]

  for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

    #Selecting events using truth
    if trueCutNC(eventTree) == False:
      continue

    if trueCutFiducials(eventTree, fiducialData) == False:
      continue
    
    if trueCutCosmic(eventTree) == False:
      continue
    
    if trueCutPionProton(eventTree) == False:
      continue
  
    truePhotonIDs = truePhotonList(eventTree, fiducialData)

    if len(truePhotonIDs) == 0:
      continue

  
    #EFFICIENCY - GRAPHING BASED ON RECO

    if recoNoVertex(eventTree) == False:
      countIncrease(truePhotonIDs, backgroundList)
      continue

    #See if the event is neutral current                                                                                                  
    if recoNeutralCurrent(eventTree) == False:
      countIncrease(truePhotonIDs, backgroundList)
      continue

    #Cut events with vertexes outside the fiducial
    if recoFiducials(eventTree, fiducialData) == False:
      countIncrease(truePhotonIDs, backgroundList)
      continue

    #Cut events with suitably energetic protons or charged pions
    if recoProton(eventTree) == False:
      countIncrease(truePhotonIDs, backgroundList)
      continue

    if recoPion(eventTree) == False:
      countIncrease(truePhotonIDs, backgroundList)
      continue

    #See if there are any photons in the event - if so, list them
    recoList = recoPhotonList(eventTree)
    if len(recoList) == len(truePhotonIDs):
      countIncrease(truePhotonIDs, signalList)
    else:
      countIncrease(truePhotonIDs, backgroundList)

  #Fill Efficiency Histograms
  for n in range(3):
    effSignalHistList[n].SetBinContent(z, signalList[n])
    effBackgroundHistList[n].SetBinContent(z, backgroundList[n])

#CALCULATING EFFICIENCY
for x in range(1, 50):
  signal = signal1.GetBinContent(x)
  total = signal1.GetBinContent(x) + background1.GetBinContent(x) + cosmic1.GetBinContent(x)
  if total > 0:
    purity = signal/total
  else:
    purity = 0
  purityHist1.Fill(x, purity)

  signal = effSignal1.GetBinContent(x)
  total = effSignal1.GetBinContent(x) + effBackground1.GetBinContent(x)
  if total > 0:
    efficiency = signal/total
  else:
    efficiency = 0
  efficiencyHist1.Fill(x, efficiency)

#LOOPS OVER - HISTOGRAM ORGANIZING NOW
Canvas1, Stack1, Legend1, Int1 = purityStack("Purity of Single-Photon Events over Fiducial", histList1, cosmic1, ntuplePOTsum, cosmicPOTsum)
Canvas2, Stack2, Legend2, Int2 = purityStack("Purity of Two-Photon Events over Fiducial", histList2, cosmic2, ntuplePOTsum, cosmicPOTsum)
Canvas3, Stack3, Legend3, Int3 = purityStack("Purity of 3+ Photon Events over Fiducial", histList3, cosmic3, ntuplePOTsum, cosmicPOTsum)

effCanvas1, effStack1, effLegend1, effInt1 = histStack("Efficiency of Single-Photon Events over Fiducial", effList1, ntuplePOTsum)
effCanvas2, effStack2, effLegend2, effInt2 = histStack("Efficiency of Two-Photon Events over Fiducial", effList2, ntuplePOTsum)
effCanvas3, effStack3, effLegend3, effInt3 = histStack("Efficiency of 3+ Photon Events over Fiducial", effList3, ntuplePOTsum)

canvasList = [Canvas1, Canvas2, Canvas3, effCanvas1, effCanvas2, effCanvas3]
stackList = [Stack1, Stack2, Stack3, effStack1, effStack2, effStack3]

for stack in stackList:
  stack.GetXaxis().SetTitle("Fiducial Width")

#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in canvasList:
  canvas.Write()
efficiencyHist1.Write()
purityHist1.Write()