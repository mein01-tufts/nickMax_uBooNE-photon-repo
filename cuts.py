#NOTES:
#All of these functions are meant to operate within a larger function that loops through every event (or the chosen events) within the eventTree
import sys, argparse
import numpy as np
import ROOT as rt

from helpers.larflowreco_ana_funcs import getCosThetaGravVector


def trueCutNC(eventTree):
#Filter for neutral current using truth - False if CC, True if NC
  if eventTree.trueNuCCNC != 1:
    return False
  else:
    return True

def trueCutPionProton(eventTree):
#Filter for pions and photons using truth - False if either (or both) are present, True otherwise
  pionPresent = False
  protonPresent = False                                    
  for x in range(len(eventTree.truePrimPartPDG)):
    if abs(eventTree.truePrimPartPDG[x]) == 211:
      pionEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
      if pionEnergy > 0.03:
        pionPresent = True
        break
    elif eventTree.truePrimPartPDG[x] == 2212:
      protonEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2)+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2)
      if protonEnergy > 0.06:
        protonPresent = True
        break
  if pionPresent == True or protonPresent == True:
    return False
  else:
    return True
  
def trueCutFiducials(eventTree, fiducialData):
#Filter by determining if the event vertex falls within the fiducial width using truth  - True if it's not within the radius, false if it is
  if eventTree.trueVtxX <= (fiducialData["xMin"] + fiducialData["width"]) or eventTree.trueVtxX >= (fiducialData["xMax"] - fiducialData["width"]) or eventTree.trueVtxY <= (fiducialData["yMin"] + fiducialData["width"]) or eventTree.trueVtxY >= (fiducialData["yMax"] - fiducialData["width"]) or eventTree.trueVtxZ <= (fiducialData["zMin"] + fiducialData["width"]) or eventTree.trueVtxZ >= (fiducialData["zMax"] - fiducialData["width"]):
    return False
  else:
    return True

def trueCutCosmic(eventTree):
  #skip events where all hits overlap with tagged cosmic rays
  if eventTree.vtxFracHitsOnCosmic >= 1.:
    return False
  else:
    return True
  
def trueCutPhotons(eventTree):
#Checks for the presence of at least one photon using truth
  if 22 in eventTree.trueSimPartPDG:
    return True
  else:
    return False

    
def trueCheckPionKaon(eventTree):
#Use truth to determine if a pion or kaon is present - True if either is present (or both), False if neither
  foundOne = False
  if 111 in eventTree.truePrimPartPDG or 311 in eventTree.truePrimPartPDG:
    return True
  else:
    return False


def trueCheckParentTracker(eventTree):
#Creates a list of primary particle TIDs, then checks the photons to see if they have them as parents. If so, returns true, otherwise returns false
  photonInSecondary = False
  #Create a list of prime particle Track IDs
  for x in range(len(eventTree.trueSimPartTID)):
    if eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
      primList.append(eventTree.trueSimPartTID[x])
  #Iterate through to find photons
  for x in range(len(eventTree.trueSimPartPDG)):
    if eventTree.trueSimPartPDG[x] == 22:
      #Check for parent particle in the primary list
      if eventTree.trueSimPartMID[x] in primList:
        photonInSecondary = True
  if photonInSecondary == False:
    return False
  else:
    return True

def trueCutVertex(eventTree):
#Iterates through each photon in the event using truth to determine if they originate within 15mm of the vertex; returns True if one does (or more), false if none do
  secondaryPhoton = False
  for x in range(len(eventTree.trueSimPartPDG)):
    if eventTree.trueSimPartPDG[x] == 22:
      if abs(eventTree.trueSimPartX[x] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[x] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[x] -eventTree.trueVtxZ) <= 0.15:
        secondaryPhoton = True
        break
  if secondaryPhoton == True:
    return True
  else:
    return False

def trueCutEDep(eventTree):
#Uses truth to determine whether the event fires photons out of the fiducial; if so, get rid of them
#FUNCTION CURRENTLY SUSPECT
  badPhoton = False
  for x in photonList:
    if eventTree.trueSimPartEDepX[x] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[x] >= (xMax - fiducialWidth) or eventTree.trueSimPartEDepY[x] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[x] >= (yMax - fiducialWidth) or eventTree.trueSimPartEDepZ[x] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[x] >= (zMax - fiducialWidth):
      badPhoton = True
      break
  if badPhoton == True:
    return False
  else:
    return True


def truePhotonList(eventTree, list1, fiducial):
#Uses truth to create a list of photons that pass the vertex and deposit tests
  list1 = []
  secondaryList = []
  primList = []
  for x in range(len(eventTree.trueSimPartTID)):
    if eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
      primList.append(eventTree.trueSimPartTID[x])
  for x in range(eventTree.nTrueSimParts):
    if eventTree.trueSimPartPDG[x] == 22:
      if eventTree.trueSimPartMID[x] in primList:
        secondaryList.append(x)
      elif abs(eventTree.trueSimPartX[x] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[x] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[x] -eventTree.trueVtxZ) <= 0.15:
        secondaryList.append(x)
  for x in secondaryList:
    if eventTree.trueSimPartEDepX[x] > (fiducial["xMin"] + fiducial["width"]) and eventTree.trueSimPartEDepX[x] < (fiducial["xMax"] - fiducial["width"]) and eventTree.trueSimPartEDepY[x] > (fiducial["yMin"] + fiducial["width"]) and eventTree.trueSimPartEDepY[x] < (fiducial["yMax"] - fiducial["width"]) and eventTree.trueSimPartEDepZ[x] > (fiducial["zMin"] + fiducial["width"]) and eventTree.trueSimPartEDepZ[x] < (fiducial["zMax"] - fiducial["width"]):
      list1.append(x)  
  return list1

def trueCutMaxProtons(eventTree):
#Filter for pions and photons using truth - False if either (or both) are present, True otherwise                      
  pionPresent = False
  protonPresent = False
  for x in range(len(eventTree.truePrimPartPDG)):
    if abs(eventTree.truePrimPartPDG[x]) == 211:
      pionEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
      if pionEnergy > 0.03:
        pionPresent = True
        break
    elif eventTree.truePrimPartPDG[x] == 2212:
      protonEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2)+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2)
      if protonEnergy > 0.06:
        protonPresent = True
        break
  if pionPresent == False and protonPresent == True:
    return True
  else:
    return False


#HISTOGRAM FUNCTIONS
def histStack(title, histList):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt. kBlue, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kBlack, rt.kYellow, rt.kGreen, rt. kOrange+1]
  targetPOT = 6.67e+20
  ntuplePOTSum = 4.675690535431973e+20
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    hist.Scale(targetPOT/ntuplePOTSum)
    hist.SetLineColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
    stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas() 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle("Photon Energy (GeV)")
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

def histStackFill(title, histList, legendTitle, xTitle, yTitle):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt. kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kBlack, rt.kYellow, rt.kGreen]
  targetPOT = 6.67e+20
  ntuplePOTSum = 4.675690535431973e+20
  integralSum = 0
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    hist.Scale(targetPOT/ntuplePOTSum)
    hist.SetFillColor(colors[x%7])
    hist.SetMarkerStyle(21)
    hist.SetMarkerColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    integralSum += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 2))+" events per 6.67e+20 POT", "f")
    stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas(str(title)) 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle(str(xTitle))
  stack.GetYaxis().SetTitle(str(yTitle))
  legendHeaderString = str(str(legendTitle) + str(round((integralSum),1)) + " Events per 6.67e+20 POT)") 
  legend.SetHeader(str(legendHeaderString), "C")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

def sStackFillS(title, hist, kColor, canvasTitle ):
#Forms a filled stacked histogram based on only one. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.45, 0.8, 0.9, 0.9)
  targetPOT = 6.67e+20
  ntuplePOTSum = 4.675690535431973e+20
  bins = hist.GetNbinsX()
  hist.Scale(targetPOT/ntuplePOTSum)
  hist.SetFillColor(kColor)
  hist.SetMarkerStyle(21)
  hist.SetMarkerColor(kColor)
  histInt = hist.Integral(1, int(bins))
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1))+" events per 6.67e+20 POT", "f")
  stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas(str(canvasTitle)) 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle("Average Photon Energy (MeV)")
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

def sStackFillNS(title, hist, kColor, canvasTitle ):
#Same as above, except doesn't scale (designed to be used after the input histogram has already been scaled in the program)
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.45, 0.8, 0.9, 0.9)
  bins = hist.GetNbinsX()
  hist.SetFillColor(kColor)
  hist.SetMarkerStyle(21)
  hist.SetMarkerColor(kColor)
  histInt = hist.Integral(1, int(bins))
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1))+" events per 6.67e+20 POT", "f")
  stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas(str(canvasTitle)) 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle("Average Photon Energy (MeV)")
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  stack.GetYaxis().SetRange(0,150)
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt


def scaleRecoEnergy(eventTree, recoIDs):
  #Uses reconstructed variables to return scaled energy and invariant mass (if the photon count is not two, the invariant mass defaults to -1)
  scaledEnergy = []
  invariantMass = -1
  for x in recoIDs:
    energyGeV = eventTree.showerRecoE[x]/1000
    scaledEnergy.append(energyGeV)

  leadingPhoton = scaledEnergy[0]
  for x in range(len(scaledEnergy)):
    if scaledEnergy[x] > leadingPhoton:
      leadingPhoton = scaledEnergy[x]

  if len(recoIDs) == 2:
    a = 0
    b = 1
    invariantMass = np.sqrt((scaledEnergy[a]*scaledEnergy[b]) - (scaledEnergy[a]*scaledEnergy[b]*(eventTree.showerStartDirX[a]*eventTree.showerStartDirX[b] + eventTree.showerStartDirY[a]*eventTree.showerStartDirY[b] + eventTree.showerStartDirZ[a]*eventTree.showerStartDirZ[b])))
  
  return leadingPhoton, invariantMass


def scaleTrueEnergy(eventTree, truePhotonIDs):
  scaledEnergy = []
  invariantMass = -1
  for x in truePhotonIDs:
    energyGeV = eventTree.trueSimPartE[x]/1000
    scaledEnergy.append(energyGeV)

  leadingPhoton = scaledEnergy[0]
  for x in range(len(scaledEnergy)):
    if scaledEnergy[x] > leadingPhoton:
      leadingPhoton = scaledEnergy[x]

    if len(truePhotonIDs) == 2:
      a = 0
      b = 1
      lengthA = np.sqrt(eventTree.trueSimPartPx[a]**2 + eventTree.trueSimPartPy[a]**2 + eventTree.trueSimPartPz[a]**2)
      lengthB = np.sqrt(eventTree.trueSimPartPx[b]**2 + eventTree.trueSimPartPy[b]**2 + eventTree.trueSimPartPz[b]**2)
      aDotB = eventTree.trueSimPartPx[a]*eventTree.trueSimPartPx[b] + eventTree.trueSimPartPy[a]*eventTree.trueSimPartPy[b] + eventTree.trueSimPartPz[a]*eventTree.trueSimPartPz[b]
      invariantMass = np.sqrt((scaledEnergy[a]*scaledEnergy[b]) - (scaledEnergy[a]*scaledEnergy[b]*(aDotB)/(lengthA*lengthB)))
  return leadingPhoton, invariantMass

#RECO FUNCTIONS
def recoNoVertex(eventTree):
  #Checks to see if the event has a reconstructed vertex; returns True if so, False if not
  if eventTree.foundVertex != 1:
    return False
  else:
    return True


def recoFiducials(eventTree, fiducial):
  #Checks to see if the reconstructed event vertex is within the fiducial volume
  if eventTree.foundVertex == 1:
    if eventTree.vtxX > (fiducial["xMin"] + fiducial["width"]) and eventTree.vtxX < (fiducial["xMax"] - fiducial["width"]) and eventTree.vtxY > (fiducial["yMin"] + fiducial["width"]) and eventTree.vtxY < (fiducial["yMax"] - fiducial["width"]) and eventTree.vtxZ > (fiducial["zMin"] + fiducial["width"]) and eventTree.vtxZ < (fiducial["zMax"] - fiducial["width"]):
      return True
    else:
      return False


def recoPhotonList(eventTree):
  #Creates a list of photons based on the showers in the event
  recoList = []
  for x in range(eventTree.nShowers):
    if eventTree.showerClassified[x] == 1:
      if eventTree.showerPID[x] == 22:
        recoList.append(x)
  return recoList

def recoPionProton(eventTree):
  #Checks for sufficiently energetic protons and charged pions in the event
  chargedPionFound = False
  protonFound = False
  for x in range(eventTree.nTracks):
    if eventTree.trackPID[x] == 2212:
      if eventTree.trackRecoE[x] >= 60:
        protonFound = True
        break
    elif abs(eventTree.trackPID[x]) == 211:
      if eventTree.trackRecoE[x] >= 30:
        chargedPionFound = True
        break
  if protonFound == True or chargedPionFound == True:
    return False
  else:
    return True

def recoNeutralCurrent(eventTree):
  #Checks for signs of a neutral current event; returns True if it thinks the event is NC, False if CC
  chargeCurrent = False
  for x in range(eventTree.nTracks):
    if eventTree.trackClassified[x] == 1:
      if eventTree.trackPID[x] == 13 and eventTree.trackRecoE[x] > 20:
        chargeCurrent = True
        break
  for x in range(eventTree.nShowers):
    if eventTree.showerClassified[x] == 1:
      if eventTree.showerPID[x] == 11 and eventTree.showerRecoE[x] > 10:
        chargeCurrent = True
        break
  if chargeCurrent == True:
    return False
  else:
    return True


def trickyPionProtonCuts(eventTree):
#Filter for pions and photons using truth - False if either (or both) are present, True otherwise
#NOTE - This is deliberately checking for photons and charged pions the WRONG WAY, in order to inspect the specific set of events with protons that fall between this threshold and the actual kinetic energy threshold. If you're looking for the right way to do it, trueCutPionProton is the way to go (or you might want trueCutMaxProton if you're looking for protons but not charged pions) 
  pionPresent = False
  protonPresent = False
  for x in range(len(eventTree.truePrimPartPDG)):
    if abs(eventTree.truePrimPartPDG[x]) == 211:
      if eventTree.truePrimPartE[x] > 0.03:
        pionPresent = True
        break
    elif eventTree.truePrimPartPDG[x] == 2212:
      if eventTree.truePrimPartE[x] > 0.06:
        protonPresent = True
        break
  if pionPresent == True or protonPresent == True:
    return True
  else:
    return False
