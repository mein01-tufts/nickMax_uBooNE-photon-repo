#NOTES:
#All of these functions are meant to operate within a larger function that loops through every event (or the chosen events) within the eventTree
import sys, argparse
import numpy as np
import ROOT as rt

from helpers.larflowreco_ana_funcs import getCosThetaGravVector


def trueCutNC(ntuple):
#Filter for neutral current using truth - False if CC, True if NC
  if ntuple.trueNuCCNC != 1:
    return False
  else:
    return True

def trueCutPionProton(ntuple):
#Filter for pions and photons using truth - False if either (or both) are present, True otherwise
  pionPresent = False
  protonPresent = False                                    
  for x in range(len(ntuple.truePrimPartPDG)):
    if abs(ntuple.truePrimPartPDG[x]) == 211:
      pionEnergy = ntuple.truePrimPartE[x] - np.sqrt(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2))
      if pionEnergy > 0.03:
        pionPresent = True
        break
    elif ntuple.truePrimPartPDG[x] == 2212:
      protonEnergy = ntuple.truePrimPartE[x] - np.sqrt(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2)+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2)
      if protonEnergy > 0.06:
        protonPresent = True
        break
  if pionPresent == True or protonPresent == True:
    return False
  else:
    return True
  
def trueCutFiducials(ntuple, fiducialData):
#Filter by determining if the event vertex falls within the fiducial width using truth  - True if it's not within the radius, false if it is
  if ntuple.trueVtxX <= (fiducialData["xMin"] + fiducialData["width"]) or ntuple.trueVtxX >= (fiducialData["xMax"] - fiducialData["width"]) or ntuple.trueVtxY <= (fiducialData["yMin"] + fiducialData["width"]) or ntuple.trueVtxY >= (fiducialData["yMax"] - fiducialData["width"]) or ntuple.trueVtxZ <= (fiducialData["zMin"] + fiducialData["width"]) or ntuple.trueVtxZ >= (fiducialData["zMax"] - fiducialData["width"]):
    return False
  else:
    return True

def trueCutCosmic(ntuple):
  #skip events where all hits overlap with tagged cosmic rays
  if ntuple.vtxFracHitsOnCosmic >= 1.:
    return False
  else:
    return True
  
def trueCutPhotons(ntuple):
#Checks for the presence of at least one photon using truth
  if 22 in ntuple.trueSimPartPDG:
    return True
  else:
    return False

    
def trueCheckPionKaon(ntuple):
#Use truth to determine if a pion or kaon is present - True if either is present (or both), False if neither
  foundOne = False
  if 111 in ntuple.truePrimPartPDG or 311 in ntuple.truePrimPartPDG:
    return True
  else:
    return False


def trueCheckParentTracker(ntuple):
#Creates a list of primary particle TIDs, then checks the photons to see if they have them as parents. If so, returns true, otherwise returns false
  photonInSecondary = False
  #Create a list of prime particle Track IDs
  for x in range(len(ntuple.trueSimPartTID)):
    if ntuple.trueSimPartTID[x] == ntuple.trueSimPartMID[x]:
      primList.append(ntuple.trueSimPartTID[x])
  #Iterate through to find photons
  for x in range(len(ntuple.trueSimPartPDG)):
    if ntuple.trueSimPartPDG[x] == 22:
      #Check for parent particle in the primary list
      if ntuple.trueSimPartMID[x] in primList:
        photonInSecondary = True
  if photonInSecondary == False:
    return False
  else:
    return True

def trueCutVertex(ntuple):
#Iterates through each photon in the event using truth to determine if they originate within 15mm of the vertex; returns True if one does (or more), false if none do
  secondaryPhoton = False
  for x in range(len(ntuple.trueSimPartPDG)):
    if ntuple.trueSimPartPDG[x] == 22:
      if abs(ntuple.trueSimPartX[x] - ntuple.trueVtxX) <= 0.15 and abs(ntuple.trueSimPartY[x] - ntuple.trueVtxY) <= 0.15 and abs(ntuple.trueSimPartZ[x] - ntuple.trueVtxZ) <= 0.15:
        secondaryPhoton = True
        break
  if secondaryPhoton == True:
    return True
  else:
    return False

def trueCutEDep(ntuple, fiducialInfo):
#Uses truth to determine whether the event fires photons out of the fiducial; if so, get rid of them
#FUNCTION CURRENTLY SUSPECT
  badPhoton = False
  for x in photonList:
    if ntuple.trueSimPartEDepX[x] <= (fiducialInfo["xMin"] + fiducialInfo["width"]) or ntuple.trueSimPartEDepX[x] >= (fiducialInfo["xMax"] - fiducialInfo["width"]) or ntuple.trueSimPartEDepY[x] <= (fiducialInfo["yMin"] + fiducialInfo["width"]) or ntuple.trueSimPartEDepY[x] >= (fiducialInfo["yMax"] - fiducialInfo["fiducialWidth"]) or ntuple.trueSimPartEDepZ[x] <= (fiducialInfo["zMin"] + fiducialInfo["width"]) or ntuple.trueSimPartEDepZ[x] >= (fiducialInfo["zMax"] - fiducialInfo["width"]):
      badPhoton = True
      break
  if badPhoton == True:
    return False
  else:
    return True


def truePhotonList(ntuple, fiducial):
#Uses truth to create a list of photons that pass the vertex and deposit tests
  list1 = []
  secondaryList = []
  primList = []
  for x in range(len(ntuple.trueSimPartTID)):
    if ntuple.trueSimPartTID[x] == ntuple.trueSimPartMID[x]:
      primList.append(ntuple.trueSimPartTID[x])
  for x in range(ntuple.nTrueSimParts):
    if ntuple.trueSimPartPDG[x] == 22:
      if ntuple.trueSimPartMID[x] in primList:
        secondaryList.append(x)
      elif abs(ntuple.trueSimPartX[x] - ntuple.trueVtxX) <= 0.15 and abs(ntuple.trueSimPartY[x] - ntuple.trueVtxY) <= 0.15 and abs(ntuple.trueSimPartZ[x] - ntuple.trueVtxZ) <= 0.15:
        secondaryList.append(x)
  for x in secondaryList:
    if ntuple.trueSimPartEDepX[x] > (fiducial["xMin"] + fiducial["width"]) and ntuple.trueSimPartEDepX[x] < (fiducial["xMax"] - fiducial["width"]) and ntuple.trueSimPartEDepY[x] > (fiducial["yMin"] + fiducial["width"]) and ntuple.trueSimPartEDepY[x] < (fiducial["yMax"] - fiducial["width"]) and ntuple.trueSimPartEDepZ[x] > (fiducial["zMin"] + fiducial["width"]) and ntuple.trueSimPartEDepZ[x] < (fiducial["zMax"] - fiducial["width"]):
      list1.append(x)  
  return list1

def trueCutMaxProtons(ntuple):
#Filter for pions and photons using truth - False if either (or both) are present, True otherwise                      
  pionPresent = False
  protonPresent = False
  for x in range(len(ntuple.truePrimPartPDG)):
    if abs(ntuple.truePrimPartPDG[x]) == 211:
      pionEnergy = ntuple.truePrimPartE[x] - np.sqrt(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2))
      if pionEnergy > 0.03:
        pionPresent = True
        break
    elif ntuple.truePrimPartPDG[x] == 2212:
      protonEnergy = ntuple.truePrimPartE[x] - np.sqrt(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2)+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2)
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


def scaleRecoEnergy(ntuple, recoIDs):
  #Uses reconstructed variables to return scaled energy and invariant mass (if the photon count is not two, the invariant mass defaults to -1)
  scaledEnergy = []
  invariantMass = -1
  for x in recoIDs:
    energyGeV = ntuple.showerRecoE[x]/1000
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


def scaleTrueEnergy(ntuple, trueIDs):
  scaledEnergy = []
  invariantMass = -1
  for x in trueIDs:
    energyGeV = ntuple.trueSimPartE[x]/1000
    scaledEnergy.append(energyGeV)

  leadingPhoton = scaledEnergy[0]
  for x in range(len(scaledEnergy)):
    if scaledEnergy[x] > leadingPhoton:
      leadingPhoton = scaledEnergy[x]

    if len(trueIDs) == 2:
      a = 0
      b = 1
      lengthA = np.sqrt(ntuple.trueSimPartPx[a]**2 + ntuple.trueSimPartPy[a]**2 + ntuple.trueSimPartPz[a]**2)
      lengthB = np.sqrt(ntuple.trueSimPartPx[b]**2 + ntuple.trueSimPartPy[b]**2 + ntuple.trueSimPartPz[b]**2)
      aDotB = ntuple.trueSimPartPx[a]*ntuple.trueSimPartPx[b] + ntuple.trueSimPartPy[a]*ntuple.trueSimPartPy[b] + ntuple.trueSimPartPz[a]*ntuple.trueSimPartPz[b]
      invariantMass = np.sqrt((scaledEnergy[a]*scaledEnergy[b]) - (scaledEnergy[a]*scaledEnergy[b]*(aDotB)/(lengthA*lengthB)))
  return leadingPhoton, invariantMass

#RECO FUNCTIONS
def recoNoVertex(ntuple):
  #Checks to see if the event has a reconstructed vertex; returns True if so, False if not
  if ntuple.foundVertex != 1:
    return False
  else:
    return True


def recoFiducials(ntuple, fiducial):
  #Checks to see if the reconstructed event vertex is within the fiducial volume
  if ntuple.foundVertex == 1:
    if ntuple.vtxX > (fiducial["xMin"] + fiducial["width"]) and ntuple.vtxX < (fiducial["xMax"] - fiducial["width"]) and ntuple.vtxY > (fiducial["yMin"] + fiducial["width"]) and ntuple.vtxY < (fiducial["yMax"] - fiducial["width"]) and ntuple.vtxZ > (fiducial["zMin"] + fiducial["width"]) and ntuple.vtxZ < (fiducial["zMax"] - fiducial["width"]):
      return True
    else:
      return False


def recoPhotonList(ntuple):
  #Creates a list of photons based on the showers in the event
  recoIDs = []
  for x in range(ntuple.nShowers):
    if ntuple.showerClassified[x] == 1:
      if ntuple.showerPID[x] == 22:
        recoList.append(x)
  return recoIDs

def recoPionProton(ntuple):
  #Checks for sufficiently energetic protons and charged pions in the event
  chargedPionFound = False
  protonFound = False
  for x in range(ntuple.nTracks):
    if ntuple.trackPID[x] == 2212:
      if ntuple.trackRecoE[x] >= 60:
        protonFound = True
        break
    elif abs(ntuple.trackPID[x]) == 211:
      if ntuple.trackRecoE[x] >= 30:
        chargedPionFound = True
        break
  if protonFound == True or chargedPionFound == True:
    return False
  else:
    return True

def recoNeutralCurrent(ntuple):
  #Checks for signs of a neutral current event; returns True if it thinks the event is NC, False if CC
  chargeCurrent = False
  for x in range(ntuple.nTracks):
    if ntuple.trackClassified[x] == 1:
      if ntuple.trackPID[x] == 13 and ntuple.trackRecoE[x] > 20:
        chargeCurrent = True
        break
  for x in range(ntuple.nShowers):
    if ntuple.showerClassified[x] == 1:
      if ntuple.showerPID[x] == 11 and ntuple.showerRecoE[x] > 10:
        chargeCurrent = True
        break
  if chargeCurrent == True:
    return False
  else:
    return True


def trickyPionProtonCuts(ntuple):
#Filter for pions and photons using truth - False if either (or both) are present, True otherwise
#NOTE - This is deliberately checking for photons and charged pions the WRONG WAY, in order to inspect the specific set of events with protons that fall between this threshold and the actual kinetic energy threshold. If you're looking for the right way to do it, trueCutPionProton is the way to go (or you might want trueCutMaxProton if you're looking for protons but not charged pions) 
  pionPresent = False
  protonPresent = False
  for x in range(len(ntuple.truePrimPartPDG)):
    if abs(ntuple.truePrimPartPDG[x]) == 211:
      if ntuple.truePrimPartE[x] > 0.03:
        pionPresent = True
        break
    elif ntuple.truePrimPartPDG[x] == 2212:
      if ntuple.truePrimPartE[x] > 0.06:
        protonPresent = True
        break
  if pionPresent == True or protonPresent == True:
    return True
  else:
    return False

# This function uses a particle's momentum vector and total energy to calculate its kinetic energy, then returns its kinetic energy in GeV
def kineticEnergyCalculator(ntuple, i):
  momentumVector = np.square(ntuple.truePrimPartPx[i]) + np.square(ntuple.truePrimPartPy[i]) + np.square(ntuple.truePrimPartPz[i])
  kineticGeV = ntuple.truePrimPartE[i] - np.sqrt((np.square(ntuple.truePrimPartE[i])) - momentumVector)
  return kineticGeV
