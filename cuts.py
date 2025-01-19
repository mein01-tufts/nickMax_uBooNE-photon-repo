#NOTES:
#All of these functions are meant to operate within a larger function that loops through every event (or the chosen events) within the eventTree
import sys, argparse
import numpy as np
import ROOT as rt
import math

from helpers.larflowreco_ana_funcs import getCosThetaGravVector


def trueParticleTallies(ntuple):
    electronIndex, photonIndex, muonIndex, pionIndex, protonIndex, otherIndex = [], [], [], [], [], [] 
    for x in range(ntuple.nTrueSimParts):
      if ntuple.trueSimPartPDG[x] == 13:
        muonIndex.append(x)
      elif abs(ntuple.trueSimPartPDG[x]) == 11:
        electronIndex.append(x)
      elif abs(ntuple.trueSimPartPDG[x]) == 2212:
        protonIndex.append(x)
      elif ntuple.trueSimPartPDG[x] == 22:
        photonIndex.append(x)
      elif abs(ntuple.trueSimPartPDG[x]) == 211:
        pionIndex.append(x)
      else:
        otherIndex.append(x)
    return electronIndex, photonIndex, muonIndex, pionIndex, protonIndex, otherIndex

def trueCutNC(ntuple):
#Filter for neutral current using truth - False if CC, True if NC
  chargePartPresent = False
  for x in range(ntuple.nTruePrimParts):
    if ntuple.truePrimPartPDG[x] == 13:
      particleEnergy = ntuple.truePrimPartE[x] - np.sqrt(abs(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2)))
      if particleEnergy >= 0.1:
        chargePartPresent = True
        break
    elif ntuple.truePrimPartPDG[x] == 11:
      particleEnergy = ntuple.truePrimPartE[x] - np.sqrt(abs(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2)))
      if particleEnergy >= 0.01:
        chargePartPresent = True
        break
  if chargePartPresent == True:
    return False
  else:
    return True

def trueCutMuons(ntuple):
  #Cut events with over-threshold muons
  muonPresent = False
  for x in range(ntuple.nTruePrimParts):
    if ntuple.truePrimPartPDG[x] == 13:
      particleEnergy = ntuple.truePrimPartE[x] - np.sqrt(abs(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2)))
      if particleEnergy >= 0.1:
        muonPresent = True
        break
  if muonPresent == True:
    return False
  else:
    return True

def trueCutElectrons(ntuple):
  #Cut events with over-threshold electrons
  electronPresent = False
  for x in range(ntuple.nTruePrimParts):
    if ntuple.truePrimPartPDG[x] == 11:
      particleEnergy = ntuple.truePrimPartE[x] - np.sqrt(abs(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2)))
      if particleEnergy >= 0.01:
        electronPresent = True
        break
  if electronPresent == True:
    return False
  else:
    return True


def trueCutPionProton(ntuple):
#Filter for pions and photons using truth - False if either (or both) are present, True otherwise
  pionsPresent = False
  protonsPresent = 0                                    
  for x in range(len(ntuple.truePrimPartPDG)):
    if abs(ntuple.truePrimPartPDG[x]) == 211:
      pionEnergy = ntuple.truePrimPartE[x] - np.sqrt(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2))
      if pionEnergy > 0.05:
        pionsPresent += 1
    elif ntuple.truePrimPartPDG[x] == 2212:
      protonEnergy = ntuple.truePrimPartE[x] - np.sqrt(ntuple.truePrimPartE[x]**2 - ((ntuple.truePrimPartPx[x]**2)+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2))
      if protonEnergy > 0.1:
        protonsPresent += 1
  return pionsPresent, protonsPresent


def trueCutProtonInclusive(ntuple):
  #Filter for pions and photons using truth - False if either (or both) are present, True otherwise
  pionPresent = False
  protonPresent = 0                                   
  for x in range(len(ntuple.truePrimPartPDG)):
    if abs(ntuple.truePrimPartPDG[x]) == 211:
      pionEnergy = ntuple.truePrimPartE[x] - np.sqrt(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2))
      if pionEnergy > 0.05:
        pionPresent = True
        break
    elif ntuple.truePrimPartPDG[x] == 2212:
      protonEnergy = ntuple.truePrimPartE[x] - np.sqrt(ntuple.truePrimPartE[x]**2 - (ntuple.truePrimPartPx[x]**2)+ntuple.truePrimPartPy[x]**2+ntuple.truePrimPartPz[x]**2)
      if protonEnergy > 0.1:
        protonPresent += 1
  if pionPresent == True:
    return False
  else:
    return protonPresent

def trueCutFiducials(ntuple, fiducialData):
#Filter by determining if the event vertex falls within the fiducial width using truth  - True if it's not within the radius, false if it is
  if ntuple.trueVtxX <= (fiducialData["xMin"] + fiducialData["width"]) or ntuple.trueVtxX >= (fiducialData["xMax"] - fiducialData["width"]) or ntuple.trueVtxY <= (fiducialData["yMin"] + fiducialData["width"]) or ntuple.trueVtxY >= (fiducialData["yMax"] - fiducialData["width"]) or ntuple.trueVtxZ <= (fiducialData["zMin"] + fiducialData["width"]) or ntuple.trueVtxZ >= (fiducialData["zMax"] - fiducialData["width"]):
    return False
  else:
    return True

def trueCutBottomlessFiducial(ntuple, fiducialData):
#Filter by determining if the event vertex falls within the fiducial width using truth, excluding the lower y-axis cut  - True if it's not within the radius, false if it is
  if ntuple.trueVtxX <= (fiducialData["xMin"] + fiducialData["width"]) or ntuple.trueVtxX >= (fiducialData["xMax"] - fiducialData["width"]) or ntuple.trueVtxY >= (fiducialData["yMax"] - fiducialData["width"]) or ntuple.trueVtxZ <= (fiducialData["zMin"] + fiducialData["width"]) or ntuple.trueVtxZ >= (fiducialData["zMax"] - fiducialData["width"]):
    return False
  else:
    return True

def trueCutCosmic(ntuple):
  #skip events where all hits overlap with tagged cosmic rays
  if ntuple.vtxFracHitsOnCosmic < 1:
    return True
  else:
    return False
  
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
    if ntuple.trueSimPartEDepX[x] > (fiducial["xMin"] + fiducial["photonWidth"]) and ntuple.trueSimPartEDepX[x] < (fiducial["xMax"] - fiducial["photonWidth"]) and ntuple.trueSimPartEDepY[x] > (fiducial["yMin"] + fiducial["photonWidth"]) and ntuple.trueSimPartEDepY[x] < (fiducial["yMax"] - fiducial["photonWidth"]) and ntuple.trueSimPartEDepZ[x] > (fiducial["zMin"] + fiducial["photonWidth"]) and ntuple.trueSimPartEDepZ[x] < (fiducial["zMax"] - fiducial["photonWidth"]):
      pixelList = [ntuple.trueSimPartPixelSumUplane[x], ntuple.trueSimPartPixelSumVplane[x], ntuple.trueSimPartPixelSumYplane[x]]
      pixelEnergy = max(pixelList)*0.0126
      if pixelEnergy > 20:
        list1.append(x)
  return list1

def trueBottomlessPhotonList(ntuple, fiducial):
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
    if ntuple.trueSimPartEDepX[x] > (fiducial["xMin"] + fiducial["width"]) and ntuple.trueSimPartEDepX[x] < (fiducial["xMax"] - fiducial["width"]) and ntuple.trueSimPartEDepY[x] < (fiducial["yMax"] - fiducial["width"]) and ntuple.trueSimPartEDepZ[x] > (fiducial["zMin"] + fiducial["width"]) and ntuple.trueSimPartEDepZ[x] < (fiducial["zMax"] - fiducial["width"]):
      list1.append(x)
  return list1

def trueCutOverlapPhotonList(ntuple, fiducial):
  list1 = []
  secondaryList = []
  primList = []
  #Check directly to see if the photon is secondary (won't catch pi0 photons)
  for x in range(len(ntuple.trueSimPartTID)):
    if ntuple.trueSimPartTID[x] == ntuple.trueSimPartMID[x]:
      primList.append(ntuple.trueSimPartTID[x])
  #Check somewhat less directly
  for x in range(ntuple.nTrueSimParts):
    if ntuple.trueSimPartPDG[x] == 22:
      if ntuple.trueSimPartMID[x] in primList:
        secondaryList.append(x)
      elif abs(ntuple.trueSimPartX[x] - ntuple.trueVtxX) <= 0.15 and abs(ntuple.trueSimPartY[x] - ntuple.trueVtxY) <= 0.15 and abs(ntuple.trueSimPartZ[x] - ntuple.trueVtxZ) <= 0.15:
        secondaryList.append(x)
  #Check for fiducial
  for x in secondaryList:
    if ntuple.trueSimPartEDepX[x] > (fiducial["xMin"] + fiducial["width"]) and ntuple.trueSimPartEDepX[x] < (fiducial["xMax"] - fiducial["width"]) and ntuple.trueSimPartEDepY[x] > (fiducial["yMin"] + fiducial["width"]) and ntuple.trueSimPartEDepY[x] < (fiducial["yMax"] - fiducial["width"]) and ntuple.trueSimPartEDepZ[x] > (fiducial["zMin"] + fiducial["width"]) and ntuple.trueSimPartEDepZ[x] < (fiducial["zMax"] - fiducial["width"]):
      #Now check to see if this photon could have overlapped with another (for our purposes, if they fell within 13 degrees of each other)
      goodPhoton = True
      for y in list1:
        dotProduct = ntuple.trackStartDirX[x]*ntuple.trackStartDirX[y] + ntuple.trackStartDirY[x]*ntuple.trackStartDirY[y] + ntuple.trackStartDirZ[x]*ntuple.trackStartDirZ[y]
        #Since cos(theta) = (dot product)/(length1*length2), and both vectors have length 1, we check to see if the dot product is greater than cos(13)
        if dotProduct < 0.9743700648:
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
#def massHistMake(histDict, bincount, lowx, highx):
#  for key in histDisct.keys():


def histStack(histName, title, histList, POTSum, axisLabel="Photon Energy (GeV)"):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack(str("Histname"), str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt. kBlue,  rt. kOrange+1, rt.kRed, rt.kCyan, rt.kMagenta, rt.kYellow+1, rt.kBlack, rt.kViolet]
  POTTarget = 4.4e+19
  histIntTotal = 0
  #Adds the histograms to the stack
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    hist.Scale(POTTarget/POTSum)
    #Make sure signal is the only one with green, for easy identification
    if x == 0:
      hist.SetLineColor(rt.kGreen+2)
    #The rest get random colors
    else:
      hist.SetLineColor(colors[x%7])
    #Now we add to the stack
    stack.Add(hist)
  #Making the legend entry for the signal first, so it goes on top
  hist = histList[0]
  histInt = hist.Integral(1, int(bins))
  histIntTotal += histInt
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  #Adding the rest of the legend entries (has to be backwards here, or the legend order will be reversed on the graph)
  listLength = (len(histList) - 1)
  for x in range(listLength, 0, -1):
    hist = histList[x]
    histInt = hist.Integral(1, int(bins))
    histIntTotal += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  legendHeaderString = "Total: " + str(round((histIntTotal),1)) 
  legend.SetHeader(str(legendHeaderString), "C")
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas() 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle(axisLabel)
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

def histStackDark(title, histList, POTSum, axisLabel="Photon Energy (GeV)"):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt.kGreen, rt. kBlue,  rt. kOrange+1, rt.kViolet+3, rt.kRed, rt.kCyan, rt.kMagenta, rt.kYellow+1, rt.kBlack, rt.kViolet]
  POTTarget = 6.67e+20
  histIntTotal = 0
  #Adds the histograms to the stack
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    hist.SetLineColor(colors[x%9])
    #Now we add to the stack
    stack.Add(hist)
  #Writing up the legend
  listLength = (len(histList) - 1)
  for x in range(listLength, -1, -1):
    hist = histList[x]
    histInt = hist.Integral(1, int(bins))
    histIntTotal += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  legendHeaderString = "Total: " + str(round((histIntTotal),1)) 
  legend.SetHeader(str(legendHeaderString), "C")
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas() 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle(axisLabel)
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt


def histStackTwoSignal(title, histList, POTSum):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt. kBlue, rt.kRed, rt.kCyan, rt.kMagenta, rt.kYellow+2, rt.kBlack, rt.kYellow, rt.kViolet, rt. kOrange+1]
  POTTarget = 4.4e+19
  histIntTotal = 0
  #Adds the histograms to the stack
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    hist.Scale(POTTarget/POTSum)
    #Make sure the signals are the only ones with green, for easy identification
    if x == 0:
      hist.SetLineColor(rt.kGreen)
    elif x == 1:
      hist.SetLineColor(rt.kGreen + 4)
    #The rest get random colors
    else:
      hist.SetLineColor(colors[x%9])
    #Now we add to the stack
    stack.Add(hist)
  #Making the legend entry for the signal first, so it goes on top
  hist = histList[0]
  histInt = hist.Integral(1, int(bins))
  histIntTotal += histInt
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  #Do it again for the second signal
  hist = histList[1]
  histInt = hist.Integral(1, int(bins))
  histIntTotal += histInt
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  #Adding the rest of the legend entries (has to be backwards here, or the legend order will be reversed on the graph)
  listLength = (len(histList) - 1)
  for x in range(listLength, 1, -1):
    hist = histList[x]
    histInt = hist.Integral(1, int(bins))
    histIntTotal += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  legendHeaderString = "Total: " + str(round((histIntTotal),1)) 
  legend.SetHeader(str(legendHeaderString), "C")
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas() 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle("Photon Energy (GeV)")
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

def histStackNoScale(title, histList, POTSum):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written. 
  #This one is designed for 1 Gamma + 0 and 1 Gamma + 1P data, so no scaling is implimented
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt. kBlue, rt.kRed, rt.kCyan, rt.kMagenta, rt.kYellow+2, rt.kBlack, rt.kYellow, rt.kViolet, rt. kOrange+1]
  POTTarget = 6.67e+20
  histIntTotal = 0
  #Adds the histograms to the stack
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    #Make sure the signals are the only ones with green, for easy identification
    if x == 0:
      hist.SetLineColor(rt.kGreen)
    elif x == 1:
      hist.SetLineColor(rt.kGreen + 4)
    #The rest get random colors
    else:
      hist.SetLineColor(colors[x%9])
    #Now we add to the stack
    stack.Add(hist)
  #Making the legend entry for the signal first, so it goes on top
  hist = histList[0]
  histInt = hist.Integral(1, int(bins))
  histIntTotal += histInt
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  #Do it again for the second signal
  hist = histList[1]
  histInt = hist.Integral(1, int(bins))
  histIntTotal += histInt
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  #Adding the rest of the legend entries (has to be backwards here, or the legend order will be reversed on the graph)
  listLength = (len(histList) - 1)
  for x in range(listLength, 1, -1):
    hist = histList[x]
    histInt = hist.Integral(1, int(bins))
    histIntTotal += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  legendHeaderString = "Total: " + str(round((histIntTotal),1)) 
  legend.SetHeader(str(legendHeaderString), "C")
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas() 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle("Photon Energy (GeV)")
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt

def histStackFill(title, histList, legendTitle, xTitle, yTitle, ntuplePOTSum):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.35, 0.5, 0.9, 0.9)
  colors = [rt.kRed, rt.kOrange, rt.kYellow+2, rt.kCyan, rt. kBlue, rt.kMagenta, rt.kViolet, rt.kBlack]
  targetPOT = 4.4e+19
  integralSum = 0
  sum = 0
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    hist.Scale(targetPOT/ntuplePOTSum)
    histInt = hist.Integral(1, int(bins))
    sum += histInt
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    if x == 0:
        hist.SetFillColor(rt.kGreen+2)
        hist.SetMarkerStyle(21)
        hist.SetMarkerColor(rt.kGreen+2)
    else:
        hist.SetFillColor(colors[x%7])
        hist.SetMarkerStyle(21)
        hist.SetMarkerColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    integralSum += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round((histInt/sum*100), 2))+" percent of events", "f")
    stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas(str(title)) 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle(str(xTitle))
  stack.GetYaxis().SetTitle(str(yTitle))
  legendHeaderString = str(str(legendTitle) + str(round((integralSum),1)) + " per 4.4e+19 POT)") 
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
  stack.GetXaxis().SetTitle("Leading Photon Energy (MeV)")
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

#takes a total event space histogram, signal histogram, blank hist to be filled with the ratio of signal to total
#returns an efficiency plot of the event reconstruction
def efficiencyPlot(totalHist, signalHist, ratioHist, title, xTitle):
  totalHist.Scale(6.67e+20/4.675690535431973e+20)
  signalHist.Scale(6.67e+20/4.675690535431973e+20)
  ratioHist.Scale(6.67e+20/4.675690535431973e+20)
  ratioHist = (signalHist)/(totalHist)
  efficiencyStack = rt.THStack("EfficiencyStack", str(title))
  ratioHist.SetFillColor(rt.kGreen+2)
  efficiencyStack.Add(ratioHist)
  efficiencyCanvas = rt.TCanvas("Efficiency Plot") 
  efficiencyStack.Draw("HIST")
  efficiencyStack.GetXaxis().SetTitle(str(xTitle))
  efficiencyStack.GetYaxis().SetTitle("Percent Efficiency")
  efficiencyCanvas.Update()

  return efficiencyCanvas, efficiencyStack  

def scaleRecoEnergy(ntuple, recoIDs, recoIDs2):
  #Uses reconstructed variables to return scaled energy and invariant mass (if the photon count is not two, the invariant mass defaults to -1)
  scaledEnergy = []
#  invariantMass = -1
  for x in recoIDs:
    energyGeV = ntuple.showerRecoE[x]/1000
    scaledEnergy.append(energyGeV)

  for x in recoIDs2:
    energyGeV = ntuple.trackRecoE[x]/1000
    scaledEnergy.append(energyGeV)

  highestPhoton = scaledEnergy[0]
  for x in range(len(scaledEnergy)):
    if scaledEnergy[x] > highestPhoton:
      highestPhoton = scaledEnergy[x]

  return highestPhoton


def scaleTrueEnergy(ntuple, trueIDs):
  scaledEnergy = []
#  invariantMass = -1
  for x in trueIDs:
    energyGeV = ntuple.trueSimPartE[x]/1000
    scaledEnergy.append(energyGeV)

  leadingPhoton = scaledEnergy[0]
  for x in range(len(scaledEnergy)):
    if scaledEnergy[x] > leadingPhoton:
      leadingPhoton = scaledEnergy[x]

#    if len(trueIDs) == 2:
#      a = 0
#      b = 1
#      lengthA = np.sqrt(ntuple.trueSimPartPx[a]**2 + ntuple.trueSimPartPy[a]**2 + ntuple.trueSimPartPz[a]**2)
#      lengthB = np.sqrt(ntuple.trueSimPartPx[b]**2 + ntuple.trueSimPartPy[b]**2 + ntuple.trueSimPartPz[b]**2)
#      aDotB = ntuple.trueSimPartPx[a]*ntuple.trueSimPartPx[b] + ntuple.trueSimPartPy[a]*ntuple.trueSimPartPy[b] + ntuple.trueSimPartPz[a]*ntuple.trueSimPartPz[b]
#      invariantMass = np.sqrt((scaledEnergy[a]*scaledEnergy[b]) - (scaledEnergy[a]*scaledEnergy[b]*(aDotB)/(lengthA*lengthB)))
  return leadingPhoton#, invariantMass

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


def recoBottomlessFiducials(ntuple, fiducial):
  #Checks to see if the reconstructed event vertex is within the fiducial volume, excluding the lower side
  if ntuple.foundVertex == 1:
    if ntuple.vtxX > (fiducial["xMin"] + fiducial["width"]) and ntuple.vtxX < (fiducial["xMax"] - fiducial["width"]) and ntuple.vtxY < (fiducial["yMax"] - fiducial["width"]) and ntuple.vtxZ > (fiducial["zMin"] + fiducial["width"]) and ntuple.vtxZ < (fiducial["zMax"] - fiducial["width"]):
      return True
    else:
      return False


def recoPhotonList(ntuple, threshold = 0):
  #Creates a list of photons based on the showers in the event
  recoIDs = []
  for x in range(ntuple.nShowers):
    if threshold == 0:
      if ntuple.showerClassified[x] == 1:
       if ntuple.showerPID[x] == 22:
        recoIDs.append(x) 
    elif ntuple.showerSize[x] > threshold:
      if ntuple.showerPID[x] == 22:
        recoIDs.append(x)
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
    if abs(ntuple.trackPID[x]) == 211:
      if ntuple.trackRecoE[x] >= 30:
        chargedPionFound = True
        break
  if protonFound == True or chargedPionFound == True:
    return False
  else:
    return True


def recoProton(ntuple, threshold = 0):
  protonsFound = 0
  #Go through tracks to see if we can find sufficiently energetic protons
  for x in range(ntuple.nTracks):
    if ntuple.trackPID[x] == 2212:
      if threshold != 0:
        if ntuple.trackSize[x] < threshold:
          continue
      distxyz = np.sqrt((ntuple.trackStartPosX[x] - ntuple.trackEndPosX[x])**2 + (ntuple.trackStartPosY[x] - ntuple.trackEndPosY[x])**2 + (ntuple.trackStartPosZ[x] - ntuple.trackEndPosZ[x])**2)
      #A very high percentage of photons reconstructed below threshold but with tracks longer than 7.5 cm are actually above threshold
      if ntuple.trackRecoE[x] > 100 or distxyz > 7.5:
        protonsFound += 1
  for x in range(ntuple.nShowers):
    if threshold != 0:
        if ntuple.showerSize[x] < threshold:
          continue
    if ntuple.showerPID[x] == 2212:
      if ntuple.showerRecoE[x] > 100:
        protonsFound += 1
  return protonsFound

def recoPion(ntuple, threshold = 0):
  pionFound = False
  for x in range(ntuple.nTracks):
    if abs(ntuple.trackPID[x]) == 211:
      if threshold != 0:
        if ntuple.trackSize[x] < threshold:
          continue
      if ntuple.trackRecoE[x] >= 50:
        pionFound = True
        break
  for x in range(ntuple.nShowers):
    if abs(ntuple.showerPID[x]) == 211:
      if threshold != 0:
        if ntuple.showerSize[x] < threshold:
          continue
      if ntuple.showerRecoE[x] >= 50:
        pionFound = True
        break
  if pionFound == True:
    return False
  else: return True
  
def recoNeutralCurrent(ntuple, threshold = 0):
  #Checks for signs of a neutral current event; returns True if it thinks the event is NC, False if CC
  chargeCurrent = False
  for x in range(ntuple.nTracks):
    if threshold == 0:
      if ntuple.trackClassified[x] == 1 and ntuple.trackProcess[x] == 0:
        chargeCurrent = True
        break
    elif ntuple.trackSize[x] > threshold and ntuple.trackProcess[x] == 0:
      if ntuple.trackPID[x] == 13 and ntuple.trackRecoE[x] > 100:
        chargeCurrent = True
        break
      if ntuple.trackPID[x] == 11 and ntuple.trackRecoE[x] > 10:
        chargeCurrent = True
        break

  for x in range(ntuple.nShowers):
    if threshold == 0:
      if ntuple.showerClassified[x] == 1 and ntuple.showerProcess[x] == 0:
        chargeCurrent = True
        break
    elif ntuple.showerSize[x] > threshold and ntuple.showerProcess[x] == 0:
      if ntuple.showerPID[x] == 11 and ntuple.showerRecoE[x] > 10:
        chargeCurrent = True
        break
      elif ntuple.showerPID[x] == 13 and ntuple.showerRecoE[x] >  100:
        chargeCurrent = True
        break
  if chargeCurrent == True:
    return False
  else:
    return True

def recoCutMuons(ntuple, threshold = 0):
  #Cuts out events with any tracks or showers identified as muons
  muonPresent = False
  for x in range(ntuple.nTracks):
    #Check to see if we have a custom track classification threshold
    if threshold == 0:
      #If not, look at all classified tracks that have the right predicted origin point
      if ntuple.trackClassified[x] == 1 and ntuple.trackProcess[x] == 0:
        if ntuple.trackPID[x] == 13 and ntuple.trackRecoE[x] > 100:
          muonPresent = True
          break
    #If so, look at any tracks exceeding the threshold that have the right predicted origin point
    elif ntuple.trackSize[x] > threshold and ntuple.trackProcess[x] == 0:
      if ntuple.trackPID[x] == 13 and ntuple.trackRecoE[x] > 100:
        muonPresent = True
        break
  for x in range(ntuple.nShowers):
    #Now we do the same thing for showers, because muons sometimes get classified as showers for some reason
    if threshold == 0:
      if ntuple.showerClassified[x] == 1 and ntuple.showerProcess[x] == 0:
        if ntuple.showerPID[x] == 13 and ntuple.showerRecoE[x] > 100:
          muonPresent = True
          break
    elif ntuple.showerSize[x] > threshold and ntuple.showerProcess[x] == 0:
      if ntuple.showerPID[x] == 13 and ntuple.showerRecoE[x] > 100:
        muonPresent = True
        break
  if muonPresent == True:
    return False
  else:
    return True

def recoCutElectrons(ntuple, threshold = 0):
  #Cut any events containing primary electrons 
  electronPresent = False
  for x in range(ntuple.nTracks):
    if threshold == 0:
      if ntuple.trackClassified[x] == 1 and ntuple.trackProcess[x] == 0:
        if ntuple.trackPID[x] == 11 and ntuple.trackRecoE[x] > 10:
          electronPresent = True
          break
    elif ntuple.trackSize[x] >threshold and ntuple.trackProcess[x] == 0:
      if ntuple.trackPID[x] == 11 and ntuple.trackRecoE[x] > 10:
        electronPresent = True
        break
  for x in range(ntuple.nShowers):
    if threshold == 0:
      if ntuple.showerClassified[x] == 1 and ntuple.showerProcess[x] == 0:
        if ntuple.showerPID[x] == 11 and ntuple.showerRecoE[x] > 10:
          electronPresent = True
          break
    elif ntuple.showerSize[x] > threshold and ntuple.showerProcess[x] == 0:
      if ntuple.showerPID[x] == 11 and ntuple.showerRecoE[x] > 10:
        electronPresent = True
        break
  if electronPresent == True:
    return False
  else:
    return True

def recoCutLowEnergy(recoList, ntuple):
#Seeks to sharply reduce background in single-photon events by slicing out events in which the photon has less than 0.1 GeV, as a considerable number of cosmic events fit that bill
  if len(recoList) == 1:
    if ntuple.showerRecoE[recoList[0]] < 70:
      return False
    else:
      return True
  else:
    return True

def CCSeeker(ntuple, recoPhotons):
  chargeParent = False
  for x in recoPhotons:
    if ntuple.showerProcess[x] == 2:
      chargeParent = True
  if chargeParent == True:
    return False
  else:
    return True
    
def recoCutElectronScore(ntuple, recoPhotons, recoPhotons2):
  #NOTE: THIS THING IS NOT WORTH IT. Cuts single photon events with high enough electron scores. Setting the threshold to 4 eliminates most cosmics at the cost of a great deal of signal
  if len(recoPhotons) + len(recoPhotons2) == 1:
    if len(recoPhotons) == 1:
      if abs(ntuple.showerElScore[recoPhotons[0]]) < 2:
        return False
      else:
        return True
    elif len(recoPhotons2) == 1:
      if abs(ntuple.trackElScore[recoPhotons2[0]]) < 2:
        return False
      else:
        return True

def recoCutShowerFromChargeScore(ntuple, recoPhotons, recoPhotons2):
  #Cuts any single-photon event that has a Shower-From-Charged Score between -5 and 0. Targets cosmics very efficiently
  if len(recoPhotons) + len(recoPhotons2) == 1:
    if len(recoPhotons) == 1 and abs(ntuple.showerFromChargedScore[recoPhotons[0]]) < 5:
      return False
    elif len(recoPhotons2) == 1 and abs(ntuple.trackFromChargedScore[recoPhotons2[0]]) < 5:
      return False
    else:
      return True

def recoCutShowerfromNeutralScore(ntuple, recoPhotons, recoPhotons2):
  if len(recoPhotons) + len(recoPhotons2) == 1:
    if len(recoPhotons) == 1 and abs(ntuple.showerFromNeutralScore[recoPhotons[0]]) > 0.4:
      return False
    elif len(recoPhotons2) == 1 and abs(ntuple.trackFromNeutralScore[recoPhotons2[0]]) > 0.4:
      return False
    else:
      return True

def recoCutMuScore(ntuple, recoPhotons):
  allGood = True
  if ntuple.nTracks > 0:
    for x in range(ntuple.nTracks):
      if abs(ntuple.trackMuScore[x]) < 1:
        allGood = False
        break 
  if allGood == False:
    return False
  else:
    return True


def recoCutCompleteness(ntuple, recoPhotons, recoPhotons2):
  #Cuts any single-photon event that has a Shower-From-Charged Score between -5 and 0. Targets cosmics very efficiently
  if len(recoPhotons) + len(recoPhotons2) == 1:
    if len(recoPhotons) == 1 and ntuple.showerComp[recoPhotons[0]] < 0.3:
      return False
    elif len(recoPhotons2) and ntuple.trackComp[recoPhotons2[0]] < 0.3:
      return False
    else:
      return True

def recoCutLongTracks(ntuple, fiducial):
  #Cut any event that has any reconstructed track longer than 20 cm
  acceptable = True
  for x in range(ntuple.nTracks):
    if ntuple.trackPID[x] != 2212:
      distxyz = np.sqrt((ntuple.trackStartPosX[x] - ntuple.trackEndPosX[x])**2 + (ntuple.trackStartPosY[x] - ntuple.trackEndPosY[x])**2 + (ntuple.trackStartPosZ[x] - ntuple.trackEndPosZ[x])**2)
      if distxyz > 25:
        acceptable = False
      elif ntuple.trackEndPosX[x] > (fiducial["xMax"] - 5) or ntuple.trackEndPosX[x] < (fiducial["xMin"] + 5) or ntuple.trackEndPosY[x] > (fiducial["yMax"] - 5) or ntuple.trackEndPosY[x] < (fiducial["yMin"] + 5) or ntuple.trackEndPosX[x] > (fiducial["zMax"] - 5) or ntuple.trackEndPosZ[x] < (fiducial["zMin"] + 5):
        acceptable = False
  if acceptable == False:
    return False
  else:
    return True 


def recoCutShortTracks(ntuple, threshold):
  #Cut any event that has any reconstructed track longer than 20 cm - NOT YET UPDATED TO ACCOMODATE THRESHOLD CHANGES
  acceptable = True
  for x in range(ntuple.nTracks):
    if threshold == 0:
      if ntuple.trackClassified[x] == 0:
        distxyz = np.sqrt((ntuple.trackStartPosX[x] - ntuple.trackEndPosX[x])**2 + (ntuple.trackStartPosY[x] - ntuple.trackEndPosY[x])**2 + (ntuple.trackStartPosZ[x] - ntuple.trackEndPosZ[x])**2)
        if distxyz < 10 and distxyz > 4:
          acceptable = False
    elif ntuple.trackSize[x] < threshold:
      distxyz = np.sqrt((ntuple.trackStartPosX[x] - ntuple.trackEndPosX[x])**2 + (ntuple.trackStartPosY[x] - ntuple.trackEndPosY[x])**2 + (ntuple.trackStartPosZ[x] - ntuple.trackEndPosZ[x])**2)
      if distxyz < 10 and distxyz > 4:
        acceptable = False
  if acceptable == False:
    return False
  else:
    return True

def recoCutManyTracks(ntuple, threshold = 0):
  #If an event has more than four tracks, it's almost certainly up to no good. Maybe.
  noTracks = 0
  noUnclassified = 0
  for track in range(ntuple.nTracks):
    noTracks += 1
    if threshold == 0:
      if ntuple.trackClassified[track] == 0:
        noClassified += 1
    elif ntuple.trackSize[x] < threshold:
      noUnclassified += 1
  if noTracks > 4:
    return False
  elif noTracks > 0:
    if noUnclassified/noTracks > 0.55:
      return False
    else:
      return True
  else:
    return True


def recoPhotonListFiducial(fiducial, ntuple, threshold = 0):
  #Creates a list of photons based on the showers in the event
  recoIDs = []
  for x in range(ntuple.nShowers):
    if threshold == 0:
      if ntuple.showerClassified[x] == 1:
        if ntuple.showerPID[x] == 22:
          if ntuple.showerStartPosX[x] > (fiducial["xMin"] + fiducial["width"]) and ntuple.showerStartPosX[x] < (fiducial["xMax"] - fiducial["width"]) and ntuple.showerStartPosY[x] > (fiducial["yMin"] + fiducial["width"]) and ntuple.showerStartPosY[x] < (fiducial["yMax"] - fiducial["width"]) and ntuple.showerStartPosZ[x] > (fiducial["zMin"] + fiducial["width"]) and ntuple.showerStartPosZ[x] < (fiducial["zMax"] - fiducial["width"]):
            recoIDs.append(x)
    elif ntuple.showerSize[x] > threshold:
      if ntuple.showerPID[x] == 22:
        #Extra check to ensure photons deposit in fiducial volume
        if ntuple.showerStartPosX[x] > (fiducial["xMin"] + fiducial["photonWidth"]) and ntuple.showerStartPosX[x] < (fiducial["xMax"] - fiducial["photonWidth"]) and ntuple.showerStartPosY[x] > (fiducial["yMin"] + fiducial["photonWidth"]) and ntuple.showerStartPosY[x] < (fiducial["yMax"] - fiducial["photonWidth"]) and ntuple.showerStartPosZ[x] > (fiducial["zMin"] + fiducial["photonWidth"]) and ntuple.showerStartPosZ[x] < (fiducial["zMax"] - fiducial["photonWidth"]):
          recoIDs.append(x)
  return recoIDs

def recoCutMaxInTime(ntuple, protonNo):
  #Cuts based on Taritree's new variable to cut down one events where the vertex was misplaced
  goodEvent = True
  if protonNo == 0 and ntuple.vtxMaxIntimePixelSum > 18000:
    goodEvent = False

  elif protonNo == 1 and ntuple.vtxMaxIntimePixelSum > 40000:
    goodEvent = False

  return goodEvent

def recoPhotonListTracks(fiducial, ntuple, threshold = 0):
  #Creates a list of photons based on the tracks in the event
  recoIDs = []
  for x in range(ntuple.nTracks):
    if threshold == 0:
      if ntuple.trackClassified[x] == 1:
        if ntuple.trackPID[x] == 22:
        #Extra check to ensure photons deposit in fiducial volume (with a 5 cm margin of error)
          if ntuple.trackStartPosX[x] > (fiducial["xMin"] + fiducial["width"]) and ntuple.trackStartPosX[x] < (fiducial["xMax"] - fiducial["width"]) and ntuple.trackStartPosY[x] > (fiducial["yMin"] + fiducial["width"]) and ntuple.trackStartPosY[x] < (fiducial["yMax"] - fiducial["width"]) and ntuple.trackStartPosZ[x] > (fiducial["zMin"] + fiducial["width"]) and ntuple.trackStartPosZ[x] < (fiducial["zMax"] - fiducial["width"]):
            recoIDs.append(x)
    elif ntuple.trackSize[x] > threshold:
      if ntuple.trackPID[x] == 22:
        #Extra check to ensure photons deposit in fiducial volume (with a 5 cm margin of error)
        if ntuple.trackStartPosX[x] > (fiducial["xMin"] + fiducial["width"]) and ntuple.trackStartPosX[x] < (fiducial["xMax"] - fiducial["width"]) and ntuple.trackStartPosY[x] > (fiducial["yMin"] + fiducial["width"]) and ntuple.trackStartPosY[x] < (fiducial["yMax"] - fiducial["width"]) and ntuple.trackStartPosZ[x] > (fiducial["zMin"] + fiducial["width"]) and ntuple.trackStartPosZ[x] < (fiducial["zMax"] - fiducial["width"]):
          recoIDs.append(x)
  return recoIDs

def recoCutMuonCompleteness(eventTree):
  goodEvent = True
  for x in range(eventTree.nTracks):
    if eventTree.trackPID[x] == 13:
      if eventTree.trackComp[x] < 0.5:
        goodEvent = False
  for x in range(eventTree.nShowers):
    if eventTree.showerPID[x] == 13:
      if eventTree.showerComp[x] < 0.5:
        goodEvent = False
  if goodEvent == False:
    return False
  else:
    return True

def recoCutPrimary(ntuple, recoPhotons, recoPhotons2):
  if len(recoPhotons) + len(recoPhotons2) == 1:
    if len(recoPhotons) == 1 and abs(ntuple.showerPrimaryScore[recoPhotons[0]]) < 1.4:
      return False
    elif len(recoPhotons2) ==1 and abs(ntuple.trackPrimaryScore[recoPhotons2[0]]) < 1.4:
      return False
    else: 
      return True
  else:
    return True
  
def recoCutTrackEnd(ntuple, recoPhotons, recoPhotons2):
  #Cuts single-photon events where the shower appears closer to the end of a muon track than to the vertex
  badEvent = False
  if len(recoPhotons) == 1 and len(recoPhotons2) == 0:
    x = recoPhotons[0]
    vertexDist = np.sqrt((ntuple.showerStartPosX[x] - ntuple.vtxX)**2 + (ntuple.showerStartPosY[x] - ntuple.vtxY)**2 +(ntuple.showerStartPosZ[x] - ntuple.vtxZ)**2)
    for y in range(ntuple.nTracks):
      if ntuple.trackPID[y] == 13:
        trackDist = np.sqrt((ntuple.showerStartPosX[x] - ntuple.trackEndPosX[y])**2 + (ntuple.showerStartPosY[x] - ntuple.trackEndPosY[y])**2 +(ntuple.showerStartPosZ[x] - ntuple.trackEndPosZ[y])**2)
        if vertexDist > trackDist:
          badEvent = True
          break
  elif len(recoPhotons) == 0 and len(recoPhotons2) == 1:
    x = recoPhotons2[0]
    vertexDist = np.sqrt((ntuple.trackStartPosX[x] - ntuple.vtxX)**2 + (ntuple.trackStartPosY[x] - ntuple.vtxY)**2 +(ntuple.trackStartPosZ[x] - ntuple.vtxZ)**2)
    for y in range(ntuple.nTracks):
      if ntuple.trackPID[y] == 13:
        trackDist = np.sqrt((ntuple.trackStartPosX[x] - ntuple.trackEndPosX[y])**2 + (ntuple.trackStartPosY[x] - ntuple.trackEndPosY[y])**2 +(ntuple.trackStartPosZ[x] - ntuple.trackEndPosZ[y])**2)
        if vertexDist > trackDist:
          badEvent = True
          break
  if badEvent == True:
    return False
  else:
    return True

def recoCutFarShowers(ntuple):
  tooLong = False
  for x in range(ntuple.nShowers):
    if ntuple.showerPID[x] != 22:
      distxyz = np.sqrt((ntuple.showerStartPosX[x] - ntuple.vtxX)**2 + (ntuple.showerStartPosY[x] - ntuple.vtxY)**2 + (ntuple.showerStartPosZ[x] - ntuple.vtxZ)**2)
      if distxyz < 10:
        tooLong = True
  if tooLong == True:
    return False
  else:
    return True

def recoCutOneProton(ntuple):
  protonCount = 0
  for x in range(ntuple.nTracks):
    if ntuple.trackPID[x] == 2212 and ntuple.trackRecoE[x] >= 100: 
      protonCount += 1
  for x in range(ntuple.nShowers):
    if ntuple.showerPID[x] == 2212 and ntuple.showerRecoE[x] >= 100:
      protonCount += 1
  return protonCount


#OTHER MISCELLANIOUS FUNCTIONS
  
# This function uses a particle's momentum vector and total energy to calculate its kinetic energy, then returns its kinetic energy in GeV
def kineticEnergyCalculator(ntuple, i):
  momentumVector = np.square(ntuple.truePrimPartPx[i]) + np.square(ntuple.truePrimPartPy[i]) + np.square(ntuple.truePrimPartPz[i])
  kineticGeV = ntuple.truePrimPartE[i] - np.sqrt((np.square(ntuple.truePrimPartE[i])) - momentumVector)
  return kineticGeV

# This function takes the event ntuple and trueSimPart index number of a particle 
# and returns the particle's track length and end distance from the detector wall
def particleDistancesCalculator(eventTree, i):
  endX, endY, endZ = eventTree.trueSimPartEndX[i], eventTree.trueSimPartEndY[i], eventTree.trueSimPartEndZ[i]
  if eventTree.trueSimPartEndX[i] >= 256:
    endX = 256
  elif eventTree.trueSimPartEndX[i] <= 0:
    endX = 0
  if eventTree.trueSimPartEndY[i] >= 116.5:
    endY = 116.5
  elif eventTree.trueSimPartEndY[i] <= -116.5:
    endY = -116.5
  if eventTree.trueSimPartEndZ[i] >= 1036:
    endZ = 1036
  elif eventTree.trueSimPartEndZ[i] <= 0:
    endZ = 0        
  deltaX = endX - eventTree.trueSimPartX[i]
  deltaY = endY - eventTree.trueSimPartY[i]
  deltaZ = endZ - eventTree.trueSimPartZ[i]
  particleDistance = np.sqrt(np.square(deltaX) + np.square(deltaY) + np.square(deltaZ))

  distToDetectorWallX, distToDetectorWallY, distToDetectorWallZ = endX, endY, endZ
  if endX >= 128:
    distToDetectorWallX = 256 - endX
  if endY >= 0:
    distToDetectorWallY = 116.5 - endY
  if endZ >= 518:
    distToDetectorWallZ = 1036 - endZ
  distanceToDetectorWall = min(distToDetectorWallX, abs(distToDetectorWallY), distToDetectorWallZ)

  return particleDistance, distanceToDetectorWall

# true cc cut
def trueCCCut(eventTree):
  if eventTree.trueNuCCNC == 0:
    return True
  else:
    return False
# reco cc cut
def recoCCCut(eventTree):
  recoPrimaryMuonTrackFound = False
  recoPrimaryMuonShowerFound = False
  recoPrimaryElectronTrackFound = False
  recoPrimaryElectronShowerFound = False
  cc = False
  unclassifiedTracks = 0
  for i in range(eventTree.nTracks):
    if eventTree.trackClassified[i] == 0:
      unclassifiedTracks += 1
    if eventTree.trackIsSecondary[i] == 0:
      if abs(eventTree.trackPID[i]) == 13:
        recoPrimaryMuonTrackFound = True
      elif abs(eventTree.trackPID[i]) == 11:
        recoPrimaryElectronTrackFound = True
  for i in range(eventTree.nShowers):
    if eventTree.showerIsSecondary[i] == 0:
      if abs(eventTree.showerPID[i]) == 11:
        recoPrimaryElectronShowerFound == True
      elif abs(eventTree.showerPID[i]) == 13:
        recoPrimaryMuonShowerFound == True
  if eventTree.nTracks != 0:
    if unclassifiedTracks/eventTree.nTracks >= 0.55:
      cc = True
  if recoPrimaryMuonTrackFound or recoPrimaryMuonShowerFound or recoPrimaryElectronTrackFound or recoPrimaryElectronShowerFound:   
    cc = True
  if eventTree.nTracks >= 4:
    cc = True

  if cc == True:
    return True
  else:
    return False 

# Takes ntuple and fiducial width, returns true if vtx is outside of fiducial volume
def trueFiducialCut(eventTree, fiducialWidth):
  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  outFiducial = False
  if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or \
    eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or \
      eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
    outFiducial = True
  return outFiducial

# Takes ntuple and fiducial width, returns true if vtx is outside of fiducial volume
def recoFiducialCut(eventTree, fiducialWidth):
  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  outFiducial = False
  if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
    eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
      eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
      outFiducial = True
  return outFiducial

# true pi+ cut: exclude any pi+ above 30MeV
def truePiPlusCut(eventTree):
  truePiPlusPresent = False
  for i in range(len(eventTree.truePrimPartPDG)):
    if abs(eventTree.truePrimPartPDG[i]) == 211:
      kE = kineticEnergyCalculator(eventTree, i)
      if kE >= 0.03:
        truePiPlusPresent = True
        break
  return truePiPlusPresent
  
# reco pi+ cut: exclude any suitably energetic charged pions  
def recoPiPlusCut(eventTree):  
  recoPiPlusPresent = False
  for i in range(eventTree.nTracks):
    if abs(eventTree.trackPID[i]) == 211:
      if eventTree.trackRecoE[i] >= 30:
        recoPiPlusPresent = True
        break
  return recoPiPlusPresent

# True proton selection: takes ntuple, returns count of all primary protons and float of TID (will match correctly if nprotons = 1)
def trueProtonSelection(eventTree):
  nTrueProtons = 0
  trueProtonTID = 0
  trueProtonIndex = 0
  for i in range(eventTree.nTrueSimParts):
    if eventTree.trueSimPartProcess[i] == 0 and eventTree.trueSimPartPDG[i] == 2212:
      momentumVector = np.square(eventTree.trueSimPartPx[i]) + np.square(eventTree.trueSimPartPy[i]) + np.square(eventTree.trueSimPartPz[i])
      kineticMeV = eventTree.trueSimPartE[i] - np.sqrt((np.square(eventTree.trueSimPartE[i])) - momentumVector)
      if kineticMeV >= 60:
        nTrueProtons += 1
        trueProtonTID = eventTree.trueSimPartTID[i]
        trueProtonIndex = i
  return nTrueProtons, trueProtonTID, trueProtonIndex

# reco proton selection: returns tid as well for matching later
def recoProtonSelection(eventTree):
  nRecoProtons = 0
  recoProtonTID = 0
  recoProtonIndex = 0
  for i in range(eventTree.nTracks):
    if eventTree.trackPID[i] == 2212 and eventTree.trackRecoE[i] >= 60:
      nRecoProtons += 1
      recoProtonTID = eventTree.trackTrueTID[i]
      recoProtonIndex = i
  return nRecoProtons, recoProtonTID, recoProtonIndex

# true photon process
# Find all secondary photons using Edep Sum
# then only count those that edep in detector
# then compute leading photon energy
def truePhotonSelection(eventTree, fiducialWidth):
  photonInSecondary = False
  photonIndexList = []
  truePhotonTIDList = []
  photonEDepOutsideFiducial = 0
  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  for i in range(eventTree.nTrueSimParts):
    if eventTree.trueSimPartPDG[i] == 22 and eventTree.trueSimPartProcess[i] == 1:
      if abs(eventTree.trueSimPartX[i] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[i] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[i] -eventTree.trueVtxZ) <= 0.15:
        pixelEnergy = eventTree.trueSimPartPixelSumYplane[i]*0.0126
        if pixelEnergy >= 20:
          if eventTree.trueSimPartEDepX[i] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[i] >= (xMax - fiducialWidth)\
            or eventTree.trueSimPartEDepY[i] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[i] >= (yMax - fiducialWidth)\
              or eventTree.trueSimPartEDepZ[i] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[i] >= (zMax - fiducialWidth):
            photonEDepOutsideFiducial += 1
          else:
            truePhotonTIDList.append(eventTree.trueSimPartTID[i])
            photonIndexList.append(i)

  photonEnergyList = [0]
  trueLeadingPhotonEnergy = 0.
  for i in range(len(photonIndexList)):
    photonEnergyMeV = eventTree.trueSimPartE[photonIndexList[i]] 
    photonEnergyList.append(photonEnergyMeV)
  trueLeadingPhotonEnergy = max(photonEnergyList)

  return truePhotonTIDList, trueLeadingPhotonEnergy, photonIndexList

def truePhotonSelectionPiZero(eventTree, fiducialWidth):
  photonInSecondary = False
  photonIndexList = []
  truePhotonTIDList = []
  photonEDepOutsideFiducial = 0
  for i in range(eventTree.nTrueSimParts):
    if eventTree.trueSimPartPDG[i] == 22 and eventTree.trueSimPartProcess[i] == 1:
      if eventTree.trueSimPartMID[i] not in eventTree.trueSimPartTID:
        if abs(eventTree.trueSimPartX[i] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[i] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[i] -eventTree.trueVtxZ) <= 0.15:
          pixelEnergy = eventTree.trueSimPartPixelSumYplane[i]*0.0126
          if pixelEnergy >= 20:
            photonIndexList.append(i)
            photonInSecondary = True

  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  if photonInSecondary == True:
    for i in range(len(photonIndexList)):
      if eventTree.trueSimPartEDepX[photonIndexList[i]] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[photonIndexList[i]] >= (xMax - fiducialWidth)\
        or eventTree.trueSimPartEDepY[photonIndexList[i]] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[photonIndexList[i]] >= (yMax - fiducialWidth)\
        or eventTree.trueSimPartEDepZ[photonIndexList[i]] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[photonIndexList[i]] >= (zMax - fiducialWidth):
          photonEDepOutsideFiducial += 1
      else:
        truePhotonTIDList.append(eventTree.trueSimPartTID[i])

  photonEnergyList = [0]
  trueLeadingPhotonEnergy = 0.
  for i in range(len(photonIndexList)):
    photonEnergyMeV = eventTree.trueSimPartE[photonIndexList[i]] 
    photonEnergyList.append(photonEnergyMeV)
  trueLeadingPhotonEnergy = max(photonEnergyList)

  return truePhotonTIDList, trueLeadingPhotonEnergy


# true photon process:
# finds all, checks edep within fiducial, calcs leading photon energy, returns both
# but no edep pixel sum
def truePhotonSelectionOldNtuple(eventTree, fiducialWidth):
  photonInSecondary = False
  primList = []
  photonIndexList = []
  truePhotonTIDList = []
  photonEDepOutsideFiducial = 0
  for i in range(len(eventTree.trueSimPartTID)):
    if eventTree.trueSimPartTID[i] == eventTree.trueSimPartMID[i]:
      primList.append(eventTree.trueSimPartTID[i])
  for i in range(len(eventTree.trueSimPartPDG)):
    if eventTree.trueSimPartPDG[i] == 22:
      if eventTree.trueSimPartMID[i] in primList:
        photonIndexList.append(i)
        photonInSecondary = True
      elif abs(eventTree.trueSimPartX[i] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[i] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[i] -eventTree.trueVtxZ) <= 0.15:
        photonIndexList.append(i)
        photonInSecondary = True

  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  if photonInSecondary == True:
    for i in range(len(photonIndexList)):
      if eventTree.trueSimPartEDepX[photonIndexList[i]] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[photonIndexList[i]] >= (xMax - fiducialWidth)\
        or eventTree.trueSimPartEDepY[photonIndexList[i]] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[photonIndexList[i]] >= (yMax - fiducialWidth)\
        or eventTree.trueSimPartEDepZ[photonIndexList[i]] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[photonIndexList[i]] >= (zMax - fiducialWidth):
          photonEDepOutsideFiducial += 1
      else:
        truePhotonTIDList.append(eventTree.trueSimPartTID[i])

  photonEnergyList = [0]
  trueLeadingPhotonEnergy = 0.
  for i in range(len(photonIndexList)):
    photonEnergyMeV = eventTree.trueSimPartE[photonIndexList[i]] 
    photonEnergyList.append(photonEnergyMeV)
  trueLeadingPhotonEnergy = max(photonEnergyList)

  return truePhotonTIDList, trueLeadingPhotonEnergy

# reco photon process:
# finds all, checks edep within fiducial, calcs leading photon energy, returns both 
def recoPhotonSelection(eventTree, fiducialWidth):
  reco = 0
  recoPhotonTIDList = []
  recoPhotonIndexList = []
  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  for i in range(eventTree.nShowers):
    if eventTree.showerPID[i] == 22:
      if eventTree.showerStartPosX[i] <= (xMin + fiducialWidth) or eventTree.showerStartPosX[i] >= (xMax - fiducialWidth) or \
        eventTree.showerStartPosY[i] <= (yMin + fiducialWidth) or eventTree.showerStartPosY[i] >= (yMax - fiducialWidth) or \
          eventTree.showerStartPosZ[i] <= (zMin + fiducialWidth) or eventTree.showerStartPosZ[i] >= (zMax - fiducialWidth):
            reco += 1
      else:
        recoPhotonTIDList.append(eventTree.showerTrueTID[i])
        recoPhotonIndexList.append(i)

  recoPhotonEnergyList = [0]
  recoLeadingPhotonEnergy = 0.
  
  for i in range(len(recoPhotonIndexList)):
    recoPhotonEnergyMeV = eventTree.showerRecoE[recoPhotonIndexList[i]] 
    recoPhotonEnergyList.append(recoPhotonEnergyMeV)
  recoLeadingPhotonEnergy = max(recoPhotonEnergyList)

  return recoPhotonTIDList, recoLeadingPhotonEnergy, recoPhotonIndexList

def recoPhotonSelectionInvMass(eventTree, fiducialWidth):
  reco = 0
  recoPhotonTIDList = []
  recoPhotonIndexList = []
  for i in range(eventTree.nShowers):
    if eventTree.showerPID[i] == 22:
      recoPhotonIndexList.append(i)
  
  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  for i in range(len(recoPhotonIndexList)):
    if eventTree.showerStartPosX[i] <= (xMin + fiducialWidth) or eventTree.showerStartPosX[i] >= (xMax - fiducialWidth) or \
    eventTree.showerStartPosY[i] <= (yMin + fiducialWidth) or eventTree.showerStartPosY[i] >= (yMax - fiducialWidth) or \
      eventTree.showerStartPosZ[i] <= (zMin + fiducialWidth) or eventTree.showerStartPosZ[i] >= (zMax - fiducialWidth):
          reco += 1
    else:
      recoPhotonTIDList.append(eventTree.showerTrueTID[i])

  recoInvMass = 0
  if len(recoPhotonTIDList) == 2:
    for i in range(eventTree.nShowers):
      if recoPhotonTIDList[0] == eventTree.showerTrueTID[i]:
        a = i
      if recoPhotonTIDList[1] == eventTree.showerTrueTID[i]:
        b = i
    print("e1: " + str(eventTree.showerRecoE[a]) + ", e2: " + str(eventTree.showerRecoE[b]) + ", pdg1: " + str(eventTree.showerTruePID[a]) + ", pdg2: " + str(eventTree.showerTruePID[b]))
    aDotB = eventTree.showerStartDirX[a]*eventTree.showerStartDirX[b] + eventTree.showerStartDirY[a]*eventTree.showerStartDirY[b] + eventTree.showerStartDirZ[a]*eventTree.showerStartDirZ[b]
    recoInvMass = np.sqrt((2*(eventTree.showerRecoE[a])*(eventTree.showerRecoE[b]))*(1-aDotB))

  return recoPhotonTIDList, recoInvMass

# true CC cut but only 100mev+ primaries
def trueCCCutLoose(eventTree):
  primaryEMu = False
  for i in range(eventTree.nTrueSimParts):
    if eventTree.trueSimPartProcess[i] == 0: 
      if abs(eventTree.trueSimPartPDG[i]) == 13 or abs(eventTree.trueSimPartPDG[i]) == 11:
        momentumVector = np.square(eventTree.trueSimPartPx[i]) + np.square(eventTree.trueSimPartPy[i]) + np.square(eventTree.trueSimPartPz[i])
        kineticMeV = eventTree.trueSimPartE[i] - np.sqrt((np.square(eventTree.trueSimPartE[i])) - momentumVector)
        if kineticMeV >= 100:
          primaryEMu = True
          break
  return primaryEMu

# reco CC cut but only 100MeV+ primaries
def recoCCCutLoose(eventTree):
  recoPrimaryMuonTrackFound = False
  recoPrimaryMuonShowerFound = False
  recoPrimaryElectronTrackFound = False
  recoPrimaryElectronShowerFound = False
  for i in range(eventTree.nTracks):
    if eventTree.trackIsSecondary[i] == 0:
      if abs(eventTree.trackPID[i]) == 13 and eventTree.trackRecoE[i] >= 100:
        recoPrimaryMuonTrackFound = True
      elif abs(eventTree.trackPID[i]) == 11 and eventTree.trackRecoE[i] >= 100:
        recoPrimaryElectronTrackFound = True
  for i in range(eventTree.nShowers):
    if eventTree.showerIsSecondary[i] == 0:
      if abs(eventTree.showerPID[i]) == 11 and eventTree.showerRecoE[i] >= 100:
        recoPrimaryElectronShowerFound == True
      elif abs(eventTree.showerPID[i]) == 13 and eventTree.showerRecoE[i] >= 100:
        recoPrimaryMuonShowerFound == True
  if recoPrimaryMuonTrackFound or recoPrimaryMuonShowerFound or recoPrimaryElectronTrackFound or recoPrimaryElectronShowerFound:   
    return True
  else:
    return False 

# function takes the ntuple, true proton index, and list of photon indices, and returns invariant mass calculations
def trueInvariantMassCalculations(eventTree, pi, truePhotonIndexList):
  # Create variables for photon indices
  gi1, gi2 = truePhotonIndexList[0], truePhotonIndexList[1]
  # Obtain proton, photon energies
  pE = eventTree.trueSimPartE[pi]
  gE1, gE2 = eventTree.trueSimPartE[gi1], eventTree.trueSimPartE[gi2]

  # Obtain momentum component variables 
  pPx, pPy, pPz = eventTree.trueSimPartPx[pi], eventTree.trueSimPartPy[pi], eventTree.trueSimPartPz[pi]
  gPx1, gPy1, gPz1 = eventTree.trueSimPartPx[gi1], eventTree.trueSimPartPy[gi1], eventTree.trueSimPartPz[gi1]
  gPx2, gPy2, gPz2 = eventTree.trueSimPartPx[gi2], eventTree.trueSimPartPy[gi2], eventTree.trueSimPartPz[gi2]

  # compute final quantities
  trueProtonInv = np.sqrt( np.square(pE) - ( np.square(pPx) + np.square(pPy) + np.square(pPz)))
  truePiZeroInv = np.sqrt( np.square(gE1 + gE2) - ( np.square(gPx1 + gPx2) + np.square(gPy1 + gPy2) + np.square(gPz1 + gPz2)))
  trueDeltaInv = np.sqrt( np.square(pE + gE1 + gE2) - ( np.square(pPx + gPx1 + gPx2) + np.square(pPy + gPy1 + gPy2) + np.square(pPz + gPz1 + gPz2)))

  return trueProtonInv, truePiZeroInv, trueDeltaInv

def recoInvariantMassCalculations(eventTree, recoProtonIndex, recoPhotonIndexList):
  pMass = 938.272089
  # obtain indices 
  pi = recoProtonIndex
  gi1, gi2 = recoPhotonIndexList[0], recoPhotonIndexList[1]
  # obtain energies
  pEtot = eventTree.trackRecoE[pi] + pMass
  gE1, gE2 = eventTree.showerRecoE[gi1], eventTree.showerRecoE[gi2]
  # proton: obtain start direction unit vector and momentum norm to create momentum vector
  nx, ny, nz = eventTree.trackStartDirX[pi], eventTree.trackStartDirY[pi], eventTree.trackStartDirZ[pi]
  pNorm = np.sqrt(np.square(pEtot) - np.square(pMass))
  pPx, pPy, pPz = pNorm * nx, pNorm * ny, pNorm * nz
  # obtain photon momentum vectors
  gPx1, gPy1, gPz1 = gE1 * eventTree.showerStartDirX[gi1], gE1 * eventTree.showerStartDirY[gi1], gE1 * eventTree.showerStartDirZ[gi1]
  gPx2, gPy2, gPz2 = gE2 * eventTree.showerStartDirX[gi2], gE2 * eventTree.showerStartDirY[gi2], gE2 * eventTree.showerStartDirZ[gi2]
  # final computations
  recoProtonInv = np.sqrt(np.square(pEtot) - (np.square(pPx) + np.square(pPy) + np.square(pPz)))
  recoPiZeroInv = np.sqrt(np.square(gE1 + gE2) - (np.square(gPx1 + gPx2)+(np.square(gPy1 + gPy2))+np.square(gPz1 + gPz2)))
  recoDeltaInv = np.sqrt( np.square(pEtot + gE1 + gE2) - ( np.square(pPx + gPx1 + gPx2) + np.square(pPy + gPy1 + gPy2) + np.square(pPz + gPz1 + gPz2)))

  return recoProtonInv, recoPiZeroInv, recoDeltaInv
