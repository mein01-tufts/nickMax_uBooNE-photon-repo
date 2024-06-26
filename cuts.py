#NOTES:
#All of these functions are meant to operate within a larger function that loops through every event (or the chosen events) within the eventTree
import sys, argparse
import numpy as np
import ROOT as rt

from helpers.larflowreco_ana_funcs import getCosThetaGravVector


def trueCutNCC(eventTree):
#Filter for neutral current using truth - False if CC, True if NC
  if eventTree.trueNuCCNC != 1:
    return False
  else:
    return True

def trueCutPionPhoton(eventTree):
#Filter for pions and photons using truth - False if either (or both) are present, True otherwise
  pionPresent = False
  protonPresent = False                                    
  for x in range(len(eventTree.truePrimPartPDG)):
    if eventTree.truePrimPartPDG[x] == 211 and eventTree.truePrimPartE[x] >= 0.03:
      pionPresent = True
      break
    elif eventTree.truePrimPartPDG[x] == 2212 and eventTree.truePrimPartE[x] >= 0.06:
      protonPresent = True
      break
  if pionPresent == True or protonPresent == True:
    return False
  else:
    return True
  
def trueCutFiducials(eventTree, fiducialWidth):
#Filter by determining if the event vertex falls within the fiducial width using truth  - True if it's not within the radius, false if it is
  if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
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
    #Check if there is actually a detectable photon
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







  
def trueSignalFinder(eventTree):
#Uses a series of truth-based checks to determine signal as thoroughly as possible
  xMin, xMax = 0, 256
  yMin, yMax = -116.5, 116.5
  zMin, zMax = 0, 1036
  fiducialWidth = 10

  badEvent = False

  if eventTree.trueNuCCNC != 1:
    badEvent = True

#fiducial cut removes everything within fiducial

  if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or \
    eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or \
    eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
    badEvent = True
        
  #skip events where all hits overlap with tagged cosmic rays
  if eventTree.vtxFracHitsOnCosmic >= 1.:
    badEvent = True

  photonInSecondary = False
  primList = []
  photonIDList = []
  protonPresent = False
  PionPresent = False

  #Determining presence of suitably energetic protons and charged pions - if they're present, the event is beyond the scope of our investigation
  for x in range(len(eventTree.truePrimPartPDG)):
    if eventTree.truePrimPartPDG[x] == 211:
      pionEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
      if pionEnergy >= 0.03:
        pionPresent = True
        badEvent = True
    elif eventTree.truePrimPartPDG[x] == 2212:
      protonEnergy = eventTree.truePrimPartE[x] - np.sqrt(eventTree.truePrimPartE[x]**2 - (eventTree.truePrimPartPx[x]**2+eventTree.truePrimPartPy[x]**2+eventTree.truePrimPartPz[x]**2))
      if protonEnergy >= 0.06:
        protonPresent = True
        badEvent = True
  
  #Check if Neutral Pion and Kaon in primaries
#  if 111 in eventTree.truePrimPartPDG or 311 in eventTree.truePrimPartPDG:
    #Check if there is actually a detectable photon
#    for x in range(eventTree.nTrueSimParts):
#      if eventTree.TrueSimPartPDG[x] == 22:
        photonInSecondary = True
#  else:
  #Create a list of prime particle Track IDs
  for x in range(len(eventTree.trueSimPartTID)):
    if eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
      primList.append(eventTree.trueSimPartTID[x])
    #Iterate through to find photons
  for x in range(len(eventTree.trueSimPartPDG)):
    if eventTree.trueSimPartPDG[x] == 22:
      photonIDList.append(x)
        #Check for parent particle in the primary list
      if eventTree.trueSimPartMID[x] in primList:
        photonInSecondary = True
      #Failing that, check if the photon has coordinates within 15 mm of those of the event vertex
      elif abs(eventTree.trueSimPartX[x] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[x] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[x] -eventTree.trueVtxZ) <= 0.15:
        photonInSecondary = True
  #Discard event unless a secondary photon is found
  if photonInSecondary == False:
    badEvent = True  

    #Discard event unless the photon begins to deposit energy within the fiducial volume
  for x in range(len(photonIDList)):
    if eventTree.trueSimPartEDepX[x] <= (xMin + fiducialWidth) or eventTree.trueSimPartEDepX[x] >= (xMax - fiducialWidth) or eventTree.trueSimPartEDepY[x] <= (yMin + fiducialWidth) or eventTree.trueSimPartEDepY[x] >= (yMax - fiducialWidth) or eventTree.trueSimPartEDepZ[x] <= (zMin + fiducialWidth) or eventTree.trueSimPartEDepZ[x] >= (zMax - fiducialWidth):
      badEvent = True


  if badEvent == True:
    return False
  
  else:
    return True
