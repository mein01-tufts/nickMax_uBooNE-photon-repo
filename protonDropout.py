import sys, argparse
import numpy as np
import ROOT as rt

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
args = parser.parse_args()

# needed for proper scaling of error bars:
#rt.TH1.SetDefaultSumw2(rt.kTRUE)

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles

#calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
kEabove60 = 0
kEbelow60 = 0
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT


xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

protonCounter = 0
kEabove60 = 0
kEbelow60 = 0
noProtonEvents = 0
ccCut = 0
fiducialCut = 0
cosmicCut = 0
piPlusCut = 0
noPhotonCut = 0
totalEntries = 0
protonEvents = 0

for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
    totalEntries += 1

    protonPresent = False
    for x in range(len(eventTree.truePrimPartPDG)):
        if eventTree.truePrimPartPDG[x] == 2212:
            protonPresent = True
            protonCounter += 1
            pVector = np.square(eventTree.truePrimPartPx[x]) + np.square(eventTree.truePrimPartPy[x]) + np.square(eventTree.truePrimPartPz[x])
            kE = eventTree.truePrimPartE[x] - np.sqrt((np.square(eventTree.truePrimPartE[x])) - pVector)
            if kE >= 0.06:
                kEabove60 += 1
            else:
                kEbelow60 += 1
    if protonPresent == True:
        protonEvents += 1
    else:
        noProtonEvents += 1
        continue

    if eventTree.trueNuCCNC != 1:
        ccCut += 1
        continue

    if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or \
        eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or \
        eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
        fiducialCut += 1
        continue

    if eventTree.vtxFracHitsOnCosmic >= 1.:
        cosmicCut += 1
        continue

    piPlusPresent = False    
    if abs(211) in eventTree.truePrimPartPDG: 
        piPlusPresent = True
        piPlusCut += 1
        continue
    
    photonInSecondary = False
    primList = []

    #check if Neutral Pion and Kaon in primaries - must be a photon if so
    if 111 in eventTree.truePrimPartPDG or 311 in eventTree.truePrimPartPDG:
        photonInSecondary = True
    else:
    #Create a list of primary particle Track IDs: 
    #Checks that a track id is equal to its mother id, this will return true for all primary particles
        for x in range(len(eventTree.trueSimPartTID)):
            if eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
                primList.append(eventTree.trueSimPartTID[x])
    #Iterate through trueSimParts to find photons
        for x in range(len(eventTree.trueSimPartPDG)):
            if eventTree.trueSimPartPDG[x] == 22:
    #Check for photon's parent particle in the primary list
                if eventTree.trueSimPartMID[x] in primList:
                    photonInSecondary = True
    #Check if the photon has starting coordinates within 1.5mm of the event vertex
                elif abs(eventTree.trueSimPartX[x] - eventTree.trueVtxX) <= 0.15 and abs(eventTree.trueSimPartY[x] - eventTree.trueVtxY) <= 0.15 and abs(eventTree.trueSimPartZ[x] -eventTree.trueVtxZ) <= 0.15:
                    photonInSecondary = True                
    #Discard event if no secondary photon is found
    if photonInSecondary == False:
        noPhotonCut += 1
        continue

print("there are " + str(totalEntries) + " total entries")
print("there are " + str(protonEvents) + " total proton events")
print("there are " + str(protonCounter) + " total protons found")
print("there are " + str(kEabove60) + " protons above 60MeV")
print("there are " + str(kEbelow60) + " protons below 60MeV")
print("there are " + str(noProtonEvents) + " events without a proton")
print("there are " + str(ccCut) + " proton events cut for CC")
print("there are " + str(fiducialCut) + " proton events outside fiducial volume")
print("there are " + str(cosmicCut) + " proton events tagged as cosmic")
print("there are " + str(piPlusCut) + " proton events cut for pi+ presence")
print("there are " + str(noPhotonCut) + " proton events cut for photon absence")
