import sys, argparse
import numpy as np
import ROOT as rt

from cuts import histStackFill

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gRecoOutputTrue.root", help="output root file name")
args = parser.parse_args()

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

#scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles

#calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

#define histograms to fill
trueSignalHist = rt.TH1F("trueSignalHist", "True NC 1 proton 2 gamma Events",60,0,3)

noVtxFoundHist = rt.TH1F("No Vertex Found", "Reco couldn't find a vertex",60,0,3)
outFiducialHist = rt.TH1F("Futside Fiducial", "Reco placed vertex outside fiducial volume",60,0,3)
chargedCurrentHist = rt.TH1F("Charged Current", "Reco identified as charged-current",60,0,3)
piPlusHist = rt.TH1F("pi+", "Reco found a pi+",60,0,3)
noProtonHist = rt.TH1F("No Protons", "Reco found no protons",60,0,3)
pluralProtonHist = rt.TH1F("Multiple Protons", "Reco found multiple protons",60,0,3)
noPhotonHist = rt.TH1F("No Photons", "Reco found no photons",60,0,3)
onePhotonHist = rt.TH1F("One Photon", "Reco found 1 photon",60,0,3)
manyPhotonHist = rt.TH1F("Many Photons", "Reco found 3+ photons",60,0,3)
recoSignalHist = rt.TH1F("Reco = True", "Correctly identified by Reco",60,0,3)

#set histogram axis titles and increase line width
def configureHist(h):
    h.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
    h.GetXaxis().SetTitle("photon energy (GeV)")
    h.SetLineWidth(2)
    return h

#scale the histograms based on total good POT
trueSignalHist = configureHist(trueSignalHist)

noVtxFoundHist = configureHist(noVtxFoundHist)
outFiducialHist = configureHist(outFiducialHist)
chargedCurrentHist = configureHist(chargedCurrentHist)
piPlusHist = configureHist(piPlusHist)
noProtonHist = configureHist(noProtonHist)
pluralProtonHist = configureHist(pluralProtonHist)
noPhotonHist = configureHist(noPhotonHist)
onePhotonHist = configureHist(onePhotonHist)
manyPhotonHist = configureHist(manyPhotonHist)
recoSignalHist = configureHist(recoSignalHist)

#set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

totalEvents = 0
recoNoPhotons = 0
reco1Photon = 0
reco2Photon = 0
reco3Photon = 0

#begin loop over events in ntuple file:
# if we start with truth then reco-match, we can ID reco's false negatives
# if we start with reco then truth-match, we can ID reco's false positives

#start by using truth to cut all background, then do the same checks in reco
#after reco-matching, sort events by how reco identifies the true signal, 
#then create a stacked histogram w/ x-axis of photon energy
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

#true charged-current cut
    if eventTree.trueNuCCNC != 1:
        continue

#true fiducial volume cut
    if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or \
        eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or \
        eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
        continue

#true cosmic background cut
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        continue

#true proton check: does its kinetic energy exceed 60mev? - count it if so
    nProtons = 0
    for i in range(len(eventTree.truePrimPartPDG)):
        if eventTree.truePrimPartPDG[i] == 2212:
            pVector = np.square(eventTree.truePrimPartPx[i]) + np.square(eventTree.truePrimPartPy[i]) + np.square(eventTree.truePrimPartPz[i])
            kE = eventTree.truePrimPartE[i] - np.sqrt((np.square(eventTree.truePrimPartE[i])) - pVector)
            if kE >= 0.06:
                nProtons += 1
    if nProtons != 1:
        continue

#true pi+ check: - we want to exclude any above 30MeV
    piPlusPresent = False
    for i in range(len(eventTree.truePrimPartPDG)):
        if abs(eventTree.truePrimPartPDG[i] == 211):
            pVector = np.square(eventTree.truePrimPartPx[i]) + np.square(eventTree.truePrimPartPy[i]) + np.square(eventTree.truePrimPartPz[i])
            kE = eventTree.truePrimPartE[i] - np.sqrt((np.square(eventTree.truePrimPartE[i])) - pVector)
            if kE >= 0.03:
                piPlusPresent = True
                break
    if piPlusPresent == True:
        continue

#true photon check: look through primary particles to count the photon daughters
    photonInSecondary = False
    primList = []
    photonIndexList = []
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
    if photonInSecondary == False:
        continue
    if len(photonIndexList) != 2:
        continue

#find leading photon energy
    photonEnergyList = []
    for i in range(len(photonIndexList)):
        photonEnergyGeV = eventTree.trueSimPartE[photonIndexList[i]] / 1000
        photonEnergyList.append(photonEnergyGeV)
    leadingPhotonEnergy = max(photonEnergyList)

    trueSignalHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)

    totalEvents += 1

    #-------- end of truth selection --------#
    #-------- start of reco matching --------#

#reco vertex check: does reco even find an event
    if eventTree.foundVertex != 1:
        noVtxFoundHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue

#reco fiducial histogram fill
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
        eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
        eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
        outFiducialHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue

# don't need reco cosmic cut - it's the same as the true

#reco charged-current cut:        
#iterate through all tracks/showers in event,
#look for non-secondary tracks identified as muons or electrons
#and for non-secondary showers identified as electrons
    recoPrimaryMuonFound = False
    recoPrimaryElectronTrackFound = False
    recoPrimaryElectronFound = False
    for i in range(eventTree.nTracks):
        if eventTree.trackIsSecondary[i] == 0:
            if abs(eventTree.trackPID[i]) == 13:
                recoPrimaryMuonFound = True
            if abs(eventTree.trackPID[i]) == 11:
                recoPrimaryElectronTrackFound = True
    for i in range(eventTree.nShowers):
        if eventTree.showerIsSecondary[i] == 0:
            if abs(eventTree.showerPID[i]) == 11:
                recoPrimaryElectronFound == True
#fill hist w/ events that have primary electrons/muons
    if recoPrimaryMuonFound or recoPrimaryElectronTrackFound or recoPrimaryElectronFound:    
        chargedCurrentHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue

#reco protons: find all protons of KE >= 60MeV
    recoNumProtons = 0
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 2212:
            if eventTree.trackRecoE[i] >= 60:
                recoNumProtons += 1
    if recoNumProtons == 0:
        noProtonHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if recoNumProtons >= 2:
        pluralProtonHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue
    
#reco pi+: fill all events w/ pi+ of KE >= 30MeV
    recoPiPlusPresent = False
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 211:
                if eventTree.trackRecoE[i] >= 30:
                    recoPiPlusPresent = True
                    break
    if recoPiPlusPresent:
        piPlusHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue

#reco photons: find and tally all photons
    recoNumPhotons = 0
    for i in range(eventTree.nShowers):
        if eventTree.showerPID[i] == 22:
            recoNumPhotons += 1
    if recoNumPhotons == 0:
        recoNoPhotons += 1
        noPhotonHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if recoNumPhotons == 1:
        reco1Photon += 1
        onePhotonHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if recoNumPhotons >= 3:
        reco3Photon += 1
        manyPhotonHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
        continue
    if recoNumPhotons == 2:
        reco2Photon += 1
        recoSignalHist.Fill(leadingPhotonEnergy, eventTree.xsecWeight)
#----- end of event loop ---------------------------------------------#

histList = [recoSignalHist, noVtxFoundHist, outFiducialHist, chargedCurrentHist, piPlusHist, noProtonHist, pluralProtonHist, noPhotonHist, onePhotonHist, manyPhotonHist]

histCanvas, stack, legend, histInt = histStackFill("Reco IDs of True NC 1 Proton, 2 Gamma Events", histList, "Reco Identification (events per 6.67e+20 POT)")

outFile = rt.TFile(args.outfile, "RECREATE")
histCanvas.Write()


def histStackFill(title, histList, legendTitle):
  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt. kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kBlack, rt.kYellow, rt.kGreen]
  targetPOT = 6.67e+20
  ntuplePOTSum = 4.675690535431973e+20
  for x in range(len(histList)):
    hist = histList[x]
    bins = hist.GetNbinsX()
    hist.Scale(targetPOT/ntuplePOTSum)
    hist.SetFillColor(colors[x%7])
    hist.SetMarkerStyle(21)
    hist.SetMarkerColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1))+" events", "f")
    stack.Add(hist)
  #Make the canvas and draw everything to it (NOTE - This component is only designed for events using 6.67e+20 scaling
  histCanvas = rt.TCanvas() 
  stack.Draw("HIST")
  stack.GetXaxis().SetTitle("Leading Photon Energy (GeV)")
  stack.GetYaxis().SetTitle("Events per 6.67e+20 POT")
  legend.SetHeader(str(legendTitle), "C")
  legend.Draw()
  histCanvas.Update()

  return histCanvas, stack, legend, histInt