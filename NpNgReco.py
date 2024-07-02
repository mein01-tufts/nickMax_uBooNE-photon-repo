import sys, argparse
import numpy as np
import ROOT as rt

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="NpNgRecoLogged.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
args = parser.parse_args()

#needed for proper scaling of error bars:
#rt.TH1.SetDefaultSumw2(rt.kTRUE)

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
#we will write histograms to output file for:
onePhotonHist = rt.TH1F("1Photon_nProtonHist", "Energy of NC events with N proton(s) and 1 photon",24,0,6)
twoPhotonHist = rt.TH1F("2Photon_nProtonHist", "Energy of NC events with N proton(s) and 2 photons",24,0,6)
threePhotonHist = rt.TH1F("3Photon_nProtonHist", "Energy of NC events with N proton(s) and 3+ photons",24,0,6)

#set histogram axis titles and increase line width
def configureHist(h):
    h.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
    h.GetXaxis().SetTitle("energy (GeV)")
    h.SetLineWidth(2)
    return h

#scale the histograms based on total good POT
onePhotonHist = configureHist(onePhotonHist)
twoPhotonHist = configureHist(twoPhotonHist)
threePhotonHist = configureHist(threePhotonHist)

#set detector min/max and fiducial width (cm)
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

recoFiducialCut = 0
recoCosmicCut = 0
recoCCcut = 0
recoCCMuon = 0
recoCCElectron = 0
NCcounter = 0
vtxID = 0
foundVTX = 0
#begin loop over events in ntuple file
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

    if eventTree.foundVertex != 1:
        foundVTX += 1
        continue

#reco fiducial cut
    if eventTree.vtxX <= (xMin + fiducialWidth) or eventTree.vtxX >= (xMax - fiducialWidth) or \
        eventTree.vtxY <= (yMin + fiducialWidth) or eventTree.vtxY >= (yMax - fiducialWidth) or \
        eventTree.vtxZ <= (zMin + fiducialWidth) or eventTree.vtxZ >= (zMax - fiducialWidth):
        #if abs(eventTree.vtxDistToTrue) <= 3:
        recoFiducialCut += 1
        continue
       
#reco cosmic cut - same as true
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        recoCosmicCut += 1
        continue

#reco charged-current cut:        
#iterate through all tracks in event,
#look for non-secondary tracks identified as muons
    primaryMuonFound = False
    primaryElectronTrackFound = False
    for i in range(eventTree.nTracks):
        if eventTree.trackIsSecondary[i] == 0:
            if abs(eventTree.trackPID[i]) == 13:
                primaryMuonFound = True
            if abs(eventTree.trackPID[i]) == 11:
                primaryElectronTrackFound = True
#cut events w/ primary muons
    if primaryMuonFound or primaryElectronTrackFound:    
        recoCCMuon += 1
        continue
    
#reco charged-current cut phase 2:
#look for non-secondary showers identified as electrions
    primaryElectronFound = False
    for i in range(eventTree.nShowers):
        if eventTree.showerIsSecondary[i] == 0:
            if abs(eventTree.showerPID[i]) == 11:
                primaryElectronFound == True
#cut events w/ primary electrons
    if primaryElectronFound:
        recoCCElectron += 1
        continue

#reco pi+ finder:
    piPlusPresent = False
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 211:
            piPlusPresent = True
            break
    if piPlusPresent:
        continue

#reco protons: find and tally
    nProtons = 0
    for i in range(eventTree.nTracks):
        if abs(eventTree.trackPID[i]) == 2212: 
            nProtons += 1
    if nProtons == 0:
        continue

#reco photons: find and tally
    nPhotonsPresent = 0
    for i in range(eventTree.nShowers):
        if eventTree.showerPID[i] == 22:
            nPhotonsPresent += 1
    if nPhotonsPresent == 0:
        continue

#fill histograms based on number of photons
    if nPhotonsPresent == 1:
        onePhotonHist.Fill(nProtons, eventTree.xsecWeight)
    elif nPhotonsPresent == 2:
        twoPhotonHist.Fill(nProtons, eventTree.xsecWeight)
    elif nPhotonsPresent >= 3:
        threePhotonHist.Fill(nProtons, eventTree.xsecWeight)

#----- end of event loop ---------------------------------------------#


onePhotonHist.Scale(targetPOT/ntuplePOTsum)
twoPhotonHist.Scale(targetPOT/ntuplePOTsum)
threePhotonHist.Scale(targetPOT/ntuplePOTsum)

#create empty stacked histogram, add to it below
histStack = rt.THStack("stackedHist", "NC Events, N Gamma + N Proton")

#this function iteratively scales and fills each histogram to the stack
def histScalerStacker(hist, kColor):
    hist.SetFillColor(kColor)
    hist.SetMarkerStyle(21)
    hist.SetMarkerColor(kColor)
    histStack.Add(hist)

#add each histogram to stack, choose their color
histScalerStacker(onePhotonHist, rt.kRed)
histScalerStacker(twoPhotonHist, rt.kCyan)
histScalerStacker(threePhotonHist, rt.kGreen)

#integrate all events in the stack
oneInt = round(onePhotonHist.Integral(0,60), 2)
twoInt = round(twoPhotonHist.Integral(0,60), 2)
threeInt = round(threePhotonHist.Integral(0,60), 2)

legend = rt.TLegend(0.5, 0.65, 0.9, 0.9)  # (x1, y1, x2, y2) in NDC coordinates

#add color key to legend
legend.AddEntry(onePhotonHist, "#splitline{1 secondary photon,}" + "{" + str(oneInt) + " events per 6.67e+20 POT}", "f")
legend.AddEntry(twoPhotonHist, "#splitline{2 secondary photons,}" + "{" + str(twoInt) + " events per 6.67e+20 POT}", "f")
legend.AddEntry(threePhotonHist, "#splitline{3+ secondary photons,}" + "{" + str(threeInt) + " events per 6.67e+20 POT}", "f")

histCanvas = rt.TCanvas()
histStack.Draw("HIST")
histStack.GetXaxis().SetTitle("Number of Protons")
histStack.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
legend.Draw()
rt.gPad.SetLogy(1)
rt.gPad.Update()

#create output root file and write histograms to file
outFile = rt.TFile(args.outfile, "RECREATE")
histCanvas.Write()