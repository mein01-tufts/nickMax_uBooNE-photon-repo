import sys, argparse
import numpy as np
import ROOT as rt

from cuts import trueCCCut, recoCCCut, trueFiducialCut, recoFiducialCut, truePiPlusCut, recoPiPlusCut, trueProtonSelection, recoProtonSelection, truePhotonSelection, recoPhotonSelection, trueCCCutLoose, recoCCCutLoose

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="1p2gNuTrackMu.root", help="output root file name")
args = parser.parse_args()

# Grab ntuple file, eventTree, and potTree for reference
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

# Scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles

# Calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
    potTree.GetEntry(i)
    ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

# Define hists to be filled:
trueNCTrackHist = rt.TH1F("nctracks", "True NC Events",60,0,50)
trueCCTrackHist = rt.TH1F("cctracks", "True CC Events",60,0,50)
# Define fiducial width
fiducialWidth = 10

vtxList = []

ccTrackMuList = []
ncTrackMuList = []
recoNCtrueNC = 0
recoNCtrueCC = 0
trueNCunclassified = 0
trueCCunclassified = 0
trueNCCclassified = 0
trueCCclassified = 0
trueNCtotalTracks = 0
trueCCtotalTracks = 0

zeroTrackNC = 0
zeroTrackCC = 0

## Event loop: Select reco signal, then of those events which are CC, 
## Plot each non-proton track in terms of its energy and length(stop length at det edge)
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
# Reco CC cut
    if recoCCCut(eventTree):
        continue
# Reco loose CC cut
#    if recoCCCutLoose(eventTree):
#        continue
# Reco fiducial cut
    if recoFiducialCut(eventTree, fiducialWidth):
        continue
# Reco cosmic cut
    if eventTree.vtxFracHitsOnCosmic >= 1.:
        continue
# Reco pi+ cut
    if recoPiPlusCut(eventTree):
        continue
# Reco proton selection
    nRecoProtons, recoProtonTID = recoProtonSelection(eventTree)
    if nRecoProtons != 1:
        continue
# Reco photon selection
    recoPhotonTIDList, recoLeadingPhotonEnergy = recoPhotonSelection(eventTree, fiducialWidth) 
    if len(recoPhotonTIDList) != 2:
        continue

# true CC: cut all non-cc events
    if trueCCCut(eventTree) == False:
        trackMuScore = []
        for i in range(eventTree.nTracks):
            if eventTree.trackClassified[i] == 1:
                trackMuScore.append(abs(eventTree.trackMuScore[i]))

        trackMuAvg = np.mean(trackMuScore)
        trueNCTrackHist.Fill(trackMuAvg, eventTree.xsecWeight)
        continue

    recoNCtrueCC += 1

# Calculate track lengths, store track energies for all 
    trackMuScore = []
    for i in range(eventTree.nTracks):
        ccTrackMuList.append(eventTree.trackMuScore[i])
        if eventTree.trackClassified[i] == 1:
            trackMuScore.append(abs(eventTree.trackMuScore[i]))
    trackMuAvg = np.mean(trackMuScore)
    trueCCTrackHist.Fill(trackMuAvg, eventTree.xsecWeight)


def histStackFill(title, histList, legendTitle, xTitle, yTitle, ntuplePOTSum):
#def histStackFill(title, histList, xTitle, yTitle, ntuplePOTSum):

  #Takes a list of histograms and converts them into one properly formatted stacked histogram. Returns the canvas on which the histogram is written
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.35, 0.6, 0.9, 0.9)
  colors = [rt.kGreen+2, rt.kRed, rt. kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kBlack, rt.kOrange]
  targetPOT = 3.2974607516739725e+20
  colors = [rt.kGreen+2, rt.kRed, rt. kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kBlack, rt.kYellow, rt.kGreen]
  targetPOT = 6.67e+20
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
    hist.SetFillColor(colors[x%7])
    hist.SetMarkerStyle(21)
    hist.SetMarkerColor(colors[x%7])
    histInt = hist.Integral(1, int(bins))
    integralSum += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round((histInt/sum*100), 2))+" percent of reconstructed signal events", "f")
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
#  return histCanvas, stack, histInt

trackCanvas, trackStack, trackLegend, trackInt = histStackFill("Neutrino Energy in Reconstructed Signal Events", [trueNCTrackHist, trueCCTrackHist], "Reconstructed Signal: (", "Average Muon Score of Classified Tracks", "events per 6.67e+20POT", ntuplePOTsum)
#trackCanvas, trackStack, trackInt = histStackFill("Neutrino Energy in Reconstructed Signal Events", [trueNCTrackHist, trueCCTrackHist], "True Neutrino Energy (GeV)", "events per 6.67e+20POT", ntuplePOTsum)

outFile = rt.TFile(args.outfile, "RECREATE")
trackCanvas.Write("Nu Energy Hist")
