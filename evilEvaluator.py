import sys, argparse
import numpy as np
import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

from cuts import trueCutNC, trueCutFiducials,trueCutCosmic, truePhotonList, trueCutPionProton, histStack, recoNoVertex, recoFiducials, recoPhotonList, recoPionProton, recoNeutralCurrent, scaleRecoEnergy, scaleTrueEnergy, recoPion, recoProton, recoCutElectronScore, recoCutShowerFromChargeScore, recoCutLongTracks, recoPhotonListFiducial, recoCutPrimary, recoCutShortTracks, recoPhotonListTracks, recoCutFarShowers, recoCutMuons, trueCutMuons, recoCutElectrons, trueCutElectrons, histStackTwoSignal

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
#PURITY HISTOGRAMS
purityTotal1 = rt.TH1F("PTotal1", "One Photon",60,0,2)
puritySignal1 =  rt.TH1F("PSignal1", "Signal",60,0,2)
purityOffSignal1 = rt.TH1F("POffSignal", "Signal (1 Proton, Undetected)",60,0,2)
#purityCC1 = rt.TH1F("PCC1", "Actually Charged Current",60,0,2)
purityMuon1 = rt.TH1F("PMuon1", "Over-threshold Muon",60,0,2)
purityElectron1 = rt.TH1F("PElectron1", "Over-threshold Electron",60,0,2)
purityFiducials1 = rt.TH1F("PFiducial1", "Out of Fiducial",60,0,2)
purityPion1 = rt.TH1F("PPionPion1", "Charged Pion",60,0,2)
purityProton1 = rt.TH1F("PPionProton1", "2+ Protons",60,0,2)
purityNoPhotons1 = rt.TH1F("PNoPhoton1", "No Real Photons",60,0,2)
purityTwoPhotons1 = rt.TH1F("PTwoPhoton1", "2 Real Photons",60,0,2)
purityManyPhotons1 = rt.TH1F("PMorePhoton1", "3+ Real Photons",60,0,2)

purityTotalProton1 = rt.TH1F("PTotal1Proton", "One Photon",60,0,2)
puritySignalProton1 =  rt.TH1F("PSignal1Proton", "Signal",60,0,2)
purityOffSignalProton1 = rt.TH1F("PoffSignalProton1", "Signal (0 Protons, 1 Detected)",60,0,2)
#purityCCProton1 = rt.TH1F("PCC1Proton", "Actually Charged Current (With Proton)",60,0,2)
purityMuonProton1 = rt.TH1F("PMuonProton1", "Over-threshold Muon",60,0,2)
purityElectronProton1 = rt.TH1F("PElectronProton1", "Over-threshold Electron",60,0,2)
purityFiducialsProton1 = rt.TH1F("PFiducialProton1", "Out of Fiducial",60,0,2)
purityPionProton1 = rt.TH1F("PPionPionProton1", "Charged Pion",60,0,2)
purityProtonProton1 = rt.TH1F("PPionProtonProton1", "2+ Protons",60,0,2)
purityNoPhotonsProton1 = rt.TH1F("PNoPhotonProton1", "No Real Photons",60,0,2)
purityTwoPhotonsProton1 = rt.TH1F("PTwoPhotonProton1", "2 Real Photons",60,0,2)
purityManyPhotonsProton1 = rt.TH1F("PMorePhotonProton1", "3+ Real Photons",60,0,2)

purityTotal2 = rt.TH1F("PTotal2", "One Photon",60,0,2)
puritySignal2 =	 rt.TH1F("PSignal2", "Signal",60,0,2)
#purityCC2 = rt.TH1F("PCC2", "Actually Charged Current",60,0,2)
purityMuon2 = rt.TH1F("PMuon2", "Over-threshold Muon",60,0,2)
purityElectron2 = rt.TH1F("PElectron2", "Over-threshold Electron",60,0,2)
purityFiducials2 = rt.TH1F("PFiducial2", "Out of Fiducial",60,0,2)
purityPion2 = rt.TH1F("PPionPion2", "Charged Pion",60,0,2)
purityProton2 = rt.TH1F("PPionProton2", "2+ Protons",60,0,2)
purityNoPhotons2 = rt.TH1F("PNoPhoton2", "No Real Photons",60,0,2)
purityOnePhoton2 = rt.TH1F("POnePhoton2", "1 Real Photon",60,0,2)
purityManyPhotons2 = rt.TH1F("PManyPhoton2", "3+ Real Photons",60,0,2)

purityTotalProton2 = rt.TH1F("PTotalProton2", "Two Photons",60,0,2)
puritySignalProton2 =	 rt.TH1F("PSignalProton2", "Signal",60,0,2)
#purityCCProton2 = rt.TH1F("PCCProton2", "Actually Charged Current",60,0,2)
purityMuonProton2 = rt.TH1F("PMuonProton2", "Over-threshold Muon",60,0,2)
purityElectronProton2 = rt.TH1F("PElectronProton2", "Over-threshold Electron",60,0,2)
purityFiducialsProton2 = rt.TH1F("PFiducialProton2", "Out of Fiducial",60,0,2)
purityPionProton2 = rt.TH1F("PPionPionProton2", "Charged Pion",60,0,2)
purityProtonProton2 = rt.TH1F("PPionProtonProton2", "2+ Protons",60,0,2)
purityNoPhotonsProton2 = rt.TH1F("PNoPhotonProton2", "No Real Photons",60,0,2)
purityOnePhotonProton2 = rt.TH1F("POnePhotonProton2", "1 Real Photon",60,0,2)
purityManyPhotonsProton2 = rt.TH1F("PManyPhotonProton2", "3+ Real Photons",60,0,2)

purityTotal3 = rt.TH1F("PTotal3", "One Photon",60,0,2)
puritySignal3 =	 rt.TH1F("PSignal3", "Signal",60,0,2)
#purityCC3 = rt.TH1F("PCC3", "Actually Charged Current",60,0,2)
purityMuon3 = rt.TH1F("PMuon3", "Over-threshold Muon",60,0,2)
purityElectron3 = rt.TH1F("PElectron3", "Over-threshold Electron",60,0,2)
purityFiducials3 = rt.TH1F("PFiducial3", "Out of Fiducial",60,0,2)
purityPion3 = rt.TH1F("PPion3", "Charged Pion",60,0,2)
purityProton3 = rt.TH1F("PProton3", "2+ Protons",60,0,2)
purityNoPhotons3 = rt.TH1F("PNoPhoton3", "No Real Photons",60,0,2)
purityOnePhoton3 = rt.TH1F("POnePhoton3", "2 Real Photons",60,0,2)
purityTwoPhotons3 = rt.TH1F("PTwoPhotons3", "3+ Real Photons",60,0,2)

purityTotalProton3 = rt.TH1F("PTotalProton3", "Three Photons",60,0,2)
puritySignalProton3 =	 rt.TH1F("PSignalProton3", "Signal",60,0,2)
#purityCCProton3 = rt.TH1F("PCCProton3", "Actually Charged Current (With Proton)",60,0,2)
purityMuonProton3 = rt.TH1F("PMuonProton3", "Over-threshold Muon",60,0,2)
purityElectronProton3 = rt.TH1F("PElectronProton3", "Over-threshold Electron",60,0,2)
purityFiducialsProton3 = rt.TH1F("PFiducialProton3", "Out of Fiducial",60,0,2)
purityPionProton3 = rt.TH1F("PPionProton3", "Charged Pion",60,0,2)
purityProtonProton3 = rt.TH1F("PPionProtonProton3", "2+ Protons",60,0,2)
purityNoPhotonsProton3 = rt.TH1F("PNoPhotonProton3", "No Real Photons",60,0,2)
purityOnePhotonProton3 = rt.TH1F("POnePhotonProton3", "2 Real Photons",60,0,2)
purityTwoPhotonsProton3 = rt.TH1F("PTwoPhotonsProton3", "3+ Real Photons",60,0,2)

#PURITY HISTLISTS
#Non-Proton Purity Histlists
totalPHists = [purityTotal1, purityTotal2, purityTotal3]
signalPHists = [puritySignal1, puritySignal2, puritySignal3]
#CCPHists = [purityCC1, purityCC2, purityCC3]
muonPHists = [purityMuon1, purityMuon2, purityMuon3]
electronPHists = [purityElectron1, purityElectron2, purityElectron3]
fiducialPHists = [purityFiducials1, purityFiducials2, purityFiducials3]
pionPHists = [purityPion1, purityPion2, purityPion3]
protonPHists = [purityProton1, purityProton2, purityProton3]
noPhotonPHists = [purityNoPhotons1, purityNoPhotons2, purityNoPhotons3]
onePhotonPHists = [puritySignal1, purityOnePhoton2, purityOnePhoton3]
twoPhotonPHists = [purityTwoPhotons1, puritySignal2, purityTwoPhotons3]
manyPhotonPHists = [purityManyPhotons1, purityManyPhotons2, puritySignal3]

#Proton Purity Histlists
totalProtonPHists = [purityTotalProton1, purityTotalProton2, purityTotalProton3]
signalProtonPHists = [puritySignalProton1, puritySignalProton2, puritySignalProton3]
#CCProtonPHists = [purityCCProton1, purityCCProton2, purityCCProton3]
muonProtonPHists = [purityMuonProton1, purityMuonProton2, purityMuonProton3]
electronProtonPHists = [purityElectronProton1, purityElectronProton2, purityElectronProton3]
fiducialProtonPHists = [purityFiducialsProton1, purityFiducialsProton2, purityFiducialsProton3]
pionProtonPHists = [purityPionProton1, purityPionProton2, purityPionProton3]
protonProtonPHists = [purityProtonProton1, purityProtonProton2, purityProtonProton3]
noPhotonProtonPHists = [purityNoPhotonsProton1, purityNoPhotonsProton2, purityNoPhotonsProton3]
onePhotonProtonPHists = [puritySignalProton1, purityOnePhotonProton2, purityOnePhotonProton3]
twoPhotonProtonPHists = [purityTwoPhotonsProton1, puritySignalProton2, purityTwoPhotonsProton3]
manyPhotonProtonPHists = [purityManyPhotonsProton1, purityManyPhotonsProton2, puritySignalProton3]

#Cosmics go here, so we can put them on purity
cosmicOnePhoton = rt.TH1F("cBackground1", "One Photon Cosmic",60,0,2)
cosmicTwoPhotons = rt.TH1F("cBackground2", "Two Photon Cosmic",60,0,2)
cosmicThreePhotons = rt.TH1F("cBackground3", "3+ Photon Cosmic",60,0,2)
cosmicList = [cosmicOnePhoton, cosmicTwoPhotons, cosmicThreePhotons]

cosmicProtonOnePhoton = rt.TH1F("cBackgroundProton1", "One Photon Cosmic",60,0,2)
cosmicProtonTwoPhotons = rt.TH1F("cBackgroundProton2", "Two Photon Cosmic",60,0,2)
cosmicProtonThreePhotons = rt.TH1F("cBackgroundProton3", "3+ Photon Cosmic",60,0,2)
cosmicProtonList = [cosmicProtonOnePhoton, cosmicProtonTwoPhotons, cosmicProtonThreePhotons]

#Big Lists, for Big Plots
pList1 = [puritySignal1, purityOffSignal1, purityMuon1, purityElectron1, purityFiducials1, purityPion1, purityProton1, purityNoPhotons1, purityTwoPhotons1, purityManyPhotons1, cosmicOnePhoton]
pList2 = [puritySignal2, purityMuon2, purityElectron2, purityFiducials2, purityPion2, purityProton2, purityNoPhotons2, purityOnePhoton2, purityManyPhotons2, cosmicTwoPhotons] 
pList3 = [puritySignal3, purityMuon3, purityElectron3, purityFiducials3, purityPion3, purityProton3,	purityNoPhotons3, purityOnePhoton3, purityTwoPhotons3, cosmicThreePhotons]

pProtonList1 = [puritySignalProton1, purityOffSignalProton1, purityMuonProton1, purityElectronProton1, purityFiducialsProton1, purityPionProton1, purityProtonProton1, purityNoPhotonsProton1, purityTwoPhotonsProton1, purityManyPhotonsProton1, cosmicProtonOnePhoton]
pProtonList2 = [puritySignalProton2, purityMuonProton2, purityElectronProton2, purityFiducialsProton2, purityPionProton2, purityProtonProton2, purityNoPhotonsProton2, purityOnePhotonProton2, purityManyPhotonsProton2, cosmicProtonTwoPhotons] 
pProtonList3 = [puritySignalProton3, purityMuonProton3, purityElectronProton2, purityFiducialsProton3, purityPionProton3, purityProtonProton3, purityNoPhotonsProton3, purityOnePhotonProton3, purityTwoPhotonsProton3, cosmicProtonThreePhotons]


#EFFICIENCY HISTOGRAMS

#Non-Proton Efficiency Hists
effTotal1 = rt.TH1F("effTotal1", "One Photon",60,0,2)
effNoVertex1 = rt.TH1F("effNoVertex1", "No Vertex Found",60,0,2)
#effCC1 = rt.TH1F("effCC1", "CC False Positive",60,0,2)
effMuon1 = rt.TH1F("effMuon1", "Muon False Positive",60,0,2)
effElectron1 = rt.TH1F("effElectron1", "Electron False Positive",60,0,2)
effFiducial1 = rt.TH1F("effFiducial1", "Placed out of Fiducial",60,0,2)
effPion1 = rt.TH1F("effPion1", "Pion False Positive",60,0,2)
effProton1 =  rt.TH1F("effProton1", "Proton False Positive",60,0,2)
effNoPhotons1 = rt.TH1F("effNoPhotons1", "No Photons Found",60,0,2)
effSignal1 = rt.TH1F("effSignal 1", "Signal",60,0,2)
effTwoPhotons1 = rt.TH1F("effTwoPhotons1", "Two Photons Found",60,0,2)
effManyPhotons1 = rt.TH1F("effManyPhotons1", "Many Photons Found",60,0,2)
effShowerCharge1 = rt.TH1F("effShowerCharge1", "Shower from Charged Cut",60,0,2)
effPrimary1 = rt.TH1F("effPrimary1", "Primary Score Cut",60,0,2)
effLongTracks1 = rt.TH1F("effLongTracks1", "Tracks with length > 20 cm",60,0,2)

effTotal2 = rt.TH1F("effTotal2", "Two Photons",60,0,2)
effNoVertex2 = rt.TH1F("effNoVertex2", "No Vertex Found",60,0,2)
#effCC2 = rt.TH1F("effCC2", "CC False Positive",60,0,2)
effMuon2 = rt.TH1F("effMuon2", "Muon False Positive",60,0,2)
effElectron2 = rt.TH1F("effElectron2", "Electron False Positive",60,0,2)
effFiducial2 = rt.TH1F("effFiducial2", "Placed out of Fiducial",60,0,2)
effPion2 = rt.TH1F("effPion2", "Pion False Positive",60,0,2)
effProton2 =  rt.TH1F("effProton2", "Proton False Positive",60,0,2)
effNoPhotons2 =	rt.TH1F("effNoPhotons2", "No Photons Found",60,0,2)
effSignal2 = rt.TH1F("effSignal2", "Signal",60,0,2)
effOnePhoton2 = rt.TH1F("effOnePhoton2", "One Photon Found",60,0,2)
effManyPhotons2 = rt.TH1F("effManyPhotons2", "Many Photons Found",60,0,2)
effShowerCharge2 = rt.TH1F("effShowerCharge2", "Shower from Charged Cut",60,0,2)
effPrimary2 = rt.TH1F("effPrimary2", "Primary Score Cut",60,0,2)
effLongTracks2 = rt.TH1F("effLongTracks2", "Tracks with length > 20 cm",60,0,2)

effTotal3 = rt.TH1F("effTotal3", "3+ Photons",60,0,2)
effNoVertex3 = rt.TH1F("effNoVertex3", "No Vertex Found",60,0,2)
#effCC3 = rt.TH1F("effCC3", "CC False Positive",60,0,2)
effMuon3 = rt.TH1F("effMuon3", "Muon False Positive",60,0,2)
effElectron3 = rt.TH1F("effElectron3", "Electron False Positive",60,0,2)
effFiducial3 = rt.TH1F("effFiducial3", "Placed out of Fiducial",60,0,2)
effPion3 = rt.TH1F("effPion3", "Pion False Positive",60,0,2)
effProton3 =  rt.TH1F("effProton3", "Proton False Positive",60,0,2)
effNoPhotons3 =	rt.TH1F("effNoPhotons3", "No Photons Found",60,0,2)
effSignal3 = rt.TH1F("effSignal 3", "Signal",60,0,2)
effOnePhoton3 = rt.TH1F("effManyPhotons3", "Many Photons Found",60,0,2)
effTwoPhotons3 = rt.TH1F("effTwoPhotons3", "Two Photons Found",60,0,2)
effShowerCharge3 = rt.TH1F("effShowerCharge3", "Shower from Charged Cut",60,0,2)
effPrimary3 = rt.TH1F("effPrimary3", "Primary Score Cut",60,0,2)
effLongTracks3 = rt.TH1F("effLongTracks3", "Tracks with length > 20 cm",60,0,2)

#Proton Efficiency Hists
effTotalProton1 = rt.TH1F("effTotalProton1", "One Photon (with proton)",60,0,2)
effNoVertexProton1 = rt.TH1F("effNoVertexProton1", "No Vertex Found (with proton)",60,0,2)
#effCCProton1 = rt.TH1F("effCCProton1", "CC False Positive (with proton)",60,0,2)
effMuonProton1 = rt.TH1F("effMuonProton1", "Muon False Positive",60,0,2)
effElectronProton1 = rt.TH1F("effElectronProton1", "Electron False Positive",60,0,2)
effFiducialProton1 = rt.TH1F("effFiducialProton1", "Placed out of Fiducial (with proton)",60,0,2)
effPionProton1 = rt.TH1F("effPionProton1", "Pion False Positive (with proton)",60,0,2)
effProtonProton1 =  rt.TH1F("effProtonProton1", "Proton False Positive (with proton)",60,0,2)
effNoPhotonsProton1 = rt.TH1F("effNoPhotonsProton1", "No Photons Found (with proton)",60,0,2)
effSignalProton1 = rt.TH1F("effSignalProton1", "Signal (with proton)",60,0,2)
effTwoPhotonsProton1 = rt.TH1F("effTwoPhotonsProton1", "Two Photons Found (with proton)",60,0,2)
effManyPhotonsProton1 = rt.TH1F("effManyPhotonsProton1", "Many Photons Found (with proton)",60,0,2)
effShowerChargeProton1 = rt.TH1F("effShowerChargeProton1", "Shower from Charged Cut (with proton)",60,0,2)
effPrimaryProton1 = rt.TH1F("effPrimaryProton1", "Primary Score Cut (with proton)",60,0,2)
effLongTracksProton1 = rt.TH1F("effLongTracksProton1", "Tracks with length > 20 cm (with proton)",60,0,2)

effTotalProton2 = rt.TH1F("effTotalProton2", "Two Photons (with proton)",60,0,2)
effNoVertexProton2 = rt.TH1F("effNoVertexProton2", "No Vertex Found (with proton)",60,0,2)
#effCCProton2 = rt.TH1F("effCCProton2", "CC False Positive (with proton)",60,0,2)
effMuonProton2 = rt.TH1F("effMuonProton2", "Muon False Positive",60,0,2)
effElectronProton2 = rt.TH1F("effElectronProton2", "Electron False Positive",60,0,2)
effFiducialProton2 = rt.TH1F("effFiducialProton2", "Placed out of Fiducial (with proton)",60,0,2)
effPionProton2 = rt.TH1F("effPionProton2", "Pion False Positive (with proton)",60,0,2)
effProtonProton2 =  rt.TH1F("effProtonProton2", "Proton False Positive (with proton)",60,0,2)
effNoPhotonsProton2 =	rt.TH1F("effNoPhotonProton2", "No Photons Found (with proton)",60,0,2)
effSignalProton2 = rt.TH1F("effSignalProton2", "Signal (with proton)",60,0,2)
effOnePhotonProton2 = rt.TH1F("effOnePhotonProton2", "One Photon Found (with proton)",60,0,2)
effManyPhotonsProton2 = rt.TH1F("effManyPhotonsProton2", "Many Photons Found (with proton)",60,0,2)
effShowerChargeProton2 = rt.TH1F("effShowerChargeProton2", "Shower from Charged Cut (with proton)",60,0,2)
effPrimaryProton2 = rt.TH1F("effPrimaryProton2", "Primary Score Cut (with proton)",60,0,2)
effLongTracksProton2 = rt.TH1F("effLongTracksProton2", "Tracks with length > 20 cm (with proton)",60,0,2)

effTotalProton3 = rt.TH1F("effTotalProton3", "3+ Photon (with proton)",60,0,2)
effNoVertexProton3 = rt.TH1F("effNoVertexProton3", "No Vertex Found (with proton)",60,0,2)
#effCCProton3 = rt.TH1F("effCCProton3", "CC False Positive (with proton)",60,0,2)
effMuonProton3 = rt.TH1F("effMuonProton3", "Muon False Positive",60,0,2)
effElectronProton3 = rt.TH1F("effElectronProton3", "Electron False Positive",60,0,2)
effFiducialProton3 = rt.TH1F("effFiducialProton3", "Placed out of Fiducial (with proton)",60,0,2)
effPionProton3 = rt.TH1F("effPionProton3", "Pion False Positive (with proton)",60,0,2)
effProtonProton3 =  rt.TH1F("effProtonProton3", "Proton False Positive (with proton)",60,0,2)
effNoPhotonsProton3 =	rt.TH1F("effNoPhotonsProton3", "No Photons Found (with proton)",60,0,2)
effSignalProton3 = rt.TH1F("effSignalProton3", "Signal (with proton)",60,0,2)
effOnePhotonProton3 = rt.TH1F("effManyPhotonsProton3", "Many Photons Found (with proton)",60,0,2)
effTwoPhotonsProton3 = rt.TH1F("effTwoPhotonsProton3", "Two Photons Found (with proton)",60,0,2)
effShowerChargeProton3 = rt.TH1F("effShowerChargeProton3", "Shower from Charged Cut (with proton)",60,0,2)
effPrimaryProton3 = rt.TH1F("effPrimaryProton3", "Primary Score Cut (with proton)",60,0,2)
effLongTracksProton3 = rt.TH1F("effLongTracksProton3", "Tracks with length > 20 cm (with proton)",60,0,2)

#Histogram Lists!
#Non-Proton Lists
effTotalList = [effTotal1, effTotal2, effTotal3]
effNoVertexHists = [effNoVertex1, effNoVertex2, effNoVertex3]
#effCCHists = [effCC1, effCC2, effCC3]
effMuonHists = [effMuon1, effMuon2, effMuon3]
effElectronHists = [effElectron1, effElectron2, effElectron3]
effFiducialHists = [effFiducial1, effFiducial2, effFiducial3]
effPionHists = [effPion1, effPion2, effPion3]
effProtonHists = [effProton1, effProton2, effProton3]
effNoPhotonHists = [effNoPhotons1, effNoPhotons2, effNoPhotons3]
effOnePhotonHists = [effSignal1, effOnePhoton2, effOnePhoton3] 
effTwoPhotonHists = [effTwoPhotons1, effSignal2, effTwoPhotons3]
effManyPhotonHists = [effManyPhotons1, effManyPhotons2, effSignal3]
effShowerChargeHists = [effShowerCharge1, effShowerCharge2, effShowerCharge3]
effLongTrackHists = [effLongTracks1, effLongTracks2, effLongTracks3]
effPrimaryHists = [effPrimary1, effPrimary2, effPrimary3]

#Proton Lists
effTotalProtonList = [effTotalProton1, effTotalProton2, effTotalProton3]
effNoVertexProtonHists = [effNoVertexProton1, effNoVertexProton2, effNoVertexProton3]
#effCCProtonHists = [effCCProton1, effCCProton2, effCCProton3]
effMuonProtonHists = [effMuonProton1, effMuonProton2, effMuonProton3]
effElectronProtonHists = [effElectronProton1, effElectronProton2, effElectronProton3]
effFiducialProtonHists = [effFiducialProton1, effFiducialProton2, effFiducialProton3]
effPionProtonHists = [effPionProton1, effPionProton2, effPionProton3]
effProtonProtonHists = [effProtonProton1, effProtonProton2, effProtonProton3]
effNoPhotonProtonHists = [effNoPhotonsProton1, effNoPhotonsProton2, effNoPhotonsProton3]
effOnePhotonProtonHists = [effSignalProton1, effOnePhotonProton2, effOnePhotonProton3] 
effTwoPhotonProtonHists = [effTwoPhotonsProton1, effSignalProton2, effTwoPhotonsProton3]
effManyPhotonProtonHists = [effManyPhotonsProton1, effManyPhotonsProton2, effSignalProton3]
effShowerChargeProtonHists = [effShowerChargeProton1, effShowerChargeProton2, effShowerChargeProton3]
effLongTrackProtonHists = [effLongTracksProton1, effLongTracksProton2, effLongTracksProton3]
effPrimaryProtonHists = [effPrimaryProton1, effPrimaryProton2, effPrimaryProton3]

#Big Lists (no protons)
effList1 = [effSignal1, effNoVertex1, effMuon1, effElectron1, effPion1, effProton1, effShowerCharge1, effNoPhotons1, effTwoPhotons1, effManyPhotons1, effPrimary1, effLongTracks1]
effList2 = [effSignal2, effNoVertex2, effMuon2, effElectron2, effPion2, effProton2, effShowerCharge1, effNoPhotons2, effOnePhoton2, effManyPhotons2, effPrimary2, effLongTracks2]
effList3 = [effSignal3, effNoVertex3, effMuon3, effElectron3, effPion3, effProton3, effShowerCharge1, effNoPhotons3, effOnePhoton3, effTwoPhotons3, effPrimary3, effLongTracks3]

#Big Lists (one proton)
effListProton1 = [effSignalProton1, effNoVertexProton1, effMuonProton1, effElectronProton1, effPionProton1, effProtonProton1, effShowerChargeProton1, effNoPhotonsProton1, effTwoPhotonsProton1, effManyPhotonsProton1, effPrimaryProton1, effLongTracksProton1]
effListProton2 = [effSignalProton2, effNoVertexProton2, effMuonProton2, effElectronProton2, effPionProton2, effProtonProton2, effShowerChargeProton1, effNoPhotonsProton2, effOnePhotonProton2, effManyPhotonsProton2, effPrimaryProton2, effLongTracksProton2]
effListProton3 = [effSignalProton3, effNoVertexProton3, effMuonProton3, effElectronProton3, effPionProton3, effProtonProton3, effShowerChargeProton1, effNoPhotonsProton3, effOnePhotonProton3, effTwoPhotonsProton3, effPrimaryProton3, effLongTracksProton3]


#Built-in functions here
def addHist(eventTree, photonList, photonList2, histList, histListProton, nProtons, variable, weight):
  #Fill non-proton histograms based on photon count
  if nProtons == 0:
    if len(photonList) + len(photonList2) == 1:
      histList[0].Fill(variable, weight)
    elif len(photonList) + len(photonList2) == 2:
      histList[1].Fill(variable, weight)
    else:
      histList[2].Fill(variable, weight)
  #Fill proton histograms based on photon count
  elif nProtons == 1:
    if len(photonList) + len(photonList2) == 1:
      histListProton[0].Fill(variable, weight)
    elif len(photonList) + len(photonList2) == 2:
      histListProton[1].Fill(variable, weight)
    else:
      histListProton[2].Fill(variable, weight)

#Built-in functions here
def addSignalHist(ntuple, photonList, photonList2, offHist, offHistProton, histList, histListProton, nProtons, nTrueProtons, variable, weight):
  #Fill non-proton histograms based on photon count
  if nProtons == 0:
    if nTrueProtons == 1 and len(photonList) + len(photonList2) == 1:
      offHist.Fill(variable, weight)
    else:
      if len(photonList) + len(photonList2) == 1:
        histList[0].Fill(variable, weight)
      elif len(photonList) + len(photonList2) == 2:
        histList[1].Fill(variable, weight)
      else:
        histList[2].Fill(variable, weight)
  #Fill proton histograms based on photon count
  elif nProtons == 1:
    if nTrueProtons == 0 and len(photonList) + len(photonList2) == 1:
      offHistProton.Fill(variable, weight)
    else:
      if len(photonList) + len(photonList2) == 1:
        histListProton[0].Fill(variable, weight)
      elif len(photonList) + len(photonList2) == 2:
        histListProton[1].Fill(variable, weight)
      else:
        histListProton[2].Fill(variable, weight)

def purityStack(title, purityList, cosmicHist, POTSum, cosmicSum):
  #Create Component Variables
  stack = rt.THStack("PhotonStack", str(title))
  legend = rt.TLegend(0.5, 0.5, 0.9, 0.9)
  colors = [rt.kBlue, rt.kRed, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow+2, rt.kGreen+4, rt.kOrange+1]
  POTTarget = 6.67e+20
  histIntTotal = 0
  #Organize other histograms
  for x in range(len(purityList)):
    hist = purityList[x]
    bins = hist.GetNbinsX()
    hist.Scale(POTTarget/POTSum)
    if x == 0:
      hist.SetLineColor(rt.kGreen)
    else:
      hist.SetLineColor(colors[x%7])    
    histInt = hist.Integral(1, int(bins))
    histIntTotal += histInt
    legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
    stack.Add(hist)

  #Format and add the cosmic histogram
  hist = cosmicHist
  bins = hist.GetNbinsX()
  hist.Scale(POTTarget/cosmicSum)
  hist.SetLineColor(rt.kBlack)
  histInt = hist.Integral(1, int(bins))
  histIntTotal += histInt
  legend.AddEntry(hist, str(hist.GetTitle())+": "+str(round(histInt, 1)), "l")
  legendHeaderString = "Total: " + str(round((histIntTotal),1)) 
  legend.SetHeader(str(legendHeaderString), "C")
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

#Variables for program review
initialCount = 0
vertexCount = 0
NCCount = 0
fiducialCount = 0
noPionCount = 0
recoCount = 0
onePhoton = 0
twoPhotons = 0
threePhotons = 0
count = 0

totalCosmics = 0
vertexCosmics = 0 
NCCosmics = 0
MattProofCosmics = 0
fiducialCosmics = 0
pionlessCosmics = 0
photonCosmics = 0
uncutCosmics = 0
untrackedCosmics = 0

passTruth = 0
hasVertex = 0
hasNC = 0
inFiducial = 0
pionProtonFine = 0
survivesCuts = 0
noEffPhotons = 0
oneEffPhoton = 0
twoEffPhotons = 0
manyEffPhotons = 0

#We put this into the addHist function for truth-based graphs
emptyList = []

#Variables for program function
fiducialData = {"xMin":0, "xMax":256, "yMin":-116.5, "yMax":116.5, "zMin":0, "zMax":1036, "width":30}


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

  #Use Matt's Cosmic Cut
  if trueCutCosmic(eventTree) == False:
    continue 

  #Make sure the event is within the fiducial volume
  if recoFiducials(eventTree, fiducialData) == False:
    continue

  #Cut events with suitably energetic protons or charged pions 
  if recoPion(eventTree) == False:
    continue

  #Cut events with far-travelling protons
  recoProtonCount = recoProton(eventTree)
  if recoProtonCount > 1:
    continue

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, eventTree)

  #List all photons classified as tracks
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree)

  if len(recoList) + len(recoTrackList) == 0:
    continue

  #Try cutting based on data for Shower from Charged
  if recoCutShowerFromChargeScore(eventTree, recoList, recoTrackList) == False:
    continue
  
  #Cut based on primary score
  if recoCutPrimary(eventTree, recoList, recoTrackList) == False:
    continue

  #Cut based on the presence of tracks over 20 cm
  if recoCutLongTracks(eventTree, fiducialData) == False:
    continue

  #Cut unclassified tracks too short to be trusted
  #if recoCutShortTracks(eventTree) == False:
  #  continue


  #PURITY - GRAPHING BASED ON TRUTH
  
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(eventTree, recoList, recoTrackList)

  #Fill totals
  addHist(eventTree, recoList, recoTrackList, totalPHists, totalProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
  initialCount += 1
  #Neutral current!
#  if trueCutNC(eventTree) == False:
#    addHist(eventTree, recoList, recoTrackList, CCPHists, CCProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
#    continue

  if trueCutMuons(eventTree) == False:
    addHist(eventTree, recoList, recoTrackList, muonPHists, muonProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
    continue

  if trueCutElectrons(eventTree) == False:
    addHist(eventTree, recoList, recoTrackList, electronPHists, electronProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  NCCount += 1

  if trueCutFiducials(eventTree, fiducialData) == False:
    addHist(eventTree, recoList, recoTrackList, fiducialPHists, fiducialProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  fiducialCount += 1

  #pions and protons!
  pionCount, protonCount = trueCutPionProton(eventTree)
  if pionCount > 0:
    addHist(eventTree, recoList, recoTrackList, pionPHists, pionProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
    continue

  if protonCount > 1:
    addHist(eventTree, recoList, recoTrackList, protonPHists, protonProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  
  noPionCount += 1

  #Now we make a list of the actual photons!
  truePhotonIDs = truePhotonList(eventTree, fiducialData)

  #Are there actually any photons?
  if len(truePhotonIDs) == 0:
    addHist(eventTree, recoList, recoTrackList, noPhotonPHists, noPhotonProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  recoCount += 1
  #Is there one Photon?
  if len(truePhotonIDs) == 1:
    addSignalHist(eventTree, recoList, recoTrackList, purityOffSignal1, purityOffSignalProton1, onePhotonPHists, onePhotonProtonPHists, recoProtonCount, protonCount, leadingPhoton, eventTree.xsecWeight)
    onePhoton += 1
  #Are there two?
  elif len(truePhotonIDs) == 2:
    addHist(eventTree, recoList, recoTrackList, twoPhotonPHists, twoPhotonProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
    twoPhotons += 1
  
  #In that case, there should be at least three
  else:
    addHist(eventTree, recoList, recoTrackList, manyPhotonPHists, manyPhotonProtonPHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
    threePhotons += 1

#BEGINNING EVENT LOOP FOR COSMICS
for i in range(cosmicTree.GetEntries()):
  cosmicTree.GetEntry(i)

  #COSMICS - SELECTING EVENTS BASED ON RECO
  #See if the event has a vertex
  totalCosmics += 1

  if recoNoVertex(cosmicTree) == False:
    continue
  vertexCosmics += 1
  #See if the event is neutral current                                                                                                  
  if recoNeutralCurrent(cosmicTree) == False:
    continue
  NCCosmics += 1
  #Use Matt's Cosmic Cut
  if trueCutCosmic(cosmicTree) == False:
    continue 
  MattProofCosmics += 1
  #Make sure the event is within the fiducial volume
  if recoFiducials(cosmicTree, fiducialData) == False:
    continue
  fiducialCosmics += 1
  #Cut events with suitably energetic protons or charged pions
  if recoPion(cosmicTree) == False:
    continue
  pionlessCosmics += 1
  #Cut events with far-travelling protons
  recoProtonCount = recoProton(cosmicTree)
  if recoProtonCount > 1:
    continue

  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, cosmicTree)

  #List all photons classified as tracks
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree)
  
  if len(recoList) + len(recoTrackList) == 0:
    continue
  photonCosmics += 1
  
  #Try cutting based on data for Shower from Charged
  if recoCutShowerFromChargeScore(cosmicTree, recoList, recoTrackList) == False:
    continue

  #Cut based on primary score
  if recoCutPrimary(cosmicTree, recoList, recoTrackList) == False:
    continue
  uncutCosmics += 1

  #Longtracks cut
  if recoCutLongTracks(cosmicTree, fiducialData) == False:
    continue

  #Cut unclassified tracks too short to be trusted
  #if recoCutShortTracks(eventTree) == False:
  #  continue

  untrackedCosmics += 1

  #graphing based on photon count
  #Calculating graphing values
  leadingPhoton = scaleRecoEnergy(cosmicTree, recoList, recoTrackList)
  addHist(cosmicTree, recoList, recoTrackList, cosmicList, cosmicProtonList, recoProtonCount, leadingPhoton, 1)

#BEGINNING EVENT LOOP FOR EFFICIENCY
for i in range(eventTree.GetEntries()):
  eventTree.GetEntry(i)

  #Selecting events using truth
  if trueCutMuons(eventTree) == False:
    continue 

  if trueCutElectrons(eventTree) == False:
    continue

  if trueCutFiducials(eventTree, fiducialData) == False:
    continue
  
  pionCount, protonCount = trueCutPionProton(eventTree)
  if pionCount > 0 or protonCount > 1:
    continue
  
  truePhotonIDs = truePhotonList(eventTree, fiducialData)

  if len(truePhotonIDs) == 0:
    continue

  
  #EFFICIENCY - GRAPHING BASED ON RECO
  passTruth += 1
  leadingPhoton = scaleTrueEnergy(eventTree, truePhotonIDs)
  recoProtonCount = recoProton(eventTree)

  if recoNoVertex(eventTree) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effNoVertexHists, effNoVertexProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  hasVertex += 1
  #See if the event is neutral current                                                                                                  
  #if recoNeutralCurrent(eventTree) == False:
  #  addHist(eventTree, truePhotonIDs, emptyList, effCCHists, effCCProtonHists, recoProtonCount, leadingPhoton, eventTree.xsecWeight)
  #  continue

  #See if we reconstruct muons or electrons over threshold
  if recoCutMuons(eventTree) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effMuonHists, effMuonProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  if recoCutElectrons(eventTree) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effElectronHists, effElectronProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue

  hasNC += 1

  #Cut events with vertexes outside the fiducial
  if recoFiducials(eventTree, fiducialData) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effFiducialHists, effFiducialProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  inFiducial += 1
  #Cut events with too many protons
  recoProtonCount = recoProton(eventTree)
  if recoProtonCount > 1:
    addHist(eventTree, truePhotonIDs, emptyList, effProtonHists, effProtonProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut events with pions present
  if recoPion(eventTree) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effPionHists, effPionProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  pionProtonFine += 1
  #See if there are any photons in the event - if so, list them
  recoList = recoPhotonListFiducial(fiducialData, eventTree)
  recoTrackList = recoPhotonListTracks(fiducialData, eventTree)

  #Try cutting for Shower from Charged Score 
  if recoCutShowerFromChargeScore(eventTree, recoList, recoTrackList) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effShowerChargeHists, effShowerChargeProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue

  #Cut based on Electron Score
  if recoCutPrimary(eventTree, recoList, recoTrackList) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effPrimaryHists, effPrimaryProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue

  #Try cutting based on data for Track Lengths
  if recoCutLongTracks(eventTree, fiducialData) == False:
    addHist(eventTree, truePhotonIDs, emptyList, effLongTrackHists, effLongTrackProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    continue
  survivesCuts += 1
  #Cut unclassified tracks too short to be trusted
  #if recoCutShortTracks(eventTree) == False:
  #  addHist(eventTree, truePhotonIDs, emptyList, effShortTrackHists, effShortTrackProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
  #  continue

  #Now we're pretty sure the event is legitimate, so we go ahead and graph based on the number of photons
  if len(recoList) + len(recoTrackList) == 0:
    addHist(eventTree, truePhotonIDs, emptyList, effNoPhotonHists, effNoPhotonProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    noEffPhotons += 1
  elif len(recoList) + len(recoTrackList) == 1:
    addHist(eventTree, truePhotonIDs, emptyList, effOnePhotonHists, effOnePhotonProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    oneEffPhoton += 1
  elif len(recoList) + len(recoTrackList) == 2:
    addHist(eventTree, truePhotonIDs, emptyList, effTwoPhotonHists, effTwoPhotonProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    twoEffPhotons += 1
  else:
    addHist(eventTree, truePhotonIDs, emptyList, effManyPhotonHists, effManyPhotonProtonHists, protonCount, leadingPhoton, eventTree.xsecWeight)
    manyEffPhotons += 1
#LOOPS OVER - HISTOGRAM ORGANIZING TIME

for hist in cosmicList:
  hist.Scale(ntuplePOTsum/cosmicPOTsum)

#Stacking non-proton histograms
purityCanvas1, purityStack1, purityLegend1, purityInt1 = histStackTwoSignal("1 Gamma + 0 Sample", pList1, ntuplePOTsum)
purityCanvas2, purityStack2, purityLegend2, purityInt2 = histStack("2 Gamma + 0 Sample", pList2, ntuplePOTsum)
purityCanvas3, purityStack3, purityLegend3, purityInt3 = histStack("3+ Gamma + 0 Sample", pList3, ntuplePOTsum)
effCanvas1, effStack1, effLegend1, effInt1 = histStack("True 1 Gamma + 0 Outcomes", effList1, ntuplePOTsum)
effCanvas2, effStack2, effLegend2, effInt2 = histStack("True 2 Gamma + 0 Outcomes", effList2, ntuplePOTsum)
effCanvas3, effStack3, effLegend3, effInt3 = histStack("True 3+ Gamma + 0 Outcomes", effList3, ntuplePOTsum)

cosmicCanvas, cosmicStack, cosmicLegend, cosmicInt = histStack("Cosmic Background, Reconstructed with 0 Protons", cosmicList, cosmicPOTsum)

#Stacking proton histograms
purityCanvasProton1, purityStackProton1, purityLegendProton1, purityIntProton1 = histStackTwoSignal("Single-Photon + Proton Purity", pProtonList1, ntuplePOTsum)
purityCanvasProton2, purityStackProton2, purityLegendProton2, purityIntProton2 = histStack("Two-Photon + Proton Purity", pProtonList2, ntuplePOTsum)
purityCanvasProton3, purityStackProton3, purityLegendProton3, purityIntProton3 = histStack("3+ Photon + Proton Purity", pProtonList3, ntuplePOTsum)
effCanvasProton1, effStackProton1, effLegendProton1, effIntProton1 = histStack("Single-Photon + Proton Efficiency", effListProton1, ntuplePOTsum)
effCanvasProton2, effStackProton2, effLegendProton2, effIntProton2 = histStack("Two-Photon + Proton Efficiency", effListProton2, ntuplePOTsum)
effCanvasProton3, effStackProton3, effLegendProton3, effIntProton3 = histStack("3+ Photon + Proton Efficiency", effListProton3, ntuplePOTsum)

cosmicProtonCanvas, cosmicProtonStack, cosmicProtonLegend, cosmicProtonInt = histStack("Cosmic Background, Reconstructed with 1 Proton", cosmicProtonList, cosmicPOTsum)

writeList = [purityCanvas1, purityCanvasProton1, purityCanvas2, purityCanvasProton2, purityCanvas3, purityCanvasProton3, effCanvas1, effCanvasProton1, effCanvas2, effCanvasProton2, effCanvas3, effCanvasProton3, cosmicCanvas, cosmicProtonCanvas]


#Now all that's left to do is write the canvases to the file
outFile = rt.TFile(args.outfile, "RECREATE")
for canvas in writeList:
  canvas.Write()

print("PURITY STATS")
print("Vertex reconstructed:", initialCount)
print("Neutral Current:", NCCount)
print("In Fiducial:", fiducialCount)
print("No pions, protons:", noPionCount)
print("Fully reconstructed:", recoCount)
print(onePhoton, "events had one photon")
print(twoPhotons, "events had two photons")
print(threePhotons, "events had 3+ photons")

print("EFFICIENCY STATS")
print("Total signal space:", passTruth)
print("Total with vertex:", hasVertex)
print("Total NC:", hasNC)
print("Total reconstructed in Fiducial:", inFiducial)
print("Total with acceptable pion/protons:", pionProtonFine)
print("Total that survive our cuts:", survivesCuts)
print("Total with no photons:", noEffPhotons)
print("Total with one photon:", oneEffPhoton)
print("Total with two photons:", twoEffPhotons)
print("Total with three photons:", manyEffPhotons)


print("POTsum:", ntuplePOTsum)