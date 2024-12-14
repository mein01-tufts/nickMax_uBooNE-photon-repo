import sys, argparse
import numpy as np
import ROOT as rt

from darkCuts import particleTallies, protonCut, muonCut, pionCut, electronCut, cleverPhotonCut, trueParticleTallies, trueProtonCut, trueMuonCut, truePionCut, truePhotonCut

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-d", "--darkFile", type=str, required=True, help="darkNu input ntuple filepath")
parser.add_argument("-s", "--stdFile", type=str, required=True, help="standard model input ntuple filepath")
parser.add_argument("-c", "--cosmicFile", type=str, required=True, help="cosmic data input ntuple filepath")
parser.add_argument("-o", "--outFile", type=str, default="MonthDayDarkMatter.root", help="output root file name")
args = parser.parse_args()

# Create root TFile classes for darkNu, standard, and cosmic files
dark_file = rt.TFile(args.darkFile)
std_file = rt.TFile(args.stdFile)
cosmic_file = rt.TFile(args.cosmicFile)


# Grab eventTrees and potTrees for each input for reference
darkTree = dark_file.Get("EventTree")
darkPotTree = dark_file.Get("potTree")

stdTree = std_file.Get("EventTree")
stdPotTree = std_file.Get("potTree")

cosmicTree = cosmic_file.Get("EventTree")

# State target POT for scaling, create POT sums for each input file
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20"
darkPOTsum, stdPOTsum = 0, 0
cosmicPOTsum = 3.2974607516739725e+20

# Find darkNu POT sum
for i in range(darkPotTree.GetEntries()):
    darkPotTree.GetEntry(i)
    darkPOTsum = darkPOTsum + darkPotTree.totGoodPOT

# Find standard POT sum
for i in range(stdPotTree.GetEntries()):
    stdPotTree.GetEntry(i)
    stdPOTsum = stdPOTsum + stdPotTree.totGoodPOT

potSumList = [darkPOTsum, stdPOTsum, cosmicPOTsum]
print(str(potSumList))


# 2 loops over dark file (efficiency, purity)
#LOOP 1 - PURITY
for i in range(darkTree.GetEntries()):
    darkTree.GetEntry(i)
    #Make lists of every track and shower that contain each kind of particle
    electronShowers, photonShowers, muonShowers, pionShowers, protonShowers, electronTracks, protonTracks, muonTracks, pionTracks, protonTracks = particleTallies(darkTree)

    #Now we cut based on the particles we don't want
    if protonCut(protonShowers, protonTracks) == False:
        continue
    if muonCut(muonShowers, muonTracks) == False:
        continue
    if pionCut(pionShowers, pionTracks) == False:
        continue

    #Try to identify events that fit our stipulations
    if electronCut == False and cleverPhotonCut == False:
        continue

    #Load up the histogram, baby!

#LOOP 2 - EFFICIENCY
    electronList, photonList, muonList, pionList, protonList = trueParticleTallies(darkTree)

    if trueProtonCut(protonShowers, protonTracks) == False:
        continue
    if trueMuonCut(muonShowers, muonTracks) == False:
        continue
    if truePionCut(pionShowers, pionTracks) == False:
        continue

