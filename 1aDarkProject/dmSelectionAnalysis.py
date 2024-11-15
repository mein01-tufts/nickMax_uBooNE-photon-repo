import sys, argparse
import numpy as np
import ROOT as rt

from darkCuts import cut

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="MonthDayDarkMatter.root", help="output root file name")
args = parser.parse_args()

# Grab ntuple file, eventTree, and potTree for reference
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")


# Event Loop: We want to use this to sort events by final states
# IE: Tallying final state particles?
eventsTotal = 0
noElec = 0
oneElec = 0
twoElec = 0
threeElec = 0

#---------------- Begin Tally Loop ----------------#
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
    electrons = 0

    if electrons == 0:
        noElec += 1
    elif electrons == 1:
        oneElec += 1
    elif electrons == 2:
        twoElec += 1
    else:
        threeElec += 1
    print("Hello World")


#---------------- End Tally Loop ----------------#

