import sys, argparse
import numpy as np
import ROOT as rt

parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
args = parser.parse_args()

#open input file and get event tree
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")

# set x/y/z min/max values (cm) of the detector, then set the amount by which to inset the fiducial volume
xMin, xMax = 0, 256
yMin, yMax = -116.5, 116.5
zMin, zMax = 0, 1036
fiducialWidth = 10

fiducialTally = 0
# This loop checks whether each entry's interaction vertex falls outside the fiducial volume, and tallies if it does
for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)

    if eventTree.trueVtxX <= (xMin + fiducialWidth) or eventTree.trueVtxX >= (xMax - fiducialWidth) or \
    eventTree.trueVtxY <= (yMin + fiducialWidth) or eventTree.trueVtxY >= (yMax - fiducialWidth) or \
    eventTree.trueVtxZ <= (zMin + fiducialWidth) or eventTree.trueVtxZ >= (zMax - fiducialWidth):
        fiducialTally += 1

print("There are " +str(fiducialTally)+ " events within " +str(fiducialWidth)+"cm of the edge of the detector.")