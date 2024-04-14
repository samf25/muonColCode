import ROOT
import sys
import fastjet as fj
import numpy as np

jetdef = fj.JetDefinition(fj.antikt_algorithm, 0.4)

# Get data
file = ROOT.TFile.Open(sys.argv[1], "READ")
tree = file.Get("MyLCTuple;42")

# Counters
count = 0
alignedBs = 0
zeroTrackpT = 0

# Constants
mmu = .511
etaPhiThresh = .2
jetpTThresh = 1000

# Data Storage Arrays
pT = []
bpT = []
jetpT = []
trackHits = []
smallDR = []
ratios = []
missedpT = []

# Function to determine the distance between two particles in eta-phi
def etaPhiDist(p1, p2):
    return np.sqrt((p1.pseudorapidity()-p2.pseudorapidity())**2+(p1.phi()-p2.phi())**2)

zeros = []
ones = []
totals = []
codes = []

for event in tree:
    sets = [0,0]
    for i in range(len(event.mcpdg)):
        if abs(event.mcpdg[i]) == 14:
            sets[event.mcgst[i]] += 1
            codes.append(abs(event.mcsst[i]))
    zeros.append(sets[0])
    ones.append(sets[1])
    totals.append(sets[0]+sets[1])

# A function to make a histogram
def histgraph(data, title, file, bins = 40):
    canvas = ROOT.TCanvas("c"+title)
    g = ROOT.TH1F(title, title, bins, 0, np.max(data))
    for p in data:
        g.Fill(p)
    g.Draw()

    canvas.Print(file)
    canvas.Close()

histgraph(totals, "Total Neutrino Count", "totalnu.pdf", int(np.max(totals)))
histgraph(codes, "Neutrino Status Codes", "nucodes.pdf")
histgraph(ones, "One Code Neutrino Count", "onesnu.pdf", int(np.max(ones)))

# Create Histograms
#histgraph(pT, "Total Track pT", "total_pT.pdf")
#histgraph(bpT, "Higgs Decay b pT", "bquark_pT.pdf")
#histgraph(trackHits, "Number of Tracks", "trackcount.pdf", int(np.max(trackHits))+1)
#histgraph(smallDR, "Smallest dR Between Jet and b", "smallDR.pdf")
#histgraph(jetpT, "Jet pT", "jetpt.pdf")
#histgraph(ratios, "Jet pT / Close b pT", "jetbpTratio.pdf")
#histgraph(missedpT, "unassigned b pT", "unassignedbpT.pdf")
