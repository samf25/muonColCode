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

for event in tree: 
    # Collect all the UNIQUE b quarks
    bqs = []
    for i in range(len(event.mcpdg)):
        # Collect the first 2 b quarks and add to list
        if abs(event.mcpdg[i]) == 5:
            check = True
            part = fj.PseudoJet(event.mcmox[i], event.mcmoy[i], event.mcmoz[i], event.mcene[i])
            bqs.append(part)
            bpT.append(part.pt())
            if len(bqs) == 2:
                break
        if event.mcpdg[i] == 25:
            higgs = fj.PseudoJet(event.mcmox[i], event.mcmoy[i], event.mcmoz[i], event.mcene[i])

    if len(bqs) == 2:
        pH = [higgs.px(), higgs.py(), higgs.pz(), higgs.E()]
        p1 = [bqs[0].px(),bqs[0].py(),bqs[0].pz(),bqs[0].E()]
        p2 = [bqs[1].px(),bqs[1].py(),bqs[1].pz(),bqs[1].E()]
        diff = [(pH[i]-p1[i]-p2[i])/pH[i] for i in range(len(pH))]
        print(diff)

# A function to make a histogram
def histgraph(data, title, file, bins = 40):
    canvas = ROOT.TCanvas("c"+title)
    g = ROOT.TH1F(title, title, bins, 0, np.max(data))
    for p in data:
        g.Fill(p)
    g.Draw()

    canvas.Print(file)
    canvas.Close()

# Create Histograms
#histgraph(pT, "Total Track pT", "total_pT.pdf")
#histgraph(bpT, "Higgs Decay b pT", "bquark_pT.pdf")
#histgraph(trackHits, "Number of Tracks", "trackcount.pdf", int(np.max(trackHits))+1)
#histgraph(smallDR, "Smallest dR Between Jet and b", "smallDR.pdf")
#histgraph(jetpT, "Jet pT", "jetpt.pdf")
#histgraph(ratios, "Jet pT / Close b pT", "jetbpTratio.pdf")
#histgraph(missedpT, "unassigned b pT", "unassignedbpT.pdf")
