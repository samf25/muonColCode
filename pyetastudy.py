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
trackEta = []
trackTheta = []
beta = []


# Function to determine the distance between two particles in eta-phi
def etaPhiDist(p1, p2):
    return np.sqrt((p1.pseudorapidity()-p2.pseudorapidity())**2+(p1.phi()-p2.phi())**2)

for event in tree:
    # Collect all the Track Particles
    vect4 = []
    totalpt = 0
    hits = 0
    # Loop through each entry. Create particle and add to list if pT is above cutoff
    for i in range(len(event.tsrpx)):
        momentum = [event.tsrpx[i]/1000, event.tsrpy[i]/1000, event.tsrpz[i]/1000]
        ene = np.sqrt(momentum[0]**2 +momentum[1]**2 +momentum[2]**2 +mmu**2)
        part = fj.PseudoJet(momentum[0], momentum[1], momentum[2], ene)
        if part.pt() > 0.5:
            vect4.append(part)
            trackEta.append(part.pseudorapidity())
            trackTheta.append(np.arctan(part.py()/part.px()))
        totalpt += part.pt()
        hits += 1
    if count % 5 == 0:
        print("EVENT "+str(count), end="\r") 
        #print("\tTotal pT: "+str(totalpt))
    pT.append(totalpt)
    trackHits.append(hits)

    count += 1
    
    # Collect all the UNIQUE b quarks
    bqs = []
    for i in range(len(event.mcpdg)):
        # Collect the first 2 b quarks and add to list
        if abs(event.mcpdg[i]) == 5:
            check = True
            part = fj.PseudoJet(event.mcmox[i], event.mcmoy[i], event.mcmoz[i], event.mcene[i])
            bqs.append(part)
            bpT.append(part.pt())
            beta.append(part.pseudorapidity())
            if len(bqs) == 2:
                break
   
# A function to make a histogram
def histgraph(data, title, file, bins = 40):
    canvas = ROOT.TCanvas("c"+title)
    mi = 0
    if np.min(data) < mi:
        mi = np.min(data)
    g = ROOT.TH1F(title, title, bins, mi, np.max(data))
    for p in data:
        g.Fill(p)
    g.Draw()

    canvas.Print(file)
    canvas.Close()

histgraph(trackEta, "Track Etas", "tracketa.pdf")
histgraph(trackTheta, "Track Thetas", "tracktheta.pdf")
histgraph(beta, "b Etas", "beta.pdf")

# Create Histograms
#histgraph(pT, "Total Track pT", "total_pT.pdf")
#histgraph(bpT, "Higgs Decay b pT", "bquark_pT.pdf")
#histgraph(trackHits, "Number of Tracks", "trackcount.pdf", int(np.max(trackHits))+1)
#histgraph(smallDR, "Smallest dR Between Jet and b", "smallDR.pdf")
#histgraph(jetpT, "Jet pT", "jetpt.pdf")
#histgraph(ratios, "Jet pT / Close b pT", "jetbpTratio.pdf")
#histgraph(missedpT, "unassigned b pT", "unassignedbpT.pdf")
