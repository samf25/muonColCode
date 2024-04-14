import ROOT
import sys
import fastjet as fj
import numpy as np

# Get Data
file = ROOT.TFile.Open(sys.argv[1], "READ")
tree = file.Get("MyLCTuple;42")

# Constants
mpi = .140
jetThresh = 5
rThresh = 0.2
jetdef = fj.JetDefinition(fj.antikt_algorithm, 0.4)

# Save Arrays
deltaPt = []
deltaR = []

# A function to make a histogram
def histgraph(data, title, file, bins = 40, minimum = 0, maximum = -1):
    canvas = ROOT.TCanvas("c"+title)
    if maximum == -1:
        maximum = np.max(data)
    if np.min(data) < 0:
        minimum = np.min(data)
    g = ROOT.TH1F(title, title, bins, minimum, maximum)
    for p in data:
        g.Fill(p)
    g.Draw()

    canvas.Print(file)
    canvas.Close()

# Loop through each event
count = 0
for event in tree:
    if len(event.mcpdg) == 0 or len(event.tsrpx) == 0:
        pass
    else:
        count += 1
        endStates = []
        tracks = []
        bs = [None, None]
        for i in range(len(event.mcpdg)):
            if event.mcgst[i] == 1:
                endStates.append(fj.PseudoJet(event.mcmox[i],event.mcmoy[i],event.mcmoz[i],event.mcene[i]))
            elif event.mcpdg[i] == 5 and bs[0] == None:
                bs[0] = fj.PseudoJet(event.mcmox[i],event.mcmoy[i],event.mcmoz[i],event.mcene[i])
            elif event.mcpdg[i] == -5 and bs[1] == None:
                bs[1] = fj.PseudoJet(event.mcmox[i],event.mcmoy[i],event.mcmoz[i],event.mcene[i])
        for i in range(len(event.tsrpx)):
            momentum = [event.tsrpx[i]/1000, event.tsrpy[i]/1000, event.tsrpz[i]/1000]
            ene = np.sqrt(momentum[0]**2 +momentum[1]**2 +momentum[2]**2 +mpi**2)
            part = fj.PseudoJet(momentum[0], momentum[1], momentum[2], ene)
            if part.pt() > 0.5:
                tracks.append(part)
        
        cluster = fj.ClusterSequence(endStates, jetdef)
        jets = cluster.inclusive_jets(1)
        
        bigJets = []  
        for jet in jets:
            if jet.pt() > jetThresh:
                bigJets.append(jet)
        for b in bs:
            distances = []
            for jet in bigJets:
                distances.append(b.delta_R(jet))
            mi = np.min(distances)
            if mi < rThresh:
                deltaPt.append(-(b.pt()-bigJets[distances.index(mi)].pt())/b.pt())
                deltaR.append(mi)
        if count % 5 == 0:
            print("Event "+ str(count), end = "\r")

print("There were "+str(len(deltaR))+" 'close' bs.")
histgraph(deltaR, "Delta R", "Graphs/deltar.pdf")
histgraph(deltaPt, "Delta pT", "Graphs/deltapt.pdf")

