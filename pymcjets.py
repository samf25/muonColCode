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
dPtBottom = []
dRBottom = []
dPtTrue = []
dRTrue = []
dPtCharged = []
dRCharged = []

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

# A function to find the closest distance between a jet/particle and a set of jets
def findClose(j1, jets, thresh):
    distances = []
    for jet in jets:
        distances.append(j1.delta_R(jet))
    if np.min(distances) < thresh:
        return np.min(distances), (jets[distances.index(np.min(distances))].pt() / j1.pt())
    else:
        return None, None 

# Loop through each event
count = 0
for event in tree:
    # Only work with non-empty events
    if len(event.mcpdg) == 0 or len(event.tsrpx) == 0:
        pass
    else:
        count += 1
        endStates = []
        chargedEndStates = []
        tracks = []
        bs = [None, None]
        # Loop through truth table and get endstates and b quarks
        for i in range(len(event.mcpdg)):
            if event.mcgst[i] == 1:
                part = fj.PseudoJet(event.mcmox[i],event.mcmoy[i],event.mcmoz[i],event.mcene[i]))
                endStates.append(part)
                if event.mccha[i] != 0:
                    chargedEndStates.append(part)
            elif event.mcpdg[i] == 5 and bs[0] == None:
                bs[0] = fj.PseudoJet(event.mcmox[i],event.mcmoy[i],event.mcmoz[i],event.mcene[i])
            elif event.mcpdg[i] == -5 and bs[1] == None:
                bs[1] = fj.PseudoJet(event.mcmox[i],event.mcmoy[i],event.mcmoz[i],event.mcene[i])
        # Loop through track data and get particles above a threshold pt
        for i in range(len(event.tsrpx)):
            momentum = [event.tsrpx[i]/1000, event.tsrpy[i]/1000, event.tsrpz[i]/1000]
            ene = np.sqrt(momentum[0]**2 +momentum[1]**2 +momentum[2]**2 +mpi**2)
            part = fj.PseudoJet(momentum[0], momentum[1], momentum[2], ene)
            if part.pt() > 0.05:
                tracks.append(part)
        
        # Run Anti-kt algorithm for out 3 types of jets
        clusterT = fj.ClusterSequence(endStates, jetdef)
        clusterC = fj.ClusterSequence(chargedEndStates, jetdef)
        clusterR = fj.ClusterSequence(tracks, jetdef)
        trueJets = clusterT.inclusive_jets(1)
        charJets = clusterC.inclusive_jets(1)
        recoJets = clusterR.inclusive_jets(1)

        # Loop through the bs and if they are close to a jet, record
        for b in bs:
            mi, pTRatio = findClose(b, recoJets, rThresh)
            if mi != None:
                dPtBottom.append(pTRatio)
                dRBottom.append(mi)
        # Loop through truth jets and if they are close to a jet, record
        for t in trueJets:
            mi, pTRatio = findClose(t, recoJets, rThresh)
            if mi != None:
                dPtTrue.append(pTRatio)
                dRTrue.append(mi)
        # Loop through charged jets and if they are close to a jet, record
        for c in charJets:
            mi, pTRatio = findClose(c, recoJets, rThresh)
            if mi != None:
                dPtCharged.append(pTRatio)
                dRCharged.append(mi)
        
        # User Output to track progress
        if count % 5 == 0:
            print("Event "+ str(count), end = "\r")

# Output Histograms
histgraph(dPtBottom, "Jet/b Responce Plot", "Graphs/bResponce.pdf")
histgraph(dPtTrue, "Jet/truth Responce Plot", "Graphs/truthResponce.pdf")
histgraph(dPtCharged, "Jet/charged Responce Plot", "Graphs/chargedResponce.pdf")
histgraph(dRBottom, "Jet-b Delta R", "Graphs/bDeltaR.pdf")
histgraph(dRTrue, "Jet-truth Delta R", "Graphs/truthDeltaR.pdf")
histgraph(dRCharged, "Jet-charged Delta R", "Graphs/chargedDeltaR.pdf")
