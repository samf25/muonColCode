import ROOT
import sys
import fastjet as fj
import numpy as np

# Get Data
file = ROOT.TFile.Open(sys.argv[1], "READ")
tree = file.Get("MyLCTuple;42")

# Constants
mmu = .511

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
            ene = np.sqrt(momentum[0]**2 +momentum[1]**2 +momentum[2]**2 +mmu**2)
            part = fj.PseudoJet(momentum[0], momentum[1], momentum[2], ene)
            if part.pt() > 0.0005:
                tracks.append(part)
        
        endEtas = [part.pseudorapidity() for part in endStates]
        endPhis = [part.phi() for part in endStates]
        bEtas = [part.pseudorapidity() for part in bs]
        bPhis = [part.phi() for part in bs]
        trackEtas = [part.pseudorapidity() for part in tracks]
        trackPhis = [part.phi() for part in tracks]
        
        canvas = ROOT.TCanvas("canvas")
        canvas.Divide(2,1)
        canvas.cd(1)
        g1 = ROOT.TGraph(len(endStates), np.asarray(endEtas), np.asarray(endPhis))
        g1.SetTitle("Event "+str(count)+": "+str(len(endStates))+" Entries (MC)"+";Pseudorapidity;Phi")
        g1.GetXaxis().SetLimits(-5, 5)
        g1.GetYaxis().SetRangeUser(0.0, 6.283)
        g1.Draw("APT")
        g4 = ROOT.TGraph(len(bs), np.asarray(bEtas), np.asarray(bPhis))
        g4.SetMarkerStyle(8)
        g4.Draw("P")
        canvas.cd(2)
        g2 = ROOT.TGraph(len(tracks), np.asarray(trackEtas), np.asarray(trackPhis))
        g2.SetTitle("Event "+str(count)+": "+str(len(tracks))+" Entries (Tracks)"+";Pseudorapidity;Phi")
        g2.GetXaxis().SetLimits(-5, 5)
        g2.GetYaxis().SetRangeUser(0.0,6.283)
        g2.Draw("APT")
        g3 = ROOT.TGraph(len(bs), np.asarray(bEtas), np.asarray(bPhis))
        g3.SetMarkerStyle(8)
        g3.Draw("P")
        
        canvas.Print("Plots/etaPhiEvent"+str(count)+".pdf")
        if count > 9:
            break

