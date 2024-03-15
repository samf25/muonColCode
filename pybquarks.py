import ROOT
import sys
import fastjet as fj
import numpy as np

jetdef = fj.JetDefinition(fj.antikt_algorithm, 0.4)

file = ROOT.TFile.Open(sys.argv[1], "READ")
tree = file.Get("MyLCTuple;42")

count = 0
alignedBs = 0
zeroTrackpT = 0

mmu = .511
etaPhiThresh = .2
jetpTThresh = 1000

pT = []
bpT = []
jetpT = []
trackHits = []
smallDR = []
ratios = []
missedpT = []

def etaPhiDist(p1, p2):
    return np.sqrt((p1.pseudorapidity()-p2.pseudorapidity())**2+(p1.phi()-p2.phi())**2)

for event in tree:
    # Collect all the Track Particles
    vect4 = []
    totalpt = 0
    hits = 0
    for i in range(len(event.tsrpx)):
        ene = np.sqrt(event.tsrpx[i]**2 +event.tsrpy[i]**2 +event.tsrpz[i]**2 +mmu**2)
        part = fj.PseudoJet(event.tsrpx[i], event.tsrpy[i], event.tsrpz[i], ene)
        if np.sqrt(part.kt2()) > 10:
            vect4.append(part)
        totalpt += np.sqrt(part.kt2())
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
        if abs(event.mcpdg[i]) == 5:
            check = True
            part = fj.PseudoJet(event.mcmox[i], event.mcmoy[i], event.mcmoz[i], event.mcene[i])
            bqs.append(part)
            bpT.append(np.sqrt(part.kt2()))
            if len(bqs) == 2:
                break
    
    # Makes jets
    cluster = fj.ClusterSequence(vect4, jetdef)
    jets = cluster.inclusive_jets(1)

    '''
    bigjets = []
    for jet in jets:
        if np.sqrt(jet.kt2()) > jetpTThresh:
            bigjets.append(jet)
    '''
    
    dr = []
    check = True
    for b in bqs:
        ds = []
        for jet in jets:
            if check:
                jetpT.append(np.sqrt(jet.kt2()))
                check = False
            ds.append(etaPhiDist(b, jet))
        if len(ds) > 0:
            d = np.min(ds)
            if d < etaPhiThresh:
                alignedBs += 1
                dr.append(d)
                ratios.append(np.sqrt(jets[ds.index(d)].kt2()) / np.sqrt(b.kt2()))
                '''
                print()
                print("\tBottom Quark: ("+str(b.pseudorapidity())+", "+str(b.phi())+")")
                print("\tJet pT:\t\t"+str(np.sqrt(jet.kt2())))
                print("\tBottom Q pT:\t"+str(np.sqrt(b.kt2())))
                '''
            else:
                missedpT.append(np.sqrt(b.kt2()))
    
    if len(dr) > 0:
        smallDR.append(np.min(dr))

print()
print("ALIGNED Bs FOUND: " + str(alignedBs))

def histgraph(data, title, file, bins = 40):
    canvas = ROOT.TCanvas("c"+title)
    g = ROOT.TH1F(title, title, bins, 0, np.max(data))
    for p in data:
        g.Fill(p)
    g.Draw()

    canvas.Print(file)
    canvas.Close()

histgraph(pT, "Total Track pT", "total_pT.pdf")
histgraph(bpT, "Higgs Decay b pT", "bquark_pT.pdf")
histgraph(trackHits, "Number of Tracks", "trackcount.pdf", int(np.max(trackHits))+1)
histgraph(smallDR, "Smallest dR Between Jet and b", "smallDR.pdf")
histgraph(jetpT, "Jet pT", "jetpt.pdf")
histgraph(ratios, "Jet pT / Close b pT", "jetbpTratio.pdf")
histgraph(missedpT, "unassigned b pT", "unassignedbpT.pdf")

'''
cleanbs = [[],[],[],[],[]]
threshhold = .1
for i in range(len(bqs)):
    for k in bqs[i]:
        check = True
        for j in range(len(cleanbs[i])):
            if phietadist(k,cleanbs[i][j]) < threshhold:
                check = False
                cleanbs[i][j] = k
                break
        if check:
            cleanbs[i].append(k)
            
for i in range(len(vect4)):
    cluster = fj.ClusterSequence(vect4[i], jetdef)
    jets = cluster.inclusive_jets(1)

    bigjets = []
    print()
    print("RECONSTRUCTED")
    for jet in jets:
        if jet.pt() > 1000:
            bigjets.append(jet)
    
    count = 1
    totalpt = 0
    for jet in bigjets:
        print("Jet "+str(count)+": ("+str(jet.pseudorapidity())+", "+str(jet.phi())+")")
        for bot in cleanbs[i]:
            if phietadist(bot, jet) < threshhold:
                print("\tBottom Quark: ("+str(bot.pseudorapidity())+", "+str(bot.phi())+")")
                print("\tJet pT:\t\t"+str(np.sqrt(jet.kt2())))
                print("\tBottom Q pT:\t"+str(np.sqrt(bot.kt2())))
        totalpt += np.sqrt(jet.kt2())
        count += 1
    print("pT: "+str(totalpt))
    print()
''' 
    
    
