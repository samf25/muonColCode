import ROOT
import sys
import fastjet as fj
import numpy as np

# Get Data
file = ROOT.TFile.Open(sys.argv[1], "READ")
tree = file.Get("MyLCTuple;42")

# Storage Arrays
remMom = [[], [], [], []]
sumMass = []
masses = []
betas = []

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
    count += 1
    higgs = None
    b1 = None
    b2 = None
    for i in range(len(event.mcpdg)):
        if event.mcpdg[i] == 25:
            higgs = fj.PseudoJet(event.mcmox[i],event.mcmoy[i], event.mcmoz[i], event.mcene[i])
        elif event.mcpdg[i] == 5 and b1 == None:
            b1 = fj.PseudoJet(event.mcmox[i],event.mcmoy[i], event.mcmoz[i], event.mcene[i])
        elif event.mcpdg[i] == -5 and b2 == None:
            b2 = fj.PseudoJet(event.mcmox[i],event.mcmoy[i], event.mcmoz[i], event.mcene[i])
        if higgs != None and b1 != None and b2 != None:
            break
    if higgs == None or b1 == None or b2 == None:
        print("EVENT "+str(count)+ " is Missing an element")
        print("\t"+str(len(event.mcpdg)))
    else:
        masses.append(b1.m())
        masses.append(b2.m())
        betas.append(abs(b1.pseudorapidity()))
        betas.append(abs(b2.pseudorapidity()))
        sumMass.append(b1.m()+b2.m())
        temp = higgs - b1
        temp = temp - b2
        remMom[0].append(temp.E())
        remMom[1].append(temp.px())
        remMom[2].append(temp.py())
        remMom[3].append(temp.pz())

        if count % 5 == 0:
            print("EVENT "+str(count), end="\r")

histgraph(remMom[0], "Remaining Energy", "Graphs/RemEne.pdf")
histgraph(remMom[1], "Remaining px", "Graphs/RemPx.pdf")
histgraph(remMom[2], "Remaining py", "Graphs/RemPy.pdf")
histgraph(remMom[3], "Remaining pz", "Graphs/RemPz.pdf")

histgraph(sumMass, "b bbar Mass sum", "Graphs/bMasses.pdf")
histgraph(masses, "b masses", "Graphs/bMassIndv.pdf")

betas = np.array(betas)
print("Fraction of b etas btwn -1.6 and 1.6: " + str(len(betas[betas<1.6])/len(betas)))
print("Fraction of b etas btwn -2 and 2: " + str(len(betas[betas<2])/len(betas)))

