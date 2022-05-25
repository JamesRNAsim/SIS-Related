#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import copy
import numpy as np
from Dcd import DcdFile
import math
from turtle import color
import toml

#Acknowledgenments
#This file was created through the work of both Dr. Naoto Hori and James Robins
#The file was created for analysis of the output files from the SIS model of RNA dynamics with a focus on Long Stranded RNA


tomldata = toml.load(sys.argv[1])
#bpein = open('5SR20Run1.bpe', 'r')
#FASTAsequence = open('1C2X_Sequence_DummyNec.txt', 'r')
#dcdin = open('5SR21Run1.dcd', 'r')
bpein = open(tomldata['files']['in']['bpe'], 'r')
FASTAsequence = open(tomldata['files']['in']['sequence'], 'r')
dcdin = tomldata['files']['in']['dcd']
#Files Out
'''CorrectedBPsout = open('5SR20Run1correctedbps.txt', 'w+')
totalBPsout = open('5SR20Run1TotalBPs.txt', 'w+')
ReducedBPsout = open('5SR20Run1ReducedBPs.txt', 'w+')
bpseq = open('5SR20Run1BpSeq.bpseq', 'w')
RGandEEout = open('5SR20Run1RGEE.txt', 'a+')
allbasepairsout = open('5SR20Run1BPM.txt', 'w+')
basepairsprobabilities = open('5SR20Run1BPP.txt', 'w+')
MountValsOut = open('5SR20Run1MountainVals.txt', 'w+')'''
CorrectedBPsout = open(tomldata['files']['out']['CorrectedBPs'], 'w+')
totalBPsout = open(tomldata['files']['out']['totalBPs'], 'w+')
ReducedBPsout = open(tomldata['files']['out']['ReducedBPs'], 'w+')
bpseq = open(tomldata['files']['out']['bpseq'], 'w')
RGandEEout = open(tomldata['files']['out']['RGandEE'], 'a+')
allbasepairsout = open(tomldata['files']['out']['allpairs'], 'w+')
basepairsprobabilities = open(tomldata['files']['out']['pairprobabilities'], 'w+')
MountValsOut = open(tomldata['files']['out']['MountainVals'], 'w+')

#Variables bpseq
Nwant = 15000
NNT = 122

# Correct BPs - Remove any pairs below Cutoff
stepcount = 0 
step = 0

#base Pair Probabilities 
BPPwant = 11000
BPPfinal = 15000
probabilitysteps = BPPfinal - BPPwant

#Correct BPs
for line in bpein:
    #b = len(line)
    #print(b)
    
    lsp = line.split()

    n = 0
    step = step + 1
    if len(lsp) % 3 != 0:
        #print('Error len(lsp) % 3 != 0')
        sys.exit(2)

    nbp = len(lsp) // 3

    if stepcount <1000:
        kbt = -1.987204259e-3 * 423.15
    elif stepcount >=1000 and stepcount <2000:
        kbt = -1.987204259e-3 * 398.15
    elif stepcount >=2000 and stepcount <3000:
        kbt = -1.987204259e-3 * 373.15
    elif stepcount >=3000 and stepcount <4000:
        kbt = -1.987204259e-3 * 353.15
    elif stepcount >=4000 and stepcount <5000:  
        kbt = -1.987204259e-3 * 343.15
    elif stepcount >=5000 and stepcount <6000:  
        kbt = -1.987204259e-3 * 333.15
    elif stepcount >=6000 and stepcount <7000:  
        kbt = -1.987204259e-3 * 323.15
    elif stepcount >=7000 and stepcount <8000:  
        kbt = -1.987204259e-3 * 318.15
    elif stepcount >=8000 and stepcount <9000:  
        kbt = -1.987204259e-3 * 313.15
    elif stepcount >=9000:
        kbt = -1.987204259e-3 * 310.15 

    Cutoff = kbt

    for ipair in range(nbp):
        i = int(lsp[ipair*3])
        j = int(lsp[ipair*3+1])
        e = float(lsp[ipair*3+2])

        if e <= Cutoff:
            n += 1
         # write the bps to file if they are greater than cutoff
            CorrectedBPsout.write(str(i) + ' ' + str(j) + ' ' + str(e) + ' ')
    stepcount = stepcount+1        
    CorrectedBPsout.write('\n')    
    totalBPsout.write(str(step) + ' ' + str(n) + '\n')
#print('Corrected Base Pairs Completed')

# Reduce Multiple Base Pairs
CorrectedBPsout.seek(0)
pairs = []
energies = []

def basepair_frequency(Nnt, pairs):
    frequency = [0] * (Nnt+1)  # Starts from index 1

    for p in pairs:
        frequency[p[0]] += 1
        frequency[p[1]] += 1

    return frequency

def reduce_multiple_basepairs(Nnt, pairs_in, energies_in):
    
    pairs = copy.deepcopy(pairs_in)
    energies = copy.deepcopy(energies_in)

    frequency = basepair_frequency(Nnt, pairs)
    max_f = max(frequency)
        
    while(max_f > 1):

        highest_e = -9999.0
        highest_idx = -1
        for idx, (p, e) in enumerate(zip(pairs, energies)):
            imp1 = p[0]
            imp2 = p[1]
            if frequency[imp1] == max_f or frequency[imp2] == max_f:
                if e > highest_e:
                    highest_e = e
                    highest_idx = idx

        if highest_idx == -1:
            raise Exception("Error: highest_idx = 0")

        del pairs[highest_idx]
        del energies[highest_idx]
        
        frequency = basepair_frequency(Nnt, pairs)
        max_f = max(frequency)
    
    return pairs, energies
f = 0
for line in CorrectedBPsout:
    lsp = line.split()
    f = f + 1
    #print(lsp)
    #for bps
    nbp = len(lsp) // 3
    #print(bonding)
    #print(nbp)
    #print("nbp = ", nbp)
    pairs = []
    for ipair in range(nbp):
        #for bps
        i = int(lsp[ipair*3])
        j = int(lsp[ipair*3+1])
        e = float(lsp[ipair*3+2])
        pairs.append((i, j))
        energies.append(e)
    
    a, v = reduce_multiple_basepairs(NNT, pairs, energies)
    #print(a)
    m = 0
    for x in a:
        i, j = a[m]
        e = v[m]
        m = m + 1
   
        ReducedBPsout.write(str(i) + ' ' + str(j) + ' ' + str(e) + ' ')
     
    ReducedBPsout.write('\n')
    #print("Step", f)
#print('reduced Base Pairs complete')   

#BpSeq
ReducedBPsout.seek(0)
t = 0
ilist = []
jlist = []
for lines in FASTAsequence:
    if lines.startswith('>'):
        continue
    seq = list(lines)

#print(seq)
for line in ReducedBPsout:
    t = t + 1
    if t == Nwant:
        lsp = line.split()
        length = len(lsp)
        nbp = length / 3
        for ipair in range(int(nbp)):
            i = int(lsp[ipair*3])
            j = int(lsp[ipair*3+1]) 
            e = float(lsp[ipair*3+2])
            ilist.append(i)
            jlist.append(j)
#print(ilist)
#print(jlist)
x = len(seq)
y = len(lsp)

#print(x, y)
m = 0
for p in range(1, x):
    #print(p)
    
    if p in ilist:
        index = ilist.index(p)
        bpseq.write(str(ilist[index]) + ' ' + str(seq[p-1]) + ' ' + str(jlist[index]) + '\n')
    if p in jlist:
        index = jlist.index(p)
        bpseq.write(str(jlist[index]) + ' ' + str(seq[p-1]) + ' ' + str(ilist[index]) + '\n')   
    if p not in ilist and p not in jlist:
        bpseq.write(str(p) + ' ' + str(seq[p-1]) + ' ' + str(0) + '\n')


#Radius of Gyration and E-E Distance
if len(sys.argv) != 2:
    print('Usage: SCRIPT [dcd]')
    print('')
    print('Returns: %5i %6.3f %6.3f %6.3f % (N, Rg, D, S)')
    print('    N  = number of particles (atoms)')
    print('    Rg = radius of gyration')
    print('    D  = Sphericity (0 <= D <= 1)')
    print('          D = 0 --> perfect sphere')
    print('          D > 0 --> anisotropic')
    print('    S  = Spheroidal shape (-1/4 <= S <= 2)')
    print('          S < 0 --> oblate')
    print('          S > 0 --> prolate')
    sys.exit(2)

dcd = DcdFile(dcdin)
dcd.open_to_read()
dcd.read_header()

'''
xyzs for N-atom PDB
 [ [ 1x, 1y, 1z],
   [ 2x, 2y, 2z],
   [ 3x, 3y, 3z],
   [ 4x,  .
          .
          .
   [ Nx, Ny, Nz] ]
'''

N = dcd._header.nmp_real
fN = float(N)

def calc_shape(xyzs):
    # Inertia tensor
    T = np.zeros((3,3)) 
    
    for i in range(N):
        for j in range(N):
            for a in range(3):
                for b in range(3):
                    T[a,b] += (xyzs[i][a] - xyzs[j][a]) * (xyzs[i][b] - xyzs[j][b])
    
    T = T / (2.0 * fN**2)
    
    w, v = np.linalg.eig( T )
    
    trT = sum(w)
    #print 'w', w
    #print 'v',v
    #print 'trT = ', trT
    Rg = math.sqrt(trT)
    
    w_avg = trT / 3.0
    
    D = 1.5 * ( (w[0]-w_avg)**2 + (w[1]-w_avg)**2 + (w[2]-w_avg)**2 ) / (trT**2)
    
    S = 27.0 * ( (w[0]-w_avg) * (w[1]-w_avg) * (w[2]-w_avg) ) / (trT**3)
    
    return Rg, D, S

def EndDistance(dcd):
    x1, y1, z1 = dcd[0]
    x2, y2, z2 = dcd[-1]
    distance = np.sqrt(((x1-x2)**2) + ((y1 - y2)**2) + ((z1 - z2)**2))
    return distance

#print("#N: ", N)
z = 0
while dcd.has_more_data():
    #print('step', z)
    z = z + 1
    data = dcd.read_onestep()
    Rg, D, S = calc_shape(data)
    distance = EndDistance(data)

    #print(('%6.3f %6.3f %6.3f %6.3f' % (Rg, D, S, distance)))
    RGandEEout.write(('%6.3f %6.3f %6.3f %6.3f' % (Rg, D, S, distance)))
    RGandEEout.write('\n')
#print('Radius of Gyration and E-E Distance complete')

#Base Pair Contacts 
ReducedBPsout.seek(0)
z = 0
probabilitybasepairs = []
BPchance = []
pairsz = []
totalbonds = 0
for line in ReducedBPsout:
    z = z + 1
    k = 0
    
    if z >= BPPwant and z <=BPPfinal:
        #print('i', i)
        BPList = line.split()
        
        length = len(BPList)
        nbp = length / 3
        for ipair in range(int(nbp)):
            i = int(BPList[ipair*3])
            j = int(BPList[ipair*3+1]) 
            e = float(BPList[ipair*3+2])
        #print('Base pairs', BPList)
            pairsz.append((i, j))
            
#print('file 1 pairs found')
probabilitybasepairs = [ii for n,ii in enumerate(pairsz) if ii not in pairsz[:n]]

lengthz = len(probabilitybasepairs)
mbp = 0
for p in range(int(lengthz)):
    count = pairsz.count(probabilitybasepairs[k])
    #print('k', k)
    probability = count / probabilitysteps
    if probability >= 1:
        print('probability error', probabilitybasepairs[k], totalbonds, probability)
    BPchance.append(probability)
    basepairsprobabilities.write(str(probability) + ' ')
    k = k + 1
for p in range(int(lengthz)):
    i, j = probabilitybasepairs[mbp]
    mbp = mbp + 1
    allbasepairsout.write(str(i) + ' ' + str(j) + ' ')

#scaled Mountain Plot 
ReducedBPsout.seek(0)
mountline = 0
mountpairs = []

for Mount in ReducedBPsout:
    mountline = mountline + 1
    if mountline == Nwant:
        Vallist = Mount.split()

        numofpair = len(Vallist) //3
        #print(Vallist)
        #print('pairs', numofpair)
        for allpair in range(int(numofpair)):
            #print('allpair',allpair)
            i = int(Vallist[allpair*3])
            j = int(Vallist[allpair*3+1]) 
            e = float(Vallist[allpair*3+2])

            mountpairs.append((i, j))
#print('Pairs',mountpairs)
def pairs2Mountain(mountpairs, nucleo, scaled=True):

    m = [0.]*nucleo
    for (nt1, nt2) in mountpairs:
        if nt2 < nt1:
            nt1, nt2 = nt2, nt1
        if scaled:
            l = nt2 - nt1
            m[nt1-1] += 1./float(l)
            m[nt2-1] += -1./float(l)
        else:
            m[nt1-1] += 1
            m[nt2-1] += -1

    s = 0.
    for i in range(nucleo):
        s += m[i]
        m[i] = s

    return m
mountainvals = pairs2Mountain(mountpairs, nucleo=NNT, scaled=True)
numofval = len(mountainvals)
for mvtw in range(numofval):
    write = mountainvals[mvtw]
    MountValsOut.write(str(write) + ' ')

#Graphing Results - Graph 1 is Base Pairs, Graph 2 is Scaled Mountain Plot, Graph 3 is Radius of Gyration and Graph 4 is End to End distance
#The final and seperate figure (fig 2) is the base pair probability matrix
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2) 
fig2 = plt.figure()
ax5 = fig2.add_subplot(111,)

#total BPs
totalBPsout.seek(0)
BPstime = []
BPslist = []
for Bps in totalBPsout:
    Bp = Bps.split()
    #print(lsp)
    Bpa = int(Bp[0])
    Bpb = int(Bp[1])
    BPstime.append(Bpa)
    BPslist.append(Bpb)

ax1.plot(BPstime, BPslist)
ax1.set_title("Base pairs")
ax1.set_xlabel("Step")
ax1.set_ylabel("Base pairs")

#Mountain Plot
MountValsOut.seek(0)
mountplotvals = []
for mounts in MountValsOut:

    valsplit = mounts.split()
    #print(valsplit)
    for valstoadd in range(NNT):
        appendval = valsplit[valstoadd]
        mountplotvals.append(float(appendval))

ax2.plot(range(1, NNT+1), mountplotvals)
ax2.set_title("Mountain Plot with PKs")
ax2.set_xlabel("Nucleotide")
ax2.set_ylabel("Base Pair Value")

#Radius of Gyration
RGandEEout.seek(0)
Rg = []
Rgtime = []
Rgn = 0
for RG in RGandEEout:
    #print('RG', RG)
    if RG.startswith(' '):
        continue
    Rgs = RG.split()
    #print('Rgs', Rgs)
    Rga = float(Rgs[0])
    Rg.append(Rga)
    Rgtime.append(Rgn)
    Rgn = Rgn + 1

ax3.plot(Rgtime, Rg)
ax3.set_title("Radius of Gyration")
ax3.set_xlabel("Step")
ax3.set_ylabel("Radius of Gyration (Å)")

#End to End Distance
RGandEEout.seek(0)
EE = []
EEtime = []
EEn = 0
for EEdis in RGandEEout:
    if EEdis.startswith('frame'):
        continue
    
    EEsplit = EEdis.split()
    EEa = float(EEsplit[3])
    EE.append(EEa)
    EEtime.append(EEn)
    EEn = EEn + 1
ax4.plot(EEtime, EE)
ax4.set_title("End to End Distance")
ax4.set_xlabel("Step")
ax4.set_ylabel("End to End Distance (Å)")

#BP Probability Matrix
allbasepairsout.seek(0)
basepairsprobabilities.seek(0)
pairs = []
x9 = []
y9= []
x8 = []
y8= []
x7 = []
y7= []
x6 = []
y6= []
x5 = []
y5= []
xother = []
yother= []
xtotal = []
ytotal = []
for BPM in allbasepairsout:
    BPMs = BPM.split()
for BPP in basepairsprobabilities:
    BPPs = BPP.split()
#print('bps', bp)
#print('probabilities', pro)

if len(BPPs) != (len(BPMs)/2):
    print('error - Probabilities and Base Pairs do not match', len(BPPs), len(BPMs))
    
long = len(BPMs)
NofN = long / 2
for BPMBPP in range(int(NofN)):
    i = int(BPMs[BPMBPP*2])
    j = int(BPMs[BPMBPP*2+1]) 
#print('Base pairs', BPList)
    pairs.append((i, j))

for npair in range(0, len(BPPs)):
    if float(BPPs[npair]) >= 0.8:
        i, j = pairs[npair]
       
        x9.append(i)
        y9.append(j)
    elif float(BPPs[npair]) >= 0.6 and float(BPPs[npair]) <= 0.8:
        i, j = pairs[npair]
        
        x8.append(i)
        y8.append(j)
    elif float(BPPs[npair]) >= 0.4 and float(BPPs[npair]) <= 0.6:
        i, j = pairs[npair]
        
        x7.append(i)
        y7.append(j)
    elif float(BPPs[npair]) >= 0.2 and float(BPPs[npair]) <= 0.4:
        i, j = pairs[npair]
        
        x6.append(i)
        y6.append(j)
    elif float(BPPs[npair]) >= 0.05 and float(BPPs[npair]) <= 0.2:
        i, j = pairs[npair]
        
        x5.append(i)
        y5.append(j)

for p in range(int(NofN)):
    i, j = pairs[p] 
    xtotal.append(i)
    ytotal.append(j)

ax5.scatter(ytotal, xtotal, c='b', marker = '.')
ax5.scatter(x5, y5, c="c", marker = '.')
ax5.scatter(x6, y6, c="g", marker = '.')
ax5.scatter(x7, y7, c="m", marker = '.')
ax5.scatter(x8, y8, c="y", marker = '.')
ax5.scatter(x9, y9, c="r", marker = '.')

ax5.set_aspect('equal')
ax5.set_title("Base Pair matrix")
ax5.set_xlabel("Base i")
ax5.set_ylabel("Base j")

#Graph Saves
fig.tight_layout()
fig.savefig('4MainGraphs.png')
fig2.tight_layout()
fig2.savefig('BasePairMatrix.png')


#Close Files
bpein.close()
FASTAsequence.close()
CorrectedBPsout.close()
totalBPsout.close()
ReducedBPsout.close()
bpseq.close()
RGandEEout.close()
basepairsprobabilities.close()
allbasepairsout.close()
MountValsOut.close()
#print('Analysis Complete')