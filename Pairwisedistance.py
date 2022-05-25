
import numpy as np
import sys
model = '2ZJRH'
structure = open('cg_2ZJRH.pdb', 'r')
pairwise = open('pairwiseout.out', 'a')
t = 0
for lines in structure:
    t = t+1
    if t == 1:
        firstnucleotide = int(lines[22:26].strip())
    finalnucleotide = int(lines[22:26].strip())
#print(firstnucleotide, finalnucleotide)
#pairwise.write('Model'+'    '+'First_Bond_Length?'+'    '+'i1'+'    '+'i2'+'    '+'Distance (Angstrom)'+'\n')

for i in range(firstnucleotide, finalnucleotide+1):
    #if i == 894:
    #
    #    print(i)
    #    sys.exit
    structure.seek(0)
    for line in structure:
        if int(line[22:26].strip()) == i:            
            testline = line
            atomx = float(testline[30:38].strip())
            atomy = float(testline[38:46].strip())
            atomz = float(testline[46:54].strip())
            previd = int(line[22:26].strip())
        else:
            continue
    if testline ==' ':
        continue
    
    structure.seek(0)
    t = 2
    for line in structure:   

        if int(line[22:26].strip()) == i+1:
            atom2x = float(line[30:38].strip())
            atom2y = float(line[38:46].strip())
            atom2z = float(line[46:54].strip())
            
            firstdistance = np.sqrt(((atomx-atom2x)**2) + ((atomy - atom2y)**2) + ((atomz - atom2z)**2))
            pairwise.write(model +'    '+ 'True' +'    '+ str(testline[22:26].strip())+'    '+str(i+1)+'    '+str(firstdistance)+'\n')
        
        if int(line[22:26].strip()) >= i+3:
            t = t+1
            atomNx = float(line[30:38].strip())
            atomNy = float(line[38:46].strip())
            atomNz = float(line[46:54].strip())

            pairwisedist = np.sqrt(((atomx-atomNx)**2) + ((atomy - atomNy)**2) + ((atomz - atomNz)**2))
            pairwise.write(model +'    '+ 'False' +'    '+ str(testline[22:26].strip())+'    '+str(line[22:26].strip())+'    '+str(pairwisedist)+'\n')
    testline = ' '