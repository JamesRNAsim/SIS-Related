#!/usr/bin/env python3

import numpy as np
import sys
import math
bonds = open('437dx3DNABonds.txt', 'r')
coordinates = open('cg_437d.xyz', 'r')
molecule = '437d'
dihedrals = open('debugdihedral.out', 'w')
#debug = open('debug.txt', 'w')
dihedrals.write('File'+'    '+'Pair'+'    '+ 'Base1'+'    '+'Base2'+'    '+'Distance(Angstrom)'+'    '+'i_pairing'+'    '+'j_pairing'+'    '+'Phi1(degrees)'+'    '+'Phi2(degrees)')
dihedrals.write('\n')
Modifiedbases = ['M2G','1MA','7MG','5MC','SAM','YYG','OMG','P5P','OMC','2MG','QUO','4SU','H2U','5MU','G7M','PSU','5BU','A2M','CCC','OMU','FHU','AET','5BU','GTP', 'UD5']
Modifedbases2 = ['YG']
base1 = []
base2 = []
allbases1 = []
allbases2 = []

for line in bonds:
    #print(line)
    if line.startswith('List of'):
        continue
    elif 'nt1' in line:
        continue
    
    #print(data)
    #pairs = data
    #print(data)
   
    if line[39:51].strip() == 'WC' or line[39:51].strip() =='Wobble':
        if all("".join(line[7:10]) not in m for m in Modifiedbases) and all("".join(line[7:10]) not in m for m in Modifedbases2):
            base1filter = filter(str.isdigit, line[7:17])
            basejoin = "".join(base1filter)
            firstbase = basejoin
            x = int(firstbase)
        elif all("".join(line[7:9]) not in m for m in Modifedbases2):
            base1filter = filter(str.isdigit, line[10:17])
            basejoin = "".join(base1filter)
            firstbase = basejoin
            x = int(firstbase)
        else:
            base1filter = filter(str.isdigit, line[9:17])
            basejoin = "".join(base1filter)
            firstbase = basejoin
            x = int(firstbase)
        if x in base1:
            continue
        else:
            base1.append(x)
            allbases1.append(x)
        
        #print("".join(line[22:25]))
        if all("".join(line[22:25]) not in m for m in Modifiedbases) and all("".join(line[22:25]) not in m for m in Modifedbases2):
            base2filter = filter(str.isdigit, line[22:31])
            base2join = "".join(base2filter)
            #print(base2join)
            secondbase = base2join
            y = int(secondbase)
        elif all("".join(line[22:25]) not in m for m in Modifedbases2):
            base2filter = filter(str.isdigit, line[25:31])
            base2join = "".join(base2filter)
            secondbase = base2join
            y = int(secondbase)
        else:
            base2filter = filter(str.isdigit, line[24:31])
            base2join = "".join(base2filter)
            secondbasebase = base2join
            y = int(secondbase)
        if y in base2:
            continue
        else:
            base2.append(y) 
            allbases2.append(y)  
        
    elif line[39:51].strip() !='Wobble' and line[39:51].strip() != 'WC':
        if all("".join(line[7:10]) not in m for m in Modifiedbases) and all("".join(line[7:10]) not in m for m in Modifedbases2):
            base1filter = filter(str.isdigit, line[7:17])
            basejoin = "".join(base1filter)
            firstbase = basejoin
            x = int(firstbase)
        elif all("".join(line[7:9]) not in m for m in Modifedbases2):
            base1filter = filter(str.isdigit, line[10:17])
            basejoin = "".join(base1filter)
            firstbase = basejoin
            x = int(firstbase)
        else:
            base1filter = filter(str.isdigit, line[9:17])
            basejoin = "".join(base1filter)
            firstbase = basejoin
            x = int(firstbase)
        allbases1.append(x)
        if all("".join(line[22:25]) not in m for m in Modifiedbases) and all("".join(line[22:25]) not in m for m in Modifedbases2):
            base2filter = filter(str.isdigit, line[22:31])
            base2join = "".join(base2filter)
            secondbase = base2join
            y = int(secondbase)
        elif all("".join(line[22:25]) not in m for m in Modifedbases2):
            base2filter = filter(str.isdigit, line[25:31])
            base2join = "".join(base2filter)
            secondbase = base2join
            y = int(secondbase)
        else:
            base2filter = filter(str.isdigit, line[24:31])
            base2join = "".join(base2filter)
            secondbasebase = base2join
            y = int(secondbase)
        allbases2.append(y)

    if len(base1) != len(base2):
        print(base1, base2)
        print(allbases1, allbases2)
        print('ERROR - Base Pairs do not match')     
        sys.exit()       
    else:
        continue

#print('First bases', base1,)
#print('second bases', base2)
#print('All First', allbases1)
#print('All second', allbases2)
#print(len(base1), len(allbases1))

base1normal = []
base2normal = []
allbases1normal = []
allbases2normal = []

for p in base1:
    q = p - 1
    base1normal.append(q)

for p2 in base2:
    q2 = p2 - 1
    base2normal.append(q2)

for p3 in allbases1:
    #print(p3, allbases1[0])
    q3 = p3 - 1
    allbases1normal.append(q3)

for p4 in allbases2:
    q4 = p4 - 1
    allbases2normal.append(q4)

if len(base2normal) != len(base1normal):
    print('Normalisation Failed')
    sys.exit()

#print('base1normal = ', base1normal)
#print('base2normal = ', base2normal)


for pair in base1normal:
    if pair == base1normal[0]:
        continue
    index = base1normal.index(pair)
    #print('Pair = ',pair, pair+3)
    partner = base2normal[index]
    #print('Partner = ', partner, partner+3)
    
    
    #firstdihedral = i-1, i, j, j-1
    t = 0
    coordinates.seek(0)
    for lines in coordinates:
        #print('lines', lines)
        t = t+1
        if t == pair+2:
            coord = lines.split()
            i1x = float(coord[1])
            #print(i1x)
            i1y = float(coord[2])
            i1z = float(coord[3])
        else:
            continue
        i1 = np.array([i1x, i1y, i1z])
        
    t = 0
    coordinates.seek(0)
    for lines in coordinates:
        t = t+1
        
        if t == pair+3:
            #print(pair, pair+3)
            #print('t', t, 'line = ', lines)
            #print('line picked')
            #print('pair+3 = ', pair+3)
            coord = lines.split()
           # print('coord', coord)
            j1x = float(coord[1])
            #print(i1x)
            j1y = float(coord[2])
            j1z = float(coord[3])
            nucleoside1 = coord[0]
            #print('base = ', coord[0])
            
        else:
            continue
        j1 = np.array([j1x, j1y, j1z])
        
    t = 0
    coordinates.seek(0)
    for lines in coordinates:
        t = t+1
        if t == partner+3:
            #print(partner, partner+3)
            #print('t', t, 'line = ', lines)
            
            coord = lines.split()
            #print('coord', coord)
            #print('coord', coord)
            k1x = float(coord[1])
            #print(i1x)
            k1y = float(coord[2])
            k1z = float(coord[3])
            nucleoside2 = coord[0]
            #print('base = ', coord[0])
        else:
            continue
        k1 = np.array([k1x, k1y, k1z])
       
    t = 0
    coordinates.seek(0)
    for lines in coordinates:
        t = t+1
        if t == partner+2:
            coord = lines.split()
            #print('coord', coord)
            l1x = float(coord[1])
            #print(i1x)
            l1y = float(coord[2])
            l1z = float(coord[3])
        else:
            continue
        l1 = np.array([l1x, l1y, l1z])
        
    
    #Dihedral Calculation
    F = i1 - j1
    G = j1 - k1
    H = l1 - k1
    A = np.cross(F, G)
    B = np.cross(H, G)
    x = np.dot(A, B)   
    y = np.linalg.norm(A)
    z = np.linalg.norm(B)
    angle = np.arccos((x)/(y*z))
    Q = np.cross(B, A)
    #print(Q)
    S = np.dot(Q, G) 
    #print('sign = ', S)
    distance = np.sqrt(((j1x-k1x)**2) + ((j1y - k1y)**2) + ((j1z - k1z)**2))
    dis = round(distance, 5)
    if dis >= 14.1 or dis <=12.9:
        print('WANRING - Distance = '+str(dis)+' - This may be due to an error, please check carefully')
    dihedrals.write(molecule + '    '+nucleoside1+nucleoside2+'   '+str(pair+1) + '    '+ str(partner+1) + '    '+str(dis)+'    ') 
    if pair-1 not in allbases1normal and pair+1 not in allbases1normal:
        dihedrals.write('i+-1_unpaired'+'    ')     
    elif pair-1 not in allbases1normal and pair+1 in allbases1normal:
        dihedrals.write('i-1_unpaired'+'    ')     
    elif pair+1 not in allbases1normal and pair-1 in allbases1normal:
        dihedrals.write('i+1_unpaired'+'    ')       
    else:
        dihedrals.write('i+-1_paired'+'    ')         
    if partner-1 not in allbases2normal and partner+1 not in allbases2normal:
        dihedrals.write('j+-1_unpaired'+'    ')      
    elif partner-1 not in allbases2normal and partner+1 in allbases2normal:
        dihedrals.write('j-1_unpaired'+'    ')     
    elif partner+1 not in allbases2normal and partner-1 in allbases2normal:
        dihedrals.write('j+1_unpaired'+'    ')
    else:
        dihedrals.write('j+-1_paired'+'    ')
       

    theta = angle*(180/np.pi)
    dihedral = round(theta, 5)
    if S <=0: 
        dihedrals.write('-'+str(dihedral))
    elif S >=0:
        dihedrals.write(str(dihedral))
    dihedrals.write('    ')
    #seconddhiedral = i+1, i, j, j+1
    t = 0
    coordinates.seek(0)
    for lines in coordinates:
        t = t+1
        #print('t', t)
        #print('pair+1', pair+1)
        if t == pair+4:
            coord = lines.split()
            #print('coord', coord)
            i2x = float(coord[1])
            #print(i1x)
            i2y = float(coord[2])
            i2z = float(coord[3])
        else:
            continue
        i2 = np.array([i2x, i2y, i2z])
    t = 0
    coordinates.seek(0)
    for lines in coordinates:
        t = t+1
        #print('t', t)
        #print('pair+1', pair+1)
        if t == pair+3:
            coord = lines.split()
            #print('coord', coord)
            j2x = float(coord[1])
            #print(i1x)
            j2y = float(coord[2])
            j2z = float(coord[3])
        else:
            continue
        j2 = np.array([j2x, j2y, j2z])
    t = 0
    coordinates.seek(0)
    for lines in coordinates:
        t = t+1
        if t == partner+3:
            coord = lines.split()
            k2x = float(coord[1])
            k2y = float(coord[2])
            k2z = float(coord[3])
        else:
            continue
        k2 = np.array([k2x, k2y, k2z])
    
    t = 0
    coordinates.seek(0)
    for lines in coordinates:
        t = t+1
        
        if t == partner+4:
            coord = lines.split()
            l2x = float(coord[1])
            l2y = float(coord[2])
            l2z = float(coord[3])
        else:
            continue
        l2 = np.array([l2x, l2y, l2z])
    
    F2 = i2 - j2
    G2 = j2 - k2
    H2 = l2 - k2
    A2 = np.cross(F2, G2)
    B2 = np.cross(H2, G2)
    x2 = np.dot(A2, B2)   
    y2 = np.linalg.norm(A2)
    z2 = np.linalg.norm(B2)
    angle2 = np.arccos((x2)/(y2*z2))
    Q2 = np.cross(B2, A2)
    #print(Q)
    W2 = np.dot(Q2, G2)
    #print(W)
    
    theta2 = angle2*(180/np.pi)
    dihedral2 = round(theta2,5)
    if W2 <=0: 
        dihedrals.write('-'+str(dihedral2))
    elif W2 >=0:
        dihedrals.write(str(dihedral2))
    dihedrals.write('\n')

    #print('dihedral 1 pair = ', pair, 'atoms = ', i1, )
    #print('dihedral 2 pair = ', pair, 'Atoms = ', i2, j2, k2, l2)
    