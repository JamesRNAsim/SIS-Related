#!/usr/bin/env python3
import sys

pdb = open('1c2x.pdb', 'r')
output = open('cg_1c2x.xyz', 'w')
x = 0.0
y = 0.0
z = 0.0
nt_pre = None
n = 0
abases = ['A', '1MA', 'A2M', 'AET', 'P5P', 'SAM']
cbases = ['C','5MC', 'OMC', 'CCC']
ubases = ['U', 'PSU', '5MU', 'H2U', '4SU', '5BU', 'OMU', 'FHU', '5BU', 'UD5',]
gbases = ['G', 'G7M', 'QUO', '2MG', 'OMG', 'YYG', '7MG', 'M2G', 'YG', 'GTP']
unknown = ['N']
prevline=''
output.write('\n')
for line in pdb:
    
    if line.startswith('ATOM') or line.startswith('HETATM'):

        Atom = line[12:16]
        nt = int(line[22:26])
        nucleotide = line[17:21]
        #print('Nucleotide =', nucleotide)
        base = ''.join(nucleotide.strip())
        #print('base = ', base)

        if line[15] == "'":
            if Atom.strip()[0] == 'H':
                continue        
            #print(Atom)
            x += float(line[30:38])
            y += float(line[38:46])
            z += float(line[46:54])
            n += 1
    
        
        if nt != nt_pre:
            
            if nt_pre is not None:
                avgx = x/float(n)
                avgy = y/float(n)
                avgz = z/float(n)

                avgx = '%.3f' %avgx
                avgy = '%.3f' %avgy
                avgz = '%.3f' %avgz           
                
                #print("{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
                x = 0.0
                y = 0.0
                z = 0.0
                n = 0

                
                if prevbase in abases:
                    output.write('\n' + 'A'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
                elif prevbase in ubases:
                    output.write('\n' + 'U'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
                elif prevbase in cbases:
                    output.write('\n' + 'C'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
                elif prevbase in gbases:
                    output.write('\n' + 'G'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
                elif prevbase in unknown:
                    output.write('\n' + 'N'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
                else:
                    print("error - Base = ", prevbase) 

            nt_pre = nt
            prevbase = base 

    if line.startswith('TER'):
        print('Finished Reading Chain'+' '+ line[21])
        
        #continue
    #if ''.join(line[17:21]).strip() =='HOH':
    #    print('WATER FOUND CHECK FILE')
    #    break
    
        
if n != 0:
    avgx = x/float(n)
    avgy = y/float(n)
    avgz = z/float(n)

    avgx = '%.3f' %avgx
    avgy = '%.3f' %avgy
    avgz = '%.3f' %avgz  

    #print("{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
    if prevbase in abases:
        output.write('\n' + 'A'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
    elif prevbase in ubases:
        output.write('\n' + 'U'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
    elif prevbase in cbases:
        output.write('\n' + 'C'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
    elif prevbase in gbases:
        output.write('\n' + 'G'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
    elif prevbase in unknown:
        output.write('\n' + 'N'+"{: >10} {: >10} {: >10}".format(avgx, avgy, avgz))
    else:
        print("error - Base = ", prevbase) 
