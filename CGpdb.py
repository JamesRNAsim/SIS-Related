import sys
import os

pdb = open('2ZJRH.pdb', 'r')
output = open('intermediate.txt', 'w+')
nospaces = open('cg_2ZJRH.pdb', 'w')
x = 0.0
y = 0.0
z = 0.0
nt_pre = None
n = 0
abases = ['A', '1MA', 'A2M', 'AET', 'P5P', 'SAM']
cbases = ['C','5MC', 'OMC', 'CCC', 'CBV']
ubases = ['U', 'PSU', '5MU', 'H2U', '4SU', '5BU', 'OMU', 'FHU', '5BU', 'UD5','UR3', 'DU']
gbases = ['G', 'G7M', 'QUO', '2MG', 'OMG', 'YYG', '7MG', 'M2G', 'YG', 'GTP']
unknown = ['N']
prevline=''
index = 0
for line in pdb:  
    if line.startswith('ATOM') or line.startswith('HETATM'):
        Atom = line[12:16]
        nt = int(line[22:26])
        nucleotide = line[17:21]
       
        base = ''.join(nucleotide.strip())
       
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

                x = 0.0
                y = 0.0
                z = 0.0
                n = 0

                if prevbase in abases:                  
                    output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:>1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
                        int(prevline[6:11]),str(line[12:16]),str(prevline[16:17]),str('A'),str(prevline[21:22]),int(prevline[22:26]),str(prevline[26:27]),
                        float(avgx),float(avgy),float(avgz),float(prevline[54:60]),float(prevline[60:66]),str(prevline[76:78]),str(prevline[78:80]))+'\n')
                    
                elif prevbase in ubases:          
                    output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
                        int(prevline[6:11]),str(line[12:16]),str(prevline[16:17]),str('U'),str(prevline[21:22]),int(prevline[22:26]),str(prevline[26:27]),
                        float(avgx),float(avgy),float(avgz),float(prevline[54:60]),float(prevline[60:66]),str(prevline[76:78]),str(prevline[78:80]))+'\n')
                    
                elif prevbase in cbases:                    
                    output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
                        int(prevline[6:11]), str(line[12:16]), str(prevline[16:17]),str('C'), str(prevline[21:22]), int(prevline[22:26]), str(prevline[26:27]),
                        float(avgx),float(avgy), float(avgz), float(prevline[54:60]), float(prevline[60:66]), str(prevline[76:78]), str(prevline[78:80]))+'\n')
                    
                elif prevbase in gbases:
                    output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
                        int(prevline[6:11]), str(line[12:16]), str(prevline[16:17]),str('G'), str(prevline[21:22]), int(prevline[22:26]), str(prevline[26:27]),
                        float(avgx),float(avgy), float(avgz), float(prevline[54:60]), float(prevline[60:66]), str(prevline[76:78]), str(prevline[78:80]))+'\n')
                   
                elif prevbase in unknown:
                    output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
                        int(prevline[6:11]), str(line[12:16]), str(prevline[16:17]),str('N'), str(prevline[21:22]), int(prevline[22:26]), str(prevline[26:27]),
                        float(avgx),float(avgy), float(avgz), float(prevline[54:60]), float(prevline[60:66]), str(prevline[76:78]), str(prevline[78:80])))
                else:
                    print("error - Base = ", prevbase) 

            nt_pre = nt
            prevbase = base 
            prevline = line

    if line.startswith('TER'):
        print('Finished Reading Chain'+' '+line[21])
        break
    #if line.startswith('END'):
    #    break


if n != 0:
    avgx = x/float(n)
    avgy = y/float(n)
    avgz = z/float(n)
    
    if prevbase in abases:                  
        output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:>1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
            int(prevline[6:11]),str(line[12:16]),str(prevline[16:17]),str('A'),str(prevline[21:22]),int(prevline[22:26]),str(prevline[26:27]),
            float(avgx),float(avgy),float(avgz),float(prevline[54:60]),float(prevline[60:66]),str(prevline[76:78]),str(prevline[78:80]))+'\n')
        
    elif prevbase in ubases: 
                      
        output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
            int(prevline[6:11]),str(line[12:16]),str(prevline[16:17]),str('U'),str(prevline[21:22]),int(prevline[22:26]),str(prevline[26:27]),
            float(avgx),float(avgy),float(avgz),float(prevline[54:60]),float(prevline[60:66]),str(prevline[76:78]),str(prevline[78:80]))+'\n')
        
    elif prevbase in cbases:                    
        output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
            int(prevline[6:11]), str(line[12:16]), str(prevline[16:17]),str('C'), str(prevline[21:22]), int(prevline[22:26]), str(prevline[26:27]),
            float(avgx),float(avgy), float(avgz), float(prevline[54:60]), float(prevline[60:66]), str(prevline[76:78]), str(prevline[78:80]))+'\n')
        
    elif prevbase in gbases:
        output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
            int(prevline[6:11]), str(line[12:16]), str(prevline[16:17]),str('G'), str(prevline[21:22]), int(prevline[22:26]), str(prevline[26:27]),
            float(avgx),float(avgy), float(avgz), float(prevline[54:60]), float(prevline[60:66]), str(prevline[76:78]), str(prevline[78:80]))+'\n')
        
    elif prevbase in unknown:
        output.write("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(str(prevline[0:6]),
            int(prevline[6:11]), str(line[12:16]), str(prevline[16:17]),str('N'), str(prevline[21:22]), int(prevline[22:26]), str(prevline[26:27]),
            float(avgx),float(avgy), float(avgz), float(prevline[54:60]), float(prevline[60:66]), str(prevline[76:78]), str(prevline[78:80])))
    else:
        print("error - Base = ", prevbase) 
        
output.seek(0)
for lines in output:
    if lines.startswith('ATOM') or lines.startswith('HETATM'):
        nospaces.write(lines)
    else:
        continue