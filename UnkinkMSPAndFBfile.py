
# coding: utf-8

#flip catching/correcting script that detects strand flips and corrects them back in RFmix output
#this version fixes the RFmix output, .msp file with LAI calls for 2-way admixed individuals and outputs a file documenting the switch locations for the RFmix LAI call windows and whether that window was switched or not in each individual
__author__ = 'egatkinson, modified by ThiagoPL'

import argparse
import re
import numpy as np

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--msp', help='path stem to RFmix msp file, not including .msp.tsv', required=True)
  parser.add_argument('--num', help='number of parental', required=True)
  args = parser.parse_args()
  return(args)    


args = parse_args()
numAncestry = int(args.num)
mspfile = open(args.msp + '.msp.tsv')

mspfile.readline()
header = mspfile.readline()
chrom, startbp, stopbp, startcm, stopcm, snps, samples = header.split('\t', 6)
state = {}
switched = {}

for i in range(int(len(samples)/2)):
    state[i] = (float(0), [-1,-1])
    switched[i] = False

out = open(args.msp + '.Unkinked.msp.tsv', 'w') #output will be input filename suffixed with "Unkinked"
out1 = open(args.msp + '.switches.txt', 'w') #also save the points where the switches happen. Export all windows, with whether switched = true or false.
out.write(header) #currently only writing the second header line with the individual information

for line in mspfile:
    chrom, startbp, stopbp, startcm, stopcm, snps, calls = line.strip().split('\t', 6)
    calls = calls.split('\t')
    output = '\t'.join([chrom, startbp, stopbp, startcm, stopcm, snps])
    outswitches = '\t'.join([chrom, startbp, stopbp, startcm, stopcm, snps])
    for i in range(int(len(calls)/2)):
        hapA = calls[2*i]
        hapB = calls[2*i + 1] 
        newhapA = hapA
        newhapB = hapB
        
        if hapA != hapB:            
            if [hapA,hapB] != state[i][1] and (float(startcm) - 1) < state[i][0]:
                switched[i] = not switched[i]
                
            if switched[i]:
                newhapA = hapB
                newhapB = hapA
           
            state[i] = (float(stopcm), [hapA,hapB])
        output += '\t' + newhapA + '\t' + newhapB
        outswitches += '\t' + str(switched[i])
    output += '\n'
    outswitches += '\n'
    out.write(output)
    out1.write(outswitches)
out.close()
out1.close()


#Modification to Unkink the FB file to be used by KART
fbfile = open(args.msp + '.fb.tsv')
out = open(args.msp + '.Unkinked.fb.tsv', 'w')
switches = open(args.msp + '.switches.txt')

#keeping both headers on fb file
header = fbfile.readline()
out.write(header)
header = fbfile.readline()
out.write(header)

for switchLine in switches:
    fbSplit = fbfile.readline().strip().split()
    fbSliced = fbSplit[4:]
    switchSplit = switchLine.strip().split()
    switchSliced = switchSplit[6:]

    out.write(f'{fbSplit[0]}\t{fbSplit[1]}\t{fbSplit[2]}\t{fbSplit[3]}')
    for i in range(0, len(switchSliced)):
        base = (i) * (numAncestry) * 2
        if switchSliced[i] == "True":
            beginHap1 = base+numAncestry
            endHap1 = base+numAncestry+numAncestry
            beginHap2 = base
            endHap2 = base + numAncestry

            for j in range(beginHap1, endHap1):
                out.write(f'\t{fbSliced[j]}')
            for j in range(beginHap2, endHap2):
                out.write(f'\t{fbSliced[j]}')
        else:
            for j in range(base, base+(2*numAncestry)):
                out.write(f'\t{fbSliced[j]}')

    out.write("\n")
out.close()