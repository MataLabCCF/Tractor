# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 15:53:35 2021

@author: PEIXOTT
"""

import os
import sys
import argparse

def getchar():
    c = sys.stdin.read(1)

parser = argparse.ArgumentParser(description='Create individuals set')

required = parser.add_argument_group("Required arguments")
required.add_argument('-f', '--fam', help='Fam file from plink with all individuals', required=True)
required.add_argument('-n', '--num', help='Number of individuals per set',
                          required=True)
required.add_argument('-c', '--correspondence', help='Correspondence between parental and IDs', required=True)

args = parser.parse_args()

fam = open(args.fam)
correspondence = open(args.correspondence)

dictRef = []
for ind in correspondence:
    split = ind.split()
    dictRef.append(split[0])

file = 0
ref = ''
notRef = []
for ind in fam:
    split = ind.split()
    if split[1] not in dictRef:
        notRef.append(split[1])
    else:
        ref = f'{ref}{split[1]}\n'

file = 0
num = int(args.num)
for i in range(0,len(notRef), num):
    extract = open(f'split{file}', 'w')
    for j in range(num):
        if i+j < len(notRef):
            extract.write(f'{notRef[i+j]}\n')
    extract.write(ref)
    extract.close()
    file = file + 1
