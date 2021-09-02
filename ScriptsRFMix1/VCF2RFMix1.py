# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:41:29 2021

@author: PEIXOTT
"""
import os
import sys
import vcfpy
import argparse

def getchar():
    c = sys.stdin.read(1)


def createBin(vcf, genmap, output, chrom):
    if genmap[-3:] == '.gz':
        command = f'cp {genmap} ./GenMap{chrom}.gz | gunzip ./GenMap{chrom}.gz'
        execute(command)
        genmap = f'./GenMap{chrom}'
    command = f'plink --vcf {vcf} --cm-map {genmap} {chrom} --make-bed --out {output}_gen --double-id'
    execute(command)
    return f'{output}_gen.bim'

def execute(command):
    print(command + "\n\n")
    os.system(command)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='VCF to RFMix1')

    required = parser.add_argument_group("Required arguments")
    # required.add_argument('-s', '--steps', help='File with the order to perform the analisys', required=True)
    required.add_argument('-v', '--vcf', help='VCF file to be converted',
                          required=True)
    required.add_argument('-c', '--correspondence', help='Correspondence between parental and IDs', required=True)
    required.add_argument('-C', '--chromosome', help='It is self-explanatory ', required=True)
    required.add_argument('-o', '--output', help='Output name file', required=True)
    required.add_argument('-m', '--map', help='Genetic map', required=False)

    args = parser.parse_args()



    
    indClass = {}
    mapDict = {}
    popDict = {}
    numParental = 0
    
    correspondence = open(f'{args.correspondence}')
    
    for line in correspondence:
        split = line.strip().split()
        if split[1] not in popDict:
            numParental = numParental+1
            popDict[split[1]] = numParental
        indClass[split[0]] = popDict[split[1]]
    classPrint = False
    
    bim = createBin(args.vcf, args.map, args.output, args.chromosome)
    
    print(f'Abrindo o {bim}')
    mapFile = open(f'{bim}')
    header = False
    for line in mapFile:
        if header:
            header = False
        else:
            split = line.strip().split()
            if(float(split[2]) < 0):
                split[2] = 0.0
            mapDict[int(split[3])] = float(split[2])
    
    correspondence.close()
    mapFile.close()
    
    count = 0
    firstInd = True
    reader = vcfpy.Reader.from_path(f'{args.vcf}')
    classes = open(f'{args.output}_classes', 'w')
    alleles = open(f'{args.output}_alleles', 'w')
    location = open(f'{args.output}_location', 'w')
    for record in reader:
        location.write(f'{mapDict[record.POS]}\n')
        callList = record.call_for_sample
        for ind in callList:
            if not classPrint:
                if ind not in indClass:
                    if firstInd:
                        firstInd = False
                        classes.write('0 0')
                    else:
                        classes.write(' 0 0')
                else:
                    if firstInd:
                        firstInd = False
                        classes.write(f'{indClass[ind]} {indClass[ind]}')
                    else:
                        classes.write(f' {indClass[ind]} {indClass[ind]}')
            
            GT = callList[ind].data.get('GT').split('|')
            alleles.write(f'{GT[0]}{GT[1]}')
        
        if not classPrint:
            classPrint = True
            classes.close()
        alleles.write('\n')
    print(f'{count} unmapped (0 cM)')
    alleles.close()
    location.close()
    
    
