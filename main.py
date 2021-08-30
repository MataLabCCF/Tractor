import sys
import gzip
import os
import time
import argparse
from multiprocessing import Pool


def getchar():
    c = sys.stdin.read(1)


class MSP:
    def __init__(self):
        self.windowStart = []
        self.windowEnd = []

        self.chr = 0

        self.individuals = {}

    def addIndividual(self, id):
        if id in self.individuals:
            finish(f"[Error in MSP File]: The ID ({id}) is duplicate")
        self.individuals[id] = mosaics()

    def addAncestryHaplotype(self, id, hapID, anc):
        self.individuals[id].addAnc(hapID, anc)

    def addInfoAboutWindow(self, chr, startPos, endPos):
        self.chr = chr
        self.windowStart.append(int(startPos))
        self.windowEnd.append(int(endPos))

    def getAncPos(self, ind, pos):
        pos = int(pos)
        if pos <= self.windowStart[0]:
            return self.individuals[ind].getAncInWindow(0)
        elif pos >= self.windowStart[-1]:
            return self.individuals[ind].getAncInWindow(-1)
        else:
            for i in range(len(self.windowStart)):
                if pos >= self.windowStart[i]:
                    if pos <= self.windowEnd[i]:
                        return self.individuals[ind].getAncInWindow(i)
                else:
                    diff1 = (pos - self.windowStart[i])*(pos - self.windowStart[i])
                    diff2 = (pos - self.windowEnd[i-1]) * (pos - self.windowEnd[i-1])

                    if diff1 < diff2:
                        return self.individuals[ind].getAncInWindow(i)
                    else:
                        return self.individuals[ind].getAncInWindow(i-1)

    def printMSP(self):

        print("================== Individual information ================== ")
        for ind in self.individuals:
            print(f'\n{ind} A', end='\t')
            self.individuals[ind].printHaplotype("A")
            print(f'\n{ind} B', end='\t')
            self.individuals[ind].printHaplotype("B")
        print("\n\n==================   Window information   ================== ")
        for i in range(len(self.windowStart)):
            print(f'{self.windowStart[i]} - {self.windowEnd[i]}')


class mosaics:
    def __init__(self):
        self.haplotypeA = []
        self.haplotypeB = []

    def printHaplotype(self, hapID):
        if hapID == "A":
            for window in self.haplotypeA:
                print(f' {window}', end="")
        elif hapID == "B":
            for window in self.haplotypeB:
                print(f' {window}', end="")

    def getAncInWindow(self, index):
        return int(self.haplotypeA[index]), int(self.haplotypeB[index])

    def addAnc(self, hapID, anc):
        if hapID == "A":
            self.haplotypeA.append(anc)
        elif hapID == "B":
            self.haplotypeB.append(anc)
        else:
            finish(f"[Error in MSP File] The haplotype {hapID} was not accepted")


def finish(message):
    sys.exit(message)


def openMSP(mspFile):
    file = open(mspFile)

    obj = MSP()

    numLine = 1
    for line in file:

        if numLine == 1:
            pass
        elif numLine == 2:
            header = line.strip().split('\t')
            for i in range(6, len(header), 2):
                obj.addIndividual(header[i])
        else:
            split = line.strip().split('\t')
            for i in range(6, len(split), 2):
                obj.addAncestryHaplotype(header[i], "A", split[i])
                obj.addAncestryHaplotype(header[i], "B", split[i + 1])

            obj.addInfoAboutWindow(split[0], split[1], split[2])
        numLine = numLine + 1
    return obj


def openCovar(covarFile, idCol):
    file = open(covarFile)
    covarDict = {}

    colNum = -1

    header = True
    for line in file:
        if header:
            headerLine = line.strip().split()
            header = False

            for i in range(len(headerLine)):
                if headerLine[i] == idCol:
                    if colNum >= 0:
                        finish(f"[Error in Covar File]: There is more than one ID column ({idCol})")
                    else:
                        colNum = i

            if colNum == -1:
                finish(f"[Error in Covar File]: The ID column ({idCol}) was not found")
        else:
            if line.strip() != "\n" and line.strip() != "":
                data = line.strip().split()

                id = data[colNum]
                if id in covarDict:
                    finish(f"[Error in Covar File]: The ID ({id}) is duplicate")
                covarDict[id] = {}

                for i in range(len(data)):
                    if i != colNum:
                        headerID = headerLine[i]
                        covarDict[id][headerID] = data[i]
    return covarDict


def checkIndividuals(vcfLine, msp, covar):
    data = vcfLine.strip().split()

    missing = False
    for i in range(9, len(data)):
        if data[i] not in covar:
            print(f"[Error in Covar]: The individual {data[i]} is present on VCF file but not in covar")
            missing = True
        if data[i] not in msp.individuals:
            print(f"[Error in MSP]: The individual {data[i]} is present on VCF file but not in MSP")
            missing = True
    if missing:
        finish("There is one or more individuals missing on Covar or MSP file")

    return data

def getGT(dataRaw):
    data = dataRaw.split(":")
    for field in data:
        if '|' in field:
            return field.split('|')

def getCovarHeader(covarDict, numAncestry):
    order = []
    header = "ID"
    order.append("ID")
    for ind in covarDict:
        for field in covarDict[ind]:
            header = header +f'\t{field}'
            order.append(field)

        for i in range(numAncestry):
            header = header + f'\tdosage{i}'

        for i in range(numAncestry-1):
            header = header + f'\tnumHap{i}'

        return header+"\n", order


def createLists(numAncesty):
    dosage = []
    haplotype = []
    for i in range(numAncesty):
        dosage.append(0)
        haplotype.append(0)
    return haplotype, dosage

def resetList(haplotype, dosage):
    for i in range(len(dosage)):
        haplotype[i] = 0
        dosage[i] = 0
    return haplotype, dosage

def runRegression(VCF, id, bcftools, Rscript, outputPrefix, msp, covarDict, numAncestry):
    extractedVCF = f'{outputPrefix}_Extracted_{id}.vcf.gz'
    os.system(f"{bcftools} view -r {id} -Oz -o {extractedVCF} {VCF}")
    covarFile, SNP = prepareCovar(extractedVCF, msp, covarDict, outputPrefix, numAncestry)
    os.system(f'{Rscript} {outputPrefix}_script.R {extractedVCF} {outputPrefix} {covarFile} {SNP}')

def callRunRegression(vcfFile, bcftools, Rscript, outputPrefix, msp, covarDict, numAncestry):
    SNPs = []

    decode = False
    if vcfFile[-2:] == "gz":
        file = gzip.open(vcfFile, 'r')
        decode = True
    elif vcfFile[-3:] == "vcf":
        file = open(vcfFile)

    header = True
    for line in file:
        if decode:
            line = line.decode("utf-8")

        if header:
            if line[0:6] == "#CHROM":
                header = False
        else:
            data = line.strip().split()
            SNPs.append(f'{data[0]}:{data[1]}')

    for id in SNPs:
        print(f'Running regression for SNP {id}')
        runRegression(vcfFile, id, bcftools, Rscript, outputPrefix, msp, covarDict, numAncestry)
        finish('check')


def createRScript(outputPrefix, statisticalModel, phenotype, id, kinship, covarDict, numAncestry):
    script = open(f'{outputPrefix}_script.R', 'w')
    print(f"Creating the R script ({outputPrefix}_script.R) ... ", end = "")
    if not statisticalModel:
        first = True
        model = ""
        ind = list(covarDict.keys())[0]
        for covar in covarDict[ind]:
            if covar != phenotype:
                if first:
                    first = False
                    model = f'{covar}'
                else:
                    model = f'{model} + {covar}'
        for i in range(numAncestry):
            model = model + f' + dosage{i}'

        for i in range(numAncestry - 1):
            model = model + f' + numHap{i}'

        script.write(f'library(GMMAT)\n'
                     f'library(SeqArray)\n'
                     f'options <- commandArgs(trailingOnly = TRUE)\n'
                     f'vcfFile <- options[0]\n'
                     f'outName <- options[1]\n'
                     f'covarFile <- options[2]\n'
                     f'SNP <- options[3]\n'
                     
                     f'gdsName= paste(outName, "_GDS", sep = \"\")\n'
                     
                     f'SeqArray::seqVCF2GDS(vcfFile, gdsName)\n'
                     
                     f'covar <- read.table(covarFile, header = T, sep = "\\t")\n'
                     
                     f'kinship <- as.matrix(read.table({kinship}, check.names=FALSE))\n'
                     
                     f'model <- glmmkin({phenotype} ~ {model} , data = covar, kins = kinship, '
                     f'id = \"ID\", family = binomial(link = \"logit\")) \n'
                     
                     f'Wald <- glmm.wald(model, data = covar, kins = kinship, id = \"ID\", infile = gdsName, '
                     f'snps = SNP, is.dosage = T)\n'
                     f'write.table(Wald, paste(outName,SNP, sep = "_"), sep="\t")\n')
        script.close()
    else:
        script.write(f'library(GMMAT)\n'
                     f'library(SeqArray)\n'
                     f'options <- commandArgs(trailingOnly = TRUE)\n'
                     f'vcfFile <- options[0]\n'
                     f'outName <- options[1]\n'
                     f'covarFile <- options[2]\n'
                     f'SNP <- options[3]\n'

                     f'gdsName= paste(outName, "_GDS", sep = \"\")\n'

                     f'SeqArray::seqVCF2GDS(vcfFile, gdsName)\n'

                     f'covar <- read.table(covarFile, header = T, sep = "\\t")\n'

                     f'kinship <- as.matrix(read.table({kinship}, check.names=FALSE))\n'

                     f'model <- glmmkin({phenotype} ~ {statisticalModel} , data = covar, kins = kinship, '
                     f'id = \"ID\", family = binomial(link = \"logit\")) \n'

                     f'Wald <- glmm.wald(model, data = covar, kins = kinship, id = \"ID\", infile = gdsName, '
                     f'snps = SNP, is.dosage = T)\n'
                     f'write.table(Wald, paste(outName,SNP, sep = "_"), sep="\t")\n')
        script.close()
    print('Done')



def prepareCovar(vcfFile, msp, covar, outputPrefix, numAncestry):
    decode = False
    if vcfFile[-2:] == "gz":
        file = gzip.open(vcfFile, 'r')
        decode = True
    elif vcfFile[-3:] == "vcf":
        file = open(vcfFile)

    covarHeader, order = getCovarHeader(covar, numAncestry)
    numHaplotype, dosageByAncestry = createLists(numAncestry)

    filesCreated = 0
    header = True
    for line in file:
        if decode:
            line = line.decode("utf-8")

        if header:
            if line[0:6] == "#CHROM":
                print("Check the individuals in all files")
                headerLine = checkIndividuals(line, msp, covar)
                print("Done")
                header = False
        else:
            data = line.strip().split()
            chr = data[0]
            pos = data[1]
            SNP = data[2]

            covarFile = open(f'{outputPrefix}_covar_{chr}_{pos}.txt', 'w')
            covarFile.write(covarHeader)
            for i in range(9, len(data)):
                outLine = ''

                ind = headerLine[i]
                anc1, anc2 = msp.getAncPos(ind, pos)

                GT = getGT(data[i])

                dosageByAncestry[anc1] = dosageByAncestry[anc1] + int(GT[0])
                dosageByAncestry[anc2] = dosageByAncestry[anc2] + int(GT[1])

                numHaplotype[anc1] = numHaplotype[anc1] + 1
                numHaplotype[anc2] = numHaplotype[anc2] + 1

                for field in order:
                    if(field == "ID"):
                        outLine = ind
                    else:
                        outLine = f'{outLine}\t{covarDict[ind][field]}'

                for j in range(numAncestry):
                    outLine = f'{outLine}\t{dosageByAncestry[j]}'

                for j in range(numAncestry-1):
                    outLine = f'{outLine}\t{numHaplotype[j]}'

                covarFile.write(f'{outLine}\n')
            covarFile.close()
    return f'{outputPrefix}_covar_{chr}_{pos}.txt', SNP

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tractor without Hail')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-m', '--msp', help='MSP file from RFMix (preferred the unkinked)', required=True)
    required.add_argument('-v', '--vcf', help='VCF file (preferred the unkinked)', required=True)
    required.add_argument('-n', '--numAncestry', help='Number of ancestries in MSP file', required=True)
    required.add_argument('-c', '--covar', help='Table with the covariatives', required=True)
    required.add_argument('-i', '--id', help='Name of the column with individual ID on covar file', required=True)
    required.add_argument('-p', '--phenotype', help='Name of the column with phenotype', required=True)
    required.add_argument('-k', '--kinship', help='File with kinship coefficientx', required=True)
    required.add_argument('-o', '--output', help='Output file name', required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-s', '--statisticalModel',
                          help='Statistical Model. If not provided, we will include all covariatives', required=False,
                          default='')
    optional.add_argument('-t', '--threads', help='Number of threads to run (default = 1)', required=False, default=1)
    optional.add_argument('-b', '--bcftools', help='Path for the BCFTOOLS (default bcftools)',
                          required=False, default='bcftools')
    optional.add_argument('-R', '--Rscript', help='Path for the Rscript with SeqArray and GMMAT installed (default Rscript)',
                          required=False, default='bcftools')

    args = parser.parse_args()
    msp = openMSP(args.msp)
    covarDict = openCovar(args.covar, args.id)
    createRScript(args.output, args.statisticalModel, args.phenotype, args.id, args.kinship, covarDict, int(args.numAncestry))
    callRunRegression(args.vcf, args.bcftools, args.Rscript, args.output, msp, covarDict, int(args.numAncestry))


    