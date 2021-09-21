import sys
import gzip
import os
import time
import argparse
import multiprocessing


def init(queue):
    global idx
    idx = queue.get()

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

    def addAncestryHaplotype(self, id, hapID, anc, posteriori):
        self.individuals[id].addAnc(hapID, anc, posteriori)

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
        self.posterioriA = []
        self.haplotypeB = []
        self.posterioriB = []

    def printHaplotype(self, hapID):
        if hapID == "A":
            for i in range(len(self.haplotypeA)):
                print(f' {self.haplotypeA[i]}:{self.posterioriA[i]}', end="")
        elif hapID == "B":
            for i in range(len(self.haplotypeB)):
                print(f' {self.haplotypeB[i]}:{self.posterioriB[i]}', end="")


    def getAncInWindow(self, index):
        return int(self.haplotypeA[index]), int(self.haplotypeB[index]), float(self.posterioriA[index]), float(self.posterioriB[index])

    def addAnc(self, hapID, anc, posteriori):
        if hapID == "A":
            self.haplotypeA.append(anc)
            self.posterioriA.append(posteriori)
        elif hapID == "B":
            self.haplotypeB.append(anc)
            self.posterioriB.append(posteriori)
        else:
            finish(f"[Error in MSP File] The haplotype {hapID} was not accepted")


def finish(message):
    sys.exit(message)


def openMSP(mspFile, fowardBackward, numAnc):
    file = open(mspFile+'.msp.tsv')
    if fowardBackward:
        filePosteriori = open(mspFile + '.fb.tsv')
        #Remove both headers
        line = filePosteriori.readline()
        line = filePosteriori.readline()

    obj = MSP()

    headerRead = False
    for line in file:
        if not headerRead:
            if line[0:3] == '#ch':
                header = line.strip().split('\t')
                for i in range(6, len(header), 2):
                    obj.addIndividual(header[i])
                headerRead = True
        else:
            split = line.strip().split('\t')
            splitSlice = split[6:]
            if fowardBackward:
                fbLine = filePosteriori.readline()
                splitFB = [float(item) for item in fbLine.strip().split('\t')]
                fbSlice = splitFB[5:]

            for i in range(0, len(splitSlice), 2):
                if fowardBackward:
                    baseFB = (i)*numAnc

                    posteriori1 = max(fbSlice[baseFB:baseFB+numAnc])
                    posteriori2 = max(fbSlice[baseFB+numAnc:baseFB+numAnc*2])

                    obj.addAncestryHaplotype(header[i+6], "A", splitSlice[i], posteriori1)
                    obj.addAncestryHaplotype(header[i+6], "B", splitSlice[i+1], posteriori2)
                else:
                    obj.addAncestryHaplotype(header[i+6], "A", splitSlice[i], 1)
                    obj.addAncestryHaplotype(header[i+6], "B", splitSlice[i+1], 1)

            obj.addInfoAboutWindow(split[0], split[1], split[2])
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

def getDosage(dataRaw, description, dosageField, dosageSeparator, dosagePosition):
    position = dosagePosition.split()
    position[0] = int(position[0])
    position[1] = int(position[1])
    data = dataRaw.split(":")
    descriptionNames = description.split(":")
    for i in range(len(descriptionNames)):
        if descriptionNames[i] == dosageField:
            interest = i

    dataSplit = data[interest].split(dosageSeparator)
    return dataSplit[position[0]], dataSplit[position[1]]

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

        for i in range(numAncestry):
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

def runRegression(args):
    (VCF, chrPos, id, bcftools, Rscript, outputPrefix, msp, covarDict, numAncestry, gdsName, delete, dosageField, dosageSeparator, dosagePosition) = args
    global idx
    extractedVCF = f'{outputPrefix}_Extracted_{chrPos}_Worker{idx}.vcf.gz'
    os.system(f"{bcftools} view -r {chrPos} -Oz -o {extractedVCF} {VCF}")
    covarFile, SNP = prepareCovar(extractedVCF, msp, covarDict, outputPrefix, numAncestry, id, dosageField, dosageSeparator, dosagePosition)
    gdsName = gdsName.replace('IDTOPOOL', f'{idx}')
    outputWald = f'{outputPrefix}_Wald_{SNP.replace(":", "_")}'
    os.system(f'{Rscript} {outputPrefix}_script.R {outputWald} {covarFile} {id} {gdsName}')
    if delete:
        os.system(f'rm {extractedVCF} {covarFile}')



def convertToGDS(vcfFile, Rscript, outName, numProcess):
    script = open(f'{outName}_convertGDS.R', 'w')
    script.write(f'library(GMMAT)\n'
                 f'library(SeqArray)\n'

                 f'gdsName= \"{outName}_GDS_IDTOPOOL.gds\"\n'

                 f'SeqArray::seqVCF2GDS(\"{vcfFile}\", gdsName, parallel={numProcess})\n'
                 )
    script.close()
    os.system(f'{Rscript} {outName}_convertGDS.R')
    return (f"{outName}_GDS_IDTOPOOL.gds")

def callRunRegression(vcfFile, bcftools, Rscript, outputPrefix, msp, covarDict, numAncestry, numProcess, delete, dosageField, dosageSeparator, dosagePosition, insertHeader):
    params = []

    #print("Creating the GDS\n")
    #startTime = time.time()
    #gdsName = convertToGDS(vcfFile, Rscript, outputPrefix, numProcess)
    #endTime = time.time()
    #print(f"GDS created ({endTime-startTime} s)")

    gdsName = "Temporary"



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
            params.append([vcfFile, f'{data[0]}:{data[1]}',f'{data[2]}', bcftools, Rscript, outputPrefix, msp,
                           covarDict, numAncestry, gdsName, delete, dosageField, dosageSeparator, dosagePosition])

    #====================================================================================================
    #Pool to save time
    manager = multiprocessing.Manager()
    poolIDs = manager.Queue()
    for proc in range(numProcess):
        #os.system(f"cp {gdsName} {gdsName.replace('IDTOPOOL', f'{proc}')}")
        poolIDs.put(proc)


    pool = multiprocessing.Pool(numProcess, init, (poolIDs,))
    pool.map(runRegression, params)

    #======================================================================================================
    finalFile = open(f'{outputPrefix}_AllWald.txt', 'w')
    print(f'Merging all Wald files ({outputPrefix}_AllWald.txt)')

    if insertHeader:
        finalFile.write("SNP\tCHR\tPOS")
        for i in range(numAncestry):
            finalFile.write(f'\tBETA_anc{i}\tSE_anc{i}\tPVALUE_anc{i}')
        finalFile.write('\n')

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
            SNP = data[2]

            fileWaldName = f'{outputPrefix}_Wald_{SNP.replace(":", "_")}'
            fileWald = open(fileWaldName)

            finalFile.write(f"{data[2]}\t{data[0].replace('chr', '')}\t{data[1]}")

            for lineR in fileWald:
                split = lineR.strip().split()
                for i in range(0, len(split)):
                    withoutQuotes=split[i].strip().replace("\"", "")
                    finalFile.write(f'\t{withoutQuotes}')
                finalFile.write('\n')
            fileWald.close()
            if delete:
                os.system(f'rm {outputPrefix}_Wald_{SNP.replace(":", "_")}')
    finalFile.close()

def createAllModels(model, numAncestry, phenotype):
    modelWithDosage = model
    for i in range(numAncestry):
        modelWithDosage = modelWithDosage + f' + dosage{i}'

    for i in range(numAncestry-1):
        modelWithDosage = modelWithDosage + f' + numHap{i}'
    line = ""

    line = line + f'modelFull <- glm({phenotype} ~ {modelWithDosage}, data = covar, family = binomial)\n'

    for i in range(0, numAncestry):
        line = line+ f'beta{i} = coef(modelFull)[[\"dosage{i}\"]]\n'
        line = line +f'SE{i} = sqrt(vcov(modelFull)[\"dosage{i}\",\"dosage{i}\"])\n'
        line = line +f'pval{i} = pnorm(-abs(beta{i})/SE{i})*2\n'

    line = line + f'writeLines(paste('
    for i in range(0, numAncestry):
        line = line+f'beta{i}, SE{i}, pval{i},'
    line = line + 'sep = \" \"), fileOut)\n'
    line = line + f'close(fileOut)\n'
    return line



def createRScript(outputPrefix, statisticalModel, phenotype, id, kinship, covarDict, numAncestry, vcfFile, Rscript, tractorCovar):
    script = open(f'{outputPrefix}_script.R', 'w')
    print(f"Creating the R script ({outputPrefix}_script.R) ... ", end = "")

    script.write(f'library(GMMAT)\n'
                 f'library(SeqArray)\n'
                 f'library(lmtest)\n'
                 f'\noptions <- commandArgs(trailingOnly = TRUE)\n'
                 f'outName <- options[1]\n'
                 f'covarFile <- options[2]\n'
                 f'SNP <- options[3]\n'
                 f'GDS <- options[4]\n'
                 f'\nfileOut <- file(outName)\n'

                 f'\ncovar <- read.table(covarFile, header = T, sep = "\\t")\n'
                 f'covar$ID <- as.factor(covar$ID)\n'
                 f'covar${phenotype} <- as.factor(covar${phenotype})\n\n')

    #TODO KINSHIP
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

        lines = createAllModels(model, numAncestry, phenotype)
        script.write(lines)
    else:
        if tractorCovar:
            lines = createAllModels(statisticalModel, numAncestry, phenotype)
            script.write(lines)



    #     for i in range(numAncestry):
    #         model = model + f' + dosage{i}'
    #
    #     for i in range(numAncestry - 1):
    #         model = model + f' + numHap{i}'
    #
    #     if kinship:
    #         script.write(f'kinship <- as.matrix(read.table("{kinship}", check.names=FALSE))\n'
    #                      f'Wald <- glmm.wald({phenotype} ~ {model}, data = covar, kins = kinship, id = \"ID\", infile = GDS, '
    #                      f'snps = SNP, is.dosage = T)\n')
    #
    #     else:
    #         script.write(f'Wald <- glmm.wald({phenotype} ~ {model}, data = covar, id = \"ID\", infile = GDS, '
    #                      f'snps = SNP, is.dosage = T)\n')
    # else:
    #     if kinship:
    #         script.write(f'kinship <- as.matrix(read.table("{kinship}", check.names=FALSE))\n'
    #                      f'Wald <- glmm.wald({phenotype} ~ {statisticalModel}, data = covar, kins = kinship, id = \"ID\", infile = GDS, '
    #                      f'snps = SNP, is.dosage = T)\n')
    #     else:
    #         script.write(f'Wald <- glmm.wald({phenotype} ~ {statisticalModel}, data = covar, id = \"ID\", infile = GDS, '
    #                      f'snps = SNP, is.dosage = T)\n')
    #
    # script.write(f'write.table(Wald, outName, sep="\t")\n')
    script.close()



def prepareCovar(vcfFile, msp, covar, outputPrefix, numAncestry, targetSNP, dosageField, dosageSeparator, dosagePosition):
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
                #print("Check the individuals in all files")
                headerLine = checkIndividuals(line, msp, covar)
                #print("Done")
                header = False
        else:
            data = line.strip().split()
            chr = data[0]
            pos = data[1]
            SNP = data[2]

            if SNP == targetSNP:
                covarFile = open(f'{outputPrefix}_covar_{SNP.replace(":", "_")}.txt', 'w')
                covarFile.write(covarHeader)
                for i in range(9, len(data)):
                    outLine = ''

                    ind = headerLine[i]
                    anc1, anc2, posteriori1, posteriori2 = msp.getAncPos(ind, pos)
                    numHaplotype, dosageByAncestry = resetList(numHaplotype, dosageByAncestry)

                    dosages = getDosage(data[i], data[8], dosageField, dosageSeparator, dosagePosition)

                    dosageByAncestry[anc1] = dosageByAncestry[anc1] + float(posteriori1)*float(dosages[0])
                    dosageByAncestry[anc2] = dosageByAncestry[anc2] + float(posteriori2)*float(dosages[1])

                    numHaplotype[anc1] = numHaplotype[anc1] + 1
                    numHaplotype[anc2] = numHaplotype[anc2] + 1

                    for field in order:
                        if(field == "ID"):
                            outLine = ind
                        else:
                            outLine = f'{outLine}\t{covarDict[ind][field]}'

                    for j in range(numAncestry):
                        outLine = f'{outLine}\t{dosageByAncestry[j]}'

                    for j in range(numAncestry):
                        outLine = f'{outLine}\t{numHaplotype[j]}'

                    covarFile.write(f'{outLine}\n')
                covarFile.close()
                return f'{outputPrefix}_covar_{SNP.replace(":", "_")}.txt', SNP

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tractor without Hail')

    requiredGroup1 = parser.add_argument_group("Required arguments related to genetic files")
    requiredGroup1.add_argument('-m', '--mspPrefix', help='MSP prefix (without .msp.fb) file from RFMix (preferred the unkinked)', required=True)
    requiredGroup1.add_argument('-v', '--vcf', help='VCF file (preferred the unkinked)', required=True)
    requiredGroup1.add_argument('-n', '--numAncestry', help='Number of ancestries in MSP file', required=True)

    requiredGroup2 = parser.add_argument_group("Required arguments related to covar files")
    requiredGroup2.add_argument('-c', '--covar', help='Table with the covariatives', required=True)
    requiredGroup2.add_argument('-i', '--id', help='Name of the column with individual ID on covar file', required=True)
    requiredGroup2.add_argument('-p', '--phenotype', help='Name of the column with phenotype', required=True)

    requiredGroup3 = parser.add_argument_group("Required arguments related to output")
    requiredGroup3.add_argument('-o', '--output', help='Output prefix ', required=True)
    requiredGroup3.add_argument('-f', '--folder', help='Folder to store the results ', required=True)

    optionalGroup1 = parser.add_argument_group("Optional arguments to model")

    optionalGroup1.add_argument('-k', '--kinship', help='File with kinship coefficient', required=False, default = '')
    optionalGroup1.add_argument('-s', '--statisticalModel',
                          help='Statistical Model. If not provided, we will include all covariatives', required=False,
                          default='')
    optionalGroup1.add_argument('-T', '--Tractor',
                          help='If the user select the flag -s, this flag insert the Tractor covariates (Default True)', required=False,
                          default= True, action = 'store_false')

    optionalGroup2 = parser.add_argument_group("Optional arguments related to system")
    optionalGroup2.add_argument('-t', '--threads', help='Number of threads to run (default = 1)', required=False, default=1)
    optionalGroup2.add_argument('-b', '--bcftools', help='Path for the BCFTOOLS (default bcftools)',
                          required=False, default='bcftools')
    optionalGroup2.add_argument('-R', '--Rscript', help='Path for the Rscript with SeqArray and GMMAT installed (default Rscript)',
                          required=False, default='Rscript')
    optionalGroup2.add_argument('-e', '--exclude',
                          help='Set to exclude temporary files (default = False)',
                          required=False, default=False, action='store_true')
    optionalGroup2.add_argument('-H', '--Header',
                                help='Set to not output the header in all p-values merged (default = True, ie, output has a header)',
                                required=False, default=True, action='store_false')

    optionalGroup3 = parser.add_argument_group("Optional arguments related to dosage")
    optionalGroup3.add_argument('-D', '--dosageField',
                          help='Select the field to extract the dosage value (default: GT)',
                          required=False, default="GT")
    optionalGroup3.add_argument('-S', '--separatorDosageField',
                          help='Separator from select dosage field (default | for GT)',
                          required=False, default="|")
    optionalGroup3.add_argument('-P', '--positionDosageField',
                          help='Position form the values splited on dosage field (default \"0 1\" for GT)',
                          required=False, default="0 1")

    optionalGroup4 = parser.add_argument_group("Optional arguments related to local ancestry")
    optionalGroup4.add_argument('-F', '--fowardBackward',
                          help='Use the posterior probability from Foward backward to count the dosage per ancestry (default False)',
                          required=False, default=False, action= "store_true")

    args = parser.parse_args()
    msp = openMSP(args.mspPrefix, args.fowardBackward, int(args.numAncestry))

    os.system(f'mkdir {args.folder}')
    args.output = args.folder+'/'+args.output

    covarDict = openCovar(args.covar, args.id)
    createRScript(args.output, args.statisticalModel, args.phenotype, args.id, args.kinship, covarDict,
                  int(args.numAncestry), args.vcf, args.Rscript, args.Tractor)
    callRunRegression(args.vcf, args.bcftools, args.Rscript, args.output, msp, covarDict, int(args.numAncestry),
                      int(args.threads), args.exclude, args.dosageField, args.separatorDosageField,
                      args.positionDosageField, args.Header)

