import argparse
import sys

def getchar():
    c = sys.stdin.read(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PLINK to Table')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-i', '--input', help='Input from KING', required=True)
    required.add_argument('-o', '--output', help='Output name', required=True)

    args = parser.parse_args()

    dictInd = {}

    file = open(args.input)
    header = True
    for line in file:
        if header:
            header = False
        else:
            split = line.strip().split()

            ind1 = split[1]
            ind2 = split[3]
            kin = split[9]

            if ind1 not in dictInd:
                dictInd[ind1] = {}
            if ind2 not in dictInd:
                dictInd[ind2] = {}

            dictInd[ind1][ind2] = kin
            dictInd[ind2][ind1] = kin

    keys = list(dictInd.keys())
    fileOut = open(args.output, 'w')
    #fileOut.write('ID')
    for ind1 in keys:
        fileOut.write(f'\t{ind1}')
    fileOut.write('\n')


    for ind1 in keys:
        fileOut.write(f'{ind1}')
        for ind2 in keys:
            if ind1 == ind2:
                fileOut.write('\t1.0')
            else:
                fileOut.write(f'\t{dictInd[ind1][ind2]}')
        fileOut.write(f'\n')
    fileOut.close()