# RFScripts - Scripts to run RFMix 

In this folder I describe the scripts implemented to perform the RFMix1 from a VCF File with admixed and references (only biallelics and commons SNPs)

## Step 1: Creating sets of individuals 

The RFMix1 is time and memory demanding. To run a file with more than 1500 individuals to chromosome 1 you need more than 500 GB of RAM. With this we have two options:

- Buy a new computer
- Split the data in sets

We do not have money enough to buy a computer with more than 500 GB of RAM, so we implemented the sets. To create the set list you have to run:

```
python3 prepareData.py -f <Fam file from Plink with all individuals> -n <number of individuals per set> -c <correspondence list>
```

The Fam file can be replaced by a list of individual per line (in theory). 

In our experiments we used 100 individuals per set. 

The correspondence list is a file with just individuals from reference. It is two collumsn, the first the individual ID and the second is the ancestry

```
NA20799 EUR
NA20808 EUR
NA20813 EUR
NA20815 EUR
NA20822 EUR
NA20832 EUR
```

## Step 2: The bot 

We ran this pipeline in a HPC with Slurm and our administrator allowed to run 28 jobs per time (with 24 processor each and 90 GB of RAM each). To avoid lose time, we created a 
bot that looks the queue each 1 second and, if there is less jobs than 26, create the bash file and submit to the queue. I am providing this script because they have all steps:

- Split VCF file
- Convert VCF to RFMix1 input
- Run RFMix1

## Step 3: RFmix1To2

With all sets finished, the next step is merge and convert to RFMix2 output format (FB and MSP)

```
python3.8 RFMix1ToRFMix2.py -v <Original VCF, with all individuals> -s <prefix of the sets with the number replaced by \*> 
-m <correspondence file> -b <first set (0)> -e <last set (N)> -C <path from classes file from RFMix10 with the number replaced by \* >
-F <path from FowardBackward file from RFMix10 with the number replaced by \* > -o <output prefix> -c <chr num> -g <genetic map used> 
-S <SNP per window file from RFMix1>
```


#### Acknowledgment
  
Developer: Thiago Peixoto Leal, PhD (thpeixotol@hotmail.com)

Contributors: Victor Borda, PhD - University of Maryland

Special Thanks: Michael Weiner, the HPC Administrator
