import os
import time

for i in range(21, 0, -1):
    for j in range(0,15):
        putToRun = False
        while not putToRun:
            time.sleep(1)
            os.system('squeue > batch.txt')
            file = open('batch.txt')
            numLine = 0
            for line in file:
                numLine = numLine+1
            if numLine < 27: #I can run 28, 1 line is header and one is to this bot
                print(f'Put to run chr{i} set{j}')
                fileSlurm = open(f'slurmSet_chr{i}_{j}.pbs', 'w')
                fileSlurm.write('#!/bin/sh\n')
                fileSlurm.write('#SBATCH --mail-type=END,FAIL\n')
                fileSlurm.write('#SBATCH --mail-user=ABC@DEF.GHI\n')
                fileSlurm.write('#SBATCH -n 24\n')
                fileSlurm.write('#SBATCH --mem 90000\n')
        
                fileSlurm.write(f'#SBATCH --job-name=RF_{i}_{j}\n')
                fileSlurm.write(f'#SBATCH -o RF_{i}_{j}.out        # Standard output\n')
                fileSlurm.write(f'#SBATCH -e RF_{i}_{j}.err        # Standard error\n')
        
                fileSlurm.write('module load python/2.7.12\n')
                fileSlurm.write('module load python/3.8.6\n')
                fileSlurm.write('module load RFMix/1.5.4\n')
                fileSlurm.write('module load bcftools/1.9\n')
                fileSlurm.write('module load plink/1.90\n')
        
                fileSlurm.write(f'mkdir /home/USER/RFMix/inRF_chr{i}\n')
                fileSlurm.write(f'bcftools view -S /home/USER/RFMix/split{j} -Oz -o /home/USER/RFMix/inRF_chr{i}/Set{j}_chr{i}.vcf.gz /home/USER/RFMix/data{i}/Merged{i}.vcf.gz \n')
                fileSlurm.write(f'python3.8 /home/USER/RFMix/VCF2RFMix1.py --vcf /home/USER/RFMix/inRF_chr{i}/Set{j}_chr{i}.vcf.gz -c /home/USER/3A/TestAllRun/RFMix/ListParentalMapping.txt -m /home/USER/RFMix/chr{i}.b38.gmap -o /home/USER/RFMix/inRF_chr{i}/RFMi_chr{i}_set{j} -C {i}\n')
        
                fileSlurm.write(f'mkdir /home/USER/RFMix/outRF_chr{i}\n')
                fileSlurm.write(f'python2.7 /cm/shared/apps/RFMix/1.5.4/RunRFMix.py PopPhased /home/USER/RFMix/inRF_chr{i}/RFMi_chr{i}_set{j}_alleles /home/USER/RFMix/inRF_chr{i}/RFMi_chr{i}_set{j}_classes /home/USER/RFMix/inRF_chr{i}/RFMi_chr{i}_set{j}_location -o /home/USER/RFMix/outRF_chr{i}/output_chr{i}_set{j} --num-threads 24 -e 2 -w 0.2 --forward-backward --skip-check-input-format --succinct-output\n')
        
                fileSlurm.close()
        
                os.system(f'sbatch slurmSet_chr{i}_{j}.pbs')
                putToRun = True
