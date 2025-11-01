#!/bin/bash

#SBATCH --account=bgmp                          #REQUIRED: which account to use
#SBATCH --partition=bgmp                        #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                      #number of CPU cores to use, default is 1
#SBATCH --mem=16GB                              #memory per node, default is 1GB            
#SBATCH --job-name=dupededupe                         #optional: job name, default is the user name
#SBATCH --time=01:00:00                         #optional: time limit, default is 1 day
#SBATCH --output=/Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/unit_tests_%j.out                    #optional: output file name, default is slurm-<jobid>.out
#SBATCH --error=/Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/unit_tests_%j.err                     #optional: error file name,

#SBATCH --mail-type=END,FAIL                    #send email when job ends or fails
#SBATCH --mail-user=amdo@uoregon.edu            #email address to send notifications to 
module load python/3.8
module load samtools/1.10

test=/Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/test1.sam
sorted="${test}_sorted1.1.sam"

python ./dougherty_deduper.py -f $test -s $sorted -o /Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/deduped.sam -u /Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/STL96.txt

# python3 ./dougherty_deduper.py -f /Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/test1.sam -o deduped.sam -u /Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/STL96.txt