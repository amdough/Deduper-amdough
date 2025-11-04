#!/bin/bash

#SBATCH --account=bgmp                          #REQUIRED: which account to use
#SBATCH --partition=bgmp                        #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                      #number of CPU cores to use, default is 1
#SBATCH --mem=120GB                              #memory per node, default is 1GB            
#SBATCH --job-name=dupededupe4                        #optional: job name, default is the user name
#SBATCH --time=04:00:00                         #optional: time limit, default is 1 day
#SBATCH --output=./output/ddupe1%j.out                    #optional: output file name, default is slurm-<jobid>.out
#SBATCH --error=./output/ddupe1%j.err                     #optional: error file name,

#SBATCH --mail-type=END,FAIL                    #send email when job ends or fails
#SBATCH --mail-user=amdo@uoregon.edu            #email address to send notifications to 

# module load python/3.8
# Load modules
module load samtools/1.19
# module load python/3.8  # if needed, depending on your environment

# Define variables
input_sam=/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam
sorted_sam=/gpfs/projects/bgmp/amdo/bioinfo/Bi624/Deduper/C1_SE_uniqAlign_sorted.sam
deduped_out=/gpfs/projects/bgmp/amdo/bioinfo/Bi624/Deduper/C1_SE_deduped.sam
umis=/gpfs/projects/bgmp/amdo/bioinfo/Bi624/Deduper/STL96.txt
# tmp_dir=/gpfs/projects/bgmp/amdo/bioinfo/Bi624/Deduper/tmp

# Make sure all output directories exist
mkdir -p /gpfs/projects/bgmp/amdo/bioinfo/Bi624/Deduper/output
# mkdir -p $tmp_dir

# Run Python script with timing
/usr/bin/time python ./dougherty_deduper.py \
    -f $input_sam \
    -s $sorted_sam \
    -o $deduped_out \
    -u $umis \

# python3 ./dougherty_deduper.py -f /Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/test1.sam -o deduped.sam -u /Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/STL96.txt