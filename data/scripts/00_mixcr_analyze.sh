#!/bin/bash
# 00_mixcr_analyze.sh
# auxiliary file: SEEDFILE=sampleslist.txt

# All lines with # before means comment to bash
# However, if you put #SBATCH, the SBATCH program uses every line as parameter

# This script is launched in shell as:
# $sbatch sbatch_onejob_multitask.sh
# Usefull to launch one job to the queue system but using an arrayjob, i.e. one core per task.
# This is usefull when you want to lauch one program several times but in different files.
# You will have the same bash line but changing input and output.


#SBATCH --job-name=arrayJob
#SBATCH --output=arrayJob_%A_%a.out #archivo_%j.out
#SBATCH --error=arrayJob_%A_%a.err #archivo_%j.err
#SBATCH --array=1-22
#SBATCH --partition=long
#SBATCH --cpus-per-task 3 # -c
#SBATCH --mem 10G



# Modules
module load Java

# Input variables
SEEDFILE=sampleslist.txt
SEED1=$(cat $SEEDFILE | sed -n 1~2p | head -n $SLURM_ARRAY_TASK_ID | tail -n 1) #sed silences odds lines
SEED2=$(cat $SEEDFILE | sed -n 0~2p | head -n $SLURM_ARRAY_TASK_ID | tail -n 1) #sed silences even lines

# Data directory
DATADIR=data_RNASeq/

# Launch commands
mixcr analyze shotgun \
-t 3 \
--species hs \
--starting-material rna \
--only-productive \
${DATADIR}$SEED1 \
${DATADIR}$SEED2 \
${SEED1%read1.fastq.gz}
