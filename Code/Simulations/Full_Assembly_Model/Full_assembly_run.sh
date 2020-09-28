#!/bin/sh
#$ -cwd              # Set the working directory for the job to the current directory
#$ -V
#$ -pe smp 1         # Request 1 CPU cores
#$ -l h_rt=24:00:00  # Request 12 hour runtime
#$ -l h_vmem=2G      # Request 5GB RAM/core, i.e. 20GB total in total
#$ -t 1-50           # Request number of nodes

module load use.dev
module load singularity
module load R
shopt -s extglob
set -m
time=$(date)


INPUT_FILE1=$(sed -n "${SGE_TASK_ID}p" list_of_files_run.txt) # Run 
INPUT_FILE2=$INPUT_FILE1; INPUT_FILE2+=.cfg	
number=0
target_n=1200

while [ ${number} -le $target_n ]
do

	(cd ../Data/$INPUT_FILE1 && singularity exec ../PDMM.img NewWeb/build/NewWorld $INPUT_FILE2);

	
	number="$(ls -t *bz2 | wc -l)"

	INPUT_FILE2="$(ls -t *bz2 | head -1)"
done

time=$(date)





