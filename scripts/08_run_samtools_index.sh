#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=run_samtools_index
#SBATCH --array=1-15
#SBATCH --output=/data/users/jli/rnaseq_course/data/processed_data/bam_sorted/output_index_%A.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/processed_data/bam_sorted/error_index_%A.e

# set variable
SAMPLES_FILE="/data/users/jli/rnaseq_course/data/processed_data/samples.txt"
INPUT_DIR="/data/users/jli/rnaseq_course/data/processed_data/bam_sorted"
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${SAMPLES_FILE})

# execute the samtools index command using slurm array
 apptainer exec --bind ${INPUT_DIR} \
 /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
 samtools index -@ 4 ${INPUT_DIR}/${SAMPLE}_sorted.bam