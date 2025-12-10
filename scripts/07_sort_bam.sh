#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=sort_bam
#SBATCH --array=1-15
#SBATCH --output=/data/users/jli/rnaseq_course/data/processed_data/bam_sorted/output_%A.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/processed_data/bam_sorted/error_%A.e

# set variable
SAMPLES_FILE="/data/users/jli/rnaseq_course/data/processed_data/samples.txt"
INPUT_DIR="/data/users/jli/rnaseq_course/data/processed_data/bam"
OUTPUT_DIR="/data/users/jli/rnaseq_course/data/processed_data/bam_sorted"
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${SAMPLES_FILE})
# create output folder
mkdir -p ${OUTPUT_DIR}
# execute the samtools sort command using slurm array
 apptainer exec --bind ${INPUT_DIR} --bind ${OUTPUT_DIR} \
 /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
 samtools sort -m 6G -@ 2 -o "${OUTPUT_DIR}/${SAMPLE}_sorted.bam" \
 -T temp_${SLURM_ARRAY_TASK_ID} ${INPUT_DIR}/${SAMPLE}_aligned.bam

