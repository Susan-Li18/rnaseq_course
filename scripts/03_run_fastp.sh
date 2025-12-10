#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=run_fastp
#SBATCH --array=1-15
#SBATCH --output=/data/users/jli/rnaseq_course/data/raw_data/clean_fastp/output_%A.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/raw_data/clean_fastp/error_%A.e


# set the directory variable
OUTPUT_DIR=/data/users/jli/rnaseq_course/data/raw_data/clean_fastp
INPUT_DIR=/data/courses/rnaseq_course/toxoplasma_de/reads_Lung
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /data/users/jli/rnaseq_course/data/processed_data/samples.txt)

# create the output directory
mkdir -p $OUTPUT_DIR
# run the fastp
apptainer exec --bind $INPUT_DIR --bind $OUTPUT_DIR \
 /containers/apptainer/fastp_0.24.1.sif \
 fastp -i ${INPUT_DIR}/${SAMPLE}_1.fastq.gz \
       -I ${INPUT_DIR}/${SAMPLE}_2.fastq.gz \
       -o ${OUTPUT_DIR}/${SAMPLE}_1_clean.fastq.gz \
       -O ${OUTPUT_DIR}/${SAMPLE}_2_clean.fastq.gz \
       --detect_adapter_for_pe
 
