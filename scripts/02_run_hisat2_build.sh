#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=run_hisat2
#SBATCH --output=/data/users/jli/rnaseq_course/data/processed_data/hisat2/output_hisat2_%j.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/processed_data/hisat2/error_hisat2_%j.e

OUTPUT_DIR=/data/users/jli/rnaseq_course/data/processed_data/hisat2_index
INPUT_DIR=/data/users/jli/rnaseq_course/data/external_data

mkdir -p $OUTPUT_DIR
apptainer exec --bind $INPUT_DIR --bind $OUTPUT_DIR \
 /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
 hisat2-build -f $INPUT_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa $OUTPUT_DIR/Mus_musculus