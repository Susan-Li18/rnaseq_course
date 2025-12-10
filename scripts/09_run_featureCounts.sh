#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=run_featureCounts
#SBATCH --output=/data/users/jli/rnaseq_course/data/processed_data/featureCounts/output_%j.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/processed_data/featureCounts/error_%j.e


# set variable
INPUT_DIR="/data/users/jli/rnaseq_course/data/processed_data/bam_sorted"
OUTPUT_DIR="/data/users/jli/rnaseq_course/data/processed_data/featureCounts"
ANNOTATION_DIR="/data/users/jli/rnaseq_course/data/external_data"
# create output folder
mkdir -p ${OUTPUT_DIR}
# execute the featureCounts command
 apptainer exec --bind ${INPUT_DIR} --bind ${OUTPUT_DIR} --bind ${ANNOTATION_DIR} \
 /containers/apptainer/subread_2.0.6.sif \
 featureCounts -T ${SLURM_CPUS_PER_TASK} -s 2 -p --countReadPairs -B -C \
 -a ${ANNOTATION_DIR}/Mus_musculus.GRCm39.115.gtf \
 -t exon -g gene_id -o ${OUTPUT_DIR}/counts.txt ${INPUT_DIR}/*.bam