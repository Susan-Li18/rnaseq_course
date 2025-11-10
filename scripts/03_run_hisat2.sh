#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=run_hisat2_align
#SBATCH --array=1-15
#SBATCH --output=/data/users/jli/rnaseq_course/data/processed_data/hisat2/output_hisat2_%j.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/processed_data/hisat2/error_hisat2_%j.e
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /data/users/jli/rnaseq_course/data/processed_data/samples.txt)
# set variable of input path and output path
INPUT_DIR="/data/users/jli/rnaseq_course/data/raw_data/fastq_file/reads_Lung"
OUTPUT_DIR="/data/users/jli/rnaseq_course/data/processed_data/hisat2"
INDEX_DIR="/data/users/jli/rnaseq_course/data/processed_data/hisat2_index"
# create output directory
mkdir -p ${OUTPUT_DIR}
# run the hisat2 for alignment
apptainer exec \
  --bind /data/users/jli/rnaseq_course/data/raw_data/fastq_file/reads_Lung \
  --bind /data/users/jli/rnaseq_course/data/processed_data \
  /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
  hisat2 -x ${INDEX_DIR}/Mus_musculus \
  -1 ${INPUT_DIR}/${SAMPLE}_1.fastq.gz \
  -2 ${INPUT_DIR}/${SAMPLE}_2.fastq.gz \
  -S ${OUTPUT_DIR}/${SAMPLE}_aligned.sam \
  --rna-strandness RF \
  -p ${SLURM_CPUS_PER_TASK}







