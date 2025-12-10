#!//usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=run_fastqc_clean
#SBATCH --output=/data/users/jli/rnaseq_course/data/processed_data/fastqc_clean/output_%j.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/processed_data/fastqc_clean/error_%j.e
INPUT_DIR=/data/users/jli/rnaseq_course/data/raw_data/clean_fastp
OUTPUT_DIR=/data/users/jli/rnaseq_course/data/processed_data/fastqc_clean
mkdir -p $OUTPUT_DIR
apptainer exec --bind $INPUT_DIR --bind $OUTPUT_DIR /containers/apptainer/fastqc-0.12.1.sif fastqc $INPUT_DIR/*.fastq.gz -t 4 -o $OUTPUT_DIR  
apptainer exec --bind $OUTPUT_DIR /containers/apptainer/multiqc-1.19.sif multiqc $OUTPUT_DIR -o $OUTPUT_DIR