#!//usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=run_fastqc
#SBATCH --output=/data/users/jli/rnaseq_course/data/processed_data/fastqc/output_fastqc_%j.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/processed_data/fastqc/error_fastqc_%j.e
FASTQ_DIR=/data/courses/rnaseq_course/toxoplasma_de
OUTPUT_DIR=/data/users/jli/rnaseq_course/data/processed_data/fastqc
mkdir -p $OUTPUT_DIR
apptainer exec --bind $FASTQ_DIR --bind $OUTPUT_DIR /containers/apptainer/fastqc-0.12.1.sif fastqc $FASTQ_DIR/reads_Lung/*.fastq.gz -t 4 -o $OUTPUT_DIR  
apptainer exec --bind $OUTPUT_DIR /containers/apptainer/multiqc-1.19.sif multiqc $OUTPUT_DIR -o $OUTPUT_DIR