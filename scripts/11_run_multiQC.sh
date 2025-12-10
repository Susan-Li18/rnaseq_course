#!//usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:30:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=run_multiqc
#SBATCH --output=/data/users/jli/rnaseq_course/results_%j.o
#SBATCH --error=/data/users/jli/rnaseq_course/results_%j.e


INPUT_DIR="/data/users/jli/rnaseq_course/data/processed_data"
OUTPUT_DIR="/data/users/jli/rnaseq_course/results"

apptainer exec --bind $OUTPUT_DIR --bind $INPUT_DIR \
 /containers/apptainer/multiqc-1.19.sif multiqc $INPUT_DIR \
 --ignore "samples.txt" --ignore "counts_modified.txt" \
 --ignore "*.o" --ignore "*.e" \
 --ignore "$INPUT_DIR/hisat2_index/*" \
 -o $OUTPUT_DIR -n "rnaseq_analysis_report" --force