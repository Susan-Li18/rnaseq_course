#!//usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=0:05:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=modifiy_file
#SBATCH --output=/data/users/jli/rnaseq_course/data/processed_data/featureCounts/output_modified_%j.o
#SBATCH --error=/data/users/jli/rnaseq_course/data/processed_data/featureCounts/error_modified_%j.e
# set DIR variable
INPUT_DIR="/data/users/jli/rnaseq_course/data/processed_data/featureCounts"
cd ${INPUT_DIR}
tail -n +2 counts.txt | cut -f 2,3,4,5,6 --complement > counts_modified.txt