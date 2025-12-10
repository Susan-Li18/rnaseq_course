#!//usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=00:10:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=create_sample_list

ls /data/courses/rnaseq_course/toxoplasma_de/reads_Lung/*_1.fastq.gz | \
awk -F'/' '{print $NF}' | \
awk -F'_1.fastq.gz' '{print $1}' > /data/users/jli/rnaseq_course/data/processed_data/samples.txt