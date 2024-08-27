#!/bin/bash

#SBATCH --job-name=R26_trimming
#SBATCH --partition=general
#SBATCH --output=/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/logFiles/R26_trimming_output_%A_%a.txt
#SBATCH --error=/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/logFiles/R26_trimming_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mkramer@danforthcenter.org
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-96

##----------------------------------------------------------------------------------------
echo "$(hostname), reporting for duty."

## LOAD MODULES

module load miniconda3/4.10.3_gcc_12.3.0

source activate aio

cd "/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data"

##---------------------------------------------------------------------------------

## FILENAMES

file="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/01_fastqFiles/samples.txt"

##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES
INDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/01_fastqFiles"
OUTDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/03_trimmedFasta"

if [ ! -d "$OUTDIR" ]; then
    echo "Directory does not exist. Creating $OUTDIR..."
    mkdir "$OUTDIR"
    echo "Directory $OUTDIR created."
else
    echo "Directory $OUTDIR already exists."
fi

##----------------------------------------------------------------------------------------
## Run command

#while IFS=" " read -r sample
#do {
#	gzip -dc $INDIR/${sample}.fq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -c -Q 33 -v | fastq_to_fasta -Q 33 -v > $OUTDIR/${sample}_trimmed_unfiltered.fa
#
#} done <"$file"


sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$file")

gzip -dc "$INDIR/${sample}.fq.gz" | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -c -Q 33 -v | fastq_to_fasta -Q 33 -v > "$OUTDIR/${sample}_trimmed_unfiltered.fa"
