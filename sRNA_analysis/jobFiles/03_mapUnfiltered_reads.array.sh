#!/bin/bash

#SBATCH --job-name=R26_preFilter_mapping
#SBATCH --partition=general
#SBATCH --output=/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/logFiles/R26_preFilter_mapping_output_%A_%a.txt
#SBATCH --error=/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/logFiles/R26_preFilter_mapping_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mkramer@danforthcenter.org
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-96
#SBATCH --mem-per-cpu=16G

## LOAD MODULES

module load miniconda3/4.10.3_gcc_12.3.0

source activate bowtie

cd "/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data"

##---------------------------------------------------------------------------------

## FILENAMES

file="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/01_fastqFiles/samples.txt"

##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES
INDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/03_trimmedFasta"
OUTDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/04_mapUnfiltered"
INDEX="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/annotations/bowtie_indexes/TAIR10_Chr-all_Araport11+current_transgenes.plus_Cuscuta_Ppatens"

## Check if out directory exists, if it doesn't, make it
if [ ! -d "$OUTDIR" ]; then
    echo "Directory does not exist. Creating $OUTDIR..."
    mkdir "$OUTDIR"
    echo "Directory $OUTDIR created."
else
    echo "Directory $OUTDIR already exists."
fi

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$file")

bowtie -v 0 -p 4 -fS $INDEX "$INDIR/${sample}_trimmed_unfiltered.fa" | samtools view -bS -F 4 -@ 4 - | samtools sort -@ 4 -  -o "$OUTDIR/${sample}.unfiltered.mapped.bam"
