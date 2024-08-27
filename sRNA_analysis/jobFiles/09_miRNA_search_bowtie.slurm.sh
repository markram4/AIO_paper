#!/bin/bash

#SBATCH --job-name=R25_miRNAsearch
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --output=/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/jobFiles/logFiles/miRNAsearch_output_%A_%a.txt
#SBATCH --error=/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/jobFiles/logFiles/miRNAsearch_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-62
#SBATCH --mem-per-cpu=16G

## LOAD MODULES

module load miniconda3/4.10.3_gcc_12.3.0

source activate ShortStack4 

cd "/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data"

##---------------------------------------------------------------------------------

## FILENAMES

#file="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/01_fastqFiles/samples.dcl_RUBY.txt"
file=$1

##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES
INDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/06_rRNA_tRNA_free"
OUTDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search"
INDEX="/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/bowtie_indexes/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.fa"

## Check if out directory exists, if it doesn't, make it
if [ ! -d "$OUTDIR" ]; then
    echo "Directory does not exist. Creating $OUTDIR..."
    mkdir "$OUTDIR"
    echo "Directory $OUTDIR created."
else
    echo "Directory $OUTDIR already exists."
fi

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$file")
echo $sample

bowtie -p 4 -v 3 -a --best --strata $INDEX "$INDIR/${sample}_trimmed.mito_chloroFree.sizeFilt.rRNA_tRNA_free.fasta" -fS | samtools view -h -b -F 4 - | samtools sort > "$OUT
DIR/${sample}.miRNA_search.bam" 

samtools index "$OUTDIR/${sample}.miRNA_search.bam"
