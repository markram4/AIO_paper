#!/bin/bash

#SBATCH --job-name=R26_map_to_targets
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --output=/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/jobFiles/logFiles/R26_map_to_targets_output_%A_%a.txt
#SBATCH --error=/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/jobFiles/logFiles/R26_map_to_targets_%A_%a.err
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
OUTDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/shortStack"
INDEX="/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/TAIR10_nuclear-Chr-all_Araport11.ruby_35s.transgene.fa"
MIRNA="/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath_mature_miRNAs_seq.miRBase.fa"
LOCI="/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.loci.shortstack.RUBY.txt"


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
ShortStack --known_miRNAs $MIRNA --locifile $LOCI --mmap f --threads 4  --outdir "$OUTDIR/${sample}_mapped" --readfile "$INDIR/${sample}_trimmed.mito_chloroFree.sizeFilt.rRNA_tRNA_free.fasta" --genomefile $INDEX

