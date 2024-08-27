#!/bin/bash

##---------------------------------------------------------------------------------

## FILENAMES

PROJECT=$1
SRNA=$2
SAMPLES_FILE=$3

##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES
JOBDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/$SRNA/jobFiles" 
INDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/$SRNA/data/05_sizeFiltered"
OUTDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/$SRNA/data/06_rRNA_tRNA_free"
INDEX="/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/tRNA_and_rRNA.KP"

## Check if out directory exists, if it doesn't, make it
if [ ! -d "$OUTDIR" ]; then
    echo "Directory does not exist. Creating $OUTDIR..."
    mkdir -p "$OUTDIR"
else
    echo "Directory $OUTDIR already exists."
fi

##----------------------------------------------------------------------------------------
# Count the number of lines in the samples file
NUM_SAMPLES=$(wc -l < "$SAMPLES_FILE")

# JOB SCRIPT NAME
# Get the current date and time
CURRENT_DATE=$(date +"%Y-%m-%d")
CURRENT_TIME=$(date +"%H-%M-%S")  # Use hyphens instead of colons to avoid issues in filenames

JOBSCRIPT="$JOBDIR/rmrRNA_tRNA_${PROJECT}_${CURRENT_TIME}_${CURRENT_DATE}.slurm.sh"



##----------------------------------------------------------------------------------------
# Create the SLURM job script

cat <<EOT > $JOBSCRIPT
#!/bin/bash

#SBATCH --job-name=R26_rmrRNA_tRNA
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --output=$JOBDIR/logs/rmrRNA_tRNA_mapping%A_%a.txt
#SBATCH --error=$JOBDIR/logs/rmrRNA_tRNA_mapping%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-$NUM_SAMPLES
#SBATCH --mem-per-cpu=16G

## LOAD MODULES

module load miniconda3/4.10.3_gcc_12.3.0

source activate bowtie


sample=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

bowtie -v 0 -p 4 -fS $INDEX "$INDIR/\${sample}.mito_chloroFree.fa" | samtools fasta -f 4 - >  "$OUTDIR/\${sample}.mito_chloroFree.rRNA_tRNA_free.fasta"

EOT

echo "Job file is : $JOBSCRIPT"
sbatch $JOBSCRIPT
