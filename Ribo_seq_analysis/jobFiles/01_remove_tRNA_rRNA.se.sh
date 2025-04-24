!/bin/bash


##---------------------------------------------------------------------------------

## FILENAMES

SAMPLES_FILE=$1

##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES
JOBDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/jobFiles"
INDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/01_trimmedFastq"
OUTDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/03_rRNA_tRNA_free"

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

JOBSCRIPT="$JOBDIR/rm_tRNA_rRNA_${CURRENT_DATE}_${CURRENT_TIME}.slurm.sh"

##----------------------------------------------------------------------------------------
# Create the SLURM job script

cat <<EOT > $JOBSCRIPT
#!/bin/bash

#SBATCH --job-name=rm_tRNA_rRNA
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --output=$JOBDIR/logs/rm_tRNA_rRNA%A_%a.txt
#SBATCH --error=$JOBDIR/logs/rm_tRNA_rRNA%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-$NUM_SAMPLES
#SBATCH --mem-per-cpu=16G

## LOAD MODULES

module load miniconda3/4.10.3_gcc_12.3.0

source activate ribo_seq

sample=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")


bowtie2 --local -L 10 -p 8 -x /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/bowtie_indexes/tRNA_and_rRNA.MK -U "$INDIR/\${sample}.fastq.gz" --un "$OUTDIR/\${sample}.rRNA_tRNA_free.fastq"

gzip -f "$OUTDIR/\${sample}.rRNA_tRNA_free.fastq"

EOT

echo "Job file is : $JOBSCRIPT"
sbatch $JOBSCRIPT
