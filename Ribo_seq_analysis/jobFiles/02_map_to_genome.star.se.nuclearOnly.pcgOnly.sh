#!/bin/bash


##---------------------------------------------------------------------------------

## FILENAMES

SAMPLES_FILE=$1

##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES
JOBDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/jobFiles"
INDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/03_rRNA_tRNA_free"
OUTDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/04_map_to_genome/nuclearOnly_pcg"

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

JOBSCRIPT="$JOBDIR/map_to_genome_${CURRENT_TIME}_${CURRENT_DATE}.slurm.sh"

##----------------------------------------------------------------------------------------
# Create the SLURM job script

cat <<EOT > $JOBSCRIPT
#!/bin/bash

#SBATCH --job-name=map_to_genome
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --output=$JOBDIR/logs/map_to_genome%A_%a.txt
#SBATCH --error=$JOBDIR/logs/map_to_genome%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1-$NUM_SAMPLES
#SBATCH --mem-per-cpu=16G

## LOAD MODULES

module load miniconda3/4.10.3_gcc_12.3.0

source activate star_260c


sample=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

if [ ! -d "$OUTDIR/\${sample}" ]; then
    echo "Directory does not exist. Creating $OUTDIR..."
    mkdir -p "$OUTDIR/\${sample}"
else
    echo "Directory $OUTDIR/\${sample} already exists."
fi

STAR --runThreadN 10 --outReadsUnmapped Fastx --genomeDir /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/TAIR_RUBY_STAR_noPtMt_pcgOnly --readFilesCommand zcat --readFilesIn "$INDIR/\${sample}.rRNA_tRNA_free.fastq" --alignIntronMax 5000 --alignIntronMin 15 --outFilterMismatchNmax 2 --outFilterMultimapNmax 20 --outFilterType BySJout --alignSJoverhangMin 4 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outSAMmultNmax 1 --outMultimapperOrder Random --outFileNamePrefix "$OUTDIR/\${sample}/\${sample}" 


EOT

echo "Job file is : $JOBSCRIPT"
sbatch $JOBSCRIPT
