#!/bin/bash

##---------------------------------------------------------------------------------

## FILENAMES
PROJECT=$1
SAMPLES_FILE=$2


##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES
JOBDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/$PROJECT/jobFiles" 
INDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/$PROJECT/data/04_mapUnfiltered"
OUTDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/$PROJECT/data/05_sizeFiltered"

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

JOBSCRIPT="$JOBDIR/rmChlMito_${PROJECT}_${CURRENT_TIME}_${CURRENT_DATE}.slurm.sh"

##----------------------------------------------------------------------------------------
# Create the SLURM job script

cat <<EOT > $JOBSCRIPT
#!/bin/bash

#SBATCH --job-name=R26_rmChlMito_mapping
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --output=$JOBDIR/logs/R26_rmChlMito_mapping%A_%a.txt
#SBATCH --error=$JOBDIR/logs/R26_rmChlMito_mapping%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-$NUM_SAMPLES
#SBATCH --mem-per-cpu=16G

## LOAD MODULES

module load samtools

##---------------------------------------------------------------------------------
sample=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

samtools view -h "$INDIR/\${sample}.unfiltered.mapped.bam" | awk '\$1 ~ /^@/ || (\$3 != "ChrM" && \$3 != "ChrC" && \$3 != "BK059222.1" && \$3 != "BK016277.1"&& \$3 != "Ppat
ens_mito_AB251495.1"&& \$3 != "Ppatens_chloro_AP005672.2" )' | samtools fasta - > "$OUTDIR/\${sample}.mito_chloroFree.fa"
EOT

echo "Job file is : $JOBSCRIPT"
sbatch $JOBSCRIPT
