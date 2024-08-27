#!/bin/bash

# Get file of samples to map
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/*fastq; do echo $(basename $file .fastq) >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt ; done 


## Usage: 
#bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh\
#ruby_round2\
#/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt\
#/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.fa
##---------------------------------------------------------------------------------

## FILENAMES
PROJECT=$1
SAMPLES_FILE=$2
ANNOTATION=$3

##----------------------------------------------------------------------------------------
## DEFINE DIRECTORIES and VARIABLES
JOBDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/$PROJECT/jobFiles"
BASEDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/$PROJECT/data"
MAPDIR="6_minimap"
INDIR="$BASEDIR/$MAPDIR/1_map_to_rRNA_tRNA/output/rRNA_free_fastq"
OUTDIR="2_map_to_capture"

# Check if logs directory exists, if not, make it
create_dir() {
    local dir_path="$1"
    if [ ! -d "$dir_path" ]; then
        echo "Directory does not exist. Creating $dir_path..."
        mkdir -p "$dir_path"
        echo "Directory $dir_path created."
    else
        echo "Directory $dir_path already exists."
    fi
}

## Check if out directory exists, if it doesn't, make it
MAPPEDDIR="$BASEDIR/$MAPDIR/$OUTDIR/output/mapped_to_targets"
UNMAPPEDDIR="$BASEDIR/$MAPDIR/$OUTDIR/output/non_targets_fastq"

create_dir "$JOBDIR/logs"
create_dir "$BASEDIR/$MAPDIR/$OUTDIR"
create_dir "$MAPPEDDIR"
create_dir "$UNMAPPEDDIR"


##---------------------------------------------------------------------------------
## Separate mapped and unmapped reads

##----------------------------------------------------------------------------------------
# Count the number of lines in the samples file
NUM_SAMPLES=$(wc -l < "$SAMPLES_FILE")

# JOB SCRIPT NAME
# Get the current date and time
CURRENT_DATE=$(date +"%Y-%m-%d")
CURRENT_TIME=$(date +"%H-%M-%S")  # Use hyphens instead of colons to avoid issues in filenames

JOBSCRIPT="$JOBDIR/7_map_to_targets_${CURRENT_TIME}_${CURRENT_DATE}.slurm.sh"

##----------------------------------------------------------------------------------------
# Create the SLURM job script
cat <<EOT > $JOBSCRIPT
#!/bin/bash

#SBATCH --job-name=target_map
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --output=$JOBDIR/logs/map_targets_%A_%a.txt
#SBATCH --error=$JOBDIR/logs/map_targets_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --array=1-$NUM_SAMPLES
#SBATCH --mem-per-cpu=16G

## LOAD MODULES

module load minimap2
module load samtools

##---------------------------------------------------------------------------------
## MAP
sample=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

minimap2 -t 6 -G12k -ax splice "$ANNOTATION" "$INDIR/\${sample}.fastq" 2> "$BASEDIR/$MAPDIR/$OUTDIR/logFiles/\${sample}.minimap.log" | samtools view -u - > "$BASEDIR/$MAPDIR/$OUTDIR/output/\${sample}.bam"

##---------------------------------------------------------------------------------

echo "Separating mapped/unmapped"
samtools view -h -b -F 4 "$BASEDIR/$MAPDIR/$OUTDIR/output/\${sample}.bam" | samtools sort - > "$MAPPEDDIR/\${sample}.mapped_to_targets.bam"
samtools fastq -f 4 "$BASEDIR/$MAPDIR/$OUTDIR/output/\${sample}.bam" >  "$UNMAPPEDDIR/\${sample}.non_targets.fastq"
samtools index "$MAPPEDDIR/\${sample}.mapped_to_targets.bam"

EOT

echo "Job file is : $JOBSCRIPT"
sbatch $JOBSCRIPT
