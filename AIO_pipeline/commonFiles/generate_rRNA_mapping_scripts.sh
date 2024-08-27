#!/bin/bash
# Get list of samples to map
# for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/*.trimmed.3p.trimmed.5p.fastq; do echo $(basename $file .fastq)| sed 's/.trimmed.3p.trimmed.5p//' >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/samples.txt; done 

## Usage :\
# sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_rRNA_mapping_scripts.sh\
# ruby_round2\
# /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/samples.txt

##---------------------------------------------------------------------------------

## FILENAMES
PROJECT=$1
SAMPLES_FILE=$2

##----------------------------------------------------------------------------------------
## DEFINE DIRECTORIES and VARIABLES
JOBDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/$PROJECT/jobFiles"
BASEDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/$PROJECT/data"
MAPDIR="6_minimap"
INDIR="$BASEDIR/4_cutadapt_trim_5p_3p_adapters/5p"
OUTDIR="1_map_to_rRNA_tRNA"
ANNOTATION="/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/tRNA_and_rRNA.KP.fas"

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

# Check and create directories
create_dir "$JOBDIR/logs"
create_dir "$BASEDIR/$MAPDIR/$OUTDIR/output"
create_dir "$BASEDIR/$MAPDIR/$OUTDIR/logFiles"
create_dir "$BASEDIR/$MAPDIR/$OUTDIR/output/mapped_to_rRNA"
create_dir "$BASEDIR/$MAPDIR/$OUTDIR/output/rRNA_free_fastq"

##----------------------------------------------------------------------------------------
# Count the number of lines in the samples file
NUM_SAMPLES=$(wc -l < "$SAMPLES_FILE")

# JOB SCRIPT NAME
# Get the current date and time
CURRENT_DATE=$(date +"%Y-%m-%d")
CURRENT_TIME=$(date +"%H-%M-%S")  # Use hyphens instead of colons to avoid issues in filenames

JOBSCRIPT="$JOBDIR/6_map_to_rRNA_${CURRENT_TIME}_${CURRENT_DATE}.slurm.sh"

##----------------------------------------------------------------------------------------
# Create the SLURM job script
cat <<EOT > $JOBSCRIPT
#!/bin/bash

#SBATCH --job-name=rRNA_map
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --output=$JOBDIR/logs/map_rRNA_%A_%a.txt
#SBATCH --error=$JOBDIR/logs/map_rRNA_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-$NUM_SAMPLES
#SBATCH --mem-per-cpu=16G
#SBATCH --time=72:00:00

## LOAD MODULES

module load minimap2
module load samtools

##----------------------------------------------------------------------------------------
## MAP TO FILE
sample=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

minimap2 -t 6 -G 12k -a "$ANNOTATION" "$INDIR/\${sample}.trimmed.3p.trimmed.5p.rRNA_free.fastq" 2> "$BASEDIR/$MAPDIR/$OUTDIR/logFiles/\${sample}.minimap.log" | samtools view -u - > "$BASEDIR/$MAPDIR/$OUTDIR/output/\${sample}.bam"

##----------------------------------------------------------------------------------------
## Separate into mapped and unmapped
echo "Separating mapped/unmapped"
samtools view -h -b -F 4 "$BASEDIR/$MAPDIR/$OUTDIR/output/\${sample}.bam" | samtools sort - > "$BASEDIR/$MAPDIR/$OUTDIR/output/$MAPPEDDIR/\${sample}.mapped_to_rRNA.bam"
samtools fastq -f 4 "$BASEDIR/$MAPDIR/$OUTDIR/output/\${sample}.bam" >  "$BASEDIR/$MAPDIR/$OUTDIR/output/$UNMAPPEDDIR/\${sample}.rRNA_free.fastq"

EOT
echo "$JOBSCRIPT"
sbatch $JOBSCRIPT 
