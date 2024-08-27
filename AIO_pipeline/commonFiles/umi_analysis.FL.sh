#!/bin/bash

#SBATCH --job-name=umiFL_analysis
#SBATCH --partition=slotkinr-lab
#SBATCH --account=slotkinr-lab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=72:00:00
#SBATCH --output=/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/logs/umiFL_%A_%a.txt
#SBATCH --error=/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/logs/umiFL_%A_%a.err

######################################################
## Examine how many UMIs are present in subsets of data
######################################################
## Define project and base directory
PROJECT="$1"
ARRAY="$2"
JOBDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/$PROJECT/jobFiles"
BASEDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/$PROJECT/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles"
SAMPLES="$BASEDIR/strandedness/sampleNames.txt"
UMIDIR="$BASEDIR/full_length/umi"
INDIR="$BASEDIR/full_length"
TARGETS="/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_$ARRAY/At_array.$ARRAY.targets_class.all.txt"
FASTQDIR="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/$PROJECT/data/4_cutadapt_trim_5p_3p_adapters/3p_noTrim"
SCRIPTS="/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/umi"


##----------------------------------------------------------------------------------------------------------------
## Function to check and create directories
check_and_create_directory() {
    local directory="$1"
    if [ ! -d "$directory" ]; then
        echo "Directory does not exist. Creating $directory..."
        mkdir -p "$directory"
        echo "Directory $directory created."
    fi
}

## Check and create main UMI directory
check_and_create_directory "$UMIDIR"

## List of UMI subdirectories and sample subdirectories to check/create
directories=("01_readNames" "02_subset_seq" "03_extract" "04_filtered" "05_clustered_cd")
file_with_dirs="$SAMPLES"

## Loop through directories and create them if necessary
for dir in "${directories[@]}"; do
    full_path="$UMIDIR/$dir"
    check_and_create_directory "$full_path"
    
    while IFS= read -r subdir; do
        full_sub_dir="$full_path/$subdir"
        check_and_create_directory "$full_sub_dir"
    done < "$file_with_dirs"
done


##----------------------------------------------------------------------------------------------------------------
## Define subdirs
READNAMESDIR="$UMIDIR/01_readNames"
SUBSETDIR="$UMIDIR/02_subset_seq"
EXTRACTDIR="$UMIDIR/03_extract"
FILTERDIR="$UMIDIR/04_filtered"
CLUSTERDIR="$UMIDIR/05_clustered_cd"


#1# Read Names Pull our read names that map to each gene in the sense or antisense orientation
## Define files

echo "Step 1 : Get Read Names"

for file in "$INDIR"/*full_length.bed12; do
    # Loop over each target in TARGETS
    while IFS= read -r target; do
        
        # Create array with file name
        IFS="." read -r -a array <<< $(basename $file)
        if [[ "${array[0]}" == *MK* ]] || [[ "${array[0]}" == *Ls* ]]; then
        	samplename=$(echo "${array[0]}"".""${array[1]}"".""${array[2]}")
        	classname=$(echo "${array[5]}"".""${array[6]}"".""${array[7]}")
        else
        	samplename=$(echo "${array[0]}")
        	classname=$(echo "${array[3]}"".""${array[4]}""."${array[5]})
        fi
        # Get class of read from file name
        
        
        # Construct the output filename
        output_file="$READNAMESDIR/$samplename/$samplename.$target.$classname.txt"
        
        # Perform grep, cut, and sort operations
        grep "$target" "$file" | cut -f16 | sort -u > "$output_file"
    done < "$TARGETS"
done

##----------------------------------------------------------------------------------------------------------------
#2# subset Pull out the sequences from the oriented but not trimmed fastq file for each target

echo "Step 2 : Subseq Fastq files"

## Extract subsequences
module load seqtk
rm $SUBSETDIR/*/*fastq

for dir in $READNAMESDIR/*
	do echo "Starting $(basename $dir)"
	for file in $dir/*txt
		do
			if [ -s $file ]; then
				
				## Get genotype for subdir and subsequent naming
				IFS="." read -r -a array <<< $(basename $file)
        		if [[ "${array[0]}" == *MK* ]] || [[ "${array[0]}" == *Ls* ]]; then
        			genotype=$(echo "${array[0]}"".""${array[1]}"".""${array[2]}")
        		else
        			genotype=$(echo "${array[0]}")
        		fi

				## Create out file
				outFile=$SUBSETDIR/$genotype/$(echo $(basename $file) | sed 's/.txt/.orient.fastq/')
				## Get fastq file name
				inSeq=$genotype".trimmed.3p_noTrim.fastq"
				## Subset
				seqtk subseq $FASTQDIR/$inSeq  $file > $outFile
			fi
		done
	echo "Finished $(basename $dir)"
done

##----------------------------------------------------------------------------------------------------------------
#3# Extract UMIs
echo "Step 3 : Extracting UMIs"

module load miniconda3
source activate UMIC-seq

for dir in $SUBSETDIR/*
	do echo "Starting $(basename $dir)"
	for file in $dir/*fastq
		do
			## Get genotype for subdir and subsequent naming
			IFS="." read -r -a array <<< $(basename $file)
        	if [[ "${array[0]}" == *MK* ]] || [[ "${array[0]}" == *Ls* ]]; then
        		genotype=$(echo "${array[0]}"".""${array[1]}"".""${array[2]}")
        	else
        		genotype=$(echo "${array[0]}")
        	fi

			## Get outfile name
			outFile=$EXTRACTDIR/$genotype/$(echo $(basename $file) | sed 's/orient.fastq/UMIs.fasta/')

			## extract UMIs
			if [ ! -f $outFile ]
				then
					python $SCRIPTS/UMIC-seq.py UMIextract --input $file --probe $SCRIPTS/probe_seq.fasta --umi_loc down --umi_len 18 --output $outFile
			fi
		done
	echo "Finished $(basename $dir)"
done 


##----------------------------------------------------------------------------------------------------------------
#4# Filter UMI
echo "Step 4 : Filtering UMIs"

for dir in $EXTRACTDIR/*
do echo "Starting $(basename $dir)"
	for file in $dir/*fasta
		do
			## Get genotype for subdir and subsequent naming
			IFS="." read -r -a array <<< $(basename $file)
        	if [[ "${array[0]}" == *MK* ]] || [[ "${array[0]}" == *Ls* ]]; then
        		genotype=$(echo "${array[0]}"".""${array[1]}"".""${array[2]}")
        	else
        		genotype=$(echo "${array[0]}")
        	fi

			## Get out file name
			outFile=$FILTERDIR/$genotype/$(echo $(basename $file) | sed 's/fasta/filtered.fasta/')

			grep -B1 -E "[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}" $file | sed '/^--$/d' > $outFile
		done
	echo "Finished $(basename $dir)"
done


##----------------------------------------------------------------------------------------------------------------
#5# Cluster UMI 
echo "Step 5 : Clustering UMIs"

## CD-HIT
module load miniconda3
source activate cd-hit

for dir in $FILTERDIR/*
do echo "Starting $(basename $dir)"
	for file in $dir/*fasta
		do
			## Get genotype for subdir and subsequent naming
			IFS="." read -r -a array <<< $(basename $file)
        	if [[ "${array[0]}" == *MK* ]] || [[ "${array[0]}" == *Ls* ]]; then
        		genotype=$(echo "${array[0]}"".""${array[1]}"".""${array[2]}")
        	else
        		genotype=$(echo "${array[0]}")
        	fi

			## Get out file name
			outFile=$CLUSTERDIR/$genotype/$(echo $(basename $file) | sed 's/fasta/clustered_cd.fasta/')

			cd-hit-est -i $file -o $outFile -c 0.8888
		done
	echo "Finished $(basename $dir)"
done

##----------------------------------------------------------------------------------------------------------------
#6# Count number of unique umis mapping to each target
echo "Step 6 : Count Unique UMIs"

## Initiate files with header
while IFS= read -r line; do echo -e "Gene\tCount\tSample" > "$CLUSTERDIR"/$line".mapped_to_targets.umiCount.FL.txt"; done < $SAMPLES 

for dir in $CLUSTERDIR/*; do
    if [ -d "$dir" ]; then
        echo "Starting $(basename $dir)"
        for file in $dir/*fasta; do
            ## Get genotype for subdir and subsequent naming
            IFS="." read -r -a array <<< $(basename $file)
        	if [[ "${array[0]}" == *MK* ]] || [[ "${array[0]}" == *Ls* ]]; then
        		genotype=$(echo "${array[0]}"".""${array[1]}"".""${array[2]}")
        		TARGET="${array[3]}"
            	CLASS="${array[4]}"
            	STRAND="${array[5]}"
        	else
        		genotype=$(echo "${array[0]}")
        		TARGET="${array[1]}"
            	CLASS="${array[2]}"
            	STRAND="${array[3]}"
        	fi

            ## Get out file
            outFile="$CLUSTERDIR/$genotype.mapped_to_targets.umiCount.FL.txt"
    
    
            GENE="$TARGET.$CLASS.$STRAND"
            line_count=$(wc -l < "$file")
            COUNT=$((line_count / 2))
            
            echo -e "$GENE\t$COUNT\t$genotype" >> "$outFile"
        done
        echo "Finished $(basename $dir)"
    fi
done


## Combine all Genotypes
echo "Combining all genotypes"
echo -e "Gene\tCount\tSample" > $CLUSTERDIR/"all.combined_numReads_mapped_to_targets.strandedness.perpA.$PROJECT.umi.FL.txt"
for file in $CLUSTERDIR/*.umiCount.FL.txt; do tail -n +2 $file >> $CLUSTERDIR/"all.combined_numReads_mapped_to_targets.strandedness.perpA.$PROJECT.umi.FL.txt"; done

echo "UMI script complete"
