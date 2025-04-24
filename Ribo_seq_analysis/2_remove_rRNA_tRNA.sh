################
## 02 Remove rRNA/tRNA sequences
################
## Get file with sample names
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/01_trimmedFastq
for file in *fastq.gz; do echo $(basename $file .fastq.gz) >> samples.txt ;done

## Create a Bowtie2 index for contamination sequences and remove contaminating sequences

bowtie2-build /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/tRNA_and_rRNA.MK.fas /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/bowtie_indexes/tRNA_and_rRNA.MK

# Remove rRNA/tRNA using bowtie
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/jobFiles/01_remove_tRNA_rRNA.se.sh /cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/01_trimmedFastq/samples.txt
