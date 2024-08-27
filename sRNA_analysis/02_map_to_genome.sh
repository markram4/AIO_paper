############################
## 04 Map to genome, chloroplast and mitochondria and current transgenes
############################


#---------------- RUN COMMAND -------------------------
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/03_mapUnfiltered_reads.array.sh
#------------------------------------------------------


## The number of mapped reads from these files is used for downstream normalization, so make file with #mapped reads in each sample
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/04_mapUnfiltered
echo "Sample" > tmpA.txt &&
for file in *bam; do echo $(echo $file | sed 's/.unfiltered.mapped.bam//') >> tmpA.txt; done &

echo "Num_Mapped_Reads" > tmpB.txt &&
for file in *bam; do samtools view -c $file >> tmpB.txt; done &

paste -d "\t" tmpA.txt tmpB.txt > number_mapped_reads.R26.txt
