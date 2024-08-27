############################
## 03 Trim Illumina Universal adapter
############################
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/01_fastqFiles/
ls *fq.gz | sed 's/.fq.gz//' > samples.txt

sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/02_trim_adapters.array.sh

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

############################
## 05 Filter for reads that don't map to chloroplast and mitochondria, and for size 18-28
############################
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/05_sizeFiltered
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/04_mapUnfiltered

#---------------- RUN COMMAND -------------------------
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/04_sizeFilter_rmChlMito.slurm.sh
#------------------------------------------------------

###########################
## 06 Remove sRNAs from rRNA and tRNA using bowtie
############################
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/06_rRNA_tRNA_free

#---------------- RUN COMMAND -------------------------
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/05_remove_rRNA_tRNA.slurm.array.sh
#------------------------------------------------------
