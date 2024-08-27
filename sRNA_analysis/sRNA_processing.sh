############################
## 03 Trim Illumina Universal adapter
############################
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/01_fastqFiles/
ls *fq.gz | sed 's/.fq.gz//' > samples.txt

sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/jobFiles/02_trim_adapters.array.sh

############################
## 04 Map to genome, chloroplast and mitochondria and current transgenes
############################
#---------------- RUN COMMAND -------------------------
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/jobFiles/03_mapUnfiltered_reads.array.sh
#------------------------------------------------------

## The number of mapped reads from these files is used for downstream normalization, so make file with #mapped reads in each sample
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/04_mapUnfiltered
echo "Sample" > tmpA.txt &&
for file in *bam; do echo $(echo $file | sed 's/.unfiltered.mapped.bam//') >> tmpA.txt; done &

echo "Num_Mapped_Reads" > tmpB.txt &&
for file in *bam; do samtools view -c $file >> tmpB.txt; done &

paste -d "\t" tmpA.txt tmpB.txt > number_mapped_reads.R25.txt

############################
## 05 Filter for reads that don't map to chloroplast and mitochondria, and for size 18-28
############################
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/05_sizeFiltered
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/04_mapUnfiltered

#---------------- RUN COMMAND -------------------------
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/04_rmChlMito.slurm.sh R26 /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/01_fastqFiles/samples.dcl_RUBY.wt.txt
#------------------------------------------------------

###########################
## 06 Remove sRNAs from rRNA and tRNA using bowtie
############################
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/06_rRNA_tRNA_free

#---------------- RUN COMMAND -------------------------
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/05_remove_rRNA_tRNA.slurm.array.sh ruby_round5 R26 /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/01_fastqFiles/samples.dcl_RUBY.wt.txt
#------------------------------------------------------

############################
## 07 Map specific samples to TAIR + RUBY using ShortStack
############################
## Get known miRNAs from https://www.mirbase.org/download/ ath.gff3 and mature.fa
# rename
mv /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath.gff3 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath_known_miRNAs.miRBase.gff3 
# Pull out sequences
grep -A1 "Arabidopsis thaliana" /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mature.fa > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath_mature_miRNAs_seq.miRBase.fa

## Make loci file
awk '{FS=OFS="\t"}{print $1":1-"$2,$1}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.txt > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.loci.shortstack.txt
grep RUBY /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.loci.shortstack.txt > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.loci.shortstack.RUBY.txt

## Concatenate TAIR + RUBY and build bowtie index
bowtie-build /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/TAIR10_nuclear-Chr-all_Araport11.ruby_35s.transgene.fa

mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/shortStack/
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/shortStack/logFiles


## Map to 35S:RUBY and genome
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/jobFiles/08_map_to_TAIR_plus_RUBY.shortstack.35s.slurm.sh /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/06_rRNA_tRNA_free/ruby_samples.txt
