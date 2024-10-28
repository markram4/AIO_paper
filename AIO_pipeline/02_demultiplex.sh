######################################
## 2 Demultiplex with cutadapt
######################################
### 1 Trim the sequences 5' of the barcodes before demulitplexing
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/logFiles

# require both the 5' and 3' adapters
## concatenate all fastq files from base calling - then use cutadapt to remove squencing adapters up/down stream of the barcodes
cat /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/1_guppyOUT/fc1/pass/* \
/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/1_guppyOUT/fc2/pass/* \
/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/1_guppyOUT/fc3/pass/* | \
/usr/local/bin/cutadapt -e 0.30 -j 5 -g file:/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/seq_5p_of_barcodes.fa \
-o /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/ruby.round2.{name}.trimmed.fastq \
--rc --action=trim \
--untrimmed-output /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/ruby.round2.untrimmed.fastq - 

## Count the number of reads that get trimmed
wc -l /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/ruby.round2.i5_i7.trimmed.fastq  | awk 'BEGIN {FS=" "}{OFS="\t"}{print $1/4, "trimmed"}' > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/ruby.round2.trimmed_count.txt &&
wc -l /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/ruby.round2.untrimmed.fastq | awk 'BEGIN {FS=" "}{OFS="\t"}{print $1/4, "untrimmed"}' >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/ruby.round2.trimmed_count.txt 

### 2 Demultiplex
## make file with barcode pairs used
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/logFiles
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/barcodes

## Make file with well number and sample name then pull out barcodes based on this number. Third column has product number where the barcodes are from.
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/barcodes
awk 'NR==FNR{A[$2]=$1;next}{if ($2 in A) print ">"A[$2] "\n" "^" $6 "..." $4 "$"}' /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/barcodes/ruby_round2_barcodes_used.txt /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/E6442_illumina_barcodes.txt > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/barcodes/ruby_round2_barcodes_used.anchored.fa

## Submit demultiplexing job
condor_submit /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/jobFiles/2_cutadapt_demultiplex.job 

## Count number of reads in each sample
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex

for file in *_*fastq; do awk 'BEGIN {FS=" "}{OFS="\t"}{n=split(FILENAME,a,".fastq")}END{print NR/4,a[1]}' $file >> fastq_read_counts.ruby_round2.txt; done &
for file in *fastq; do awk 'BEGIN {FS=" "}{OFS="\t"}{n=split(FILENAME,a,".fastq")}END{print NR/4,a[1]}' $file >> fastq_read_counts.plus_unassigned.txt; done  &

