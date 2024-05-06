######################################
## 2 Demultiplex with cutadapt
######################################
#####
### 1 Trim the sequences 5' of the barcodes before demulitplexing
#####
mkdir ~/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters
mkdir ~/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/logFiles

# Concatenate all fastq files that passed filters during basecalling 
cat ~/projects/target_capture/ruby_round2/data/1_guppyOUT/fc1/pass/*\
  ~/projects/target_capture/ruby_round2/data/1_guppyOUT/fc2/pass/*\
  ~/projects/target_capture/ruby_round2/data/1_guppyOUT/fc3/pass/* | \
## Use cutadapt to remove sequencing adapter sequences before the barcodes
/usr/local/bin/cutadapt\
  -e 0.30\
  -j 5\
  -g file:~/projects/target_capture/commonFiles/seq_5p_of_barcodes.fa \
  -o ~/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/ruby.round2.{name}.trimmed.fastq\
  --rc\
  --action=trim\
  --untrimmed-output ~/projects/target_capture/ruby_round2/data/2_cutadapt_trim_adapters/ruby.round2.untrimmed.fastq \
  - 

#####
### 2 Demultiplex
#####
## make file with barcode pairs used
mkdir ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex
mkdir ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/logFiles
mkdir ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/barcodes


## Make file with well number and sample name then pull out barcodes based on this number. Third column has product number where the barcodes are from.
## Use this file to pull out the appropriate barcodes from the file provided by NEB
cd /home/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/barcodes
awk 'NR==FNR{A[$2]=$1;next}{if ($2 in A) print ">"A[$2] "\n" "^" $6 "..." $4 "$"}'\
  ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/barcodes/ruby_round2_barcodes_used.txt\
  ~/projects/target_capture/commonFiles/E6442_illumina_barcodes.txt > \
  ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/barcodes/ruby_round2_barcodes_used.anchored.fa

## Submit demultiplexing job
condor_submit ~/projects/target_capture/ruby_round2/jobFiles/2_cutadapt_demultiplex.job 

/usr/local/bin/cutadapt \
  -e 0.25 \
  -j 10 \
  -g file:$(out_dir)/barcodes/ruby_round2_barcodes_used.anchored.fa\
  --info-file  $(out_dir)/cutadapt_demultplex.info.out  \
  -o $(out_dir)/{name}.fastq \
  --rc --action=lowercase \
  --untrimmed-output $(out_dir)/unassigned.fastq $(in_dir)/ruby.round2.i5_i7.trimmed.fastq

## Count number of reads in each sample
cd ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex

for file in *_*fastq; do awk 'BEGIN {FS=" "}{OFS="\t"}{n=split(FILENAME,a,".fastq")}END{print NR/4,a[1]}' $file >> fastq_read_counts.ruby_round2.txt; done &
for file in *fastq; do awk 'BEGIN {FS=" "}{OFS="\t"}{n=split(FILENAME,a,".fastq")}END{print NR/4,a[1]}' $file >> fastq_read_counts.plus_unassigned.txt; done  &

Rscript /home/mkramer/projects/target_capture/commonFiles/pieChartWClist.pctDemultiplex.R /home/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/fastq_read_counts.ruby_round2.txt /home/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/ruby_round2.demultiplex.pieChart.pdf

Rscript /home/mkramer/projects/target_capture/commonFiles/pieChartWClist.pctDemultiplex.R /home/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/fastq_read_counts.plus_unassigned.txt /home/mkramer/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/ruby_round2.demultiplex.pieChart.plus_unassigned.pdf
