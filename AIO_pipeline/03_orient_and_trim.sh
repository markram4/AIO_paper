######################################
## 3 Orient and trim adapters 
######################################
## 3 Trim 3' adapter with cutadapt
cd ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/

cut -f 2 ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/fastq_read_counts.txt > ~/projects/target_capture/ruby_round2/data/ruby_round3_fileNames.txt

mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p/logFiles

condor_submit /home/mkramer/projects/target_capture/ruby_round2/jobFiles/3_cutadapt_trim_3p_adapters.job
/usr/local/bin/cutadapt -a CTGTAGGCACCATCAAT \
  -j 6 -n 3 --rc --info-file $(basedir)/$(outdir)/$(region)/logFiles/$(sample).info.$(region).out\
  -e 0.15 -o $(basedir)/$(outdir)/$(region)/$(sample).trimmed.$(region).fastq\
  --untrimmed-output $(basedir)/$(outdir)/$(region)/$(sample).untrimmed.$(region).fastq \
  $(basedir)/$(indir)/$(sample).fastq

## 4 Trim 5' adapter with cutadapt
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/logFiles
cd ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p

ls -ltr *_*fastq | awk 'BEGIN {FS=" "}{n=split($9,a,".fastq")}{print a[1]}' > samples.txt

condor_submit ~/projects/target_capture/ruby_round2/jobFiles/4_cutadapt_trim_5p_adapters.job
/usr/local/bin/cutadapt -g GGTATCAACGCAGAGTACATGGG -e 0.15 -n 2 -j 3  \
  --info-file $(basedir)/$(outdir)/$(region)/logFiles/$(sample).info.$(region).out \
  -o $(basedir)/$(outdir)/$(region)/$(sample).trimmed.$(region).fastq \
  --untrimmed-output $(basedir)/$(outdir)/$(region)/$(sample).untrimmed.$(region).fastq \
  $(basedir)/$(outdir)/$(indir)/$(sample).fastq

## How many reads have 5p and or 3p adapter
cd ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p

echo -e "sample" > tmp1.txt &&
ls ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/*fastq | awk 'BEGIN {FS=" "}{OFS="\t"}{n=split($1,a,"/")}{n=split(a[10],b,".")}{print b[1]}' | sort -u >> tmp1.txt &&

echo -e "5p_trimmed_only" > tmp2.txt &&
for sample in ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/*.untrimmed.3p.trimmed.5p.fastq  ; do awk 'BEGIN {FS=" "}{OFS="\t"}END{print NR/4}' $sample >> tmp2.txt ; done &&

echo -e "3p_trimmed_only" > tmp3.txt &&
for sample in ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/*.trimmed.3p.untrimmed.5p.fastq ; do awk 'BEGIN {FS=" "}{OFS="\t"}END{print NR/4}' $sample >> tmp3.txt ; done &&

echo -e "Both_trimmed" > tmp4.txt &&
for sample in ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/*.trimmed.3p.trimmed.5p.fastq ; do awk 'BEGIN {FS=" "}{OFS="\t"}END{print NR/4}' $sample >> tmp4.txt ; done &&

echo -e "untrimmed" > tmp5.txt &&
for sample in ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/*.untrimmed.3p.untrimmed.5p.fastq ; do awk 'BEGIN {FS=" "}{OFS="\t"}END{print NR/4}' $sample >> tmp5.txt ; done &&

paste -d "\t" tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt > number_reads_trimmed.txt &&
rm *tmp* &


/usr/bin/Rscript ~/projects/target_capture/commonFiles/barplot.cutadapt_trimming.R /home/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/number_reads_trimmed.txt

