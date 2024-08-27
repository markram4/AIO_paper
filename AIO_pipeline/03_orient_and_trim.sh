######################################
## 3 Orient and trim adapters 
######################################
## 3 Trim 3' adapter with cutadapt
cd ~/projects/target_capture/ruby_round2/data/3_cutadapt_demultiplex/

## Make required output directories
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p/logFiles


/usr/local/bin/cutadapt -a CTGTAGGCACCATCAAT \
  -j 6\
  -n 3\
  --rc\
  -e 0.15\
  --info-file ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p/logFiles/02_MK017_rep1.35S.8green.info.3p.out\
  -o ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p/02_MK017_rep1.35S.8green.trimmed.3p.fastq\
  --untrimmed-output ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p/02_MK017_rep1.35S.8green.untrimmed.3p.fastq \
  ~/projects/target_capture/ruby_round2/data/3p/02_MK017_rep1.35S.8green.fastq

## 4 Trim 5' adapter with cutadapt
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p
mkdir ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/logFiles
cd ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p

/usr/local/bin/cutadapt -g GGTATCAACGCAGAGTACATGGG\
  -e 0.15\
  -n 2\
  -j 3 \
  --info-file ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/logFiles/02_MK017_rep1.35S.8green.trimmed.3p.info.5p.out \
  -o ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/02_MK017_rep1.35S.8green.trimmed.3p.trimmed.5p.fastq \
  --untrimmed-output ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/02_MK017_rep1.35S.8green.trimmed.3p.untrimmed.5p.fastq \
  ~/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p/02_MK017_rep1.35S.8green.trimmed.3p.fastq

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

