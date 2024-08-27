##########################################
## 4 Map to rRNA/Chloroplast/Mitochondria/tRNA
##########################################
## MAP and then separate mapped / unmapped reads into appropriate directories and files
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/*.trimmed.3p.trimmed.5p.fastq; do echo $(basename $file .fastq)| sed 's/.trimmed.3p.trimmed.5p//' >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/samples.txt; done 

sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_rRNA_mapping_scripts.sh ruby_round2 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/5p/samples.txt
