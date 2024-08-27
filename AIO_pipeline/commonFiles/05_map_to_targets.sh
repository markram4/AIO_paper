##########################################
## 5 Map to capture sequences
##########################################
cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/ruby_35s.transgene.fa > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.fa

## This file will map then split into unmapped and mapped directories and will index bam files
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/*fastq; do echo $(basename $file .fastq) >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt ; done 

bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round2 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.fa
