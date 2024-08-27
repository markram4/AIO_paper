##########################################
## 6 Map to genome
##########################################
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/non_targets_fastq/*fastq; do echo $(basename $file .fastq)| sed 's/.trimmed.3p.trimmed.5p//' >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/non_targets_fastq/samples.txt; done 


bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_genome_mapping_scripts.sh ruby_round2 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/non_targets_fastq/samples.txt
