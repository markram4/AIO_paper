## Separate polyA and non-polyA reads
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets

for file in *targets.bam; do python /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/separate_pA_nonpA.py $file ; done &
