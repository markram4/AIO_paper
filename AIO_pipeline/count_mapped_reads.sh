##########################################
## 7 Count number of reads from rRNA, targets, genome, unmapped
##########################################
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/mapped_pct
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/mapped_pct

echo -e "sample" > tmp1.txt &&
ls /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/*bam | awk 'BEGIN {FS=" "}{n=split($1,a,".bam")}{n=split(a[1],b,"/")}{print b[13]}' | sort -u >> tmp1.txt &&

echo -e "rRNA" > tmp2.txt &&
for sample in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/mapped_to_rRNA/*mapped_to_rRNA.bam ; do samtools view $sample | cut -f 1 | awk '{if($1 !~ "@") print}'| sort -u | wc -l >> tmp2.txt ; done &&

echo -e "Targets" > tmp3.txt &&
for sample in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/*.mapped_to_targets.bam ; do samtools view $sample | cut -f 1 | awk '{if($1 !~ "@") print}'| sort -u | wc -l >> tmp3.txt ; done &&

echo -e "Non-targets" > tmp4.txt &&
for sample in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/3_map_to_genome/output/mapped_to_genome/*mapped_to_genome.bam ; do samtools view $sample | cut -f 1 | awk '{if($1 !~ "@") print}'| sort -u | wc -l >> tmp4.txt ; done &&

echo -e "unmapped" > tmp5.txt &&
for sample in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/3_map_to_genome/output/unmapped/*unmapped.bam ; do samtools view $sample | cut -f 1 | awk '{if($1 !~ "@") print}'| sort -u | wc -l >> tmp5.txt ; done &&

paste tmp1.txt tmp2.txt tmp3.txt tmp4.txt tmp5.txt > number_reads_mapped.ruby_round2.txt &&
rm tmp*

module load r

Rscript /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/barplots.mapping.R number_reads_mapped.ruby_round2.txt
