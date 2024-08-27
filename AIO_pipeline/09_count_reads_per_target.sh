##########################################
## 8 Number of reads mapping to each target 
##########################################
#################################### 
## Convert target bam files to bed
#################################### 
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles

for file in *pA.bam; do bedtools bamtobed -i $file | awk '{OFS="\t"}{print $1, $2+1, $3,$4,$5,$6}' | sort -k1,1 -k2,2n > $(echo /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/$file | sed 's/.bam$/.bed/'); done &&

## Add class
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles &&
for file in *pA.bed; do awk '{OFS="\t"}NR==FNR{A[$1]=$2;next}{if ($1 in A) print $1,$2,$3,$4,$5,$6,$7,A[$1]}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.txt $file > $(echo $file | sed 's/.bed$/.class.bed/');done &

cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/strandedness

for file in *class.bed; do intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed -b $file -wb -s -nonamecheck | awk '{OFS="\t"}{n=split($4,a,".")}{print $7, $8, $9, $10, $4,a[1],$5, $6, $12,$13}' > $(echo /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/strandedness/$file | sed 's/.bed$/.intersect.sense.bed/');done &
for file in *class.bed; do intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed  -b $file -wb -S -nonamecheck| awk '{OFS="\t"}{n=split($4,a,".")}{print $7, $8, $9, $10, $4,a[1],$5, $6, $12,$13}' > $(echo /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/strandedness/$file | sed 's/.bed$/.intersect.antisense.bed/');done &

#################################### 
## Count number of reads per target
#################################### 
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/strandedness
for file in *.sense.bed ; do cut -f 4,6,7,10 $file | sort -u | cut -f 2,3,4 | sort | uniq -c | awk '{OFS="\t"}{print $1,$2,$3,$4}' > $(echo $file | sed 's/.bed$/.counts.txt/') ; done &
for file in *.antisense.bed ; do cut -f 4,6,7,10 $file | sort -u | cut -f 2,3,4 | sort | uniq -c | awk '{OFS="\t"}{print $1,$2,$3,$4}' > $(echo $file | sed 's/.bed$/.counts.txt/') ; done &

#### 
## Combine all files
#### 
## Generate file with genotypes
ls *.non_pA.class.intersect.antisense.bed | awk '{n=split($1,a,".")}{print a[1]"."a[2]"."a[3]}' > sampleNames.txt

while IFS= read -r line; do Rscript /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/combine_files.R $line; done < sampleNames.txt

## Combine all Genotypes
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/strandedness
echo -e "Gene\tName\tClass\tSense_non_pA\tAnti_non_pA\tSense_pA\tAnti_pA\tTotal\tSample" > all.combined_numReads_mapped_to_targets.strandedness.perpA.round2.txt
for file in *MK*.combined_numReads_mapped_to_targets.strandedness.perpA.txt; do tail -n +2 $file >> all.combined_numReads_mapped_to_targets.strandedness.perpA.round2.txt; done
