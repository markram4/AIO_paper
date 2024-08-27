######################################################
## Full length, polyadenylated reads
######################################################
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/full_length

## Intersect reads with annotation and report reads that start before or at the start codon AND after or at the stop codon
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles

for file in *class.bed; do intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed12 -b $file -wb -s -nonamecheck | awk '{OFS="\t"}{n=split($4,a,".")}{if ($17 >= $8 && $16 <= $7) print $1,$2,$3,$4,a[1],$5,$6,$7,$8,$9,$10,$11,$12,$21,$14,$18}' > full_length/$(echo $file | sed 's/class.bed$/sense.full_length.bed12/'); done &
for file in *class.bed; do intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed12 -b $file -wb -S -nonamecheck | awk '{OFS="\t"}{n=split($4,a,".")}{if ($17 >= $8 && $16 <= $7) print $1,$2,$3,$4,a[1],$5,$6,$7,$8,$9,$10,$11,$12,$21,$14,$18}' > full_length/$(echo $file | sed 's/class.bed$/antisense.full_length.bed12/'); done &

#### 
## Count number of reads per target
#### 
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/full_length
for file in *.bed12 ; do cut -f 5,14,15,16 $file | sort -u | cut -f 1,2,3 | sort | uniq -c | awk '{OFS="\t"}{print $1,$2,$4,$3}' > $(echo $file | sed 's/.bed12$/.counts.txt/') ; done &

#### 
## Combine all reads reads
#### 
## Generate file with genotypes
ls *.non_pA*.sense*bed12 | awk '{n=split($1,a,".")}{print a[1]"."a[2]"."a[3]}' > sampleNames.txt

while IFS= read -r line; do Rscript /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/combine_files.full_length.R $line; done < sampleNames.txt

## Combine all Genotypes
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/full_length
echo -e "Gene\tName\tClass\tSense_non_pA\tAnti_non_pA\tSense_pA\tAnti_pA\tTotal\tSample" > all.combined_numReads_mapped_to_targets.strandedness.perpA.round2.fl.txt
for file in *MK*.combined_numReads_mapped_to_targets.strandedness.perpA.txt; do tail -n +2 $file >> all.combined_numReads_mapped_to_targets.strandedness.perpA.round2.fl.txt; done
