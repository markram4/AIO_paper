####################################################
## Examine coverage at 3' ends at 5' splice sites (Supplemental Figure 6F)
######################################################
# Generate file with position of 3' splice sites
cd ~/annotations/arabidopsis/At_v1

grep "Protein-coding" At_array.v1.targets_class.txt > At_array.v1.targets_class.pcg.txt
awk 'NR==FNR{A[$1]=$1;next}{if ($1 in A) print $0}' At_array.v1.targets_class.pcg.txt At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.gtf > At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.gtf

awk '{FS="\t"}{OFS="\t"}{n=split($9,A,"\"")}  {if ($3 == "exon" && $7 == "-") print $1,$4-25,$4+25, $3"_"A[6], A[4], $7, $4-1; else if ($3 == "exon" && $7 == "+") print $1,$5-25,$5+25, $3"_"A[6], A[4], $7, $5-1}' At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.gtf | awk '{FS="\t"}{OFS="\t"}{if ($2 <0) print $1,"1",$3,$4,$5,$6,$7;else print $0}' > At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.5pss.gtf
awk '{FS="\t"}{OFS="\t"}{n=split($9,A,"\"")}  {if ($3 == "exon" && $7 == "+") print $1,$4-25,$4+25, $3"_"A[6], A[4], $7, $4-1; else if ($3 == "exon" && $7 == "-") print $1,$5-25,$5+25, $3"_"A[6], A[4], $7, $5-1}' At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.gtf | awk '{FS="\t"}{OFS="\t"}{if ($2 <0) print $1,"1",$3,$4,$5,$6,$7;else print $0}'> At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.3pss.gtf

## Add strand to bgr
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files

awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,".","-"}' /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.bgr > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.mod.bed
awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,".","+"}' /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.bgr > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.mod.bed

awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,".","-"}' /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.3p.bgr > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.3p.mod.bed
awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,".","+"}' /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.3p.bgr > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.3p.mod.bed

awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,".","-"}' /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.5p.bgr > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.5p.mod.bed
awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,".","+"}' /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.5p.bgr > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.5p.mod.bed

## merge strands
cat col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.mod.bed col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.mod.bed | sort -k1,1 -k2,2n > col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.mod.bed &

cat col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.3p.mod.bed col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.3p.mod.bed | sort -k1,1 -k2,2n > col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.3p.mod.bed &

cat col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.plus.5p.mod.bed col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.minus.5p.mod.bed | sort -k1,1 -k2,2n > col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.5p.mod.bed &

intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.3p.mod.bed -b /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v1/At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.3pss.gtf  -s -wb -nonamecheck | awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,$6,$8,$9,$10,$11,$13}' > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.pcg.3pends_at_3pss.readCount.txt &

intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.3p.mod.bed -b /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v1/At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.5pss.gtf  -s -wb -nonamecheck | awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,$6,$8,$9,$10,$11,$13}' > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files/col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.pcg.3pends_at_5pss.readCount.txt

## Across the first intron compared to all others (Supplemental Figure 6G)
# Generate file with position of 3' splice sites
cd /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v1

grep "Protein-coding" At_array.v1.targets_class.txt > At_array.v1.targets_class.pcg.txt
awk 'NR==FNR{A[$1]=$1;next}{if ($1 in A) print $0}' At_array.v1.targets_class.pcg.txt At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.gtf > At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.gtf

## Generate bed file with genes
awk '{FS="\t"}{OFS="\t"}{n=split($9,A,"\"")}  {if ($3 == "transcript") print A[4],$4,$5,$1,A[2],$7}' At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.gtf > At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.transcript.bed

## Generate bed file with only exons
awk '{FS="\t"}{OFS="\t"}{n=split($9,A,"\"")}  {if ($3 == "exon") print A[4],$4,$5,$1,$3"_"A[6],$7}' At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.gtf > At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.exons.bed

## Generate subtract exons from genes to get introns
bedtools subtract -s -a At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.transcript.bed -b At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.exons.bed | awk '{FS="\t"}{OFS="\t"}{print $1,$2+1,$3-1,$4,$5,$6}' > At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.introns.bed 


# Add intron # using R
>R
>library(dplyr)
>introns <- read.table("At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.introns.bed",col.names=c("txt","start","stop","chr","gene","strand"))
>t <- introns %>% dplyr::group_by(txt) %>% dplyr::mutate(intron = case_when(strand == "-" ~ paste("intron",row_number(desc(start)),sep="_"), strand == "+" ~ paste("intron",row_number(start),sep="_"))) %>% relocate(chr,start,stop,txt,intron,strand) %>% select(!gene)
>write.table(t, file = "At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.introns.mod.bed",row.names=F,col.names=F,quote=F,sep="\t")

## Intersect bedgraph files with intron bed file
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/8_merge_br/bgr_files

intersectBed -a col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.mod.bed -b /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v1/At_array.v1.Arabidopsis_thaliana.TAIR10.56.targetChr.all.pcg.introns.mod.bed -s -wb -nonamecheck | awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,$6,$8,$9,$10,$11}' > col0_rep123.non_pA.filtered.rRNA_free.mapped_to_targets.cat.pcg.introns.readCount.txt



## Examine read length for SimpleHat in Col0, ddm1, polv, polivpolv (Figure 4E)
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles
mkdir simpleHatLengths

for file in *targets.bed; do awk '{FS="\t"}{OFS="\t"}{n=split(FILENAME,a,".")}{if ($1 ~ "AT3TE70480") print $3-$2, $1, a[1]"."a[2]}' $file > simpleHatLengths/$(echo $file | sed 's/.bed$/.simpleHatLen.txt/'); done

cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/at_round3/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/simpleHatLengths
cat *simpleHatLen.txt > all_samples.readSize.mappedToSimpleHat.at_round3.txt

