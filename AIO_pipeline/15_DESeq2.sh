####################################################
## Create file of psuedo htseq files for input into DESeq
######################################################
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/deseq

for file in *.non_pA.bed; do cut -f1,2,3,6 $file | awk '{OFS=FS="\t"}{print $1"."$2"_"$3".nonpA_"$4}' | sort | uniq -c | awk '{OFS="\t"}{FS=" "}{if ($1 > 3) print $2,$1}' > deseq/$(echo $file | sed 's/.non_pA.bed/.htseq.txt/'); done && 
for file in *.pA.bed; do cut -f1,2,3,6 $file | awk '{OFS=FS="\t"}{print $1"."$2"_"$3".pA_"$4}' | sort | uniq -c | awk '{OFS="\t"}{FS=" "}{if ($1 > 3) print $2,$1}' >> deseq/$(echo $file | sed 's/.pA.bed/.htseq.txt/'); done & 

python /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/collapse_transcripts.py 03_MK016_rep1.35S.1red.rRNA_free.mapped_to_targets.htseq.bed | head

cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/deseq

python /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/collapse_transcripts_full.py  /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/deseq 5 

Rscript /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/deseq_aio_collapsed_transcripts.R
