####################################################
## Quantify the number of 3' ends across RUBY between phenotypes
######################################################
awk '{OFS="\t"}NR==FNR{A[$1]=$6;next}{if ($1 in A) print $1,$2,A[$1]}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.txt > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.strand.txt 

# Normalize bgr files so the nt with the highest read coverage = 1
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/7_coverage_plots/bgr_files
mkdir number_3p_ends
for file in *non_pA.plus.3p.bgr; do Rscript /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/number_3p_ends.ruby.R /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/7_coverage_plots/bgr_files/number_3p_ends $(echo $file | sed 's/.rRNA_free.mapped_to_targets.non_pA.plus.3p.bgr//') ; done &

## Combine all Genotypes
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/7_coverage_plots/bgr_files/number_3p_ends
head -n 1 01_MK041A_rep1.35S.3redParts.norm_count.txt > all.combined_normalized_3p_read_count.ruby_round2.txt &&
for file in *MK*txt; do tail -n +2 $file >> all.combined_normalized_3p_read_count.ruby_round2.txt; done
