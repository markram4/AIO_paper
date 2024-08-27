####################################################
## Generate bed12 files from merged files to generate coverage plots
######################################################
## Make relevant directories
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bgr_files
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bw_files
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bed12_files

## Convert to bgr
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bamFiles
for sample in *bam; do genomeCoverageBed -ibam $sample -bg -strand + -split | LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bgr_files/$(echo $sample | sed 's/.bam$/.plus.bgr/') ; done &
for sample in *bam; do genomeCoverageBed -ibam $sample -bg -strand - -split | LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bgr_files/$(echo $sample | sed 's/.bam$/.minus.bgr/') ; done &

for sample in *bam; do genomeCoverageBed -ibam $sample -bg -strand + -5 | LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bgr_files/$(echo $sample | sed 's/.bam$/.plus.5p.bgr/') ; done &
for sample in *bam; do genomeCoverageBed -ibam $sample -bg -strand - -5 | LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bgr_files/$(echo $sample | sed 's/.bam$/.minus.5p.bgr/') ; done &

for sample in *bam; do genomeCoverageBed -ibam $sample -bg -strand + -3 | LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bgr_files/$(echo $sample | sed 's/.bam$/.plus.3p.bgr/') ; done &
for sample in *bam; do genomeCoverageBed -ibam $sample -bg -strand - -3 | LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bgr_files/$(echo $sample | sed 's/.bam$/.minus.3p.bgr/') ; done &

## Convert to bw
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bgr_files
for sample in *bgr; do /cluster/pixstor/slotkinr-lab/mkramer/utils/bedGraphToBigWig $sample /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.txt /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bw_files/$(echo $sample | sed 's/.bgr$/.bw/');done &

## Convert bw to bed12
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bw_files
module load bwtool
for file in *plus.bw; do python2.7 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/extract_bigwig_to_bed12.elim_no_coverage.py   -plus $file -minus $(echo $file | sed 's/plus.bw$/minus.bw/') -NAval 0  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.2.bed12 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bed12_files/$(echo $file | sed 's/plus.bw$/bed12/'); done &
for file in *plus.5p.bw; do python2.7 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/extract_bigwig_to_bed12.elim_no_coverage.py -plus $file -minus $(echo $file | sed 's/plus.5p.bw$/minus.5p.bw/') -NAval 0  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.2.bed12 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bed12_files/$(echo $file | sed 's/plus.5p.bw$/5p.bed12/'); done &
for file in *plus.3p.bw; do python2.7 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/extract_bigwig_to_bed12.elim_no_coverage.py -plus $file -minus $(echo $file | sed 's/plus.3p.bw$/minus.3p.bw/') -NAval 0  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.2.bed12 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bed12_files/$(echo $file | sed 's/plus.3p.bw$/3p.bed12/'); done &

# Antisense
for file in *plus.bw; do python2.7 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/extract_bigwig_to_bed12.elim_no_coverage.py   -minus $file -plus $(echo $file | sed 's/plus.bw$/minus.bw/') -NAval 0  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.2.bed12 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bed12_files/$(echo $file | sed 's/plus.bw$/antisense.bed12/'); done &
for file in *plus.5p.bw; do python2.7 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/extract_bigwig_to_bed12.elim_no_coverage.py -minus $file -plus $(echo $file | sed 's/plus.5p.bw$/minus.5p.bw/') -NAval 0  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.2.bed12 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bed12_files/$(echo $file | sed 's/plus.5p.bw$/5p.antisense.bed12/'); done &
for file in *plus.3p.bw; do python2.7 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/extract_bigwig_to_bed12.elim_no_coverage.py -minus $file -plus $(echo $file | sed 's/plus.3p.bw$/minus.3p.bw/') -NAval 0  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.2.bed12 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/8_merge_br/bed12_files/$(echo $file | sed 's/plus.3p.bw$/3p.antisense.bed12/'); done &
