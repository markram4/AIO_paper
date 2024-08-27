################################
## Examine sRNA build up at known miRNA targets
################################
## Download names and sequences of miRNAs from miRBase for A.thaliana and P.patens - since the dcl samples were spiked with P.patens, ideally, I don't want to look at miRNAs that are conserved in P.patens to avoid any confusion

# Create file with unique miRNA names in P patens
grep ">" /cluster/pixstor/slotkinr-lab/mkramer/annotations/ppatens/mirbase_ppatens.fa | awk -F'[- ]' '{gsub(/[a-zA-Z]$/, "", $2); print $2}' | sort -u  > /cluster/pixstor/slotkinr-lab/mkramer/annotations/ppatens/mirbase_ppatens.uniq.names.txt

# Create file with unique miRNA names in Arabidopsis
grep ">" /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.fa | awk -F'[- ]' '{gsub(/[a-zA-Z]$/, "", $2); print $2}' | sort -u  > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.txt

# Pull out miRNAs that are NOT conserved in P. patens
grep -v -f /cluster/pixstor/slotkinr-lab/mkramer/annotations/ppatens/mirbase_ppatens.uniq.names.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.txt | sort > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.no_ppt.txt

## Download miRNA known targets from TarDB : http://www.biosequencing.cn/TarDB/download.html 
cd /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis
# Pull out the gene names for miRNA targets
grep -i -f mirbase_athaliana.uniq.names.no_ppt.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath/ath.deg | awk -F'[\t]' '{print $2}' | sort -u >  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.no_ppt.targets.ath_deg.txt 

## Create fasta file from weird format downloaded from TarDB
python /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath/get_fasta_from_files.py ath.deg > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath/ath.deg.fasta

## Extract cDNA sequences for genes of interest
cd /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis
module load seqtk

seqtk subseq /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_seq_20160703.fa  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.no_ppt.targets.ath_deg.txt > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.no_ppt.targets.ath_deg.Araport11_genes.201606.transcript.fasta  

## Extract cDNA sequences for protein coding genes from target capture as control loci
awk '{OFS=FS="\t"}{n=split($1,A,"_")}{print A[1]}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.pcg.txt | grep -f - /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_seq_20160703.fa | awk '{FS=" "}{n=split($1,A,">")}{print A[2]}'  > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/headers.tmp.txt &&
seqtk subseq /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_seq_20160703.fa  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/headers.tmp.txt > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_array.v1.targets_class.pcg.Araport11_genes.201606.transcript.fasta && 
rm /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/headers.tmp.txt

## Combine miRNA targets, pcg target no miRNA site control, and RUBY fasta files
cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.no_ppt.targets.ath_deg.Araport11_genes.201606.transcript.fasta  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_array.v1.targets_class.pcg.Araport11_genes.201606.transcript.fasta /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/ruby_35s.transgene.fa >  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.fa

## Map sRNA files to this new merged fasta file, allowing for some mismatches and multimappers
bowtie-build /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.fa /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/bowtie_indexes/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.fa 

sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/jobFiles/09_miRNA_search_bowtie.slurm.sh /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/06_rRNA_tRNA_free/ruby_samples.txt
#bowtie -p 8 -v 3 -a --best --strata /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/bowtie_indexes/Araport11_genes.201606.cdna.miRNA_targets.pcg_v1.35S_RUBY.fa /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/06_rRNA_tRNA_free/R1S03_MK016_BR1.1red_trimmed.mito_chloroFree.sizeFilt.rRNA_tRNA_free.fasta -fS 2> tmp.txt | samtools view -h -b -F 4 - > tmp.bam


## Map target sites to their own genes to get locations of where the miRNAs should bind
bowtie  -q -v 0 -a --best --strata /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/bowtie_indexes/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.fa  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath/ath.deg.fasta -fS 2> tmp.txt | samtools view -h -b -F4 - | bamToBed -i -  | sort -k1,1 -k2,2n | awk '{FS=OFS="\t"}{n=split($4,A,".")}{print $1,$2,$3,A[1],$5,$6}' | sort -u > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath/ath.deg.mirna_coord_in_targets.bed

## Get gtf file with relative coordinates for plotting, and invert strandedness if it's - strand for miRNA targets

 python /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/get_rel_coords_mir_targets.py /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.no_ppt.targets.ath_deg.txt  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.gtf

## Get gtf file with relative coordinates for plotting, and invert strandedness if it's - strand for targets
awk '{OFS=FS="\t"}{print $1".1"}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.pcg.txt > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.pcg.transcript.txt
 python /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/get_rel_coords_mir_targets.py /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.pcg.transcript.txt  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.gtf


## Combine miRNA targets and AIO targets rel gtf
cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/mirbase_athaliana.uniq.names.no_ppt.targets.ath_deg.Araport11_GTF_genes_transposons.current_relCoords.gtf /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.pcg.transcript.Araport11_GTF_genes_transposons.current_relCoords.gtf /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/features_from_35S_ruby_transgene.gtf > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.gtf

## Get chr lengths 
awk '{OFS="\t"}{if (substr($1,1,1)==">") if (NR>1) printf "\n%s ", substr($1,2,length($1)-1); else printf "%s ", substr($1,2,length($1)-1); else printf "%s", $0}END{printf "\n"}'  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.fa | awk 'BEGIN {FS=" "}{OFS="\t"}{print $1,length($2)}' >  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.len.txt

## Create bed12 file 
awk '{FS=OFS="\t"}{print $1,"1",$2,$1,$1,"+","1",$2,".","0,","0,","0,"}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.len.txt > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.len.bed12

## Merge bio reps
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search

samtools merge /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/MK016_rep12.35S.1red.sRNA.miRNA_search.bam /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/*MK016_*.1red*bam &
samtools merge /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/MK016_rep345.35S.3redParts.sRNA.miRNA_search.bam /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/*MK016_*.3redParts*bam &
samtools merge /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/MK016_rep345.35S.4greenParts.sRNA.miRNA_search.bam /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/*MK016_*.4greenParts*bam &
samtools merge /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/MK016_rep345.35S.5fullRed.sRNA.miRNA_search.bam /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/*MK016_*.5fullRed*bam &
samtools merge /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/MK016_rep345.35S.6fullGreen.sRNA.miRNA_search.bam /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/*MK016_*.6fullGreen*bam &

## Get values for normalization
grep MK016 /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R25/data/04_mapUnfiltered/number_mapped_reads.R25.txt

R1S03_MK016_BR1.1red	15450923	27783052
R1S04_MK016_BR2.1red	12332129	
R1S05_MK016_BR3.3redParts	17594923	38240100
R1S10_MK016_BR4.3redParts	12666801	
R2S34_MK016_BR5.3redParts	7978376	
R1S06_MK016_BR3.4greenParts	17091524	50149447
R1S11_MK016_BR4.4greenParts	11822929	
R2S35_MK016_BR5.4greenParts	21234994	
R1S07_MK016_BR3.5fullRed	15635649	51400623
R1S12_MK016_BR4.5fullRed	14186087	
R2S36_MK016_BR5.5fullRed	21578887	
R1S08_MK016_BR3.6fullGreen	11663425	37354452
R1S13_MK016_BR4.6fullGreen	12875692	
R2S37_MK016_BR5.6fullGreen	12815335	

## Create bedgraph files
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search
genomeCoverageBed -ibam MK016_rep12.35S.1red.sRNA.miRNA_search.bam -bg -strand + -scale 27.783052|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep12.35S.1red.sRNA.miRNA_search.plus.bgr &
genomeCoverageBed -ibam MK016_rep12.35S.1red.sRNA.miRNA_search.bam -bg -strand - -scale 27.783052|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep12.35S.1red.sRNA.miRNA_search.minus.bgr &

genomeCoverageBed -ibam MK016_rep345.35S.3redParts.sRNA.miRNA_search.bam -bg -strand + -scale 38.240100|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep345.35S.3redParts.sRNA.miRNA_search.plus.bgr &
genomeCoverageBed -ibam MK016_rep345.35S.3redParts.sRNA.miRNA_search.bam -bg -strand - -scale 38.240100|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep345.35S.3redParts.sRNA.miRNA_search.minus.bgr &

genomeCoverageBed -ibam MK016_rep345.35S.4greenParts.sRNA.miRNA_search.bam -bg -strand + -scale 50.149447|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep345.35S.4greenParts.sRNA.miRNA_search.plus.bgr &
genomeCoverageBed -ibam MK016_rep345.35S.4greenParts.sRNA.miRNA_search.bam -bg -strand - -scale 50.149447|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep345.35S.4greenParts.sRNA.miRNA_search.minus.bgr &

genomeCoverageBed -ibam MK016_rep345.35S.5fullRed.sRNA.miRNA_search.bam -bg -strand + -scale 51.400623|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep345.35S.5fullRed.sRNA.miRNA_search.plus.bgr &
genomeCoverageBed -ibam MK016_rep345.35S.5fullRed.sRNA.miRNA_search.bam -bg -strand - -scale 51.400623|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep345.35S.5fullRed.sRNA.miRNA_search.minus.bgr &

genomeCoverageBed -ibam MK016_rep345.35S.6fullGreen.sRNA.miRNA_search.bam -bg -strand + -scale 37.354452|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep345.35S.6fullGreen.sRNA.miRNA_search.plus.bgr &
genomeCoverageBed -ibam MK016_rep345.35S.6fullGreen.sRNA.miRNA_search.bam -bg -strand - -scale 37.354452|  LC_COLLATE=C sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files/MK016_rep345.35S.6fullGreen.sRNA.miRNA_search.minus.bgr &



## Convert to bw
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bgr_files

for sample in *bgr; do /cluster/pixstor/slotkinr-lab/mkramer/utils/bedGraphToBigWig $sample /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.len.txt /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bw_files/$(echo $sample | sed 's/.bgr$/.bw/');done &

## Make into bed12
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bw_files
module load bwtool

for file in *plus.bw; do python2.7 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/extract_bigwig_to_bed12.elim_no_coverage.py -plus $file -minus $(echo $file | sed 's/plus.bw$/minus.bw/') -NAval 0  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.len.bed12 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bed12_files/$(echo $file | sed 's/plus.bw$/plus.bed12/'); done &
for file in *plus.bw; do python2.7 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/extract_bigwig_to_bed12.elim_no_coverage.py -minus $file -plus $(echo $file | sed 's/plus.bw$/minus.bw/') -NAval 0  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/Araport11_genes.201606.transcript.miRNA_targets.pcg_v2.35S_RUBY.len.bed12 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br/bed12_files/$(echo $file | sed 's/plus.bw$/minus.bed12/'); done &


## Intersect merged bam files with miRNA sites to get the best candidates to make plots from
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/9_sRNA/miRNA_search/merge_br
for file in *bam; do intersectBed -a $file -b /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ath/ath.deg.mirna_coord_in_targets.bed -f 1 -S -bed -wb -r | cut -f 1,16 | sort | uniq -c |awk '{FS=OFS="\t"}{print $1,$2,$3}' > $(echo $file | sed 's/bam/counts.txt/'); done &
