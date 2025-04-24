###########
## 02 Map to TAIR10 Genome + RUBY
################
## filter out ChrMt and ChrPt from gtf
awk '{OFS=FS="\t"}{if ($1 != "Mt" && $1 != "Pt") print}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.gtf  > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.gtf


## Filter for only protein-coding genes in gtf file
awk '{OFS=FS="\t"}{if ($0 ~ "#") print; else if ($9 ~ /protein_coding/ || $1 == "35S_RUBY_transgene") print}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.gtf > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.gtf

## Create STAR index
conda activate star_260c
## create genome index
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/TAIR_RUBY_STAR_noPtMt_pcgOnly --genomeFastaFiles /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/TAIR10_nuclear-Chr-all_Araport11.ruby_35s.transgene.noChr.fa --sjdbGTFfile /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.gtf --sjdbOverhang 34 &

## Create slurm script and run
bash /cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/jobFiles/02_map_to_genome.star.se.nuclearOnly.pcgOnly.sh /cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/01_trimmedFastq/samples.txt
