################
## Run Kallisto for QC correlation plots
################
conda activate kallisto_env

## Create fasta file with sequence of protein-coding genes used in mapping
/cluster/pixstor/slotkinr-lab/mkramer/utils/gffread/gffread /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.gtf -g /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/TAIR10_nuclear-Chr-all_Araport11.ruby_35s.transgene.noChr.fa -w /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.transcriptome.fasta

# build an index from a FASTA formatted file of target sequences
kallisto index -i /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.transcriptome.kallisto.idx /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.transcriptome.fasta -k 19 &

# Run kallisto to count number of reads mapping to each transcript for input into something to create correlation plot
kallisto quant -i /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.transcriptome.kallisto.idx -o /cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/01_trimmedFastq/kallisto/RIBO-MK016_Rep_1_trimmed -t 10 --single -l 28 -s 2 RIBO-MK016_Rep_1_trimmed.fastq.gz &
kallisto quant -i /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.transcriptome.kallisto.idx -o /cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/01_trimmedFastq/kallisto/RIBO-MK040_Rep_1_trimmed -t 10 --single -l 28 -s 2 RIBO-MK040_Rep_1_trimmed.fastq.gz &
kallisto quant -i /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.transcriptome.kallisto.idx -o /cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/01_trimmedFastq/kallisto/RIBO-MK045_Rep_1_trimmed -t 10 --single -l 28 -s 2 RIBO-MK045_Rep_1_trimmed.fastq.gz &

## Create file with transcript and gene name to convert kallisto output for input into DESeq
awk '{FS=OFS="\t"}{n=split($9,A,"\"")}{if ($3 == "transcript") print A[4],A[2]} ' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.gtf |sort -u > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.txt_to_gene.txt
