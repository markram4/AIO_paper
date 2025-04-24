library(ggrepel)
library(Biostrings)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(devtools)
library(here)

setwd("/cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/06_ribowaltz/nuclear_pcg")

## Input annotation and make txdb. Write txdb for plotting purposes
annotation_dt<- riboWaltz::create_annotation(gtfpath ="/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.gtf")
write.table(annotation_dt, file = "/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.gtf.txdb", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## Create list of input files
reads_list <- bamtolist(bamfolder="/cluster/pixstor/slotkinr-lab/mkramer/projects/ruby_ribo_seq/data/06_ribowaltz/nuclear_pcg", annotation = annotation_dt)

## Create list of input names
input_samples <- c("RIBO-MK016_Rep_1_trimmedAligned.toTranscriptome.out", "RIBO-MK040_Rep_1_trimmedAligned.toTranscriptome.out", "RIBO-MK045_Rep_1_trimmedAligned.toTranscriptome.out")

### Generate QC Plots

## Read length distribution
# shows the distribution of reads lengths for one or multiple samples. It can be exploited for i) identifying distinct populations of read lengths, associated with different ribosome conformations and ii) exploring the contribution of each length to the final P-site determination
pdf("ribo_seq_length_distributions_facet.dodge.pdf")
rlength_distr(reads_list,sample = input_samples,multisamples = "independent",plot_style = "dodge",cl = 95, colour = c("#333f50", "#39827c","pink"))
dev.off()


## Plot Read extremity localization 
# consists of four metaheatmaps which displays the abundance of the 5' and 3' extremity of reads mapping on and close to the start and the stop codon of annotated CDSs, stratified by their length.
pdf("ribo_seq_extremity_localization.MK016.pdf")
rends_heat(reads_list, annotation_dt ,sample = "RIBO-MK016_Rep_1_trimmedAligned.toTranscriptome.out", cl = 85,utr5l = 25, cdsl = 40, utr3l = 25, colour = "#333f50")
dev.off()
pdf("ribo_seq_extremity_localization.MK040.pdf")
rends_heat(reads_list, annotation_dt ,sample = "RIBO-MK040_Rep_1_trimmedAligned.toTranscriptome.out", cl = 85,utr5l = 25, cdsl = 40, utr3l = 25, colour = "#333f50")
dev.off()
pdf("ribo_seq_extremity_localization.MK045.pdf")
rends_heat(reads_list, annotation_dt ,sample = "RIBO-MK045_Rep_1_trimmedAligned.toTranscriptome.out", cl = 85,utr5l = 25, cdsl = 40, utr3l = 25, colour = "#333f50")
dev.off()

###-------------------------- Calculate p-site with length filter--------------------------
## Filter for sizes based on QC data
filtered_list <- length_filter(data = reads_list,length_filter_mode = "custom", length_range = 31:34)

## Calculate p-sites
psite_offset <- psite(filtered_list, flanking = 6, extremity = "auto",plot=T)
reads_psite_list <- psite_info(filtered_list, psite_offset)

## Save files for each samples containing p-site information
for (i in seq_along(reads_psite_list)) {
  # Create a dynamic filename for each dataframe in the list
  file_name <- paste0(names(reads_psite_list)[[i]],"_reads_psite_list.sizeFiltered.txt")
  print(file_name)
  # Save each dataframe as a text file
  write.table(reads_psite_list[[i]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

### Generate Plots for p-sites

## P-sites per region
# Ribosome profiling data should define the CDS of transcripts as the region with the highest percentage of reads. To confirm this property the function region_psite computes the percentage of P-sites falling in the three annotated transcript regions (5' UTR, CDS and 3' UTR). To verify the accumulation of reads on the CDS the resulting plot includes an additional set of bars called "RNAs" displaying the expected read distribution from a random fragmentation of RNA. 
input_samples <- list("MK016" = c("RIBO-MK016_Rep_1_trimmedAligned.toTranscriptome.out"),"MK040" = c("RIBO-MK040_Rep_1_trimmedAligned.toTranscriptome.out"),"MK045" = c("RIBO-MK045_Rep_1_trimmedAligned.toTranscriptome.out"))

pdf("ribo_seq_psite_stacked_bar_plot.filtered.pdf")
region_psite(reads_psite_list, annotation_dt,sample = input_samples, multisamples = "average",plot_style = "stack",cl = 85, colour = c("#333f50", "gray70", "#39827c"))
dev.off()

## Trinucleotide periodicity
# A fundamental characteristic of ribosome profiling data is the trinucleotide periodicity of ribosome footprints along the coding sequences. frame_psite show if, and to which extent, the identified P-sites results in codon periodicity on the CDSs. Computes the percentage of P-sites falling in the three possible translation reading frames for the 5’ UTRs, CDSs and 3’ UTRs. 

## Frame of the P-site for the CDS, not stratified by read length.
pdf("ribo_seq_frame_psite.filtered.pdf")
frame_psite(reads_psite_list, annotation_dt, sample = input_samples, multisamples = "independent",plot_style = "facet",region = "all",colour = c("#333f50", "#39827c","pink"))
dev.off()

## Metaplots
input_samples <- list("MK016" = c("RIBO-MK016_Rep_1_trimmedAligned.toTranscriptome.out"),"MK040" = c("RIBO-MK040_Rep_1_trimmedAligned.toTranscriptome.out"),"MK045" = c("RIBO-MK045_Rep_1_trimmedAligned.toTranscriptome.out"))

# Line plots to show periodicity
# A visual representation of the trinucleotide periodicity along the coding sequences is provided by the function metaprofile_psite. It generates metaprofiles by collapsing transcript-specific profiles based on the P-sites mapping on and close to the start and the stop codon of annotated CDSs.

pdf("ribo_seq_metaprofile_psite.filtered.overlap.pdf")
metaprofile_psite(reads_psite_list, annotation_dt,sample = input_samples,multisamples = "average",plot_style = "overlap",utr5l = 20, cdsl = 40, utr3l = 20,colour = c("#333f50", "#39827c","pink"))
dev.off()

