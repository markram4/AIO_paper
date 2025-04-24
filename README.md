# AIO_paper_analyis

* This repo contains the scripts for "Identification of an aberrant RNA associated with the initiation of transgene silencing" by Marianne C Kramer et al. 
* This includes all data analysis steps for "All-in-One" RNA-seq, Random Forest machine learning, small RNA sequencing, and RMarkdown files for each figure in the manuscript. 

## Content

### AIO_pipeline
* commonFiles -- Files required for data analysis that are common among all experiments. Also in this directory are the barcodes used for one example experiment (ruby_round2). Barcodes for all other experiments in this manuscript can be found in Supplemental Table 1.
* jobFiles -- jobFiles submitted for each processing step referenced in the .sh files
* .sh files are shell scripts for each data processing step outlined in **AIO_pipeline_readme.md** This includes running job files, QC steps, and creating directory structures to fit the analysis

### Figures
* .Rmd and .md files used to generate all figures in the manuscript, that were created in R. Some panels are missing because they were not created in R (such as phenotyping images).
* There are also sub folders containing the images required to embed in .md
   
### Random-Forest Classification modeling

* Random forest binary classification modeling using sequence features to diffrentiate phenotypes 
  
### small RNA analysis
* jobFiles -- shell scripts or slurm scripts used for data analysis, references in sRNA_readsMe.md
* sRNA_readMe.sh -- shell script of commands ran for sRNA-seq data processing
* examine_known_mirna_sites.sh -- shell script of commands ran to examine sRNA coverage at known miRNA targets -- Supplemental Figure 13.

* ### Ribo-seq analysis
* ribo_seq_readMe.sh -- shell script of commands ran for sRNA-seq data processing
* .sh -- shell scripts of commands ran to process ribosome profiling sequencing for input into riboWaltz and DESeq2
* R scripts for DESeq2 and riboWaltz

