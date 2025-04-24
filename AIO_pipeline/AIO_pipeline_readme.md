## 'All-in-One' Data Analysis Pipeline
* The general directory structure for AIO data analysis is : 
* Home dir: /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2
* data
  - 2_cutadapt_trim_adapters
  - 3_cutadapt_demultiplex
  - 4_cutadapt_trim_5p_3p_adapters
  - 5_minimap
  - 6_coverage_plots
* jobFiles
* scripts

  
1. After sequencing, raw fast5 files are base called using guppy basecaller
2. Base called fastq files are concatenated and then demultiplexed using a two-step process. First the 5' and 3' sequencing adapters are removed, then the barcodes are examined and files are separated based on the barcode sequences.
3. Reads are then oriented based on the sequence of the 3' adapter, then the adapter is either trimmed (for mapping) or untrimmed (for UMI examination). The 5' adapter is then removed
4. Reads are mapped to rRNA >
5. targets >
6. genome
7. Reads mapped to targets are then separated into polyA- and polyA+ reads
8. The number / percentage of reads mapping to each class of genes (rRNA/tRNA, targets, genome, unmapped) are calculated
9. Count the number of reads that map to each target for each sample 
10. Count the number of UMIs present
11. Count the number of full-length polyA+ reads
12. Count the number of UMIs for full length polyA+ reads present
13. Create files for coverage plots. This involved bam > bgr > bw > bed12. Rscripts can be found in Figures sub folder
14. Count the number of 3' ends
15. Examine the number reads ending at 5' splice sites and number of reads covering introns
16. RIBO-seq pre-processing
17. riboWaltz
