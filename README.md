# AIO_paper_analyis

* This repo contains the AIO paper pipelines ,analysis and ML scripts


## Folder structure

* AIO_pipeline
Directory structure:
/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2
├── data
│   ├── 2_cutadapt_trim_adapters
│   │   └── logFiles
│   ├── 3_cutadapt_demultiplex
│   │   ├── barcodes
│   │   └── logFiles
│   ├── 4_cutadapt_trim_5p_3p_adapters
│   │   ├── 3p
│   │   │   └── logFiles
│   │   ├── 3p_noTrim
│   │   │   └── logFiles
│   │   └── 5p
│   │       └── logFiles
│   ├── 6_minimap
│   │   ├── 1_map_to_rRNA_tRNA
│   │   │   ├── logFiles
│   │   │   └── output
│   │   │       ├── mapped_to_rRNA
│   │   │       └── rRNA_free_fastq
│   │   ├── 2_map_to_capture
│   │   │   ├── logFiles
│   │   │   └── output
│   │   │       ├── mapped_to_targets
│   │   │       │   └── bedFiles
│   │   │       │       ├── around_3343
│   │   │       │       ├── deseq
│   │   │       │       ├── full_length
│   │   │       │       │   └── umi
│   │   │       │       │       ├── 01_readNames
│   │   │       │       │           └── [directory for each sample -- removed for clarity]
│   │   │       │       │       ├── 02_subset_seq
│   │   │       │       │           └── [directory for each sample -- removed for clarity]
│   │   │       │       │       ├── 03_extract
│   │   │       │       │           └── [directory for each sample -- removed for clarity]
│   │   │       │       │       ├── 04_filtered
│   │   │       │       │           └── [directory for each sample -- removed for clarity]
│   │   │       │       │       └── 05_clustered_cd
│   │   │       │       │           └── [directory for each sample -- removed for clarity]
│   │   │       │       ├── strandedness
│   │   │       │       └── umi
│   │   │       │           ├── 01_readNames
│   │   │       │               └── [directory for each sample -- removed for clarity]
│   │   │       │           ├── 02_subset_seq
│   │   │       │               └── [directory for each sample -- removed for clarity]
│   │   │       │           ├── 03_extract
│   │   │       │               └── [directory for each sample -- removed for clarity]
│   │   │       │           ├── 04_filtered
│   │   │       │               └── [directory for each sample -- removed for clarity]
│   │   │       │           └── 05_clustered_cd
│   │   │       │               └── [directory for each sample -- removed for clarity]
│   │   │       └── non_targets_fastq
│   │   ├── 3_map_to_genome
│   │   │   ├── logFiles
│   │   │   └── output
│   │   │       ├── mapped_to_genome
│   │   │       └── unmapped
│   │   └── mapped_pct
│   ├── 7_coverage_plots
│   │   ├── bed12_files
│   │   ├── bgr_files
│   │   │   └── number_3p_ends
│   │   └── bw_files
│   ├── 8_merge_br
│   │   ├── bamFiles
│   │   ├── bed12_files
│   │   ├── bgr_files
│   │   ├── bw_files
│   │   └── plots
├── jobFiles
│   └── logs
└── scripts

      
* sRNA_analysis
