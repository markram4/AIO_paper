Directory structure:
/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2
        data
                2_cutadapt_trim_adapters
                3_cutadapt_demultiplex
                4_cutadapt_trim_5p_3p_adapters
                        3p
                        3p_noTrim
                        5p
                6_minimap
                        1_map_to_rRNA_tRNA
                                logFiles
                                output
                                        mapped_to_rRNA
                                        rRNA_free_fastq
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
│   └── 9_sRNA
│       ├── bowtie
│       │   ├── bamFiles
│       │   │   └── bedFiles
│       │   │       └── sRNA_sequences
│       │   ├── bed12_files
│       │   │   └── *
│       │   │       └── tmp_bwscore_3103638
│       │   ├── bgr_files
│       │   ├── bw_files
│       │   ├── plots
│       │   └── sRNA_counts
│       ├── merge_br
│       │   ├── bamFiles
│       │   │   └── bedFiles
│       │   │       └── bins
│       │   ├── bed12_files
│       │   ├── bgr_files
│       │   ├── bw_files
│       │   └── plots
│       ├── miRNA_search
│       │   └── merge_br
│       │       ├── bed12_files
│       │       ├── bgr_files
│       │       └── bw_files
│       ├── shortStack
│       │   ├── bed12_files
        │   ├── bedFiles
│       │   │   └── sRNA_sequences
│       │   ├── logFiles
│       │   ├── merge_br
│       │   │   ├── bamFiles
│       │   │   └── bedFiles
│       │   │       └── sRNA_sequences
│       │   ├── R1S02_MK017R_BR1.8green_mapped
│       │   ├── R1S03_MK016_BR1.1red_mapped
│       │   ├── R1S04_MK016_BR2.1red_mapped
│       │   ├── R1S05_MK016_BR3.3redParts_mapped
│       │   ├── R1S06_MK016_BR3.4greenParts_mapped
│       │   ├── R1S07_MK016_BR3.5fullRed_mapped
│       │   ├── R1S08_MK016_BR3.6fullGreen_mapped
│       │   ├── R1S09_MK034R_BR1.2reddish_mapped
│       │   ├── R1S10_MK016_BR4.3redParts_mapped
│       │   ├── R1S11_MK016_BR4.4greenParts_mapped
│       │   ├── R1S12_MK016_BR4.5fullRed_mapped
│       │   ├── R1S13_MK016_BR4.6fullGreen_mapped
│       │   ├── R1S14_MK040_BR1.3redParts_mapped
│       │   ├── R1S15_MK040_BR1.4greenParts_mapped
│       │   ├── R1S16_MK040_BR1.5fullRed_mapped
│       │   ├── R1S17_MK040_BR1.6fullGreen_mapped
│       │   ├── R1S18_MK040_BR2.3redParts_mapped
│       │   ├── R1S19_MK040_BR2.4greenParts_mapped
│       │   ├── R1S20_MK040_BR2.5fullRed_mapped
│       │   ├── R1S21_MK040_BR2.6fullGreen_mapped
│       │   ├── R1S22_MK010R_BR1.1red_mapped
│       │   ├── R2S01_MK041A_BR1.3redParts_mapped
│       │   ├── R2S02_MK041A_BR1.4greenParts_mapped
│       │   ├── R2S03_MK041A_BR1.5fullRed_mapped
│       │   ├── R2S04_MK041A_BR2.3redParts_mapped
│       │   ├── R2S05_MK041A_BR2.4greenParts_mapped
│       │   ├── R2S06_MK041A_BR2.5fullRed_mapped
│       │   ├── R2S07_MK042_BR1.1red_mapped
│       │   ├── R2S08_MK042_BR2.1red_mapped
│       │   ├── R2S09_MK042_BR3.3redParts_mapped
│       │   ├── R2S10_MK042_BR3.4greenParts_mapped
│       │   ├── R2S11_MK042_BR3.5fullRed_mapped
│       │   ├── R2S12_MK042_BR3.6fullGreen_mapped
│       │   ├── R2S13_MK043_BR1.3redParts_mapped
│       │   ├── R2S14_MK043_BR1.4greenParts_mapped
│       │   ├── R2S15_MK043_BR1.5fullRed_mapped
│       │   ├── R2S16_MK043_BR1.6fullGreen_mapped
│       │   ├── R2S17_MK043_BR2.3redParts_mapped
│       │   ├── R2S18_MK043_BR2.4greenParts_mapped
│       │   ├── R2S19_MK043_BR2.5fullRed_mapped
│       │   ├── R2S20_MK044_BR1.1red_mapped
│       │   ├── R2S21_MK044_BR2.1red_mapped
│       │   ├── R2S22_MK044_BR3.3redParts_mapped
│       │   ├── R2S23_MK044_BR3.4greenParts_mapped
│       │   ├── R2S24_MK044_BR3.5fullRed_mapped
│       │   ├── R2S25_MK044_BR3.6fullGreen_mapped
│       │   ├── R2S26_MK045_BR1.1red_mapped
│       │   ├── R2S27_MK045_BR2.1red_mapped
│       │   ├── R2S28_MK045_BR3.3redParts_mapped
│       │   ├── R2S29_MK045_BR3.4greenParts_mapped
│       │   ├── R2S30_MK045_BR3.5fullRed_mapped
│       │   ├── R2S31_MK045_BR3.6fullGreen_mapped
│       │   ├── R2S32_MK034_BR2.2reddish_mapped
│       │   ├── R2S33_MK017_BR2.8green_mapped
│       │   ├── R2S34_MK016_BR5.3redParts_mapped
│       │   ├── R2S35_MK016_BR5.4greenParts_mapped
│       │   ├── R2S36_MK016_BR5.5fullRed_mapped
│       │   ├── R2S37_MK016_BR5.6fullGreen_mapped
│       │   ├── R2S38_MK025_BR1.3redParts_mapped
│       │   ├── R2S39_MK025_BR1.4greenParts_mapped
│       │   ├── R2S40_MK025_BR1.5fullRed_mapped
│       │   ├── sRNA_counts
│       │   └── wtCol_pSAE005-T1_B_mapped
├── jobFiles
│   └── logs
└── scripts
