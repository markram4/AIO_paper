## 'All-in-One' Data Analysis Pipeline
1. After sequencing, raw fast5 files are base called using guppy basecaller
2. Base called fastq files are concatenated and then demultiplexed using a two-step process. First the 5' and 3' sequencing adapters are removed, then the barcodes are examined and files are separated based on the barcode sequences.
4. 
