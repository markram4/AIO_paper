## Data processing for small RNA-seq data
1. Sequencing adapters are trimmed
2. All reads are mapped to the genome + any transgene used in the lab previously. The number of mapped reads is used for normalization.
3. Filter reads that don't map to chloroplast and mitochondria, and only look for reads that are 18 - 28 nt in length
4. Map to 35S:RUBY + genome
