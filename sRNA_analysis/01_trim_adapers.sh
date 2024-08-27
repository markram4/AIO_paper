############################
## 03 Trim Illumina Universal adapter
############################
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/01_fastqFiles/
ls *fq.gz | sed 's/.fq.gz//' > samples.txt

sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/02_trim_adapters.array.sh
