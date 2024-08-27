############################
## 05 Filter for reads that don't map to chloroplast and mitochondria, and for size 18-28
############################
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/05_sizeFiltered
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/04_mapUnfiltered

#---------------- RUN COMMAND -------------------------
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/04_sizeFilter_rmChlMito.slurm.sh
#------------------------------------------------------

###########################
## 06 Remove sRNAs from rRNA and tRNA using bowtie
############################
mkdir /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/data/06_rRNA_tRNA_free

#---------------- RUN COMMAND -------------------------
sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/sRNA/R26/jobFiles/05_remove_rRNA_tRNA.slurm.array.sh
#------------------------------------------------------
