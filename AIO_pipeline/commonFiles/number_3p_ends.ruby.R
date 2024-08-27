# List of package names to check, install, and load
package_names <- c("tidyverse", "ggpubr")

# Function to check if package is installed, install if not, and load
check_install_load_package <- function(pkg) {
  # Set CRAN mirror
  options(repos = "https://cran.rstudio.com/")
  
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Loop through each package in the list
for (pkg in package_names) {
  cat("Checking and loading package:", pkg, "\n")
  check_install_load_package(pkg)
}


args <- commandArgs(T)
#setwd("/Volumes/research/slotkin_lab/Marianne/RUBY/Target capture/round1/data/")
## Get file name, extract transgene and sample number
outDir = args[1] #"/Users/mariannekramer/Desktop/mkramer_server/projects/target_capture/ruby/data/7_coverage_plots/bgr_files/"
sample = args[2] #"02_MK017R.2pooled.8green"

inClass <- read.table("/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.targets_class.strand.txt",col.names=c("name","class","strand"))
senseTargets <- inClass %>% filter(strand == "+")
antiTargets <- inClass %>% filter(strand == "-")

## non_pA
data.plus <- read.table(paste(sample,".rRNA_free.mapped_to_targets.non_pA.plus.3p.bgr",sep=""),col.names = c("name","start", "stop","count"))
data.minus <- read.table(paste(sample,".rRNA_free.mapped_to_targets.non_pA.minus.3p.bgr",sep=""),col.names = c("name","start", "stop","count"))

outData.plus <- data.plus  %>% inner_join(senseTargets, by="name",keep = F) %>% group_by(name) %>%
  mutate(maximum = max(count,na.rm = TRUE), normScore = count/maximum,pheno = sample) 

outData.minus <- data.minus  %>% inner_join(antiTargets, by="name",keep = F) %>% group_by(name) %>%
  mutate(maximum = max(count,na.rm = TRUE), normScore = count/maximum,pheno = sample) 

outData <- rbind(outData.plus,outData.minus)
outFile = paste(outDir,"/",sample,".norm_count.txt",sep="")

write.table(outData, file = outFile, quote=F, row.names = F, col.names=T, sep= "\t")
