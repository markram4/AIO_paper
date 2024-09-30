library(plyr)
library(tidyverse)
library(dplyr)

#setwd("/Users/mariannekramer/Desktop/mkramer_server/projects/target_capture/at_round2/data/6_minimap/2_map_to_capture/samFiles/mapped_to
_targets/bedFiles/strandedness/")
args <- commandArgs(T)

sample <- args[1] #"1_col0"

non_pA.sense <- read.table(paste(sample,"rRNA_free.mapped_to_targets.non_pA.class.intersect.sense.counts.txt",sep="."),col.names = c("Sen
se_non_pA","Gene","Name", "Class"))
non_pA.antisense <- read.table(paste(sample,"rRNA_free.mapped_to_targets.non_pA.class.intersect.antisense.counts.txt",sep="."),col.names 
= c("Anti_non_pA","Gene","Name","Class"))
pA.sense <- read.table(paste(sample,"rRNA_free.mapped_to_targets.pA.class.intersect.sense.counts.txt",sep="."),col.names = c("Sense_pA","
Gene","Name","Class"))
pA.antisense <- read.table(paste(sample,"rRNA_free.mapped_to_targets.pA.class.intersect.antisense.counts.txt",sep="."),col.names = c("Ant
i_pA","Gene","Name","Class"))

a <- non_pA.sense %>% select(Gene, Name, Class, Sense_non_pA) %>%
  full_join(non_pA.antisense,by=c("Gene","Name","Class")) %>% 
  full_join(pA.sense,by=c("Gene","Name","Class")) %>%
  full_join(pA.antisense,by=c("Gene","Name","Class")) %>% rowwise() %>%
  dplyr::mutate(Total = sum(Sense_non_pA, Anti_non_pA, Sense_pA, Anti_pA,na.rm = TRUE)) %>% 
  mutate(Sample = sample) %>% mutate_if(is.numeric, ~replace_na(., 0))

outFile <- paste(sample,"combined_numReads_mapped_to_targets.strandedness.perpA.txt",sep=".")
write.table(a,outFile,sep="\t",row.names=FALSE,quote=FALSE)
