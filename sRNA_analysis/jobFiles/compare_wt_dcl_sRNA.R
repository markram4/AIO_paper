library(ggplot2)
library(ggVennDiagram)
library("ggvenn")
library(tidyverse)
library(dplyr)
library(eulerr)
options(bedtools.path = "/usr/local/bin")
library(bedtoolsr)
library(gggenes)
library(cowplot)
library(ggpubr)
library(scales)

setwd("/Users/mariannekramer/Google Drive/Kramer_et_al_AIO/Figures/ruby_dcl_vs_wt")

wt.red <- c(list.files(path = "/Users/mariannekramer/Google Drive/Kramer_et_al_AIO/Figures/ruby_dcl_vs_wt", 
                       pattern = "*35S.1red.trimmed.mito_chloroFree.sizeFilt.rRNA_tRNA_free.final.RUBY.bed",full.names=T),
            list.files(path = "/Users/mariannekramer/Google Drive/Kramer_et_al_AIO/Figures/ruby_dcl_vs_wt", 
                       pattern = "*35S.2reddish.trimmed.mito_chloroFree.sizeFilt.rRNA_tRNA_free.final.RUBY.bed",full.names=T))

in.wt.red <- wt.red %>% 
  map_dfr(~read.table(.,col.names = c("transgene","start","stop","read","seq","strand")))  %>% select(!read) %>%
  mutate(size=nchar(seq)) %>% tidyr::unite("seq",start:size,sep="_") 

wt.fullRed <- c(list.files(path = "/Users/mariannekramer/Google Drive/Kramer_et_al_AIO/Figures/ruby_dcl_vs_wt", 
                       pattern = "*35S.5fullRed.trimmed.mito_chloroFree.sizeFilt.rRNA_tRNA_free.final.RUBY.bed",full.names=T))
in.wt.fullRed <- wt.fullRed %>% 
  map_dfr(~read.table(.,col.names = c("transgene","start","stop","read","seq","strand")))  %>% select(!read) %>%
  mutate(size=nchar(seq)) %>% tidyr::unite("seq",start:size,sep="_")

dcl.red <- c(list.files(path = "/Users/mariannekramer/Google Drive/Kramer_et_al_AIO/Figures/ruby_dcl_vs_wt", 
                       pattern = "*35S_dcl1234.1red.sRNA.filtered.rRNA_tRNA_free.final.RUBY.bed",full.names=T))
in.dcl.red <- dcl.red %>% 
  map_dfr(~read.table(.,col.names = c("transgene","start","stop","read","seq","strand")))  %>% select(!read) %>%
  mutate(size=nchar(seq)) %>% tidyr::unite("seq",start:size,sep="_") 

## Compare dcl and Wt Always Red
# List containing the vectors
wt_red_vs_dcl <- list(
  dcl = in.dcl.red$seq,
  wt.red = in.wt.red$seq 
)


# Function to filter sequences appearing more than 2 times
filter_sequences <- function(seq_list) {
  counts <- table(unlist(seq_list))
  filtered_items <- names(counts[counts > 2])
  lapply(seq_list, function(seq) seq[seq %in% filtered_items])
}

# Apply the function to the list x
filtered_wt_red_vs_dcl <- filter_sequences(wt_red_vs_dcl)

wt_red_vs_dcl.plot <- ggvenn(filtered_wt_red_vs_dcl, fill_color=c("#853061","#B03060"),
                        fill_alpha = 0.8,auto_scale = T,text_size = 3,set_name_size = 3,
                        stroke_size = 0.5) 
ggsave("wt_red_vs_dcl_sRNA.venn.filter2.pdf",plot=wt_red_vs_dcl.plot)

##
themes <- theme(plot.title = element_text(hjust = 0.5),axis.ticks.length=unit(0.0516,"in"),
                axis.text = element_text(size=8,color="black"),
                axis.title.x = element_blank(),
                line = element_line(color = 'black',linewidth=0.3,lineend="round"),
                axis.title = element_text(color = "black"),
                strip.text = element_text(color = "black",size=8),
                legend.position = 'right')

# Get unique sequences in each filtered list
unique_dcl <- unique(filtered_wt_red_vs_dcl$dcl)
unique_wt_red <- unique(filtered_wt_red_vs_dcl$wt.red)

# Get unique sequences in dcl but not in wt.red
dcl.specific.list <- setdiff(unique_dcl, unique_wt_red)
dcl.specific <- as.data.frame(dcl.specific.list) %>% separate(dcl.specific.list,into = c("start","stop","seq","strand","length"),sep="_") %>% mutate(gene = "35S_RUBY_transgene",class = "dcl-specific") %>% relocate(gene,start,stop,seq,class,strand,length)
write.table(dcl.specific,file = "wt_red_vs_dcl.dcl_specific_sRNAs.txt",col.names = T, row.names = F,quote = F,sep="\t")

# Get unique sequences in wt but not in dcl
wt.specific.list <- setdiff(unique_wt_red,unique_dcl)
wt.specific <- as.data.frame(wt.specific.list) %>% separate(wt.specific.list,into = c("start","stop","seq","strand","length"),sep="_") %>% mutate(gene = "35S_RUBY_transgene",class = "wt-specific")%>% relocate(gene,start,stop,seq,class,strand,length)
write.table(wt.specific,file = "wt_red_vs_dcl.wt_specific_sRNAs.txt",col.names = T, row.names = F,quote = F,sep="\t")

# Get unique sequences present in both dcl and wt.red
com.specific.list <- intersect(unique_wt_red,unique_dcl)
common.srna <- as.data.frame(com.specific.list) %>% separate(com.specific.list,into = c("start","stop","seq","strand","length"),sep="_") %>% mutate(gene = "35S_RUBY_transgene",class = "common")%>% relocate(gene,start,stop,seq,class,strand,length)
write.table(common.srna,file = "wt_red_vs_dcl.shared_sRNAs.txt",col.names = T, row.names = F,quote = F,sep="\t")

wt_red_vs_dcl.toPlot.specifics <- rbind(wt.specific,dcl.specific,common.srna) %>% select(!seq) %>% group_by(strand,length,class) %>%
  dplyr::summarize(count=n()) %>% mutate(count = case_when(strand == "-" ~ count *-1, TRUE ~ count),length = as.numeric(length)*-1)
wt_red_vs_dcl.toPlot.specifics$class <- factor(wt_red_vs_dcl.toPlot.specifics$class, levels = c("dcl-specific", "common", "wt-specific"))

toLabel <- wt_red_vs_dcl.toPlot.specifics %>% mutate(count = case_when(strand == "-" ~ count *-1, TRUE ~ count)) %>%
  group_by(class) %>% dplyr::summarize(total=sum(count))

wt_red_vs_dcl.plot.specific<- ggplot(wt_red_vs_dcl.toPlot.specifics,aes(x=class,y=count,fill=as.factor(length))) +
  geom_col() + 
  theme_bw() + themes + scale_fill_brewer(palette="Spectral",direction=-1)+
  theme(legend.position = "right",legend.text = element_text(size=5),legend.key.size = unit(0.1,'in'),legend.title = element_blank() )+
  geom_text(inherit.aes = F, data=toLabel,aes(x=class,y=-1000,label = total),size=3,color="black")+
  scale_y_continuous(labels=scales::label_number(scale_cut = cut_short_scale()))
ggsave("wt_red_and_dcl_specific_sRNAs_size.filter2.pdf",plot=wt_red_vs_dcl.plot.specific,height=4,width=5.5,units="in")

wt_red_vs_dcl.venn_plus_bar <- ggarrange(wt_red_vs_dcl.plot,wt_red_vs_dcl.plot.specific,ncol=1,heights=c(1,1))
ggsave("wt_red_and_dcl_specific_sRNAs_size_plus_venn.filter2.pdf",plot=wt_red_vs_dcl.venn_plus_bar,height=3,width=3.5,units="in")



#----------------------------------
## Compare dcl and Full Red
# List containing the vectors
wt_fullRed_vs_dcl <- list(
  dcl = in.dcl.red$seq,
  wt.fullRed = in.wt.fullRed$seq 
)

# Apply the function to the list x
filtered_wt_fullRed_vs_dcl <- filter_sequences(wt_fullRed_vs_dcl)

wt_fullRed_vs_dcl.plot <- ggvenn(filtered_wt_fullRed_vs_dcl, fill_color=c("#853061","#C97795"),
                        fill_alpha = 0.8,auto_scale = T,text_size = 3,set_name_size = 3,
                        stroke_size = 0.5) 
ggsave("wt_fullRed_vs_dcl_sRNA.venn.filter2.pdf",plot=wt_fullRed_vs_dcl.plot)

# Get unique sequences in each filtered list
unique_dcl <- unique(filtered_wt_fullRed_vs_dcl$dcl)
unique_wt_fullRed <- unique(filtered_wt_fullRed_vs_dcl$wt.fullRed)

# Get unique sequences in dcl but not in wt.fullRed
dcl.specific.list <- setdiff(unique_dcl, unique_wt_fullRed)
dcl.specific <- as.data.frame(dcl.specific.list) %>% separate(dcl.specific.list,into = c("start","stop","seq","strand","length"),sep="_")%>% mutate(gene = "35S_RUBY_transgene",class = "dcl-specific")%>% relocate(gene,start,stop,seq,class,strand,length)
write.table(dcl.specific,file = "wt_fullRed_vs_dcl.dcl_specific_sRNAs.txt",col.names = T, row.names = F,quote = F,sep="\t")

# Get unique sequences in wt but not in dcl
wt.specific.list <- setdiff(unique_wt_fullRed,unique_dcl)
wt.specific <- as.data.frame(wt.specific.list) %>% separate(wt.specific.list,into = c("start","stop","seq","strand","length"),sep="_") %>% mutate(gene = "35S_RUBY_transgene",class = "wt-specific")%>% relocate(gene,start,stop,seq,class,strand,length)
write.table(wt.specific,file = "wt_fullRed_vs_dcl.wt_specific_sRNAs.txt",col.names = T, row.names = F,quote = F,sep="\t")

# Get unique sequences present in both dcl and wt.red
com.specific.list <- intersect(unique_wt_fullRed,unique_dcl)
common.srna <- as.data.frame(com.specific.list) %>% separate(com.specific.list,into = c("start","stop","seq","strand","length"),sep="_") %>% mutate(gene = "35S_RUBY_transgene",class = "common")%>% relocate(gene,start,stop,seq,class,strand,length)
write.table(common.srna,file = "wt_fullRed_vs_dcl.shared_sRNAs.txt",col.names = T, row.names = F,quote = F,sep="\t")

wt_fullRed_vs_dcl.toPlot.specifics <- rbind(wt.specific,dcl.specific,common.srna) %>% select(!seq) %>% group_by(strand,length,class) %>%
  dplyr::summarize(count=n()) %>% mutate(count = case_when(strand == "-" ~ count *-1, TRUE ~ count),length = as.numeric(length)*-1)
wt_fullRed_vs_dcl.toPlot.specifics$class <- factor(wt_fullRed_vs_dcl.toPlot.specifics$class, levels = c("dcl-specific", "common", "wt-specific"))

toLabel <- wt_fullRed_vs_dcl.toPlot.specifics %>% mutate(count = case_when(strand == "-" ~ count *-1, TRUE ~ count)) %>%
  group_by(class) %>% dplyr::summarize(total=sum(count))
wt_fullRed_vs_dcl.toPlot.specifics <- na.omit(wt_fullRed_vs_dcl.toPlot.specifics)

wt_fullRed_vs_dcl.plot.specific<- ggplot(wt_fullRed_vs_dcl.toPlot.specifics,aes(x=class,y=count,fill=as.factor(length))) +
  geom_col() + 
  theme_bw() + themes + scale_fill_brewer(palette="Spectral",direction=-1)+
  theme(legend.position = "right",legend.text = element_text(size=5),legend.key.size = unit(0.1,'in'),legend.title = element_blank() )+
  geom_text(inherit.aes = F, data=toLabel,aes(x=class,y=-1000,label = total),size=3,color="black")
ggsave("wt_fullRed_and_dcl_specific_sRNAs_size.filter2.pdf",plot=wt_fullRed_vs_dcl.plot.specific,height=4,width=5.5,units="in")

wt_fullRed_vs_dcl.venn_plus_bar <- ggarrange(wt_fullRed_vs_dcl.plot,wt_fullRed_vs_dcl.plot.specific,ncol=1,heights=c(1,1))
ggsave("wt_fullRed_and_dcl_specific_sRNAs_size_plus_venn.filter2.pdf",plot=wt_fullRed_vs_dcl.venn_plus_bar,height=3,width=3.5,units="in")

merge_dcl_vs_wt_red_fullRed <- ggarrange(wt_red_vs_dcl.plot,wt_fullRed_vs_dcl.plot,
          wt_red_vs_dcl.plot.specific+ rremove("legend"),wt_fullRed_vs_dcl.plot.specific+ rremove("legend"),
          ncol=2,nrow=2)

ggsave("wt_red_or_fullRed_and_dcl_specific_sRNAs_size_plus_venn.filter.pdf",plot=merge_dcl_vs_wt_red_fullRed,height=3,width=4,units="in")



