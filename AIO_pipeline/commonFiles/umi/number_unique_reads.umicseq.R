library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(janitor)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(plyr)
setwd("/Users/mariannekramer/Desktop/mkramer_server/projects/target_capture/at_round2/data/6_minimap/2_map_to_capture/samFiles/mapped_to_targets/bedFiles/umi/extract/filtered/clustered_umiseq//")
themes <- theme(legend.position = "top", 
                axis.ticks.length.x=unit(0.25,"cm"),
                axis.ticks.length.y=unit(0.25,"cm"), 
                plot.title=element_text(hjust=0.5, color = 'black', size = 12, face = 'bold'), 
                axis.title.x= element_blank(), 
                axis.text.x=element_text(color = 'black',size = 8, face = 'bold'), 
                axis.line.x=element_line(color = 'black'), 
                axis.text.y= element_text(color = 'black',size = 8, face = 'bold'),
                axis.line.y=element_line(color = 'black'), 
                axis.title.y= element_text(color = 'black',size = 8, face = 'bold'),
                panel.background = element_blank(),
                legend.text=element_text(color = 'black',size = 8, face = 'bold'),
                legend.key.size= unit(0.5,"cm"),
                legend.title = element_text(color = 'black',size = 8, face = 'bold'))
colName <- "1_col0.mapped_to_TEs.umiCount.umic.txt"
ddmName <- "9_ddm1.mapped_to_TEs.umiCount.umic.txt"
classes <- "/Users/mariannekramer/Desktop/mkramer_server/annotations/arabidopsis/at_round2/At_array.v1.targets_class.txt"
inCol <- read.table(colName, header = T, as.is = T)
inDDM <- read.table(ddmName, header = T, as.is = T)
inClass <- read.table(classes, col.names = c("FaName","Class"), as.is = T)

inClass <- inClass %>% separate(FaName, into = c("Name"),sep="_")

a<-inCol %>% right_join(inDDM,by=c("Gene"),suffix=c("col0","ddm")) %>% transmute(Gene,col0=Countcol0,ddm1=Countddm) %>%
  pivot_longer(cols=!Gene, names_to = "Sample") %>% separate(Gene,into=c("Name","type"), sep = "[.]") %>% dplyr::group_by(Name,Sample) %>%
  dplyr::summarize(Total = sum(value,na.rm = TRUE))

ratios <- a %>% pivot_wider(names_from=Sample,values_from =Total) %>% mutate(ratio = ddm1/col0) 

b<-inClass %>% filter(grepl("TE",Class)) %>% rowwise() %>%right_join(a,by="Name")

methPlot <- b %>% filter(Class=="Methylated_TE") %>% ggplot(aes(x=Name,y=Total,fill=Sample,label=Total)) + geom_bar(stat="identity",position=position_dodge(0.9)) + 
  scale_y_continuous(labels=scales::comma) +
  geom_text(aes(label = prettyNum(round(Total,1),big.mark=",",scientific=F)), position = position_dodge(0.9),vjust = -1,size=3)+
  themes+ ylab("Number of unique reads") + facet_wrap(~Class)
unmethPlot <- b %>% filter(Class=="Unmethylated_TE") %>% ggplot(aes(x=Name,y=Total,fill=Sample,label=Total)) + geom_bar(stat="identity",position=position_dodge(0.9)) + 
  scale_y_continuous(labels=scales::comma) +
  geom_text(aes(label = prettyNum(round(Total,1),big.mark=",",scientific=F)), position = position_dodge(0.9),vjust = -1,size=3)+
  themes+ ylab("Number of unique reads") + facet_wrap(~Class)

outPlot <- ggarrange(methPlot , unmethPlot , ncol=2,nrow=1,widths=c(1,0.5),align="v",common.legend = T)
ggsave("col_ddm1_tes.umi_collapsed.umicseq.pdf",plot=outPlot,width=12,height=5.5, units="in")

