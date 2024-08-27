# List of package names to check, install, and load
package_names <- c("tidyverse", "ggplot2", "ggrepel", "ggpubr")

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
if (length(args)!=1) {
  cat("please put input file")
  q()
}

#setwd("/Users/mariannekramer/Desktop/mkramer_server/projects/target_capture/ruby_round3/data/6_minimap/mapped_pct/")
fileName<- args[1] #"/Users/mariannekramer/Desktop/mkramer_server/projects/target_capture/ruby_round1/data/6_minimap/mapped_pct/number_reads_mapped.txt" #args[1]
inFile <- read.table(fileName,header=T, as.is=T)
themes <- theme(legend.position = "right", 
                axis.ticks.length.x=unit(0.25,"cm"),
                axis.ticks.length.y=unit(0.25,"cm"), 
                plot.title=element_text(hjust=0.5, color = 'black', size = 12, face = 'bold'), 
                axis.title.x= element_blank(), 
                axis.text.x=element_text(color = 'black',size = 10, face = 'bold',angle=90,hjust=1,vjust=0.5), 
                axis.line.x=element_line(color = 'black'), 
                axis.text.y= element_text(color = 'black',size = 10, face = 'bold'),
                axis.line.y=element_line(color = 'black'), 
                axis.title.y= element_text(color = 'black',size = 10, face = 'bold'),
                panel.background = element_blank(),
                legend.text=element_text(color = 'black',size = 8, face = 'bold'),
                legend.key.size= unit(0.5,"cm"),
                legend.title = element_text(color = 'black',size = 8, face = 'bold'))

toPlot <- inFile %>% 
  rowwise() %>%
  mutate(sum=sum(c(rRNA,Targets,Non.targets,unmapped))) %>%
  separate(sample, c("name","promoter","pheno"),"[.]") %>%
  mutate(name, 
            rRNApct =rRNA/sum*100, 
            targetpct =Targets/sum*100,
            genomepct =Non.targets/sum*100,
            unmappedpct =unmapped/sum*100) %>% 
  pivot_longer(c(rRNApct,targetpct,genomepct,unmappedpct),names_to = "Class") %>% unite("name",name:pheno,sep=".") %>%
  transmute(name, Class, value, sum,reads =value * sum/100) %>%
  mutate(across(value,round,2)) 
  
mapped_plot <- toPlot %>% 
  ggplot(aes(x=name,y=value,fill=Class)) + 
  geom_col()+
  geom_text(label = paste(toPlot$value,"%","\n",prettyNum(toPlot$reads,big.mark=",",scientific=F),sep=""), position = position_stack(vjust = 0.5),size=3)+ 
  themes + ggtitle("Percent Mapped for all reads") + ylab("Percentage")

outFile1 <- paste(strsplit(fileName,".txt")[1][[1]],"mapped.barplots.pdf",sep=".")
ggsave(outFile1,mapped_plot,height=5,width=7,unit="in")
