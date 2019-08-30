library(Biostrings) 
library(ggplot2)
library(plyr)
library(stringr)
library(gridExtra)
library(reshape2)
library(knitr)
library(vegan)
library(DESeq2)
library(phyloseq)
library(scales)
library(ggthemes)
library(xlsx)
library(rJava)

setwd("M:/2016/Workspaces/M_final_files_for_r_HG0115")

load("HG0115_subset_by_gen")

RSV_a # confirms that this is a phyloseq object
head(sample_data(RSV_a)$DatePostInfection)
sample_sums(RSV_a) # don't need to prune out any samples as all have over 10000 reads

deseq_RSV_a = phyloseq_to_deseq2(RSV_a, ~ DatePostInfection) # Converts phyloseq object to a DESeqDataSet
deseq_RSV_a = DESeq(deseq_RSV_a, test="Wald", fitType="parametric") # Does the actual testing (conbines it with above data set and includes a default multiple-infreence correction)
res_RSV_a = results(deseq_RSV_a, cooksCutoff = FALSE, contrast = c("DatePostInfection", "D7", "D0")) # Creates a table of the results of the DEseq test
alpha = 0.01 # what we want the significant p value to be
sigtab_RSV_a = res_RSV_a[which(res_RSV_a$padj < alpha), ] # orders the results by the adjusted p value
head(sigtab_RSV_a)
dim(sigtab_RSV_a)
sigtab_RSV_a = cbind(as(sigtab_RSV_a, "data.frame"), as(tax_table(RSV_a)[rownames(sigtab_RSV_a), ], "matrix")) # creates an actual table of the results


##### PLotting  D0 vs D7


theme_set(theme_bw())

x = tapply(sigtab_RSV_a$log2FoldChange, sigtab_RSV_a$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab_RSV_a$Family = factor(as.character(sigtab_RSV_a$Family), levels=names(x))

x = tapply(sigtab_RSV_a$log2FoldChange, sigtab_RSV_a$OTUID, function(x) max(x))
x = sort(x, TRUE)
sigtab_RSV_a$OTUID = factor(as.character(sigtab_RSV_a$OTUID), levels = names(x))

ggplot(sigtab_RSV_a, aes(x=OTUID, y = log2FoldChange, color=Family)) +geom_point(size=4) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  ggtitle(expression(atop(bold("RSV"), atop(italic("Changes in individual OTU abundance from day 0 to day 7"),""))))+
  geom_hline(yintercept = 0) +
  ylim(-10, 10) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(plot.title = element_text(face = "bold", size = 18)) + 
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(plot.title = element_text(margin=margin(0,0,0.000000000001,0))) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_colour_manual(values = c("aquamarine", "#F0E442", "cadetblue3", "#E68F00", "#D55E00", "purple", "#009E73", "#999999", "#CC79A7"))


