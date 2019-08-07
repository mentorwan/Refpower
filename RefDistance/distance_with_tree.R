#setwd("/data/wany/Mitch/Power/Power_Unifrac")

library(phyloseq)
library(GUniFrac)
#library(data.table)

#source("GUniFrac_single.R")
source("GUniFrac.R")

#args = commandArgs(trailingOnly=TRUE)

#row <- args[1]
#row = as.numeric(row)

tax.level <- "Rank2"

tree <- read_tree_greengenes("rep_set_v35.tre")
#tree$edge.length[which(is.na(tree$edge.length))] <- 0  #Modify branch distance in tree NA value to 0

#biomfile <- read.csv("HMP_YW_even1000_L2.txt",check.names = FALSE, sep = "\t")

biomfile <- import_biom("HMP_YW_even1000.biom")   #Rarefaction happens here.

biomfile <- tax_glom(biomfile,tax.level)  #Collapse into different taxonomy level if level2 using "Rank2"
otu.table <- otu_table(biomfile,taxa_are_rows = TRUE)

otu.table <- otu.table@.Data
otu.table <- t(otu.table)
#GUniFrac_single(otu.table,tree,alpha=c(0,0.5,1),row)


unifracs <- GUniFrac(otu.table,tree,alpha=c(0,0.05,1))

unifrac_uw <- unifracs$unifracs[,,4]  #d_UW: unweighted unifrac
unifrac_wu <- unifracs$unifracs[,,3]  #d_1: weighted unifrac

#This is based on Level 2 if different level, needs to update taxa_level

save(unifrac_uw, file="Unweighted_Unifrac_l2.phylum.RData")
save(unifrac_wu, file="Weighted_Unifrac_l2.phylum.RData")

