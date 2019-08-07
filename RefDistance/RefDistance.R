
#######################################################################
# Yunhu Wan 
# Aug 6, 2019                          
# The purpose of this script is to generate BC and Unifrac distance to 
# HMP reference stool and HMP reference nasal samples.
#######################################################################

setwd("~/Desktop/Mitch_Power/Dropbox/Code/standard_aug2019/Deposiit/") # Working folder.Change to directory it needs.
rm(list=ls())

source('microbiome-fixed-reference-functions.r')

load("taxonomy-GG13-8-and-HMP-reference-sets-L2-L6.rdata")  #Reference Stool and Nasal samples


#For Bray Curtis, it needs Relative abundance table for each taxonomy level

tax.level = "l2.phylum"
BC_file_name = "HMP_L2.txt"

rel_file <- read.table(BC_file_name,check.names=F,header = T)
hmp.dist.bc <- HMPdistance(tax.level, d.new.filename = "HMP_L2.txt", d.new.ix.col.not.rel.abu = 1:7, measure = 'bc')
rownames(hmp.dist.bc) <- rel_file[,1]


#Unweighted Unifrac 

tax.level = "Rank2"
data = "HMP_YW_even1000.biom"
treefile = "rep_set_v35.tre"
unifrac_uw = Calculate_Unifrac("HMP_YW_even1000.biom","Rank2","rep_set_v35.tre","uw")
sample_list <- rownames(unifrac_uw)
dist.to.stool <- hmp.distance.helper.new(d.ref = as.character(d.ref.stool.info$sample.id), # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "uw")

dist.to.nasal <- hmp.distance.helper.new(d.ref = as.character(d.ref.nasal.info$sample.id), # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "uw")
hmp.dist.uw <- cbind(dist.to.stool, dist.to.nasal)
colnames(hmp.dist.uw) <- c('UW dist.to.ref.stool', 'UW dist.to.ref.nasal')
rownames(hmp.dist.uw) <- as.integer(sample_list)


#Weighted Unifrac
tax.level = "Rank2"
data = "HMP_YW_even1000.biom"
treefile = "rep_set_v35.tre"
unifrac_wu = Calculate_Unifrac("HMP_YW_even1000.biom","Rank2","rep_set_v35.tre","wu")
sample_list <- rownames(unifrac_wu)
dist.to.stool <- hmp.distance.helper.new(d.ref = as.character(d.ref.stool.info$sample.id), # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "wu")

dist.to.nasal <- hmp.distance.helper.new(d.ref = as.character(d.ref.nasal.info$sample.id), # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "uw")
hmp.dist.wu <- cbind(dist.to.stool, dist.to.nasal)
colnames(hmp.dist.wu) <- c('WU dist.to.ref.stool', 'WU dist.to.ref.nasal')
rownames(hmp.dist.uw) <- as.integer(sample_list)

All_matrix <- cbind(hmp.dist.bc,hmp.dist.uw,hmp.dist.wu)

write.table(All_matrix,"output.txt",quote=FALSE,sep="\t")



