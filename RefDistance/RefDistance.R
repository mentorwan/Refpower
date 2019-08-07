
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

tax.level.bc = "l2.phylum" #"l2.phylum", "l3.class", "l4.order", "l5.family","l6.genus"
BC_file_name = "HMP_L2.txt" #phylum: HMP_L2.txt, class: HMP_L3.txt, order: HMP_L4.txt, family: HMP_L5.txt, genus: HMP_L6.txt

rel_file <- read.table(BC_file_name,check.names=F,header = T)
hmp.dist.bc <- HMPdistance(tax.level.bc, d.new.filename = BC_file_name, d.new.ix.col.not.rel.abu = 1:7, measure = 'bc')
rownames(hmp.dist.bc) <- rel_file[,1]


#Unweighted Unifrac 

tax.level.unifrac = "Rank2" #phylum: "Rank2", class: "Rank3", order: "Rank4", family: "Rank5", genus: "Rank6"
data = "HMP_YW_even1000.biom"
treefile = "rep_set_v35.tre"
unifrac_uw = Calculate_Unifrac(data,tax.level.unifrac,treefile,"uw")
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
tax.level.unifrac = "Rank2" #phylum: "Rank2", class: "Rank3", order: "Rank4", family: "Rank5", genus: "Rank6"
data = "HMP_YW_even1000.biom"
treefile = "rep_set_v35.tre"
unifrac_wu = Calculate_Unifrac(data,tax.level.unifrac,treefile,"wu")
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



