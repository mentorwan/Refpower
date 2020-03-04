
#######################################################################
# Dr. Yunhu Wan 
# Dr. Mitch Gail
# Mar 4, 2020                          
# The purpose of this script is to generate BC and Unifrac distances to 
# HMP 92 reference stool samples and HMP 74 reference nasal samples.
#######################################################################

setwd(".") # Working folder. Change to directory it needs.
rm(list=ls())

# Notes: The phyloseq package needs to be installed

source('microbiome-fixed-reference-functions.r')

load("taxonomy-GG13-8-and-HMP-reference-sets-L2-L6.rdata")  #Reference Stool and Nasal samples


args = commandArgs(trailingOnly = TRUE)


BC_input = args[1]
Unifrac_input = args[2]
level = args[3]
HMP_tag = args[4]
Output_file = args[5]

#For Bray Curtis, it needs Relative abundance table for each taxonomy level

#level is one of "l2.phylum", "l3.class", "l4.order", "l5.family","l6.genus"
if(level == 'L2' || level == "l2")
{
  tax.level.bc = "l2.phylum" 
  tax.level.unifrac="Rank2"
} else if (level == 'L3' || level == "l3")
{
  tax.level.bc = "l3.class"
  tax.level.unifrac="Rank3"
}else if (level == 'L4' || level == "l4")
{
  tax.level.bc = "l4.order"
  tax.level.unifrac="Rank4"
}else if (level == 'L5' || level == "l5")
{
  tax.level.bc = "l5.family"
  tax.level.unifrac="Rank5"
}else if (level == 'L6' || level == "l6")
{
  tax.level.bc = "l6.genus"
  tax.level.unifrac = "Rank6"
}

#BC_file_name = "HMP_L2.txt" #phylum: HMP_L2.txt, class: HMP_L3.txt, order: HMP_L4.txt, family: HMP_L5.txt, genus: HMP_L6.txt
BC_file_name = BC_input

rel_file <- read.table(BC_file_name,check.names=F,header = T)

if(HMP_tag == "Y")   #GG 2011 database
{
  hmp.dist.bc <- HMPdistance(tax.level.bc, d.new.filename = BC_file_name, d.new.ix.col.not.rel.abu = 1:7, measure = 'bc')
}else if (HMP_tag == "N") ## GG 13_8 database
{
  hmp.dist.bc <- HMPdistance(tax.level.bc, d.new.filename = BC_file_name, d.new.ix.col.not.rel.abu = 1:14, measure = 'bc')
}

rownames(hmp.dist.bc) <- rel_file[,1]


#Unweighted Unifrac 

#tax.level.unifrac = "Rank2" #phylum: "Rank2", class: "Rank3", order: "Rank4", family: "Rank5", genus: "Rank6"
#data = "HMP_YW_even1000.biom"
data = args[2]

if(HMP_tag == "Y")   #GG 2011 database
{
  treefile = "rep_set_v35.tre"
}else if (HMP_tag == "N") ## GG 13_8 database
{
  treefile = "97_otus.tree"
}

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
#tax.level.unifrac = "Rank2" #phylum: "Rank2", class: "Rank3", order: "Rank4", family: "Rank5", genus: "Rank6"
#data = "HMP_YW_even1000.biom"
#treefile = "rep_set_v35.tre"
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

write.table(All_matrix,Output_file,quote=FALSE,sep="\t")



