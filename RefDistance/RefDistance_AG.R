
#######################################################################
# Yunhu Wan 
# Aug 20, 2019                          
# The purpose of this script is to generate BC and Unifrac distance to 
# HMP reference stool and HMP reference nasal samples.
#######################################################################

setwd("~/Desktop/Mitch_Power/Dropbox/Code/standard_aug2019/Deposiit/") # Working folder.Change to directory it needs.
rm(list=ls())

source('microbiome-fixed-reference-functions.r')

load("taxonomy-GG13-8-and-HMP-reference-sets-L2-L6.rdata")  #Reference Stool and Nasal samples

tax.level.bc = "l6.genus"
tax.level.unifrac = "Rank6"
BC_file_name = "AG_1879_L6.txt"
treefile = "97_otus.tree"
resultfile = "AG_1879_Level6.txt"
sample_data = "AG_1879_final.biom"
full_data = "merge_HMP_AG_even1000.biom"

#For Bray Curtis, it needs Relative abundance table for each taxonomy level

#L2
#tax.level.bc = "l2.phylum" #"l2.phylum", "l3.class", "l4.order", "l5.family","l6.genus"
#BC_file_name = "AG_1879_L2.txt" #phylum: HMP_L2.txt, class: HMP_L3.txt, order: HMP_L4.txt, family: HMP_L5.txt, genus: HMP_L6.txt

#L3
#tax.level.bc = "l3.class" #"l2.phylum", "l3.class", "l4.order", "l5.family","l6.genus"
#BC_file_name = "AG_1879_L3.txt" #phylum: HMP_L2.txt, class: HMP_L3.txt, order: HMP_L4.txt, family: HMP_L5.txt, genus: HMP_L6.txt


rel_file <- read.table(BC_file_name,check.names=F,header = T,sep="\t")
hmp.dist.bc <- HMPdistance(tax.level.bc, d.new.filename = BC_file_name, d.new.ix.col.not.rel.abu = 1:14, measure = 'bc')
rownames(hmp.dist.bc) <- rel_file[,1]


#Unweighted Unifrac 

#tax.level.unifrac = "Rank3" #phylum: "Rank2", class: "Rank3", order: "Rank4", family: "Rank5", genus: "Rank6"
#full_data = "merge_HMP_AG_even1000.biom"
#treefile = "97_otus.tree"
unifrac_uw = Calculate_Unifrac(full_data,tax.level.unifrac,treefile,"uw")
#sample_list <- rownames(unifrac_uw)
#sample_data = "AG_1879_final.biom"
biomfile <- import_biom(sample_data)
otu.table <- otu_table(biomfile,taxa_are_rows = FALSE)
#unifrac_uw_test <- Calculate_Unifrac(data,tax.level.unifrac,treefile,"uw")
sample_list <- rownames(t(otu.table@.Data))

dist.to.stool <- hmp.distance.helper.new.AG(d.ref = as.character(d.ref.stool.info$sample.id), # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "uw")

d.ref = as.character(d.ref.nasal.info$sample.id)
d.ref = d.ref [-c(11,16,39,53,54,63,64)] # remove several missing for HMP GG even1000 sample because of rarefaction
dist.to.nasal <- hmp.distance.helper.new.AG(d.ref, # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "uw")
hmp.dist.uw <- cbind(dist.to.stool, dist.to.nasal)
colnames(hmp.dist.uw) <- c('UW dist.to.ref.stool', 'UW dist.to.ref.nasal')
rownames(hmp.dist.uw) <- as.character(sample_list)


#Weighted Unifrac
#tax.level.unifrac = "Rank2" #phylum: "Rank2", class: "Rank3", order: "Rank4", family: "Rank5", genus: "Rank6"
#full_data = "merge_HMP_AG_even1000.biom"
#treefile = "97_otus.tree"
unifrac_wu = Calculate_Unifrac(full_data,tax.level.unifrac,treefile,"wu")
#sample_list <- rownames(unifrac_wu)
dist.to.stool <- hmp.distance.helper.new.AG(d.ref = as.character(d.ref.stool.info$sample.id), # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "wu")
d.ref = as.character(d.ref.nasal.info$sample.id)
d.ref = d.ref [-c(11,16,39,53,54,63,64)] # remove several missing for HMP GG even1000 samples.

dist.to.nasal <- hmp.distance.helper.new.AG(d.ref, # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "wu")
hmp.dist.wu <- cbind(dist.to.stool, dist.to.nasal)
colnames(hmp.dist.wu) <- c('WU dist.to.ref.stool', 'WU dist.to.ref.nasal')
rownames(hmp.dist.uw) <- as.character(sample_list)

All_matrix <- cbind(hmp.dist.bc,hmp.dist.uw,hmp.dist.wu)

write.table(All_matrix,resultfile,quote=FALSE,sep="\t")



