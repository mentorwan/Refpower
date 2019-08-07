# RefDistance

This function will generate Bray Curtis and Unweighted/Weighted Unifrac distance to HMP reference stool sample and HMP reference nasal samples

Example code is in RefDistance.R file.

## Bray_Curtis Distance

Input: 
tax.level.bc #  taxonomy level, one of the following: "l2.phylum", "l3.class", "l4.order", "l5.family", "l6.genus"
BC_file_name # file name for relative abundance file in given taxnomy level

Output:
hmp.dist.bc  # n*2 matrix
```
tax.level.bc = "l2.phylum"
BC_file_name = "HMP_L2.txt"

rel_file <- read.table(BC_file_name,check.names=F,header = T)
hmp.dist.bc <- HMPdistance(tax.level.bc, d.new.filename = BC_file_name, d.new.ix.col.not.rel.abu = 1:7, measure = 'bc')
rownames(hmp.dist.bc) <- rel_file[,1]

```

## Unweighted Unifrac Distance

Input: 
tax.level.unifrac #  taxonomy level depending on biom file. In this case. "Rank2","Rank3","Rank4","Rank5","Rank6"
data  #  original biom file 
tree # phylogentic tree.

Output:
hmp.dist.uw  # n*2 matrix


```
tax.level.unifrac = "Rank2"
data = "HMP_YW_even1000.biom"
treefile = "rep_set_v35.tre"
unifrac_uw = Calculate_Unifrac(data,tax.level.unifrac,treefile,"uw")
sample_list <- rownames(unifrac_uw)
dist.to.stool <- hmp.distance.helper.new(d.ref = as.character(d.ref.stool.info$sample.id), 
                                         d.new = sample_list,
                                         measure = "uw")

dist.to.nasal <- hmp.distance.helper.new(d.ref = as.character(d.ref.nasal.info$sample.id), # matrix of rel abu
                                         d.new = sample_list,
                                         measure = "uw")
hmp.dist.uw <- cbind(dist.to.stool, dist.to.nasal)

```

## Weighted Unifrac Distance

Input: 
tax.level.unifrac #  taxonomy level depending on biom file. In this case. "Rank2","Rank3","Rank4","Rank5","Rank6"
data  #  original biom file 
tree # phylogentic tree.

Output:
hmp.dist.wu  # n*2 matrix


```
tax.level = "Rank2"
data = "HMP_YW_even1000.biom"
treefile = "rep_set_v35.tre"
unifrac_uw = Calculate_Unifrac(data,tax.level.unifrac,treefile,"uw")
sample_list <- rownames(unifrac_uw)
dist.to.stool <- hmp.distance.helper.new(d.ref = as.character(d.ref.stool.info$sample.id), 
                                         d.new = sample_list,
                                         measure = "wu")

dist.to.nasal <- hmp.distance.helper.new(d.ref = as.character(d.ref.nasal.info$sample.id), 
                                         d.new = sample_list,
                                         measure = "wu")
hmp.dist.wu <- cbind(dist.to.stool, dist.to.nasal)

```