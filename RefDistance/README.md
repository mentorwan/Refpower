# RefDistance

The purpose of this script is to compute mean Bray_Curtis, Unweighted Unifrac and Weighted Unifrac distances from test samples to 92 HMP (Human Microbiome Project) reference stool samples and 74 reference nasal samples. If the test sample is an HMP sample, the 2011 Greengene reference database is needed. If the test sample is not from HMP, the Greengene 13.8 reference database is needed.

## Examples
* Usage 1 (HMP samples): Rscript RefDistance.R HMP_L2.txt HMP_YW_even1000.biom L2 Y output.txt
* Usage 2 (non HMP samples): Rscript RefDistance.R AG_1879_L2.txt AG_1879_final.biom L2 N output.txt

## Input: 

* Relative abundance file (Examples: HMP_L2.txt or AG_1879_L2.txt)
* sample data biom file (Examples:HMP_YW_even1000.biom, AG_1879_final.biom)
* tax.level (Examples: L2,L3,L4,L5,L6)
* Y/N  : (Y for HMP test samples with HMP using 2011 Greengene database. N for test samples using Greengene 13.8)

## Output: 
(Example: output.txt)

A table with matrix n * 6 which n is number of samples. For each distance, it includes 2 distances to HMP stool samples and HMP nasal samples. The first two columns are for the Bray_Curtis distance, the next two for the Unweighted Unifrac distance and the last two columns for weighted unifrac distance.

