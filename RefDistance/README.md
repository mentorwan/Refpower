# RefDistance

The purpose of this script is to generate BC and Unifrac distances to 
HMP 92 reference stool samples and HMP 74 reference nasal samples.

## Examples
Usage 1 (HMP samples): Rscript RefDistance.R HMP_L2.txt HMP_YW_even1000.biom L2 Y output.txt
Usage 2 (non HMP samples): Rscript RefDistance.R AG_1879_L2.txt AG_1879_final.biom L2 N output.txt

## Input: 

Relative abundance file (Examples: HMP_L2.txt or AG_1879_L2.txt)
sample data biom file (Examples:HMP_YW_even1000.biom, AG_1879_final.biom)
tax.level (Examples: L2,L3,L4,L5,L6)
Y/N  : (Y for test samples with HMP using 2011 greengene database. N for test samples using Greengene 13.8)

## Output: (Example: output.txt)

A table with matrix n*6 which n is number of samples. For each distance, it includes 2 distances to HMP stool and HMP nasal.
