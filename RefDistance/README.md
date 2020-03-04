# RefDistance
#
# #######################################################################
# Dr. Yunhu Wan 
# Dr. Mitch Gail
# Dr. Jianxin Shi
# Mar 4, 2020                          
# The purpose of this script is to generate BC and Unifrac distances to 
# HMP 92 reference stool samples and HMP 74 reference nasal samples.
# ######################################################################

Usage 1 (HMP samples): Rscript RefDistance.R HMP_L2.txt HMP_YW_even1000.biom L2 Y output.txt
Usage 2 (non HMP samples): Rscript RefDistance.R AG_1879_L2.txt AG_1879_final.biom L2 N output.txt

Input: 

Relative abundance file (HMP_L2.txt)
sample data biom file (HMP_YW_even1000.biom)
tax.level (L2)
Y/N  (test samples are HMP sample or non-HMP samples)

Output: (output.txt)

A table with matrix n*6 which n is number of samples. For each distance, it includes 2 distances to HMP stool and HMP nasal.