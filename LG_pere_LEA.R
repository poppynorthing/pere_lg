#LG_pere LEA
#setting up pipeline for when I sort out issues with reads

#enter LEA_analyses directory within the pere_lg R project directory

#load libraries
library(LEA)

#convert my VCF with population SNPS to the lfmm and geno formats that LEA needs
vcf2lfmm("populations.snps.vcf")
vcf2geno("populations.snps.vcf")

pc = pca("populations.snps.vcf", scale = TRUE)
tw = tracy.widom(pc)

?vcf2lfmm()