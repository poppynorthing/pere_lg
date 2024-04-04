#LG_pere LEA
#setting up pipeline for when I sort out issues with reads

#enter LEA_analyses directory within the pere_lg R project directory
setwd("LEA_analyses")

#load libraries
library(LEA)

#convert my VCF with population SNPS to the lfmm and geno formats that LEA needs
output = vcf2lfmm("populations.snps.vcf", "pere_genotypes.lfmm")
vcf2geno("populations.snps.vcf","pere_genotypes.geno")

pc = pca("populations.snps.lfmm", scale = TRUE)
tw = tracy.widom(pc)

?vcf2lfmm()
