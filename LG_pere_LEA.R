#LG_pere LEA
#setting up pipeline for when I sort out issues with reads

#enter LEA_analyses directory within the pere_lg R project directory

#load libraries
library(LEA)


#convert my VCF with population SNPS to the lfmm and geno formats that LEA needs
vcf2lfmm("LEA/populations.snps.vcf")
vcf2geno("LEA/populations.snps.vcf")

project.missing = snmf("LEA/populations.snps.lfmm", K = 4,
                       entropy = TRUE, repetitions = 10,
                       project = "new")

best = which.min(cross.entropy(project.missing, K = 4))
impute(project.missing, "LEA/populations.snps.lfmm",
       method = 'mode', K = 4, run = best)


dat.imp = read.lfmm("LEA/populations.snps.lfmm_imputed.lfmm")
mean(lfmm[dat == 9] == dat.imp[dat == 9] )


pc <- pca("LEA/populations.snps.lfmm_imputed.lfmm", scale = TRUE)
tw = tracy.widom(pc)

# plot the percentage of variance explained by each component
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)

# main options
# K = number of ancestral populations
# entropy = TRUE computes the cross-entropy criterion,
# CPU = 4 is the number of CPU used (hidden input)
project = NULL
project = snmf("LEA/populations.snps.geno",
               K = 1:9,
               entropy = TRUE,
               repetitions = 10,
               project = "new")
plot(project, col = "blue", pch = 19, cex = 1.2)


best = which.min(cross.entropy(project, K = 9))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold", "lavender", "grey", "pink", "seagreen", "black")
barchart(project, K = 9, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)


# Genome scan for selection: opulation differentiation tests
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 9)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)
abline(a=5, b=0)
