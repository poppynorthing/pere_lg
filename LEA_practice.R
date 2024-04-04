#Practice with LEA vignette

library(LEA)

# Creation of a genotype matrix data file: "genotypes.lfmm"
# The data include 400 SNPs for 50 individuals.
data("tutorial")

# Write genotypes in the lfmm format
write.lfmm(tutorial.R, "genotypes.lfmm")

# Write genotypes in the geno format
write.geno(tutorial.R, "genotypes.geno")

# creation of an environment gradient file: gradient.env.
# The .env file contains a single ecological variable
write.env(tutorial.C, "gradients.env")


# run of pca
# Available options, K (the number of PCs),
# center and scale.
# Create files: genotypes.eigenvalues - eigenvalues,
# genotypes.eigenvectors - eigenvectors,
# genotypes.sdev - standard deviations,
# genotypes.projections - projections,
# Create a pcaProject object: pc.
pc = pca("genotypes.lfmm", scale = TRUE)
tw = tracy.widom(pc)

# plot the percentage of variance explained by each component
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)

# main options
# K = number of ancestral populations
# entropy = TRUE computes the cross-entropy criterion,
# CPU = 4 is the number of CPU used (hidden input)
project = NULL
project = snmf("genotypes.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

# plot cross-entropy criterion for all runs in the snmf project
plot(project, col = "blue", pch = 19, cex = 1.2)

# select the best run for K = 4 clusters
best = which.min(cross.entropy(project, K = 4))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
barchart(project, K = 4, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)

# Genome scan for selection: population differentiation tests
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 4)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)

# creation of a genotype matrix with missing genotypes
dat = as.numeric(tutorial.R)
dat[sample(1:length(dat), 100)] <- 9
dat <- matrix(dat, nrow = 50, ncol = 400)
write.lfmm(dat, "genoM.lfmm")

project.missing = snmf("genoM.lfmm", K = 4,
                       entropy = TRUE, repetitions = 10,
                       project = "new")

# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.missing, K = 4))

# Impute the missing genotypes
impute(project.missing, "genoM.lfmm",
       method = 'mode', K = 4, run = best)
dat.imp = read.lfmm("genoM.lfmm_imputed.lfmm")
mean( tutorial.R[dat == 9] == dat.imp[dat == 9] )


#Running environmental analysis
# Runs with K = 6 using 5 repetitions.
project = NULL
project = lfmm("genotypes.lfmm",
               "gradients.env",
               K = 6,
               repetitions = 5,
               project = "new")

# compute adjusted p-values
p = lfmm.pvalues(project, K = 6)
pvalues = p$pvalues

# GEA significance test
par(mfrow = c(2,1))
hist(pvalues, col = "lightblue")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)

for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("Expected FDR:", alpha))
  L = length(pvalues)
  # return a list of candidates with expected FDR alpha.
  # Benjamini-Hochberg's algorithm:
  w = which(sort(pvalues) < alpha * (1:L) / L)
  candidates = order(pvalues)[w]
  # estimated FDR and True Positive Rate
  Lc = length(candidates)
  estimated.FDR = sum(candidates <= 350)/Lc
  print(paste("Observed FDR:",
              round(estimated.FDR, digits = 2)))
  estimated.TPR = sum(candidates > 350)/50
  print(paste("Estimated TPR:",
              round(estimated.TPR, digits = 2)))
}
