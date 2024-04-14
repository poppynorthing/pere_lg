#adegenet

#Load libraries
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library(hierfstat)
library(ade4)

#Read in SNP data from STACKS
gind <- read.genepop("adegenet/populations.snps.gen") #creates a genind object from a genopop input file
gpop <- genind2genpop(gind) #creates a genpop object from a genind object

# test with structure data format instead
gind_test <- read.structure("adegenet/populations.stru")
gpop_test <- genind2genpop(gind_test)

gind_sum <- summary(gind_test)
gpop_sum <- summary(gpop_test)

# Lets look at some summary stats for the data: 
par(mfrow=c(2,2))

plot(gind_sum$n.by.pop, gind_sum$pop.n.all, xlab="Pop sample size",
     ylab="Number of alleles", main="Allele numbers and sample size", type="n")
text(gind_sum$n.by.pop, gind_sum$pop.n.all, lab=names(gind_sum$n.by.pop))

barplot(gind_sum$loc.n.all, ylab = "Number of alleles", 
        main = "Number of Alleles per Locus")

barplot(gind_sum$Hexp-gind_sum$Hobs, main="Heterozygosity: expected-observed", 
        ylab="Hexp-Hobs")

barplot(gind_sum$n.by.pop, main="Sample sizes per pop", 
        ylab="Number of genotypes", las=3)

# Is mean obsv H siginificantly lower than mean exp H?
bartlett.test(list(gind_sum$Hexp, gind_sum$Hobs))
t.test(gind_sum$Hexp, gind_sum$Hobs, pair=T,var.equal=TRUE,alter="greater")

  # Yes, it is!

# Calculate F-statistics
# one way is to use the wc() function from hierfstat package
wc(gind_test) #FST = 0.6976778 FIS = 0.9348866

# To generate Fst stats per locus, use Fst()
ftab <- Fst(as.loci(gind_test))
ftab
colMeans(ftab)

# Generate F-stats confidence intervals
nc <- genind2hierfstat(gind_test)
boot.vc(nc[1], nc[-1])$ci

# Generate a pairwise Fst matrix
matFst <- genet.dist(gind_test, method = "Nei87")
matFst #has missing values
is.euclid(matFst) # is not euclidian

# Now lets look at inbreeding
temp <- inbreeding(gind_test, N=100)
class(temp)

par(mfrow=c(1,1))

Fbar <- sapply(temp, mean)
hist(Fbar, col="firebrick", main="Average inbreeding in P. recurvata")


# Multivariate analysis
X_test <- tab(gind_test, freq = TRUE, NA.method="mean")
pca1 <- dudi.pca(X_test, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

s.label(pca1$li)
title("PCA of PERE dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(gind_test))
title("PCA of PERE dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li,pop(gind_test),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA of microbov dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)

col <- funky(15)
s.class(pca1$li, pop(gind_test),xax=1,yax=3, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE, addaxes = TRUE)
abline(v=0,h=0,col="grey", lty=2)

colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of PERE dataset\naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

# Generate a PCA of the SNP data
X <- tab(gind, freq = TRUE, NA.method="mean") #replace NA's with the mean allele frequency, generate a table called 'X' for PCA
pca1 <- dudi.pca(X, scannf = F, scale = F)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
temp <- as.integer(pop(gind))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(2,4)[temp]

plot(pca1$li, col=myCol, cex=5, pch=myPch)

abline(h=0,v=0,col="grey",lty=2)
s.arrow(pca1$c1*.5, add.plot=TRUE)
legend("topright", pch=c(15,17), col=transp(c("blue","red"),.7),
       leg=c("Group A","Group B"), pt.cex=2)


pop(g) = sapply(strsplit(indNames(g), "_"), function(g){g[1]})

g.pca = dudi.pca(g, scannf=FALSE)

grp <- find.clusters(g, max.n.clust=40)

dapc1 <- dapc(g, grp$grp)
scatter(dapc1)

?dudi.pca

?read.genepop

?genind
