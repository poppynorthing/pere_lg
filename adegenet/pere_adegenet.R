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
#gind <- read.genepop("adegenet/populations.snps.gen") #creates a genind object from a genopop input file
#gpop <- genind2genpop(gind) #creates a genpop object from a genind object

# Lets use Structure data format instead
gind <- read.structure("adegenet/populations_final.stru")
gpop <- genind2genpop(gind_test)


tab <- read.table("adegenet/populations_final.stru")
gind <- df2genind(tab, ploidy = 2, ncode = 1)

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


# Discriminant Analysis of PC ---------------------------------------------
grp <- find.clusters(gind_test, max.n.clust=20)
table(pop(gind_test), grp$grp)

table.value(table(pop(gind_test), grp$grp), col.lab=paste("inf", 1:6),
            row.lab=paste("ori", 1:6))

dapc1 <- dapc(gind_test, grp$grp)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)

scatter(dapc1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:14))
