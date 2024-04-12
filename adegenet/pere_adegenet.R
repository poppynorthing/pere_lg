#adegenet

#Load libraries
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

gind <- read.genepop("adegenet/populations.snps.gen") #creates a genind object from a genopop input file
gpop <- genind2genpop(gind) #creates a genpop object from a genind object

X <- tab(gind, NA.method="mean") #replace NA's, generate a table called 'X' for PCA
pca1 <- dudi.pca(X, scannf = F, scale = F)
temp <- as.integer(pop(gind))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]

plot(pca1$li, col=myCol, cex=3, pch=myPch)

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
