#adegenet

#Load libraries
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

g <- read.genepop("adegenet/populations.snps.gen")
gg <- df2genind(g, sep= "_")

pop(g) = sapply(strsplit(indNames(g), "_"), function(g){g[1]})

g.pca = dudi.pca(g, scannf=FALSE)

grp <- find.clusters(g, max.n.clust=40)

dapc1 <- dapc(g, grp$grp)
scatter(dapc1)

?dudi.pca

?read.genepop

?genind
