#adegenet

#Load libraries
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

read.structure("adegenet/populations.stru")

?read.structure

g <- read.genepop("adegenet/populations.snps.gen")

grp <- find.clusters(g, max.n.clust=40)

dapc1 <- dapc(g, grp$grp)
scatter(dapc1)