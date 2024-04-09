# Alignment troubleshooting
# Poppy Northing
# April 2024

# Load Libraries
library(ggplot2)
library(dplyr)
library(ggthemr)

# Read in data
s <- read.csv("alignment_troubleshooting.csv", header = TRUE)
View(s)

# Look at frequency of read mapping proportions 
ggplot(data = s) + aes(x = mapped_proportion) + 
  geom_histogram() + theme_classic()

d <- density(s$mapped_proportion)
plot(d)

# Plot to explore if the # of reads in a samples explains mapping proportion
plot(x = s$mapped_proportion, y = s$total_reads)

ggplot(data = s) + aes(x = total_reads, y = mapped_proportion) + 
  geom_point(size = 1) +
  theme_classic()

# Zoom in to see what's going on a low read numbers
s_lowreads <- s %>% filter(mapped_reads < 10000)

ggplot(data = s_lowreads) + aes(x = total_reads, y = mapped_proportion) + 
  geom_point() + theme_classic()