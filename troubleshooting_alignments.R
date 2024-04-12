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

phiX <- read.csv("phiX_troubleshooting.csv", header = TRUE)
View(phiX)
summary(phiX)

popmap <- read.csv("lgpopmap_294.csv", header = T)

# combine data into one dataframe

all <- merge(s, phiX, by="sample_ID")
all <- rename(all, PERE_mapped_prop = mapped_percent,
              phiX_mapped_prop = mapped_proportion,
              PERE_total_reads = total_reads,
              PERE_mapped_reads = mapped_reads,
              PERE_QC_failed_reads = QC_failed_reads)
all <- merge(all, popmap, by="sample_ID")

#write.csv(all, "summarized_mapping_data.csv")

# Visualize

  # PERE alignments
  ggplot(data = s) + aes(x = mapped_proportion) + 
    geom_histogram() + theme_classic()
  
  d <- density(s$mapped_proportion)
  plot(d)
  
  # phiX-174 alignments
  ggplot(data = phiX) + aes(x = mapped_percent) + 
    geom_histogram() + theme_classic()
  
  phi_d <- density(phiX$mapped_percent) %>% plot()

# Plot to explore if the # of reads in a samples explains mapping proportion
  # PERE alignments
  plot(x = s$mapped_proportion, y = s$total_reads)
  
  ggplot(data = s) + aes(x = total_reads, y = mapped_proportion) + 
    geom_point(size = 1) +
    theme_classic()
  
  # Zoom in to see what's going on a low read numbers
  s_lowreads <- s %>% filter(mapped_reads < 10000)
  
  ggplot(data = s_lowreads) + aes(x = total_reads, y = mapped_proportion) + 
    geom_point() + theme_classic()
  
  # phiX-174 alignments + PERE alignments
  s <- arrange(s, desc(mapped_proportion))
  
  barplot(height = phiX$mapped_percent, names.arg = phiX$sample_ID, cex.axis = 0.05)
  barplot(height = s$mapped_proportion, col = "blue")
  
  ggplot(data = all, aes(x = sample_ID)) +
    geom_bar(aes(y = phiX_mapped_prop), stat="identity", position ="identity", alpha=.8, fill='pink', color='red') +
    geom_bar(aes(y = PERE_mapped_prop), stat="identity", position ="identity", alpha=.3, fill='lightblue', color='lightblue4')
    
  all <- mutate(all, added_prop = PERE_mapped_prop + phiX_mapped_prop) %>% arrange(desc(added_prop))
  barplot(height = all$added_prop,
          xlab = "Samples",
          ylab = "phiX + PERE read mapping proportions")
  abline(h = 1, lty = 5)
  abline(h = 0.85, lty = 5)
  abline(h = 0.9, lty = 5)
  
  # PERE mapping prop against phiX mapping prop
  plot(x = all$phiX_mapped_prop, y = all$PERE_mapped_prop, 
       xlab = "phiX mapping proportion",
       ylab = "PERE mapping proportion")
  abline(a = 1, b = -1)

# What's going on with samples with unexplained reads (ie, low combined proportions of mapping to either phiX or PERE)
all_low_prop <- filter(all, added_prop < 0.85) #only 23 samples
bad_samples <- c(all_low_prop$sample_ID, all_low_prop$population)

  # Is it that there are lots of low quality reads, so QC is throwing them out?
  
  plot(x = all_low_prop$PERE_QC_failed_reads, y = all_low_prop$added_prop,
       xlab = "# of QC failed reads from PERE alignments",
       ylab = "Added mapping proportion of reads to phiX & PERE")
    # seems like that's not the answer...

  plot(x = all_low_prop$phiX_mapped_prop, y = all_low_prop$PERE_mapped_prop)
  abline(a = 0.8, b = -1)
  
  # What do the populations look like themselves
  
  ggplot(data=all, mapping=aes(x=population, y=added_prop)) + 
    stat_summary(fun.data=mean_sdl, geom="bar") + 
    scale_x_continuous(breaks=seq(1,16,1))
  
  ggplot(data=all, mapping=aes(x=population, y=PERE_mapped_prop)) + 
    stat_summary(fun.data=mean_sdl, geom="bar") +
    scale_x_continuous(breaks=seq(1,16,1))
  
  ggplot(data=all, mapping=aes(x=population, y=phiX_mapped_prop)) + 
    stat_summary(fun.data=mean_sdl, geom="bar") +
    scale_x_continuous(breaks=seq(1,16,1))

  