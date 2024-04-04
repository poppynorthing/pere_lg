#Alignment troubleshooting

library(ggplot2)
library(dplyr)

s <- read.csv("alignment_troubleshooting.csv", header = TRUE)
View(s)

ggplot(data = s) + aes(x = mapped_proportion, y = total_reads) + 
  geom_point()

ggplot(data = s) + aes(x = mapped_proportion, y = mapped_reads) + 
  geom_point()

s_lowreads <- s %>% filter(mapped_reads < 10000)

ggplot(data = s_lowreads) + aes(x = mapped_proportion, y = total_reads) + 
  geom_point()

ggplot(data = s_lowreads) + aes(x = mapped_proportion, y = mapped_reads) + 
  geom_point()

ggplot(data = s) + aes(x = mapped_proportion) + 
  geom_histogram()

d <- density(s$mapped_proportion)
plot(d)

#write.csv(x=summary(s), file = "summary.csv")