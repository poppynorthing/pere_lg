#more troubleshooting

gstacks <- read.csv("gstacks_coverage.csv", header = T)

gstacks <- gstacks[1:260,1:5]
str(gstacks)

summary(gstacks)
sd(gstacks$mean_cov)

sd(gstacks$mean_cov)/sqrt(length((gstacks$mean_cov)))

min(gstacks$mean_cov)
max(gstacks$mean_cov)

hist(gstacks$mean_cov)