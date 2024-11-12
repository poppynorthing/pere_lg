# Author: Poppy C. Northing
# Last Edited: 10NOV2024
# Title: Evaluating raw variants to determine filtering parameters

# Load libraries
library(tidyverse)

## Process stats for original VCF run (d = 250)

# Read in VCF stats files
idepth <- read_table('raw_variant_stats/pere_raw_variants.idepth') # per individual depth
imiss <- read_table('raw_variant_stats/pere_raw_variants.imiss') # per individual missing data
het <- read_table('raw_variant_stats/pere_raw_variants.het') # per individual heterozygosity
ldepth <- read_table('raw_variant_stats/pere_raw_variants.ldepth.mean') # site mean depth
lmiss <- read_table('raw_variant_stats/pere_raw_variants.lmiss') # site missing data

# Visualize depth per individual and per site

idepth %>% ggplot() + geom_histogram(aes(x=MEAN_DEPTH))
summary(idepth$MEAN_DEPTH)

ldepth %>% ggplot() + geom_histogram(aes(x=MEAN_DEPTH))
ldepth %>% ggplot() + geom_point(aes(x=POS,y=MEAN_DEPTH))
summary(ldepth$MEAN_DEPTH)

# Visualize missing data per individual and per site

imiss %>% ggplot() + geom_histogram(aes(x=F_MISS))
summary(imiss$F_MISS)

lmiss %>% ggplot() + geom_histogram(aes(x=F_MISS))
lmiss %>% ggplot() + geom_point(aes(x=POS,y=F_MISS))
summary(lmiss$F_MISS)

# Visualize heterozygosity and Fis

het %>% ggplot(aes(x=`O(HOM)`)) + geom_histogram()

  # transform into observed heterozygosity
  het %>%
    mutate(HET_O = 1-(`O(HOM)`/N_SITES)) %>%
    ggplot(aes(x=HET_O)) +
    geom_histogram()
  
het %>% ggplot(aes(x=F)) + geom_histogram()

## Look at individual depth vs individual missing data
idepth %>%
  left_join(imiss,by="INDV") %>%
  ggplot(aes(x=MEAN_DEPTH,y=F_MISS)) +
  geom_point()

# Look at heterozygosity vs individual missing data
imiss %>%
  left_join(het,by="INDV") %>%
  mutate(HET_O = 1-(`O(HOM)`/N_SITES)) %>%
  ggplot(aes(x=F_MISS,y=HET_O)) +
  geom_point()

## Do this all again but for the d=1000 data

# Read in VCF stats files
idepth_d1000 <- read_table('raw_variant_stats/pere_raw_variants.d1000.idepth') # per individual depth
imiss_d1000 <- read_table('raw_variant_stats/pere_raw_variants.d1000.imiss') # per individual missing data
het_d1000 <- read_table('raw_variant_stats/pere_raw_variants.d1000.het') # per individual heterozygosity
ldepth_d1000 <- read_table('raw_variant_stats/pere_raw_variants.d1000.ldepth.mean') # site mean depth
lmiss_d1000 <- read_table('raw_variant_stats/pere_raw_variants.d1000.lmiss') # site missing data

# Visualize depth per individual and per site

idepth_d1000 %>% ggplot() + geom_histogram(aes(x=MEAN_DEPTH))
summary(idepth_d1000$MEAN_DEPTH)

ldepth_d1000 %>% ggplot() + geom_histogram(aes(x=MEAN_DEPTH))
ldepth_d1000 %>% ggplot() + geom_point(aes(x=POS,y=MEAN_DEPTH))
summary(ldepth_d1000$MEAN_DEPTH)

# Visualize missing data per individual and per site

imiss_d1000 %>% ggplot() + geom_histogram(aes(x=F_MISS))
summary(imiss_d1000$F_MISS)

lmiss_d1000 %>% ggplot() + geom_histogram(aes(x=F_MISS))
lmiss_d1000 %>% ggplot() + geom_point(aes(x=POS,y=F_MISS))
summary(lmiss_d1000$F_MISS)

# Visualize heterozygosity and Fis

het_d1000 %>% ggplot(aes(x=`O(HOM)`)) + geom_histogram()

# transform into observed heterozygosity
het_d1000 %>%
  mutate(HET_O = 1-(`O(HOM)`/N_SITES)) %>%
  ggplot(aes(x=HET_O)) +
  geom_histogram()

het_d1000 %>% ggplot(aes(x=F)) + geom_histogram()

## Look at individual depth vs individual missing data
idepth_d1000 %>%
  left_join(imiss_d1000,by="INDV") %>%
  ggplot(aes(x=MEAN_DEPTH,y=F_MISS)) +
  geom_point()

# Look at heterozygosity vs individual missing data
imiss_d1000 %>%
  left_join(het_d1000,by="INDV") %>%
  mutate(HET_O = 1-(`O(HOM)`/N_SITES)) %>%
  ggplot(aes(x=F_MISS,y=HET_O)) +
  geom_point()

## Now, let's look at which sites and individuals have < 80/85/90% missing data

# sites 
lmiss_d1000_90 <- lmiss_d1000 %>% filter(F_MISS < 0.9)
nrow(lmiss_d1000_90) # 24,509

lmiss_d1000_85 <- lmiss_d1000 %>% filter(F_MISS < 0.85)
nrow(lmiss_d1000_85) # 16,072

lmiss_d1000_80 <- lmiss_d1000 %>% filter(F_MISS < 0.8)
nrow(lmiss_d1000_80) # 4,483

# individuals
imiss_d1000_90 <- imiss_d1000 %>% filter(F_MISS < 0.9)
nrow(imiss_d1000_90) # 69

imiss_d1000_85 <- imiss_d1000 %>% filter(F_MISS < 0.85)
nrow(imiss_d1000_85) # 62

imiss_d1000_80 <- imiss_d1000 %>% filter(F_MISS < 0.8)
nrow(imiss_d1000_80) # 55


write.csv(imiss_d1000_90, "imiss_d1000_90")
