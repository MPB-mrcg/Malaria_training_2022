
rm(list = ls())

input <- "data/cotedIvoire_minmax_DP30.vcf.gz"

library(data.table)
library(tidyverse)
library(tictoc)
library(statip)
options(scipen=999)

# Calculates the amount of missing data in a row or column where missing data is defined as a -, N, NA

missi <- function(x)
{
    sum(1*(x == "./."))/ length(x)
}

## TRIQUAD

triquad <- function(x)
{
    xx <- x[x != "./."];
    res <- 1*(sum(1*(xx=="0/0"))>0) + 1*(sum(1*(xx=="1/1"))>0) + 1*(sum(1*(xx=="0/1"))>0);
    res
}


#------------ Extracting genotypes from the raw data:
Genotypes <- 'Genotypes.txt'
expression <- '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
system(sprintf("bcftools query -f'%s' %s > %s", expression, input, Genotypes))

samples <- "sampleIDs.txt"
system(sprintf("bcftools query -l %s > %s", input, samples))

#------- Computing the missingness on SNPs and samples
#---- put the first four columns in a variable
Genotype <- fread(Genotypes, header = FALSE)
sampleIDs <- scan(samples, what = character())

first4Column <- Genotype %>% select(c(1:4))
Genotype <- Genotype %>% select(-c(1:4))

names(Genotype) <- sampleIDs
#===================
## Remove invariants
#==================
varsnp <- apply(Genotype, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, Genotype)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

#=====================================
# Remove SNPs with too much missingness
#=====================================
keepsnp <- which(misssnp < 0.8)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepsnp,]

geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#=======================================
# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

## Remove isos with too much missingness
isomissallow <- 0.8
geno <- as.data.frame(geno)
geno <- geno[,which(missiso < isomissallow)]
nrow(geno)
ncol(geno)

#===================
## Remove invariants
#==================
varsnp <- apply(geno, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

#=====================================
# Remove SNPs with too much missingness
#=====================================
keepsnp <- which(misssnp < 0.4)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepsnp,]

geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#=======================================
# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

## Remove isos with too much missingness
isomissallow <- 0.4
geno <- as.data.frame(geno)
geno <- geno[,which(missiso < isomissallow)]
nrow(geno)
ncol(geno)

#===================
## Remove invariants
#==================
varsnp <- apply(geno, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

#=====================================
# Remove SNPs with too much missingness
#=====================================
keepsnp <- which(misssnp < 0.1)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepsnp,]

geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#=======================================
# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

## Remove isos with too much missingness
isomissallow <- 0.05
geno <- as.data.frame(geno)
geno <- geno[,which(missiso < isomissallow)]
nrow(geno)
ncol(geno)

#===================
## Remove invariants
#==================
varsnp <- apply(geno, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

#=====================================
# Remove SNPs with too much missingness
#=====================================
keepsnp <- which(misssnp < 0.05)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepsnp,]

geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#=======================================
# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

## Remove isos with too much missingness
isomissallow <- 0.01
geno <- as.data.frame(geno)
geno <- geno[,which(missiso < isomissallow)]
nrow(geno)
ncol(geno)

#===================
## Remove invariants
#==================
varsnp <- apply(geno, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

#=====================================
# Remove SNPs with too much missingness
#=====================================
keepsnp <- which(misssnp == 0)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepsnp,]

geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#=======================================
# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

## Save positions and isolates to keep for downstream analysis
write.table(colnames(geno), "samplesTokeep.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(first4Column[,1:2], "snpsTokeep.txt", sep = '\t', 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

#=========================================================
#---- Removing SNPs and Isolates to discard from vcf file
#========================================================
samplesToKeep = "samplesTokeep.txt"
snpsToKeep = "snpsTokeep.txt"
maf <- 0.05

filtered.vcf <- "data/cotedIvoire_DP10_Q20_miss_maf5"

system(paste0("vcftools --gzvcf ", input,
              " --keep ", samplesToKeep, 
              " --positions ", snpsToKeep,
              " --maf ", maf,
              " --not-chr Pf3D7_API_v3", 
              " --recode --recode-INFO-all --out ", filtered.vcf))

# Remove temporally data
file.remove(Genotypes, samplesToKeep, snpsToKeep, samples)

## Rename filtered VCF file
system(paste0("mv ", filtered.vcf, ".recode.vcf", " ", filtered.vcf, ".vcf"))

## Compress and index filtered vcf file
system(paste0("bgzip ", filtered.vcf, ".vcf"))
system(paste0("tabix ", filtered.vcf, ".vcf.gz"))

