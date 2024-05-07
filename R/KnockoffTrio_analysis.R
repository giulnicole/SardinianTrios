# Project: SardinianTrios 
# Knockoff on trios considering cleaned data (426 subjects)

# setwd main directory (SardinianTrios)

load("checkpoint_knockoff/check_ordered_data.RData")
# data is already ordered and cleaned after genotype imputation
# after cleaning with PLINK for: 
# mendel 0.05, 0.05
# hwe 1e-6
# maf 0.01
# geno 0.01
# only snps (single nucleotides)
# Obtained: 2534 variants from the 2537 imputed
# Total in trios: 426 subjects cleaned

# Adjusting data
dim(ordered_geno)
order.id <- ordered_geno$IID
ordered_geno2 <- geno.sub[match(order.id, geno.sub$IID), ]

# Removing fam columns
which(ordered_geno2$IID != order.id)
ordered_trios <- ordered_geno2[,1:6]
ordered_geno2 <- ordered_geno2[,-(1:6)]

# Checks
dim(ordered_geno2)
dim(ordered_haplo2)

which(ordered_haplo2 !=1 & ordered_haplo2 !=0)
ordered_haplo2 <- as.matrix(ordered_haplo2)

# 1. Knockoof haplotypes 
library(KnockoffTrio)
pos <- bim$V4
res_knock8 <- create_knockoff(ordered_haplo2, pos, maxcor=0.8, M=8)
# res_knock10 <- create_knockoff(ordered_haplo2, pos, maxcor=0.8, M=10) # covariance matrix is not non-negative definite

sex <- as.factor(ordered_trios$SEX)
pheno <- ordered_trios$PHENOTYPE
famid <- ordered_trios$FID

pheno[pheno==1] <- 0
pheno[pheno==2] <- 1

model<- glm(pheno ~ 1 + sex, family=binomial)
residuals<- model$residuals
summary(model)

ordered_geno2 <- as.matrix(ordered_geno2)

ordered_geno2[is.na(ordered_geno2)] <- 1
ordered_geno2 <- as.matrix(ordered_geno2)

# 2. Knockoff function on genotypes
knock_trio_res8 <- KnockoffTrio(ordered_geno2,
                                dat.ko = res_knock8,
                                pos,
                                start = 30820506,
                                end = 32483270,
                                size = c(500, 1000, 2000, 5000, 10000, 15000, 20000),
                                p_value_only = FALSE,
                                adjust_for_cov = FALSE,
                                y = residuals,
                                chr = "17",
                                xchr = FALSE,
                                sex = FALSE)


# 3. Results
library(htmlTable)
library(magrittr)


# M=8 (M=10 did not converge)
res8 <- knock_trio_res8[abs(knock_trio_res8$w)>1,] 
res8 <- na.omit(res8)
res8.sign <- res8[res8$p < 0.005,]

knockoff_results<- list(res8, res8.sign)
save(knockoff_results, file= "results_from_knockoff.RData")




