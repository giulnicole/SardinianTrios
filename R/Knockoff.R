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

library(KnockoffTrio)
pos <- bim$V4
res_knock1 <- create_knockoff(ordered_haplo2, pos, maxcor=0.8, M=1)
res_knock2 <- create_knockoff(ordered_haplo2, pos, maxcor=0.8, M=2)
res_knock4 <- create_knockoff(ordered_haplo2, pos, maxcor=0.8, M=4)
res_knock6 <- create_knockoff(ordered_haplo2, pos, maxcor=0.8, M=6)
res_knock8 <- create_knockoff(ordered_haplo2, pos, maxcor=0.8, M=8)
res_knock10 <- create_knockoff(ordered_haplo2, pos, maxcor=0.8, M=10)# covariance matrix is not non-negative definite

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
knock_trio_res1 <- KnockoffTrio(ordered_geno2,
                                dat.ko = res_knock1,
                                pos,
                                start = 30820506,
                                end = 32483270,
                                size = c(500, 1000, 2000, 5000, 10000, 15000, 20000),
                                p_value_only = FALSE,
                                adjust_for_cov = FALSE,
                                y = residulas,
                                chr = "17",
                                xchr = FALSE,
                                sex = FALSE
)

 

knock_trio_res2 <- KnockoffTrio(ordered_geno2,
                                dat.ko = res_knock2,
                                pos,
                                start = 30820506,
                                end = 32483270,
                                size = c(500, 1000, 2000, 5000, 10000, 15000, 20000),
                                p_value_only = FALSE,
                                adjust_for_cov = FALSE,
                                y = residulas,
                                chr = "17",
                                xchr = FALSE,
                                sex = FALSE
)



knock_trio_res4 <- KnockoffTrio(ordered_geno2,
                                dat.ko = res_knock4,
                                pos,
                                start = 30820506,
                                end = 32483270,
                                size = c(500, 1000, 2000, 5000, 10000, 15000, 20000),
                                p_value_only = FALSE,
                                adjust_for_cov = FALSE,
                                y = residuals,
                                chr = "17",
                                xchr = FALSE,
                                sex = FALSE
)



knock_trio_res6 <- KnockoffTrio(ordered_geno2,
                                dat.ko = res_knock6,
                                pos,
                                start = 30820506,
                                end = 32483270,
                                size = c(500, 1000, 2000, 5000, 10000, 15000, 20000),
                                p_value_only = FALSE,
                                adjust_for_cov = FALSE,
                                y = residulas,
                                chr = "17",
                                xchr = FALSE,
                                sex = FALSE
)




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
                                sex = FALSE
)





knock_trio_res10 <- KnockoffTrio(as.matrix(ordered_geno2),
                                dat.ko = res_knock10,
                                pos,
                                start = 30820506,
                                end = 32483270,
                                size = c(500, 1000, 2000, 5000, 10000, 15000, 20000),
                                p_value_only = FALSE,
                                adjust_for_cov = FALSE,
                                y = residuals,
                                chr = "17",
                                xchr = FALSE,
                                sex = FALSE
)




# Results
library(htmlTable)
library(magrittr)

# M=1
res1 <- knock_trio_res1[knock_trio_res1$w >1,] 
res1 <- na.omit(res1)
res1.sign <- res1[res1$p < 0.01,]

res1.sign %>%  htmlTable

# M=2
res2 <- knock_trio_res2[knock_trio_res2$w >0,] 
res2 <- na.omit(res2)
res2.sign <- res2[res2$p < 0.01,]

# M=4
res4 <- knock_trio_res4[knock_trio_res4$w >0,] 
res4 <- na.omit(res4)
res4.sign <- res4[res4$p < 0.01,]

# M=6
res6 <- knock_trio_res6[knock_trio_res6$w >0,] 
res6 <- na.omit(res6)
res6.sign <- res6[res6$p < 0.01,]

# M=8
res8 <- knock_trio_res8[knock_trio_res8$w >0,] 
res8 <- na.omit(res8)
res8.sign <- res8[res8$p < 0.01,]

# M=10
res10 <- knock_trio_res10[knock_trio_res10$w >0,] 
res10 <- na.omit(res10)
res10.sign <- res10[res10$p < 0.01,]

final8 <- causal_loci(knock_trio_res8, fdr = 0.35, M=8)
window8 <- final8$window
window8 <- window8[window8$w >1,]
window8 <- na.omit(window8)
window8 <- window8[window8$p <0.05,] 

window8<- window8[window8$q<0.35,]

# Saving results
library(openxlsx)
setwd("C:/Users/gnbal/Desktop/GitHub/SardinianTrios/results_knockoff/knockoff_2534snps")
write.xlsx(window8, file="res_causal_knockoff.xlsx")

write.table(window8, file="res_causal_knockoff.txt", quote = F)

