# Project: SardinianTrios 
# Knockoff on trios

# setwd main directory (SardinianTrios)

library(MultiPhen)
library(snpStats)
library(admixr)

# Binary files from Plink
bed <- "imputed_genotypes/QC/chr17-corretto-snps-imputed2-001.bed"
fam <- "imputed_genotypes/QC/chr17-corretto-snps-imputed2-001.fam"
bim <- "imputed_genotypes/QC/chr17-corretto-snps-imputed2-001.bim"

chr17 <- read.plink(bed, bim, fam)
chr17$fam$pedigree

# Subjects filtered with mendel check
# Corrected for mendelian errors
bed_me <- "imputed_genotypes/QC/chr17-mendel.bed"
fam_me <- "imputed_genotypes/QC/chr17-mendel.fam"
bim_me <- "imputed_genotypes/QC/chr17-mendel.bim"

chr17_me <- read.plink(bed_me, bim_me, fam_me)
chr17_me$fam$pedigree

fam.full <- chr17$fam
fam.sub <- chr17_me$fam

related <- read.delim("imputed_genotypes/related.fam", header=FALSE)
# from 556 to 512

load("../pedigree_plots/data_for_pedigreePlots.RData")
unrelated <- read.table("C:/Users/gnbal/Desktop/GitHub/SardinianTrios/data/imputed_genotypes/unrelated.fam", quote="\"", comment.char="")
unrelated$V1

id1 <- trio2$ID
ID_TRIO2 <- read.csv("C:/Users/gnbal/Desktop/GitHub/SardinianTrios/data/imputed_genotypes/ID_TRIO2.csv")
id2<- ID_TRIO2$x
id2 <- gsub("_.*", "", id2)

id <- c(id1, id2)

unrelated.sub <- fam.sub[fam.sub$pedigree %in% unrelated$V1,]
related.sub <- fam.sub[fam.sub$member %in% related$V2,]

colnames(fam.sub.df)
fam.sub.df <- fam.sub[fam.sub$member %in% id,]
stat <- fam.sub.df %>% group_by(pedigree) %>% count()

stat2<- stat[stat$n ==3,]
cleaned_trios <- stat2$pedigree

trios.sub <- fam.sub.df[fam.sub.df$pedigree %in% cleaned_trios,]

#write.table(trios.sub, file="trios_cleaned_mendel.txt")

# 92 nuclea cleaned

# NB haplo_imputed has the same order of subjects as the input (imputation does not change the order)

# ---------------------------------------------
# Variants filtered after mendel check 

# Cleaned variants after filtering for mendelian errors 
bim <- read.delim("imputed_genotypes/chr17-mendel.bim", header=FALSE)
variants <- bim$V2

##### Building genotype table and haplotype table
haplo_imputed <- read.table("imputed_haps/chr17-corretto-snps-imputed2.haps", quote="\"", comment.char="")

haplo <- as.data.frame(t(haplo_imputed[,-c(1:5)]))
colnames(haplo) <- variants
table(as.matrix(haplo), useNA = "always")

# Correcting haplotype table structure
odd <- paste0(fam.full$member, "_1")
even <- paste0(fam.full$member, "_2")
ID3 <- cbind(odd, even)


# Labeling the haplotype matrix
vettore1 <- NULL
vettore2 <- NULL
vettore3 <- NULL

for (i in 1:754){
  
  vettore1 <-ID3[i,1]
  vettore2 <- ID3[i,2]
  vettore3 <- c(vettore3, vettore1, vettore2)
   
  
}

names(vettore3)<- NULL
rownames(haplo) <- vettore3


# On the subset 426
odd.sub <- paste0(trios.sub$member, "_1")
even.sub <- paste0(trios.sub$member, "_2")
ID3.sub <- cbind(odd.sub, even.sub)

vettore1 <- NULL
vettore2 <- NULL
vettore3 <- NULL

for (i in 1:426){
  
  vettore1 <-ID3.sub[i,1]
  vettore2 <- ID3.sub[i,2]
  vettore3 <- c(vettore3, vettore1, vettore2)
  
  
}

names(vettore3)<- NULL
ordered_haplo <- haplo[(rownames(haplo) %in% vettore3),]
ordered_haplo <- ordered_haplo[, colnames(ordered_haplo) %in% variants]




# Checkpoints for haplotypes for 536 subjects ordered:
#save.image("imputed_haps/haplotypes_corrected.RData")

geno <- read.csv("imputed_genotypes/chr17-mendel.raw", sep="")
rownames(geno) <- geno$IID
geno.sub <- geno[,7:2540]

geno.sub <- geno[rownames(geno) %in% trios.sub$member,]

aa <- geno[,1:6]
aa<- aa[rownames(aa) %in% trios.sub$member,]
aa$label <- 0
table(aa$SEX)
table(aa$PHENOTYPE)

aa$label <- ifelse(aa$PAT == 0 & aa$MAT ==0 & aa$SEX ==1, 1, 0 )
aa$label <- ifelse(aa$PAT == 0 & aa$MAT ==0 & aa$SEX ==2, 2, aa$label)
aa$label <- ifelse(aa$PAT != 0 & aa$MAT !=0, 3, aa$label)

# check 3552
a2 <- aa[aa$FID == 3552,]
which(aa$FID == 3552)

aa[398,4] <- 31843
which(aa$FID == 3552) # ok corrected MID

# check 3360
a2 <- aa[aa$FID == 3360,]
which(aa$FID == 3360)

aa[300,3] <- 32056
which(aa$FID == 3360) # ok corrected MID

# Relaunch
aa$label <- 0
aa$label <- ifelse(aa$PAT == 0 & aa$MAT ==0 & aa$SEX ==1, 1, 0 )
aa$label <- ifelse(aa$PAT == 0 & aa$MAT ==0 & aa$SEX ==2, 2, aa$label)
aa$label <- ifelse(aa$PAT != 0 & aa$MAT !=0, 3, aa$label)

# Order the dataset first by group_variable =famid and then by order_variable =label
ordered_geno <- aa %>% 
  arrange(FID, label)


b <- unique(gsub("_.*", "", rownames(ordered_haplo)))
                   
b == ordered_geno$IID

# Reordering haplo
haplo_imputed2 <- read.table("imputed_haps/chr17-corretto-snps-imputed2.haps", quote="\"", comment.char="")

haplo2 <- as.data.frame(t(haplo_imputed2[,-c(1:5)]))
colnames(haplo2) <- variants
table(as.matrix(haplo2), useNA = "always")

members <- ordered_geno$IID
odd.ord <- paste0(members, "_1")
even.ord <- paste0(members, "_2")
ID3.ord <- cbind(odd.ord, even.ord)

vettore1 <- NULL
vettore2 <- NULL
vettore3 <- NULL

for (i in 1:426){
  
  vettore1 <-ID3.ord[i,1]
  vettore2 <- ID3.ord[i,2]
  vettore3 <- c(vettore3, vettore1, vettore2)
  
  
}

names(vettore3)<- NULL
ordered_haplo2 <-  ordered_haplo[match(vettore3, rownames(ordered_haplo)), ]

b <- unique(gsub("_.*", "", rownames(ordered_haplo2)))

b == ordered_geno$IID

# ordered_haplo and ordered_geno are the two final datasets

# Save the environment with the objects to an RData file
# save(ordered_geno2, file = "geno_ordered_after_mendel.Rdata")
# save(ordered_haplo2, file = "haplo_ordered_after_mendel.Rdata")
