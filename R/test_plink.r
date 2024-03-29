library(snpStats)
geno <- read.plink(bed = "../data/test.bed",
                   bim = "../data/test.bim",
                   fam = "../data/test.fam")

my_geno <- as(geno$genotype, "numeric")
my_geno[1:10,1:10]

mmap_geno <- read.table("../bin/geno.txt", header = FALSE, sep = " ")
na_ind <- which(mmap_geno == 3L, arr.ind = TRUE)
mmap_geno[na_ind] <- NA

sum(is.na(my_geno))
sum(is.na(mmap_geno))

dim(my_geno)
dim(mmap_geno)

bim <- read.table("../data/test.bim", stringsAsFactors = FALSE)

all.equal(bim$V2, colnames(my_geno))
colnames(mmap_geno) <- bim$V2
2 - mmap_geno[1:10,1:10]
my_geno[1:10,1:10]

all.equal(2 - as.matrix(mmap_geno), my_geno, check.attributes = FALSE)
