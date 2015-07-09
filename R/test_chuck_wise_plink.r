g_files <- dir(path = "../bin", pattern = "geno[0-9]+.txt$", full.names = TRUE)

g_files <- g_files[order(as.numeric(stringr::str_replace(g_files, "../bin/geno([0-9]+).txt", "\\1")))]

all_g <- lapply(g_files, read.table, header = FALSE, sep = " ")
lapply(all_g, dim)

all_g <- do.call(cbind, all_g)
dim(all_g)
na_ind_ag <- which(all_g == 3L, arr.ind = TRUE)
all_g[na_ind_ag] <- NA




# mmap_geno <- read.table("../bin/geno.txt", header = FALSE, sep = " ")


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
sum(is.na(all_g))

dim(my_geno)
dim(mmap_geno)
dim(all_g)



all.equal(2 - as.matrix(all_g), my_geno, check.attributes = FALSE)