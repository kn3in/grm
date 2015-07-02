read_GRMBin <- function(prefix, size = 4){
  sum_i <- function(i){
    return(sum(1:i))
  }

  ## open file connections and read in data
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb")
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile <- file(NFileName, "rb")

  ## read in the number of SNPs used to calculate the GRM (does not appear to work)
  N <- readBin(NFile, n = n*(n+1)/2, what = numeric(0), size=size)
  i <- sapply(1:n, sum_i)
  
  ## clean up file connections
  close(BinFile)
  close(NFile)

  ## pull apart diagonal and lower triagular elements
  diag.elem <- grm[i]
  off.diag.elem <- grm[-i]

  ## create the full symmetric correlation matrix
  X <- diag(diag.elem)
  X[ lower.tri(X, diag = FALSE) ] <- off.diag.elem
  X <- X + t(X) - diag(diag(X)) 

  ## add sample IDs to rownames and colnames
  rownames(X) <- id$V2
  colnames(X) <- id$V2

  ##
  diag_elem_N <- N[i]
  off_diag_elem_N <- N[-i]
  Y <- diag(diag_elem_N)
  Y[ lower.tri(Y, diag = FALSE) ] <- off_diag_elem_N
  Y <- Y + t(Y) - diag(diag(Y))

  list(grm = X, non_missing = Y)
}

grm <- read_GRMBin("../bin/a")
grmt <- read_GRMBin("../data/gcta")

grm_text <- read.table("../bin/grm.txt")
grm_text <- as.matrix(grm_text)
diff <- grm_text - grmt
hist(diff[lower.tri(diff, diag = FALSE)])

which(diff == max(diff), arr.ind = TRUE)

which(diff == max(diff[lower.tri(diff, diag = FALSE)]), arr.ind = TRUE)

2553 1860

hist(grm_text[ ,1860])


mmap_geno <- read.table("../bin/geno.txt", header = FALSE, sep = " ")
mmap_geno[1:10,1:10]
na_ind <- which(mmap_geno == 3L, arr.ind = TRUE)
mmap_geno[na_ind] <- NA

bim <- read.table("../data/test.bim", stringsAsFactors = FALSE)
colnames(mmap_geno) <- bim$V2

mmap_geno[1:4,1:4]
table(is.na(mmap_geno[1860, ]))

sum(mmap_geno[2553, ] * mmap_geno[1860, ], na.rm = TRUE)




maf <- function(geno) {
  freq <- table(geno) / sum(!is.na(geno))
  fq <- freq[1] + 0.5 * freq[2]
  ifelse(fq > 0.5, 1 - fq, fq)
}

mean(apply(mmap_geno, 2, maf) < 0.05)




n <- 3925
tst <- readBin("../bin/a.grm.N.bin", n = n*(n+1)/2, what = numeric(0), size = 4)




nm <- read.table("../bin/non_missing.txt", header = FALSE)

index_nm_eq <- nm == grm$non_missing

plot(grm$grm[index_nm_eq] / grm$non_missing[index_nm_eq], grmt$grm[index_nm_eq])










