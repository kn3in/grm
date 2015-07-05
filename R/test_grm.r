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

  ## create the full symmetric matrix
  X <- diag(diag.elem)
  X[ upper.tri(X, diag = FALSE) ] <- off.diag.elem
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

nm <- read.table("../bin/non_missing.txt", header = FALSE)


all.equal(grm$non_missing, grmt$non_missing)
all.equal(grm$grm, grmt$grm)
all.equal(grm$grm, as.matrix(grm_text), check.attributes = FALSE)
