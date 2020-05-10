
# get read IDs of m6A or unm RNAs --------------------------------------------------
# RT primer 5'-TTTTTTTTTTTTTTTTTTTTGGGGGGGGGATGACGTTTCGGACGAGAAC-3' (C9A20) for m6A RNAs
# RT primer 5'-TTTTTTTTTTTTTTTTTTTTCCCCCCCCCATGACGTTTCGGACGAGAAC-3' (G9A20) for unm RNAs


library(reshape2)
library(stringr)
library(plyr)


args <- commandArgs(trailingOnly = TRUE)


# extract read ids
getReadId <- function(f, modification) {
  
  # specify barcodes
  barcode <- barcodes[modification]
  
  # reads with the barcode
  f$last <- str_sub(f$seq, start = -40, end = -1)
  d <- f[agrepl(barcode, f$last, max.distance = 1),]
  
  # read_id list
  base <- gsub(".fasta", "", args[1])
  write.table(d$id,
              paste(base, "_", barcode, ".txt", sep = ""),
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = F)
  
}




# load fasta file
f <- read.table(args[1], sep = "\t")
f <- data.frame(id = f[seq(1,nrow(f),2),], seq = f[seq(2,nrow(f),2),])
f <- cbind(seq = f$seq, colsplit(f$id, " runid=", names = c("id","c1")))
f$id <- gsub(">", "", f$id)
f$c1 <- NULL


# barcodes
barcodes <- c("unm" = "GUCAUCCC",
              "m6A" = "GUCAUGGG")
getReadId(f, "unm")
getReadId(f, "m6A")


