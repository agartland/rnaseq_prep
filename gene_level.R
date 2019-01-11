filename <- "squant_98/quant.sf"

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

library(tximport)
library(readr)
tx2gene <- read_csv(file.path("/home/agartlan/fast/LamarAndrew/UCSC_h38", "refseq2gene.csv"))
txi <- tximport(c(filename), type = "salmon", tx2gene = tx2gene)

