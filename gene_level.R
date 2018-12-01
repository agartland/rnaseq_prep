filename <- "squant_98/quant.sf"

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

library(tximport)
txi <- tximport(c(filename), type = "salmon", tx2gene = tx2gene)

