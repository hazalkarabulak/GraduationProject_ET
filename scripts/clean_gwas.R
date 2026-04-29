library(data.table)
dt <- fread("data/raw/G250_Essential_tremor_summary.txt", 
            select = c(4, 8, 7, 17, 25, 26, 27), 
            col.names = c("SNP", "A1", "A2", "freq", "OR", "se", "p"))
dt <- dt[SNP != "." & !is.na(SNP) & OR > 0 & se > 0 & p > 0]
dt[, b := log(OR)]
dt[, n := 64614]
fwrite(dt[, .(SNP, A1, A2, freq, b, se, p, n)], "data/clean/ET_gwas_clean.ma", sep = "\t", quote = FALSE)
