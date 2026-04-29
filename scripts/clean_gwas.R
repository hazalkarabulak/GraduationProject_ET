library(data.table)

dt <- fread("G250_Essential_tremor_summary.txt", na.strings = c("NA", "NaN", ".", ""))
cat("Ham satır sayısı:", nrow(dt), "\n")

setnames(dt,
  old = c("rsID", "EA", "OA", "OR", "SE", "P", "UK_FREQ"),
  new = c("SNP",  "A1", "A2", "OR", "SE", "p", "freq")
)

valid_allele <- function(x) grepl("^[ACTGactg]+$", x)
dt <- dt[valid_allele(A1) & valid_allele(A2)]
dt <- dt[!is.na(OR) & !is.na(SE) & !is.na(p) & !is.na(freq)]
dt <- dt[OR > 0 & SE > 0 & p > 0 & p <= 1]
dt <- dt[freq > 0 & freq < 1]

dt[, b := log(OR)]
dt[, z := b / SE]
dt <- dt[abs(z) < 10]

setorder(dt, p)
dt <- dt[!duplicated(SNP)]

dt[, n := 64614]
out <- dt[, .(SNP, A1, A2, freq, b, se = SE, p, n)]

fwrite(out, "ET_gwas_clean.ma", sep = "\t", quote = FALSE)
cat("Temiz satır sayısı:", nrow(out), "\n")
