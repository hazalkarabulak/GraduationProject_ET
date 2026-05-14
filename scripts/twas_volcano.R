.libPaths(c("~/R/library", .libPaths()))
library(ggplot2)
library(ggrepel)

twas_file <- "/mnt/c/Users/Lenovo/Downloads/ET_TWAS_all.txt"
df <- read.table(twas_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
cat(sprintf("Yuklenen satirlar: %d\n", nrow(df)))

df$GeneName <- sub(".*\\.([A-Za-z0-9][A-Za-z0-9._-]*)\\.wgt.*", "\\1", df$FILE)
df$TWAS.Z <- as.numeric(df$TWAS.Z)
df$TWAS.P <- as.numeric(df$TWAS.P)
df <- df[!is.na(df$TWAS.Z) & !is.na(df$TWAS.P), ]
df <- df[is.finite(df$TWAS.Z) & is.finite(df$TWAS.P), ]
cat(sprintf("Temizleme sonrasi: %d gen\n", nrow(df)))

df$neg_log10_p <- -log10(pmax(df$TWAS.P, .Machine$double.xmin))

n_genes <- nrow(df)
bonf_p <- 0.05 / n_genes
bonf_threshold <- -log10(bonf_p)
cat(sprintf("Bonferroni esigi: -log10(p) > %.1f\n", bonf_threshold))

df$significance <- ifelse(df$TWAS.P < bonf_p, "Bonferroni < 0.05", "NS")
df$significance <- factor(df$significance, levels=c("Bonferroni < 0.05", "NS"))
n_sig <- sum(df$significance == "Bonferroni < 0.05")
cat(sprintf("Anlamli gen: %d\n", n_sig))

smr_genes <- c("KRTCAP3", "IFT172", "NRBP1", "LINC00323", "AC107871.2")
top5 <- c()
 label_list <- smr_genes
df$label <- ifelse(df$GeneName %in% label_list, df$GeneName, NA)

p <- ggplot(df, aes(x=TWAS.Z, y=neg_log10_p, color=significance)) +
  geom_point(data=subset(df, significance=="NS"), size=1.0, alpha=0.35, shape=16) +
  geom_point(data=subset(df, significance=="Bonferroni < 0.05"), size=1.6, alpha=0.85, shape=16) +
  geom_hline(yintercept=bonf_threshold, linetype="dashed", color="#333333", linewidth=0.6) +
  geom_text_repel(aes(label=label), na.rm=TRUE, size=3.2, fontface="bold.italic",
    color="black", box.padding=0.5, point.padding=0.3, max.overlaps=25,
    segment.color="#555555", segment.size=0.3, min.segment.length=0.2) +
  scale_color_manual(values=c("Bonferroni < 0.05"="#D62728", "NS"="#AAAAAA"),
    name="Significance", guide=guide_legend(override.aes=list(size=3.5, alpha=1))) +
  labs(title="TWAS Volcano Plot - Essential Tremor",
    subtitle=sprintf("GTEx Brain Cortex | %d gen | %d Bonferroni anlamli", n_genes, n_sig),
    x="TWAS Z-score", y="-log10(P)") +
  theme_classic(base_size=13) +
  theme(
    panel.background=element_rect(fill="white", color=NA),
    plot.background=element_rect(fill="white", color=NA),
    panel.border=element_rect(color="#333333", fill=NA, linewidth=0.8),
    panel.grid.major=element_line(color="#EEEEEE", linewidth=0.4),
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=14, hjust=0.5),
    plot.subtitle=element_text(size=10, hjust=0.5, color="#555555"),
    axis.title=element_text(face="bold", size=12),
    axis.text=element_text(size=10, color="black"),
    legend.position="right",
    legend.background=element_rect(fill="white", color="#CCCCCC", linewidth=0.4),
    legend.key=element_rect(fill="white"),
    legend.title=element_text(face="bold", size=11),
    plot.margin=margin(15,15,15,15)
  )

out_dir <- "/mnt/c/Users/Lenovo/Downloads"
ggsave(file.path(out_dir, "twas_volcano_ET.png"), plot=p, width=10, height=7, dpi=300, bg="white")
ggsave(file.path(out_dir, "twas_volcano_ET.pdf"), plot=p, width=10, height=7, device=cairo_pdf)
cat("Kaydedildi: twas_volcano_ET.png / .pdf / .svg\n")
