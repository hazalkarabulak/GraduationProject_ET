.libPaths(c("~/R/library", .libPaths()))
library(ggplot2)
library(ggrepel)

df <- read.table("/mnt/c/Users/Lenovo/Downloads/ET_mQTL_significant_final.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

df$ProbeChr <- as.numeric(df$ProbeChr)
df$Probe_bp <- as.numeric(df$Probe_bp)
df$p_FDR    <- as.numeric(df$p_FDR)
df <- df[!is.na(df$p_FDR) & !is.na(df$ProbeChr) & !is.na(df$Probe_bp), ]
df <- df[order(df$ProbeChr, df$Probe_bp), ]

df$neg_log10_p <- -log10(pmax(df$p_FDR, .Machine$double.xmin))

chr_list    <- sort(unique(df$ProbeChr))
chr_maxpos  <- tapply(df$Probe_bp, df$ProbeChr, max)
chr_offsets <- c(0, cumsum(as.numeric(chr_maxpos[as.character(chr_list)])))
names(chr_offsets) <- c(as.character(chr_list), "end")
df$pos_cum  <- df$Probe_bp + chr_offsets[as.character(df$ProbeChr)]
chr_centers <- tapply(df$pos_cum, df$ProbeChr, mean)

chr_colors <- setNames(rep(c("#2E75B6","#C0392B"), length.out=length(chr_list)),
                       as.character(chr_list))
df$color <- chr_colors[as.character(df$ProbeChr)]

top_cpg <- head(df[order(df$p_FDR), "probeID"], 3)
df$label <- ifelse(df$probeID %in% top_cpg, df$probeID, NA)

cat(sprintf("Toplam CpG: %d\n", nrow(df)))
cat(sprintf("En guclu: %s (p_FDR=%.2e)\n", df$probeID[which.min(df$p_FDR)], min(df$p_FDR)))

p <- ggplot(df, aes(x=pos_cum, y=neg_log10_p)) +
  geom_point(aes(color=color), size=1.2, alpha=0.7) +
  scale_color_identity() +
  geom_text_repel(aes(label=label), na.rm=TRUE, size=3.0,
    fontface="bold.italic", color="black", box.padding=0.5,
    point.padding=0.3, max.overlaps=20, segment.color="#555555",
    segment.size=0.3, min.segment.length=0.2) +
  scale_x_continuous(
    breaks=as.numeric(chr_centers),
    labels=names(chr_centers),
    expand=c(0.01,0.01)) +
  scale_y_continuous(expand=c(0.02,0.02)) +
  labs(
    title="mQTL-SMR Manhattan Plot — Essential Tremor",
    subtitle=sprintf("Brain-mMeta | %d significant CpG probes | FDR < 0.05", nrow(df)),
    x="Chromosome",
    y=expression(-log[10](italic(FDR)))) +
  theme_classic(base_size=13) +
  theme(
    panel.background=element_rect(fill="white", color=NA),
    plot.background=element_rect(fill="white", color=NA),
    panel.border=element_rect(color="#333333", fill=NA, linewidth=0.8),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="#EEEEEE", linewidth=0.4),
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=14, hjust=0.5),
    plot.subtitle=element_text(size=10, hjust=0.5, color="#555555"),
    axis.title=element_text(face="bold", size=12),
    axis.text.x=element_text(size=8, color="black", angle=45, hjust=1),
    axis.text.y=element_text(size=10, color="black"),
    legend.position="none",
    plot.margin=margin(15,20,15,15))

out_dir <- "/mnt/c/Users/Lenovo/Downloads"
ggsave(file.path(out_dir, "mqtl_manhattan_ET.png"),
  plot=p, width=14, height=7, dpi=300, bg="white")
ggsave(file.path(out_dir, "mqtl_manhattan_ET.pdf"),
  plot=p, width=14, height=7, device=cairo_pdf)
cat("mQTL Manhattan kaydedildi!\n")
