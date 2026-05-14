.libPaths(c("~/R/library", .libPaths()))
library(ggplot2)
library(ggrepel)

df <- read.table("/mnt/c/Users/Lenovo/Downloads/ET_TWAS_all.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

df$GeneName <- sub(".*\\.([A-Za-z0-9][A-Za-z0-9._-]*)\\.wgt.*", "\\1", df$FILE)
df$TWAS.Z   <- as.numeric(df$TWAS.Z)
df$TWAS.P   <- as.numeric(df$TWAS.P)
df$CHR      <- as.numeric(df$CHR)
df$P0       <- as.numeric(df$P0)
df <- df[!is.na(df$TWAS.P) & !is.na(df$CHR) & !is.na(df$P0), ]
df <- df[is.finite(df$TWAS.P), ]
df <- df[order(df$CHR, df$P0), ]

n_genes        <- nrow(df)
bonf_p         <- 0.05 / n_genes
bonf_threshold <- -log10(bonf_p)
df$neg_log10_p <- -log10(pmax(df$TWAS.P, .Machine$double.xmin))
df$significant <- df$TWAS.P < bonf_p

chr_list    <- sort(unique(df$CHR))
chr_maxpos  <- tapply(df$P0, df$CHR, max)
chr_offsets <- c(0, cumsum(as.numeric(chr_maxpos[as.character(chr_list)])))
names(chr_offsets) <- c(as.character(chr_list), "end")
df$pos_cum  <- df$P0 + chr_offsets[as.character(df$CHR)]
chr_centers <- tapply(df$pos_cum, df$CHR, mean)

chr_colors <- setNames(rep(c("#2E75B6","#C0392B"), length.out=length(chr_list)),
                       as.character(chr_list))
df$color <- chr_colors[as.character(df$CHR)]

y_limit <- 50
df$clipped     <- df$neg_log10_p > y_limit
df$y_plot      <- pmin(df$neg_log10_p, y_limit)
df$point_shape <- ifelse(df$clipped, 17, 16)
df$label       <- ifelse(df$GeneName == "KRTCAP3", "KRTCAP3", NA)

n_clipped <- sum(df$clipped)

p <- ggplot(df, aes(x=pos_cum, y=y_plot)) +
  geom_point(
    data=subset(df, !clipped),
    aes(color=color),
    size=ifelse(subset(df, !clipped)$significant, 1.4, 0.8),
    alpha=ifelse(subset(df, !clipped)$significant, 0.9, 0.4),
    shape=16) +
  geom_point(
    data=subset(df, clipped),
    aes(color=color),
    size=2.0, alpha=0.9, shape=17) +
  scale_color_identity() +
  geom_hline(yintercept=bonf_threshold, linetype="dashed",
    color="black", linewidth=0.7) +
  annotate("text",
    x=max(df$pos_cum, na.rm=TRUE) * 0.01,
    y=bonf_threshold + 1.8,
    label=sprintf("Bonferroni (p=%.0e)", bonf_p),
    size=3, hjust=0, color="black") +
  annotate("text",
    x=max(df$pos_cum, na.rm=TRUE) * 0.5,
    y=y_limit - 1.5,
    label=sprintf("▲ %d gen kesildi (-log10P > %d)", n_clipped, y_limit),
    size=3, hjust=0.5, color="#555555", fontface="italic") +
  geom_text_repel(aes(label=label), na.rm=TRUE, size=3.5,
    fontface="bold.italic", color="black", box.padding=0.5,
    point.padding=0.3, max.overlaps=20, segment.color="#555555",
    segment.size=0.3, min.segment.length=0.2,
    ylim=c(0, y_limit)) +
  scale_x_continuous(
    breaks=as.numeric(chr_centers),
    labels=names(chr_centers),
    expand=c(0.01, 0.01)) +
  scale_y_continuous(
    limits=c(0, y_limit),
    expand=c(0.02, 0.02)) +
  labs(
    title="TWAS Manhattan Plot — Essential Tremor",
    subtitle=sprintf("GTEx Brain Cortex | %d gen | %d Bonferroni anlamli", n_genes, sum(df$significant)),
    x="Kromozom",
    y=expression(-log[10](italic(P)))) +
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
ggsave(file.path(out_dir, "twas_manhattan_ET.png"),
  plot=p, width=14, height=7, dpi=300, bg="white")
ggsave(file.path(out_dir, "twas_manhattan_ET.pdf"),
  plot=p, width=14, height=7, device=cairo_pdf)
cat("Manhattan kaydedildi!\n")
