.libPaths(c("~/R/library", .libPaths()))
library(ggplot2)

smr <- data.frame(
  Gene     = c("LINC00323", "AC107871.2", "IFT172", "NRBP1", "KRTCAP3"),
  Chr      = c("Chr21", "Chr15", "Chr2", "Chr2", "Chr2"),
  b_SMR    = c(-0.171876, -0.056438, 0.26471, -0.423646, 0.323677),
  se_SMR   = c(0.0190489, 0.0129924, 0.0615227, 0.0997591, 0.0762896),
  p_FDR    = c(2.47e-15, 0.04967, 0.04967, 0.04967, 0.04967),
  p_HEIDI  = c(0.249339, 0.2269533, 0.1941502, 0.3920335, 0.05572788)
)

smr$CI_low  <- smr$b_SMR - 1.96 * smr$se_SMR
smr$CI_high <- smr$b_SMR + 1.96 * smr$se_SMR
smr$Gene    <- factor(smr$Gene, levels=rev(c("LINC00323","AC107871.2","IFT172","NRBP1","KRTCAP3")))
smr$direction <- ifelse(smr$b_SMR > 0, "Positive", "Negative")

p <- ggplot(smr, aes(x=b_SMR, y=Gene, color=direction)) +
  geom_vline(xintercept=0, linetype="dashed", color="#888888", linewidth=0.6) +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high),
    height=0.25, linewidth=0.8) +
  geom_point(size=4, shape=18) +
  scale_color_manual(values=c("Positive"="#D62728", "Negative"="#2E75B6"),
    name="Effect direction") +
  scale_y_discrete(labels=c(
    "KRTCAP3"    = "KRTCAP3 (Chr2)",
    "NRBP1"      = "NRBP1 (Chr2)",
    "IFT172"     = "IFT172 (Chr2)",
    "AC107871.2" = "AC107871.2 (Chr15)",
    "LINC00323"  = "LINC00323 (Chr21)")) +
  labs(
    title    = "eQTL-SMR Significant Genes — Essential Tremor",
    subtitle = "BrainMeta eQTL | FDR < 0.05 | HEIDI > 0.05 | Error bars: 95% CI",
    x        = "SMR effect size (b_SMR)",
    y        = NULL) +
  theme_classic(base_size=13) +
  theme(
    panel.background=element_rect(fill="white", color=NA),
    plot.background=element_rect(fill="white", color=NA),
    panel.border=element_rect(color="#333333", fill=NA, linewidth=0.8),
    panel.grid.major.x=element_line(color="#EEEEEE", linewidth=0.4),
    panel.grid.major.y=element_blank(),
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=14, hjust=0.5),
    plot.subtitle=element_text(size=10, hjust=0.5, color="#555555"),
    axis.title=element_text(face="bold", size=12),
    axis.text=element_text(size=11, color="black"),
    legend.position="right",
    legend.background=element_rect(fill="white", color="#CCCCCC", linewidth=0.4),
    legend.title=element_text(face="bold", size=11),
    plot.margin=margin(15,15,15,15))

out_dir <- "/mnt/c/Users/Lenovo/Downloads"
ggsave(file.path(out_dir, "smr_forest_ET.png"),
  plot=p, width=9, height=6, dpi=300, bg="white")
ggsave(file.path(out_dir, "smr_forest_ET.pdf"),
  plot=p, width=9, height=6, device=cairo_pdf)
cat("Forest plot saved!\n")
