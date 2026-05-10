library(ggplot2)
library(ggrepel)

df <- read.table("results/ET_SMR_results_padj.txt", header=TRUE, sep="\t")
df$z_score <- df$b_SMR / df$se_SMR
df$significant <- ifelse(df$p_FDR < 0.05 & !is.na(df$p_HEIDI) & df$p_HEIDI > 0.05, "FDR < 0.05\n(HEIDI > 0.05)", "NS")

ggplot(df, aes(x=z_score, y=-log10(p_FDR), color=significant)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("FDR < 0.05\n(HEIDI > 0.05)"="red", "NS"="grey60")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  geom_text_repel(data=subset(df, p_FDR < 0.05 & !is.na(p_HEIDI) & p_HEIDI > 0.05),
    aes(label=Gene), size=4, color="black", fontface="bold",
    box.padding=0.5, max.overlaps=20) +
  labs(title="eQTL-SMR Volcano Plot - Essential Tremor",
    x="Z-score (b_SMR / se_SMR)", y="-log10(FDR p-value)", color="Significance") +
  theme_bw(base_size=13)

ggsave("results/ET_SMR_volcano_final.png", dpi=300, width=10, height=8)
