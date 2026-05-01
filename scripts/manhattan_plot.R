.libPaths("~/R_libs")
library(ggplot2)

dt <- read.table("~/ET_Yedek/results/ET_SMR_results_padj.txt", header=TRUE, sep="\t")
dt <- dt[!is.na(dt$p_FDR) & dt$p_FDR > 0, ]
dt$logp <- -log10(dt$p_FDR)
dt$ProbeChr <- as.numeric(dt$ProbeChr)
dt <- dt[order(dt$ProbeChr, dt$Probe_bp), ]
dt$index <- 1:nrow(dt)

# Kromozom renkleri
dt$color <- ifelse(dt$ProbeChr %% 2 == 0, "#4477AA", "#EE6677")

# Anlamlı genler için etiket
sig <- dt[dt$p_FDR < 0.05, ]

ggplot(dt, aes(x=index, y=logp, color=color)) +
  geom_point(size=1.5) +
  scale_color_identity() +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_text(data=sig, aes(label=Gene), vjust=-0.5, size=3, color="black") +
  theme_classic() +
  labs(x="Chromosome", y="-log10(FDR p-value)",
       title="Essential Tremor SMR Manhattan Plot (BrainMeta cis-eQTL)") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

ggsave("~/ET_Yedek/results/ET_SMR_manhattan.png", width=14, height=6, dpi=300)
cat("Plot kaydedildi!\n")
