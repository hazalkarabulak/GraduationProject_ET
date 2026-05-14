.libPaths(c("~/R/library", .libPaths()))
library(ggplot2)

png("/mnt/c/Users/Lenovo/Downloads/venn_smr_twas_ET.png",
  width=2400, height=1800, res=300)

par(bg="white", mar=c(3,1,4,1))
plot(0, type="n", xlim=c(-2,2), ylim=c(-1.5,1.5),
  axes=FALSE, xlab="", ylab="")

symbols(x=c(-0.5, 0.5), y=c(0,0), circles=c(1,1),
  inches=FALSE, add=TRUE,
  fg=c("#2E75B6","#D62728"),
  bg=adjustcolor(c("#2E75B6","#D62728"), alpha.f=0.4),
  lwd=2)

text(-1.1, 0,   "4",   cex=2.5, font=2, col="#1a4a7a")
text(0,    0,   "1",   cex=2.5, font=2, col="white")
text(1.1,  0,   "857", cex=2.5, font=2, col="#7a1a1a")

text(-0.5,  1.25, "eQTL-SMR", cex=1.3, font=2, col="#2E75B6")
text(-0.5,  1.05, "(5 genes)", cex=1.1, col="#2E75B6")
text( 0.5,  1.25, "TWAS",      cex=1.3, font=2, col="#D62728")
text( 0.5,  1.05, "(858 genes)",cex=1.1, col="#D62728")

text(0, -0.45, "KRTCAP3", cex=1.0, font=4, col="white")

title("eQTL-SMR and TWAS Overlap - Essential Tremor",
  cex.main=1.3, font.main=2)
mtext("Intersection: KRTCAP3 (Chr2p23.3) significant in both analyses",
  side=1, line=1.5, cex=1.0, col="#555555")

dev.off()

cairo_pdf("/mnt/c/Users/Lenovo/Downloads/venn_smr_twas_ET.pdf",
  width=8, height=6, bg="white")

par(bg="white", mar=c(3,1,4,1))
plot(0, type="n", xlim=c(-2,2), ylim=c(-1.5,1.5),
  axes=FALSE, xlab="", ylab="")

symbols(x=c(-0.5, 0.5), y=c(0,0), circles=c(1,1),
  inches=FALSE, add=TRUE,
  fg=c("#2E75B6","#D62728"),
  bg=adjustcolor(c("#2E75B6","#D62728"), alpha.f=0.4),
  lwd=2)

text(-1.1, 0,   "4",   cex=2.5, font=2, col="#1a4a7a")
text(0,    0,   "1",   cex=2.5, font=2, col="white")
text(1.1,  0,   "857", cex=2.5, font=2, col="#7a1a1a")

text(-0.5,  1.25, "eQTL-SMR",  cex=1.3, font=2, col="#2E75B6")
text(-0.5,  1.05, "(5 genes)", cex=1.1, col="#2E75B6")
text( 0.5,  1.25, "TWAS",      cex=1.3, font=2, col="#D62728")
text( 0.5,  1.05, "(858 genes)",cex=1.1, col="#D62728")

text(0, -0.45, "KRTCAP3", cex=1.0, font=4, col="white")

title("eQTL-SMR and TWAS Overlap - Essential Tremor",
  cex.main=1.3, font.main=2)
mtext("Intersection: KRTCAP3 (Chr2p23.3) significant in both analyses",
  side=1, line=1.5, cex=1.0, col="#555555")

dev.off()
cat("Venn saved!\n")
