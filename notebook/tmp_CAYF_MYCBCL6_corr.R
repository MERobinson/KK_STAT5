library(tidyverse)
library(plotly)

### RNA
yf <- readRDS("../processed_data/RNAseq/STAT5_CAYF_timeseries_deres.2021-04-16.rds")
myc <- readRDS("../processed_data/RNAseq/BCL6_MYC_DE.2020-03-11.rds")
myc <- myc$NRAS_GFP_vs_Cherry
ggplot(yf, aes(x=av_expr, y=L2FC)) +
  geom_point()
ggplot(myc, aes(x=av_expr, y=-L2FC)) +
  geom_point()
df <- data.frame(gene = yf$gene_symbol,
                 CA_vs_YF = -yf$L2FC,
                 MYC_vs_BCL6 = myc$L2FC[match(yf$gene_symbol, myc$gene_symbol)])
stats <- cor.test(df$CA_vs_YF, df$MYC_vs_BCL6)
stats <- data.frame(x=-1.3, y=3, label=paste0("r = ", round(stats$estimate,3),
                                                "\np = 1e-15"))
df$col <- "B"
df[which(df$CA_vs_YF > 0.5 & df$MYC_vs_BCL6 > 0.5),]$col <- "A"
df[which(df$CA_vs_YF < -0.5 & df$MYC_vs_BCL6 < -0.5),]$col <- "C"
ggplot(df, aes(x=CA_vs_YF, y=MYC_vs_BCL6, label=gene)) +
  geom_point(aes(color=col)) +
  geom_smooth(method="lm") +
  theme_bw(base_size = 14) +
  geom_text(data=stats, aes(x=x,y=y,label=label),size=4.5) +
  xlab("STAT5-CA / STAT-YF (L2FC)") +
  ylab("MYC+ / BCL6+ (L2FC)") +
  scale_color_manual(values=c("#54823544","#66666644","#88222244"))


### Metabolites
df <- readxl::read_excel("../processed_data/metabolomics/metabolite_data_KK.xlsx")
colnames(df) <- c("fc_CA_vs_YF_48h", "fc_YF_vs_CA", "pval_CA_vs_YF", "YF_metab", "match1", 
                  "ccle_metab", "other_name", "fc_bcells_vs_solid", "pval_bcell_vs_other","match2")
df$fc_CA_vs_YF_48h <- log2(df$fc_CA_vs_YF_48h)
df$fc_YF_vs_CA <- log2(df$fc_YF_vs_CA)
df$fc_bcells_vs_solid <- log2(df$fc_bcells_vs_solid)
(stats <- cor.test(df$fc_CA_vs_YF_48h, df$fc_bcells_vs_solid))
stats <- data.frame(x=-1.3, y=0.8, label=paste0("r = ", round(stats$estimate,3),
                                              "\np = 0.797"))
df$col <- "B"
df[which(df$fc_CA_vs_YF_48h > 0.25 & df$fc_bcells_vs_solid > 0.25),]$col <- "A"
df[which(df$fc_CA_vs_YF_48h < -0.25 & df$fc_bcells_vs_solid < -0.25),]$col <- "C"
ggplot(df, aes(x=fc_CA_vs_YF_48h, y=fc_bcells_vs_solid, label=YF_metab)) +
  geom_point(aes(color=col)) +
  geom_smooth(method="lm") +
  theme_bw(base_size = 14) +
  geom_text(data=stats, aes(x=x,y=y,label=label),size=4.5) +
  xlab("STAT5-CA / STAT-YF 48h (L2FC)") +
  ylab("B-cell / solid tumor (L2FC)") +
  scale_color_manual(values=c("#54823566","#66666666","#88222266"))
