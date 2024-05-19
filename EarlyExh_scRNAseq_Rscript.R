library(ggplot2)
library(ggrepel)
library(limma)
library(dplyr)
library(openxlsx)


#Figure 2B Volcano Plot
genes = read.table('~/Desktop/Tmpc_Vs_Tpex_DE.csv', header = TRUE, sep = ',')

genenames=c('IFIT1BL1', 'CCR2', 'CCR5', 'IL7R',
            'CLSPN', 'CDK1', 'FBXO5', 'PDCD1', 'CDKN2C', 'TOX', 'IKZF2',
            'ID3', 'XCL1')

sub_set = dplyr::filter(genes, gene %in% genenames)

genes$Significant <- ifelse((genes$qval > 0.05) | (abs(genes$log2fc) < 0.5), "Not Sig",
                            ifelse(genes$log2fc > 0.5, 'Up_Tpex', 'Up_Tmpc'))
ggplot(genes, aes(x = log2fc, y = -log10(qval))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("grey","#ca6702", '#005f73')) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_label_repel(
    data = sub_set,
    aes(label = gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )


#Figure 2E Gene set enrichment plot
#D28_Array2_DEGList_Chronic_Acute_raw.xlsx is one output genelist from Tox_Alfei paper

DEG = read.xlsx('~/Desktop/D28_Array2_DEGList_Chronic_Acute_raw.xlsx')
DEG = na.omit(DEG[DEG$adj.P.Val < 0.05, ][order(DEG$logFC, decreasing = T), ])

Tmpc = genes[genes$qval < 0.05 & genes$log2fc > 0.5, ]
tmp = c(na.omit(match(tolower(Tmpc$gene), tolower(DEG$GeneName))))

pdf('~/Desktop/Tpex_barcodeplot.pdf')
barcodeplot(DEG$logFC, index = tmp, 
            labels = c('Acute', 'Chronic'))
dev.off()

