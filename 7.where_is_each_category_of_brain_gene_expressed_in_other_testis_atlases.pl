=head

AFTER USAGE, RUN R:

# the following code produces Supplementary Figures 6, 7, 8, and 9: barplots of the number of brain-associated genes differentially expressed, at the transcript level, in one or more of cell clusters from each of four different single-cell expression atlases.

library(tidyverse)
library(ggpubr)
library(grid)
theme_set(theme_bw())

# Supplementary Figure 6: data from the 'Salehi 2023' single-cell atlas

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_macrocephaly_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6a <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("A. Macrocephaly/megalencephaly (n = 399)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_autism_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6b <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("B. Autism (n = 2562)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_schizo_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6c <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("C. Schizophrenia (n = 345)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_epilepsy_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6d <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("D. Epilepsy (n = 3372)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_iq_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6e <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("E. IQ/educational attainment (n = 2303)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_brain_weight_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6f <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("F. Brain weight (n = 879)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_har_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6g <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("G. Human-accelerated regions (n = 1608)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_2_or_more_explicit_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6h <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("H. 2 or more associations (n = 1987)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Salehi2023.data_for_making_barplot_of_germline_transcription_of_3_or_more_explicit_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
fig6i <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("I. 3 or more associations (n = 271)")

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 6.png",width=16,height=16,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=3,ncol=3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(fig6a, vp = define_region(row=1,col=1))
print(fig6b, vp = define_region(row=1,col=2))
print(fig6c, vp = define_region(row=1,col=3))
print(fig6d, vp = define_region(row=2,col=1))
print(fig6e, vp = define_region(row=2,col=2))
print(fig6f, vp = define_region(row=2,col=3))
print(fig6g, vp = define_region(row=3,col=1))
print(fig6h, vp = define_region(row=3,col=2))
print(fig6i, vp = define_region(row=3,col=3))
dev.off()

# Supplementary Figure 7: data from the 'Wang2025-DEGsAllCells' single-cell atlas

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_macrocephaly_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7a <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("A. Macrocephaly/megalencephaly (n = 399)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_autism_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7b <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("B. Autism (n = 2562)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_schizo_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7c <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("C. Schizophrenia (n = 345)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_epilepsy_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7d <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("D. Epilepsy (n = 3372)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_iq_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7e <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("E. IQ/educational attainment (n = 2303)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_brain_weight_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7f <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("F. Brain weight (n = 879)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_har_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7g <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("G. Human-accelerated regions (n = 1608)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_2_or_more_explicit_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7h <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("H. 2 or more associations (n = 1987)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsAllCells.data_for_making_barplot_of_germline_transcription_of_3_or_more_explicit_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
fig7i <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("I. 3 or more associations (n = 271)")

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 7.png",width=16,height=16,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=3,ncol=3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(fig7a, vp = define_region(row=1,col=1))
print(fig7b, vp = define_region(row=1,col=2))
print(fig7c, vp = define_region(row=1,col=3))
print(fig7d, vp = define_region(row=2,col=1))
print(fig7e, vp = define_region(row=2,col=2))
print(fig7f, vp = define_region(row=2,col=3))
print(fig7g, vp = define_region(row=3,col=1))
print(fig7h, vp = define_region(row=3,col=2))
print(fig7i, vp = define_region(row=3,col=3))
dev.off()

# Supplementary Figure 8: data from the 'Wang2025-DEGsGermCells' single-cell atlas

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_macrocephaly_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8a <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("A. Macrocephaly/megalencephaly (n = 399)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_autism_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8b <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("B. Autism (n = 2562)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_schizo_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8c <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("C. Schizophrenia (n = 345)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_epilepsy_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8d <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("D. Epilepsy (n = 3372)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_iq_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8e <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("E. IQ/educational attainment (n = 2303)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_brain_weight_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8f <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("F. Brain weight (n = 879)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_har_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8g <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("G. Human-accelerated regions (n = 1608)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_2_or_more_explicit_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8h <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("H. 2 or more associations (n = 1987)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells.data_for_making_barplot_of_germline_transcription_of_3_or_more_explicit_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
fig8i <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("I. 3 or more associations (n = 271)")

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 8.png",width=16,height=16,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=3,ncol=3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(fig8a, vp = define_region(row=1,col=1))
print(fig8b, vp = define_region(row=1,col=2))
print(fig8c, vp = define_region(row=1,col=3))
print(fig8d, vp = define_region(row=2,col=1))
print(fig8e, vp = define_region(row=2,col=2))
print(fig8f, vp = define_region(row=2,col=3))
print(fig8g, vp = define_region(row=3,col=1))
print(fig8h, vp = define_region(row=3,col=2))
print(fig8i, vp = define_region(row=3,col=3))
dev.off()

# Supplementary Figure 9: data from the 'Wang2025-DEGsGermCells_part1' single-cell atlas

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_macrocephaly_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9a <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("A. Macrocephaly/megalencephaly (n = 399)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_autism_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9b <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("B. Autism (n = 2562)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_schizo_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9c <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("C. Schizophrenia (n = 345)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_epilepsy_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9d <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("D. Epilepsy (n = 3372)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_iq_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9e <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("E. IQ/educational attainment (n = 2303)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_brain_weight_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9f <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("F. Brain weight (n = 879)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_har_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9g <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("G. Human-accelerated regions (n = 1608)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_2_or_more_explicit_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9h <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("H. 2 or more associations (n = 1987)")

df<-read.table('C:/Users/User/Desktop/testis_brain/whole_testis_atlas_DEGs/Wang2025-DEGsGermCells_part1.data_for_making_barplot_of_germline_transcription_of_3_or_more_explicit_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
fig9i <- ggplot(df, aes(x=as.factor(Cell.type),y=Number)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle("I. 3 or more associations (n = 271)")

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 9.png",width=16,height=16,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=3,ncol=3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(fig9a, vp = define_region(row=1,col=1))
print(fig9b, vp = define_region(row=1,col=2))
print(fig9c, vp = define_region(row=1,col=3))
print(fig9d, vp = define_region(row=2,col=1))
print(fig9e, vp = define_region(row=2,col=2))
print(fig9f, vp = define_region(row=2,col=3))
print(fig9g, vp = define_region(row=3,col=1))
print(fig9h, vp = define_region(row=3,col=2))
print(fig9i, vp = define_region(row=3,col=3))
dev.off()

=cut

use strict;
use warnings;

# PARAMETERS
my $deg_list = 'Salehi2023'; # 'Wang2025-DEGsAllCells'; # 'Wang2025-DEGsGermCells'; # 'Wang2025-DEGsGermCells_part1';
# cluster number-to-name lookup for Salehi 2023 (obtained from Table S2 of this paper):
my %cluster_names = (
'1' => 'Leydig',
'2' => 'Sertoli-1',
'3' => 'Endothelial',
'4' => 'Undifferentiated SSC',
'5' => 'Early round SPT',
'6' => 'Round SPT-2',
'7' => 'Elongated SPT',
'8' => 'Differentiating SSC',
'9' => 'Diplotene SPC',
'10' => 'Zygotene SPC',
'11' => 'Sertoli-3',
'12' => 'Macrophage',
'13' => 'Sertoli-2',
'14' => 'Round SPT-1',
'15' => 'Leptotene SPC',
'16' => 'Pachytene SPC',
'17' => 'Myoid-1',
'18' => 'Myoid-2'
);

# REQUIREMENTS
my $in_file1 = 'germline_expression_at_protein_level_of_brain_associated_genes.txt'; # from 8.summarise_germline_protein_expression_of_brain_genes.pl
my $in_file2 = "whole_testis_atlas_DEGs/$deg_list.txt"; # see Supplementary Table 4
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_file1 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_macrocephaly_associated_genes.txt";
my $out_file2 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_autism_associated_genes.txt";
my $out_file3 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_schizo_associated_genes.txt";
my $out_file4 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_epilepsy_associated_genes.txt";
my $out_file5 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_IQ_associated_genes.txt";
my $out_file6 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_brain_weight_associated_genes.txt";
my $out_file7 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_HAR_associated_genes.txt";
my $out_file8 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_2_or_more_explicit_genes.txt";
my $out_file9 = "whole_testis_atlas_DEGs/$deg_list.data_for_making_barplot_of_germline_transcription_of_3_or_more_explicit_genes.txt";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!; open(OUT4,'>',$out_file4) or die $!; open(OUT5,'>',$out_file5) or die $!; open(OUT6,'>',$out_file6) or die $!; open(OUT7,'>',$out_file7) or die $!; open(OUT8,'>',$out_file8) or die $!; open(OUT9,'>',$out_file9) or die $!;
print OUT1 "Cell type\tNumber\n"; print OUT2 "Cell type\tNumber\n"; print OUT3 "Cell type\tNumber\n"; print OUT4 "Cell type\tNumber\n"; print OUT5 "Cell type\tNumber\n"; print OUT6 "Cell type\tNumber\n"; print OUT7 "Cell type\tNumber\n"; print OUT8 "Cell type\tNumber\n"; print OUT9 "Cell type\tNumber\n";

# STORE THE SET OF BRAIN-ASSOCIATED GENES FOR EACH PHENOTYPE
my %brain_genes = (); my %category_per_gene = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if (!(defined($line[4])));
	  my $gene_id = $line[0]; my $macrocephaly = $line[20]; my $autism = $line[21]; my $schizo = $line[22]; my $epilepsy = $line[23]; my $iq = $line[24]; my $brain_weight = $line[25]; my $har = $line[26];
	  if ((defined($macrocephaly)) and ($macrocephaly eq 'yes')) { $brain_genes{'macrocephaly'}{$gene_id}++; $category_per_gene{$gene_id}{'macrocephaly'}++; }
	  if ((defined($autism)) 	   and ($autism 	  eq 'yes')) { $brain_genes{'autism'}{$gene_id}++; 		 $category_per_gene{$gene_id}{'autism'}++; 		 }
	  if ((defined($schizo)) 	   and ($schizo 	  eq 'yes')) { $brain_genes{'schizo'}{$gene_id}++; 		 $category_per_gene{$gene_id}{'schizo'}++; 		 }
	  if ((defined($epilepsy)) 	   and ($epilepsy 	  eq 'yes')) { $brain_genes{'epilepsy'}{$gene_id}++; 	 $category_per_gene{$gene_id}{'epilepsy'}++; 	 }
	  if ((defined($iq)) 		   and ($iq 		  eq 'yes')) { $brain_genes{'iq'}{$gene_id}++; 			 }
	  if ((defined($brain_weight)) and ($brain_weight eq 'yes')) { $brain_genes{'brain_weight'}{$gene_id}++; $category_per_gene{$gene_id}{'brain_weight'}++; }
	  if ((defined($har)) 		   and ($har 		  eq 'yes')) { $brain_genes{'har'}{$gene_id}++; 		 }
	}
close(IN) or die $!;

for(my $x=1;$x<=7;$x++)
	{ my $category;
	  if 	($x == 1) { $category = 'macrocephaly'; }
	  elsif ($x == 2) { $category = 'autism'; 		}
	  elsif ($x == 3) { $category = 'schizo'; 		}
	  elsif ($x == 4) { $category = 'epilepsy'; 	}
	  elsif ($x == 5) { $category = 'iq'; 			}
	  elsif ($x == 6) { $category = 'brain_weight'; }
	  elsif ($x == 7) { $category = 'har'; 			}
	  my $total_key_genes = scalar keys %{$brain_genes{$category}};
	  print "$category: $total_key_genes genes\n";
	  
	  # PARSE THE DEG RESULTS OF A PREVIOUSLY PUBLISHED WHOLE TESTIS SINGLE CELL ATLAS TO DETERMINE THE EXPRESSION PROFILES OF BRAIN-ASSOCIATED GENES
	  # We confirm the application of Seurat's default parameters here, requiring that each gene is identified in > 25% of the cells in a given cluster, that the log2 fold change of the average expression between any two groups (i.e. cluster X and all other clusters, taken together) is > 0.25 and that the difference is statistically significant (Wilcoxon rank sum test p < 0.05). We shall also act as if the only.pos parameter was 'TRUE', requiring that each gene was more highly expressed in a given cluster compared to the set of all other clusters.
	  my %de_in = ();
	  open(IN,$in_file2) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $avg_log2FC = $line[1]; my $pct1 = $line[2]; my $pct2 = $line[3]; my $adj_p = $line[4]; my $cluster_num = $line[5]; my $gene_name = $line[6];
		  my $diff = $pct1-$pct2;
		  my $cluster_name = $cluster_num;
		  if ($deg_list eq 'Salehi2023') { $cluster_name = $cluster_names{$cluster_num}; }
		  next unless ($adj_p < 0.05); # CHECKPOINT: restrict analysis only to genes which are significantly differentially expressed
		  next unless ($avg_log2FC > 0.25); # CHECKPOINT: restrict analysis only to genes where the log2 fold change of the average expression between any two groups (i.e. cluster X and all other clusters, taken together) is > 0.25 
		  next unless ($diff > 0); # CHECKPOINT: restrict analysis only to genes which are more highly expressed in a given cluster compared to the set of all other clusters
		  next if (!(exists($brain_genes{$category}{$gene_name}))); # CHECKPOINT: restrict analysis only to brain-developmental genes
		  $de_in{$cluster_name}{$gene_name}++;
		}
	  close(IN) or die $!;

	  # SUMMARISE THE NUMBER OF GENES DIFFERENTIALLY EXPRESSED IN EACH CLUSTER
	  my @cluster_names = ();
	  while((my $cluster_name,my $irrel)=each(%de_in))
		{ push(@cluster_names,$cluster_name); }
	  my @sorted_cluster_names = sort {$a cmp $b} @cluster_names;
	  foreach my $cluster_name (@sorted_cluster_names)
		{ my $num_degs = scalar keys %{$de_in{$cluster_name}};
		  if 	($x == 1) { print OUT1 "$cluster_name\t$num_degs\n"; }
		  elsif ($x == 2) { print OUT2 "$cluster_name\t$num_degs\n"; }
		  elsif ($x == 3) { print OUT3 "$cluster_name\t$num_degs\n"; }
		  elsif ($x == 4) { print OUT4 "$cluster_name\t$num_degs\n"; }
		  elsif ($x == 5) { print OUT5 "$cluster_name\t$num_degs\n"; }
		  elsif ($x == 6) { print OUT6 "$cluster_name\t$num_degs\n"; }
		  elsif ($x == 7) { print OUT7 "$cluster_name\t$num_degs\n"; }
		}
	}
close(OUT1) or die $!; close(OUT2) or die $!; close(OUT3) or die $!; close(OUT4) or die $!; close(OUT5) or die $!; close(OUT6) or die $!; close(OUT7) or die $!;

for(my $x=0;$x<=1;$x++)
	{ my $category;
	  if ($x == 0) { $category = '2 or more'; } elsif ($x == 1) { $category = '3 or more'; }
	  my %genes = ();
	  while((my $ens_id,my $irrel)=each(%category_per_gene))
		{ my $num = scalar keys %{$category_per_gene{$ens_id}};
		  if    (($x == 0) and ($num >= 2)) { $genes{$ens_id}++; }
		  elsif (($x == 1) and ($num >= 3)) { $genes{$ens_id}++; }
		}
	  my $total_key_genes = scalar keys %genes;
	  print "$category: $total_key_genes genes\n";
	  
	  # PARSE THE DEG RESULTS OF A PREVIOUSLY PUBLISHED WHOLE TESTIS SINGLE CELL ATLAS TO DETERMINE THE EXPRESSION PROFILES OF BRAIN-ASSOCIATED GENES
	  # We confirm the application of Seurat's default parameters here, requiring that each gene is identified in > 25% of the cells in a given cluster, that the log2 fold change of the average expression between any two groups (i.e. cluster X and all other clusters, taken together) is > 0.25 and that the difference is statistically significant (Wilcoxon rank sum test p < 0.05). We shall also act as if the only.pos parameter was 'TRUE', requiring that each gene was more highly expressed in a given cluster compared to the set of all other clusters.
	  my %de_in = ();
	  open(IN,$in_file2) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $avg_log2FC = $line[1]; my $pct1 = $line[2]; my $pct2 = $line[3]; my $adj_p = $line[4]; my $cluster_num = $line[5]; my $gene_name = $line[6];
		  my $diff = $pct1-$pct2;
		  my $cluster_name = $cluster_num;
		  if ($deg_list eq 'Salehi2023') { $cluster_name = $cluster_names{$cluster_num}; }
		  next unless ($adj_p < 0.05); # CHECKPOINT: restrict analysis only to genes which are significantly differentially expressed
		  next unless ($avg_log2FC > 0.25); # CHECKPOINT: restrict analysis only to genes where the log2 fold change of the average expression between any two groups (i.e. cluster X and all other clusters, taken together) is > 0.25 
		  next unless ($diff > 0); # CHECKPOINT: restrict analysis only to genes which are more highly expressed in a given cluster compared to the set of all other clusters
		  next if (!(exists($genes{$gene_name}))); # CHECKPOINT: restrict analysis only to brain-developmental genes
		  $de_in{$cluster_name}{$gene_name}++;
		}
	  close(IN) or die $!;

	  # SUMMARISE THE NUMBER OF GENES DIFFERENTIALLY EXPRESSED IN EACH CLUSTER
	  my @cluster_names = ();
	  while((my $cluster_name,my $irrel)=each(%de_in))
		{ push(@cluster_names,$cluster_name); }
	  my @sorted_cluster_names = sort {$a cmp $b} @cluster_names;
	  foreach my $cluster_name (@sorted_cluster_names)
		{ my $num_degs = scalar keys %{$de_in{$cluster_name}};
		  if 	($x == 0) { print OUT8 "$cluster_name\t$num_degs\n"; }
		  elsif ($x == 1) { print OUT9 "$cluster_name\t$num_degs\n"; }
		}
	}
close(OUT8) or die $!; close(OUT9) or die $!;
exit 1;