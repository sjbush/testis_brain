=head

AFTER USAGE:

# the following code produces Supplementary Figures 3 and 4: bar plots of the number of differentially expressed brain-associated genes per testicular cell cluster for each of four different single-cell expression atlases (using inclusive and conservative lists of brain-associated genes, respectively).

library(tidyverse)
library(ggpubr)
library(grid)
theme_set(theme_bw())

# Salehi 2023

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.Salehi2023.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
figS3A<-ggplot(df, aes(x=as.factor(Cell.type),y=Number,fill=Cell.type)) + geom_bar(stat="identity",fill="gray70") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle('A. Salehi 2023')

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.Salehi2023.conservative.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("Undifferentiated SSC","Differentiating SSC","Leptotene SPC","Zygotene SPC","Pachytene SPC","Diplotene SPC","Early round SPT","Round SPT-1","Round SPT-2","Elongated SPT","Myoid-1","Myoid-2","Leydig","Sertoli-1","Sertoli-2","Sertoli-3","Endothelial","Macrophage"))
figS4A<-ggplot(df, aes(x=as.factor(Cell.type),y=Number,fill=Cell.type)) + geom_bar(stat="identity",fill="gray70") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle('A. Salehi 2023')

# Wang2025-DEGsAllCells

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.Wang2025-DEGsAllCells.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
figS3B<-ggplot(df, aes(x=as.factor(Cell.type),y=Number,fill=Cell.type)) + geom_bar(stat="identity",fill="gray70") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle('B. Wang 2025 (AllCells)')

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.Wang2025-DEGsAllCells.conservative.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm","PTM cells","Fetal Leydig cells","Leydig-PTM Precursor","Leydig cells","Sertoli Precursor","Interstitial Precursor","Sertoli-Interstitial Progenitor","Sertoli cells","Red Blood cells","Immune cells","Endothelial cells","Epithelial cells","Smooth Muscle cells"))
figS4B<-ggplot(df, aes(x=as.factor(Cell.type),y=Number,fill=Cell.type)) + geom_bar(stat="identity",fill="gray70") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle('B. Wang 2025 (AllCells)')

# Wang2025-DEGsGermCells

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.Wang2025-DEGsGermCells.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
figS3C<-ggplot(df, aes(x=as.factor(Cell.type),y=Number,fill=Cell.type)) + geom_bar(stat="identity",fill="gray70") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle('C. Wang 2025 (GermCells)')

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.Wang2025-DEGsGermCells.conservative.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC","Spermatogonia","Spermatocyte","Spermatid","Sperm"))
figS4C<-ggplot(df, aes(x=as.factor(Cell.type),y=Number,fill=Cell.type)) + geom_bar(stat="identity",fill="gray70") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle('C. Wang 2025 (GermCells)')

# Wang2025-DEGsGermCells_part1

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.Wang2025-DEGsGermCells_part1.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
figS3D<-ggplot(df, aes(x=as.factor(Cell.type),y=Number,fill=Cell.type)) + geom_bar(stat="identity",fill="gray70") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle('D. Wang 2025 (GermCells_part1)')

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.Wang2025-DEGsGermCells_part1.conservative.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("PGC-1","PGC-2","Undiff1","Undiff2","Undiff3","Undiff4","Early-diff","Late-diff","Spermatocyte"))
figS4D<-ggplot(df, aes(x=as.factor(Cell.type),y=Number,fill=Cell.type)) + geom_bar(stat="identity",fill="gray70") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + ggtitle('D. Wang 2025 (GermCells_part1)')

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 3.png",width=16,height=10,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=2,ncol=2)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(figS3A, vp = define_region(row=1,col=1))
print(figS3B, vp = define_region(row=1,col=2))
print(figS3C, vp = define_region(row=2,col=1))
print(figS3D, vp = define_region(row=2,col=2))
dev.off()

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 4.png",width=16,height=10,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=2,ncol=2)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(figS4A, vp = define_region(row=1,col=1))
print(figS4B, vp = define_region(row=1,col=2))
print(figS4C, vp = define_region(row=2,col=1))
print(figS4D, vp = define_region(row=2,col=2))
dev.off()

=cut

use strict;
use warnings;

# PARAMETERS
my $conservative = 'yes'; # 'no'; # conservative criteria are to restrict analysis only to those genes associated with a defined brain phenotype: macro/megalencephaly, autism, schizophrenia, epilepsy, brain weight
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
my $in_file1 = ''; # from 2.summarise_germline_protein_expression_of_brain_genes.pl
if ($conservative eq 'no')
	{ $in_file1 = 'C:/Users/User/Desktop/testis_brain/germline_expression_at_protein_level_of_brain_associated_genes.txt'; } # from 2.summarise_germline_protein_expression_of_brain_genes.pl
elsif ($conservative eq 'yes')
	{ $in_file1 = 'C:/Users/User/Desktop/testis_brain/germline_expression_at_protein_level_of_brain_associated_genes.conservative.txt'; } # from 2.summarise_germline_protein_expression_of_brain_genes.pl
my $in_file2 = "whole_testis_atlas_DEGs/$deg_list.txt"; # manually obtained from publicly-available datasets; see Supplementary Table 4 for sources
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_file = '';
if ($conservative eq 'no')
	{ $out_file = "data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.$deg_list.txt"; }
elsif ($conservative eq 'yes')
	{ $out_file = "data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.$deg_list.conservative.txt"; }
open(OUT,'>',$out_file) or die $!;
print OUT "Cell type\tNumber\n";

# STORE THE SET OF BRAIN-ASSOCIATED GENES
my %key_genes = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if (!(defined($line[4])));
	  my $gene_name = $line[0]; my $gene_id = $line[1];
	  $key_genes{$gene_name}++;
	}
close(IN) or die $!;
my $total_key_genes = scalar keys %key_genes;
print "there are $total_key_genes key genes in total\n";

# PARSE THE DEG RESULTS OF A PREVIOUSLY PUBLISHED WHOLE TESTIS SINGLE CELL ATLAS TO DETERMINE THE EXPRESSION PROFILES OF BRAIN-ASSOCIATED GENES
# We confirm the application of Seurat's default parameters here, requiring that each gene is identified in > 25% of the cells in a given cluster, that the log2 fold change of the average expression between any two groups (i.e. cluster X and all other clusters, taken together) is > 0.25 and that the difference is statistically significant (Wilcoxon rank sum test p < 0.05). We shall also act as if the only.pos parameter was 'TRUE', requiring that each gene was more highly expressed in a given cluster compared to the set of all other clusters.
my %genes_seen = (); my %de_in = ();
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
	  next if (!(exists($key_genes{$gene_name}))); # CHECKPOINT: restrict analysis only to brain-developmental genes
	  $de_in{$cluster_name}{$gene_name}++;
	  $genes_seen{$gene_name}++;
	}
close(IN) or die $!;

# SUMMARISE THE NUMBER OF GENES DIFFERENTIALLY EXPRESSED IN EACH CLUSTER
my @cluster_names = ();
while((my $cluster_name,my $irrel)=each(%de_in))
	{ push(@cluster_names,$cluster_name); }
my @sorted_cluster_names = sort {$a cmp $b} @cluster_names;
foreach my $cluster_name (@sorted_cluster_names)
	{ my $num_degs = scalar keys %{$de_in{$cluster_name}};
	  print OUT "$cluster_name\t$num_degs\n";
	}
close(OUT) or die $!;

print "COMPLETE. Did we use a conservative gene list? $conservative\n";
my $num_genes_seen = scalar keys %genes_seen;
print "no. of genes seen (that is, significantly [adj. p < 0.05] differentially expressed in at least one cluster): $num_genes_seen\n";

exit 1;