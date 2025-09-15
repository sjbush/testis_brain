=head

AFTER USAGE, RUN R:

# the following code produces Supplementary Figure 5: barplots of the number of brain-associated genes differentially expressed, at the transcript level, in one or more of ten testicular cell clusters from our own single-cell expression atlas

library(tidyverse)
library(ggpubr)
library(grid)
theme_set(theme_bw())

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_macrocephaly_associated_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5a <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("A. Macrocephaly/megalencephaly (n = 399)")

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_autism_associated_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5b <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("B. Autism (n = 2562)")

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_schizo_associated_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5c <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("C. Schizophrenia (n = 345)")

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_epilepsy_associated_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5d <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("D. Epilepsy (n = 3372)")

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_iq_associated_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5e <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("E. IQ/educational attainment (n = 2303)")

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_brain_weight_associated_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5f <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("F. Brain weight (n = 879)")

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_har_associated_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5g <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("G. Human-accelerated regions (n = 1608)")

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_2_or_more_explicit_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5h <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("H. 2 or more associations (n = 1987)")

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_barplot_of_germline_transcription_of_3_or_more_explicit_genes.txt',header=T,sep='\t')
df$Gene.category <- factor(df$Gene.category, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
fig5i <- ggplot(df, aes(x=as.factor(Gene.category),y=Number,fill=Gene.category)) + geom_bar(stat="identity") + xlab('Cell cluster in which gene is differentially expressed') + ylab('Number of genes') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + scale_fill_manual(values=c("undiff SPG" = "red", "diff SPG" = "orange", "spermatocyte" = "yellow", "early spermatid 1" = "green", "early spermatid 2" = "forestgreen", "late spermatid 1" = "cornflowerblue", "late spermatid 2" = "purple", "myoid/Leydig" = "chocolate3", "Sertoli" = "gray70", "endothelia" = "coral4")) + ggtitle("I. 3 or more associations (n = 271)")

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 5.png",width=16,height=16,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=3,ncol=3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(fig5a, vp = define_region(row=1,col=1))
print(fig5b, vp = define_region(row=1,col=2))
print(fig5c, vp = define_region(row=1,col=3))
print(fig5d, vp = define_region(row=2,col=1))
print(fig5e, vp = define_region(row=2,col=2))
print(fig5f, vp = define_region(row=2,col=3))
print(fig5g, vp = define_region(row=3,col=1))
print(fig5h, vp = define_region(row=3,col=2))
print(fig5i, vp = define_region(row=3,col=3))
dev.off()

=cut

use strict;
use warnings;

# REQUIREMENTS
my $in_file1 = 'germline_expression_at_protein_level_of_brain_associated_genes.txt'; # from 2.summarise_germline_protein_expression_of_brain_genes.pl
my $in_file2 = 'results/human_adult.atlas.txt'; # produced by https://github.com/sjbush/spg_atlas/12a.create_summary_table_of_whole_testes_atlas.pl
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_file1 = 'data_for_making_barplot_of_germline_transcription_of_macrocephaly_associated_genes.txt';
my $out_file2 = 'data_for_making_barplot_of_germline_transcription_of_autism_associated_genes.txt';
my $out_file3 = 'data_for_making_barplot_of_germline_transcription_of_schizo_associated_genes.txt';
my $out_file4 = 'data_for_making_barplot_of_germline_transcription_of_epilepsy_associated_genes.txt';
my $out_file5 = 'data_for_making_barplot_of_germline_transcription_of_IQ_associated_genes.txt';
my $out_file6 = 'data_for_making_barplot_of_germline_transcription_of_brain_weight_associated_genes.txt';
my $out_file7 = 'data_for_making_barplot_of_germline_transcription_of_HAR_associated_genes.txt';
my $out_file8 = 'data_for_making_barplot_of_germline_transcription_of_2_or_more_explicit_genes.txt';
my $out_file9 = 'data_for_making_barplot_of_germline_transcription_of_3_or_more_explicit_genes.txt';
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!; open(OUT4,'>',$out_file4) or die $!; open(OUT5,'>',$out_file5) or die $!; open(OUT6,'>',$out_file6) or die $!; open(OUT7,'>',$out_file7) or die $!; open(OUT8,'>',$out_file8) or die $!; open(OUT9,'>',$out_file9) or die $!;
print OUT1 "Gene category\tNumber\n"; print OUT2 "Gene category\tNumber\n"; print OUT3 "Gene category\tNumber\n"; print OUT4 "Gene category\tNumber\n"; print OUT5 "Gene category\tNumber\n"; print OUT6 "Gene category\tNumber\n"; print OUT7 "Gene category\tNumber\n"; print OUT8 "Gene category\tNumber\n"; print OUT9 "Gene category\tNumber\n";

# STORE THE SET OF BRAIN-ASSOCIATED GENES FOR EACH PHENOTYPE
my %brain_genes = (); my %category_per_gene = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if (!(defined($line[4])));
	  my $gene_id = $line[1]; my $macrocephaly = $line[20]; my $autism = $line[21]; my $schizo = $line[22]; my $epilepsy = $line[23]; my $iq = $line[24]; my $brain_weight = $line[25]; my $har = $line[26];
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
	  
	  # PARSE THE WHOLE TESTIS SINGLE CELL ATLAS TO DETERMINE THE EXPRESSION PROFILES OF BRAIN-ASSOCIATED GENES
	  my $de_in_germline = 0; my $de_in_somatic = 0;my $de_in_myoid = 0; my $de_in_sertoli = 0; my $de_in_endothelia = 0; my $de_in_undiff_SPG = 0; my $de_in_diff_SPG = 0; my $de_in_scyte = 0; my $de_in_early_stid1 = 0; my $de_in_early_stid2 = 0; my $de_in_late_stid1 = 0; my $de_in_late_stid2 = 0;
	  open(IN,$in_file2) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $gene_name = $line[0]; my $ens_id = $line[1];
		  next if (!(exists($brain_genes{$category}{$ens_id}))); # CHECKPOINT: restrict analysis only to brain-developmental genes
		  my $clusters = $line[49];
		  if ($clusters =~ /Cluster 0 \(late spermatid 1\)/)  { $de_in_late_stid1++;  $de_in_germline++; }
		  if ($clusters =~ /Cluster 1 \(early spermatid 1\)/) { $de_in_early_stid1++; $de_in_germline++; }
		  if ($clusters =~ /Cluster 2 \(early spermatid 2\)/) { $de_in_early_stid2++; $de_in_germline++; }
		  if ($clusters =~ /Cluster 3 \(myoid\/Leydig\)/) 	  { $de_in_myoid++; 	  $de_in_somatic++;  }
		  if ($clusters =~ /Cluster 4 \(spermatocyte\)/) 	  { $de_in_scyte++; 	  $de_in_germline++; }
		  if ($clusters =~ /Cluster 5 \(late spermatid 2\)/)  { $de_in_late_stid2++;  $de_in_germline++; }
		  if ($clusters =~ /Cluster 6 \(SSC\)/) 			  { $de_in_undiff_SPG++;  $de_in_germline++; }
		  if ($clusters =~ /Cluster 7 \(endothelia\)/) 		  { $de_in_endothelia++;  $de_in_somatic++;  }
		  if ($clusters =~ /Cluster 8 \(spermatogonia\)/)	  { $de_in_diff_SPG++;    $de_in_germline++; }
		  if ($clusters =~ /Cluster 9 \(Sertoli\)/) 		  { $de_in_sertoli++;     $de_in_somatic++;  }	  
		}
	  close(IN) or die $!;

	  my $pct_in_myoid 	     = sprintf("%.2f",(($de_in_myoid/$total_key_genes)*100));
	  my $pct_in_sertoli     = sprintf("%.2f",(($de_in_sertoli/$total_key_genes)*100));
	  my $pct_in_endothelia  = sprintf("%.2f",(($de_in_endothelia/$total_key_genes)*100));
	  my $pct_in_undiff_SPG  = sprintf("%.2f",(($de_in_undiff_SPG/$total_key_genes)*100));
	  my $pct_in_diff_SPG    = sprintf("%.2f",(($de_in_diff_SPG/$total_key_genes)*100));
	  my $pct_in_scyte 	     = sprintf("%.2f",(($de_in_scyte/$total_key_genes)*100));
	  my $pct_in_early_stid1 = sprintf("%.2f",(($de_in_early_stid1/$total_key_genes)*100));
	  my $pct_in_early_stid2 = sprintf("%.2f",(($de_in_early_stid2/$total_key_genes)*100));
	  my $pct_in_late_stid1  = sprintf("%.2f",(($de_in_late_stid1/$total_key_genes)*100));
	  my $pct_in_late_stid2  = sprintf("%.2f",(($de_in_late_stid2/$total_key_genes)*100));

	  print "$category: $total_key_genes genes\n";
	  print "number of genes differentially expressed in undiff SPG = $de_in_undiff_SPG ($pct_in_undiff_SPG%)\n";
	  print "number of genes differentially expressed in diff SPG = $de_in_diff_SPG ($pct_in_diff_SPG%)\n";
	  print "number of genes differentially expressed in spermatocyte = $de_in_scyte ($pct_in_scyte%)\n";
	  print "number of genes differentially expressed in early spermatid 1 = $de_in_early_stid1 ($pct_in_early_stid1%)\n";
	  print "number of genes differentially expressed in early spermatid 2 = $de_in_early_stid2 ($pct_in_early_stid2%)\n";
	  print "number of genes differentially expressed in late spermatid 1 = $de_in_late_stid1 ($pct_in_late_stid1%)\n";
	  print "number of genes differentially expressed in late spermatid 2 = $de_in_late_stid2 ($pct_in_late_stid2%)\n";
	  print "number of genes differentially expressed in endothelia/macrophages = $de_in_endothelia ($pct_in_endothelia%)\n";
	  print "number of genes differentially expressed in myoid/Leydig cells = $de_in_myoid ($pct_in_myoid%)\n";
	  print "number of genes differentially expressed in Sertoli cells = $de_in_sertoli ($pct_in_sertoli%)\n";
	  print "~~\n";
	  
	  # SUMMARISE THE KEY FIGURES IN ORDER TO MAKE A BARPLOT SUMMARISING THE GERMLINE TRANSCRIPTION OF BRAIN-ASSOCIATED GENES, I.E. NUMBER OF GENES DIFFERENTIALLY EXPRESSED IN A GIVEN CLUSTER
	  if ($x == 1)
		{ print OUT1 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT1 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT1 "spermatocyte\t$de_in_scyte\n";
		  print OUT1 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT1 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT1 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT1 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT1 "endothelia\t$de_in_endothelia\n";
		  print OUT1 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT1 "Sertoli\t$de_in_sertoli\n";
		}
	  elsif ($x == 2)
		{ print OUT2 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT2 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT2 "spermatocyte\t$de_in_scyte\n";
		  print OUT2 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT2 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT2 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT2 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT2 "endothelia\t$de_in_endothelia\n";
		  print OUT2 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT2 "Sertoli\t$de_in_sertoli\n";
		}
	  elsif ($x == 3)
		{ print OUT3 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT3 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT3 "spermatocyte\t$de_in_scyte\n";
		  print OUT3 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT3 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT3 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT3 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT3 "endothelia\t$de_in_endothelia\n";
		  print OUT3 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT3 "Sertoli\t$de_in_sertoli\n";
		}
	  elsif ($x == 4)
		{ print OUT4 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT4 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT4 "spermatocyte\t$de_in_scyte\n";
		  print OUT4 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT4 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT4 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT4 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT4 "endothelia\t$de_in_endothelia\n";
		  print OUT4 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT4 "Sertoli\t$de_in_sertoli\n";
		}
	  elsif ($x == 5)
		{ print OUT5 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT5 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT5 "spermatocyte\t$de_in_scyte\n";
		  print OUT5 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT5 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT5 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT5 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT5 "endothelia\t$de_in_endothelia\n";
		  print OUT5 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT5 "Sertoli\t$de_in_sertoli\n";
		}
	  elsif ($x == 6)
		{ print OUT6 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT6 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT6 "spermatocyte\t$de_in_scyte\n";
		  print OUT6 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT6 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT6 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT6 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT6 "endothelia\t$de_in_endothelia\n";
		  print OUT6 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT6 "Sertoli\t$de_in_sertoli\n";
		}
	  elsif ($x == 7)
		{ print OUT7 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT7 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT7 "spermatocyte\t$de_in_scyte\n";
		  print OUT7 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT7 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT7 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT7 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT7 "endothelia\t$de_in_endothelia\n";
		  print OUT7 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT7 "Sertoli\t$de_in_sertoli\n";
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
	  
	  # PARSE THE WHOLE TESTIS SINGLE CELL ATLAS TO DETERMINE THE EXPRESSION PROFILES OF BRAIN-ASSOCIATED GENES
	  my $de_in_germline = 0; my $de_in_somatic = 0;my $de_in_myoid = 0; my $de_in_sertoli = 0; my $de_in_endothelia = 0; my $de_in_undiff_SPG = 0; my $de_in_diff_SPG = 0; my $de_in_scyte = 0; my $de_in_early_stid1 = 0; my $de_in_early_stid2 = 0; my $de_in_late_stid1 = 0; my $de_in_late_stid2 = 0;
	  open(IN,$in_file2) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $gene_name = $line[0]; my $ens_id = $line[1];
		  next if (!(exists($genes{$ens_id}))); # CHECKPOINT: restrict analysis only to brain-developmental genes
		  my $clusters = $line[49];
		  if ($clusters =~ /Cluster 0 \(late spermatid 1\)/)  { $de_in_late_stid1++;  $de_in_germline++; }
		  if ($clusters =~ /Cluster 1 \(early spermatid 1\)/) { $de_in_early_stid1++; $de_in_germline++; }
		  if ($clusters =~ /Cluster 2 \(early spermatid 2\)/) { $de_in_early_stid2++; $de_in_germline++; }
		  if ($clusters =~ /Cluster 3 \(myoid\/Leydig\)/) 	  { $de_in_myoid++; 	  $de_in_somatic++;  }
		  if ($clusters =~ /Cluster 4 \(spermatocyte\)/) 	  { $de_in_scyte++; 	  $de_in_germline++; }
		  if ($clusters =~ /Cluster 5 \(late spermatid 2\)/)  { $de_in_late_stid2++;  $de_in_germline++; }
		  if ($clusters =~ /Cluster 6 \(SSC\)/) 			  { $de_in_undiff_SPG++;  $de_in_germline++; }
		  if ($clusters =~ /Cluster 7 \(endothelia\)/) 		  { $de_in_endothelia++;  $de_in_somatic++;  }
		  if ($clusters =~ /Cluster 8 \(spermatogonia\)/)	  { $de_in_diff_SPG++;    $de_in_germline++; }
		  if ($clusters =~ /Cluster 9 \(Sertoli\)/) 		  { $de_in_sertoli++;     $de_in_somatic++;  }	  
		}
	  close(IN) or die $!;

	  my $pct_in_myoid 	     = sprintf("%.2f",(($de_in_myoid/$total_key_genes)*100));
	  my $pct_in_sertoli     = sprintf("%.2f",(($de_in_sertoli/$total_key_genes)*100));
	  my $pct_in_endothelia  = sprintf("%.2f",(($de_in_endothelia/$total_key_genes)*100));
	  my $pct_in_undiff_SPG  = sprintf("%.2f",(($de_in_undiff_SPG/$total_key_genes)*100));
	  my $pct_in_diff_SPG    = sprintf("%.2f",(($de_in_diff_SPG/$total_key_genes)*100));
	  my $pct_in_scyte 	     = sprintf("%.2f",(($de_in_scyte/$total_key_genes)*100));
	  my $pct_in_early_stid1 = sprintf("%.2f",(($de_in_early_stid1/$total_key_genes)*100));
	  my $pct_in_early_stid2 = sprintf("%.2f",(($de_in_early_stid2/$total_key_genes)*100));
	  my $pct_in_late_stid1  = sprintf("%.2f",(($de_in_late_stid1/$total_key_genes)*100));
	  my $pct_in_late_stid2  = sprintf("%.2f",(($de_in_late_stid2/$total_key_genes)*100));

	  print "$category: $total_key_genes genes\n";
	  print "number of genes differentially expressed in undiff SPG = $de_in_undiff_SPG ($pct_in_undiff_SPG%)\n";
	  print "number of genes differentially expressed in diff SPG = $de_in_diff_SPG ($pct_in_diff_SPG%)\n";
	  print "number of genes differentially expressed in spermatocyte = $de_in_scyte ($pct_in_scyte%)\n";
	  print "number of genes differentially expressed in early spermatid 1 = $de_in_early_stid1 ($pct_in_early_stid1%)\n";
	  print "number of genes differentially expressed in early spermatid 2 = $de_in_early_stid2 ($pct_in_early_stid2%)\n";
	  print "number of genes differentially expressed in late spermatid 1 = $de_in_late_stid1 ($pct_in_late_stid1%)\n";
	  print "number of genes differentially expressed in late spermatid 2 = $de_in_late_stid2 ($pct_in_late_stid2%)\n";
	  print "number of genes differentially expressed in endothelia/macrophages = $de_in_endothelia ($pct_in_endothelia%)\n";
	  print "number of genes differentially expressed in myoid/Leydig cells = $de_in_myoid ($pct_in_myoid%)\n";
	  print "number of genes differentially expressed in Sertoli cells = $de_in_sertoli ($pct_in_sertoli%)\n";
	  print "~~\n";
	  
	  if ($x == 0)
		{ print OUT8 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT8 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT8 "spermatocyte\t$de_in_scyte\n";
		  print OUT8 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT8 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT8 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT8 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT8 "endothelia\t$de_in_endothelia\n";
		  print OUT8 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT8 "Sertoli\t$de_in_sertoli\n";
		}
	  elsif ($x == 1)
		{ print OUT9 "undiff SPG\t$de_in_undiff_SPG\n";
		  print OUT9 "diff SPG\t$de_in_diff_SPG\n";
		  print OUT9 "spermatocyte\t$de_in_scyte\n";
		  print OUT9 "early spermatid 1\t$de_in_early_stid1\n";
		  print OUT9 "early spermatid 2\t$de_in_early_stid2\n";
		  print OUT9 "late spermatid 1\t$de_in_late_stid1\n";
		  print OUT9 "late spermatid 2\t$de_in_late_stid2\n";
		  print OUT9 "endothelia\t$de_in_endothelia\n";
		  print OUT9 "myoid/Leydig\t$de_in_myoid\n";
		  print OUT9 "Sertoli\t$de_in_sertoli\n";
		}
	}

close(OUT8) or die $!; close(OUT9) or die $!;
exit 1;