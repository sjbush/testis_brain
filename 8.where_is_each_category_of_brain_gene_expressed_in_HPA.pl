=head

AFTER USAGE, RUN R:

# the following code produces Supplementary Figure 10: bar plots of the number of genes per brain phenotype that are expressed at the protein level in each tissue

library(tidyverse)
library(ggpubr)
library(grid)
theme_set(theme_bw())

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_macrocephaly_associated_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10a <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("A. Macrocephaly/megalencephaly (n = 399)")

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_autism_associated_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10b <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("B. Autism (n = 2562)")

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_schizo_associated_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10c <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("C. Schizophrenia (n = 345)")

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_epilepsy_associated_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10d <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("D. Epilepsy (n = 3372)")

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_iq_associated_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10e <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("E. IQ/educational attainment (n = 2303)")

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_brain_weight_associated_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10f <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("F. Brain weight (n = 879)")

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_har_associated_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10g <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("G. Human-accelerated regions (n = 1608)")

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_2_or_more_explicit_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10h <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("H. 2 or more associations (n = 1987)")

df<-read.table('C:/Users/User/Desktop/testis_brain/protein_exp_of_brain_genes_in_HPA/data_for_making_barplot_of_protein_expression_of_3_or_more_explicit_genes.txt',header=T,sep='\t')
df$Tissue <- factor(df$Tissue, levels=rev(sort(df$Tissue)))
fig10i <- ggplot(df, aes(x=Number, y=Tissue, fill=Tissue)) + geom_bar(stat="identity") + xlab('Number of proteins') + ylab('Tissue') + theme(legend.position = "none") + scale_fill_manual(values=c("caudate" = "blue", "cerebellum" = "blue", "cerebral cortex" = "blue", "choroid plexus" = "blue", "dorsal raphe" = "blue", "hippocampus" = "blue", "substantia nigra" = "blue", "testis" = "red")) + ggtitle("I. 3 or more associations (n = 271)")

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 10.png",width=16,height=20,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=3,ncol=3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(fig10a, vp = define_region(row=1,col=1))
print(fig10b, vp = define_region(row=1,col=2))
print(fig10c, vp = define_region(row=1,col=3))
print(fig10d, vp = define_region(row=2,col=1))
print(fig10e, vp = define_region(row=2,col=2))
print(fig10f, vp = define_region(row=2,col=3))
print(fig10g, vp = define_region(row=3,col=1))
print(fig10h, vp = define_region(row=3,col=2))
print(fig10i, vp = define_region(row=3,col=3))
dev.off()

=cut

use strict;
use warnings;

# REQUIREMENTS
my $in_file1 = 'HPA/normal_tissue.tsv'; # from https://v23.proteinatlas.org/about/download # protein-level expression per cell type
my $in_file2 = 'germline_expression_at_protein_level_of_brain_associated_genes.txt'; # from 8.summarise_germline_protein_expression_of_brain_genes.pl
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_dir = 'protein_exp_of_brain_genes_in_HPA';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $out_file1 = "$out_dir/data_for_making_barplot_of_protein_expression_of_macrocephaly_associated_genes.txt";
my $out_file2 = "$out_dir/data_for_making_barplot_of_protein_expression_of_autism_associated_genes.txt";
my $out_file3 = "$out_dir/data_for_making_barplot_of_protein_expression_of_schizo_associated_genes.txt";
my $out_file4 = "$out_dir/data_for_making_barplot_of_protein_expression_of_epilepsy_associated_genes.txt";
my $out_file5 = "$out_dir/data_for_making_barplot_of_protein_expression_of_IQ_associated_genes.txt";
my $out_file6 = "$out_dir/data_for_making_barplot_of_protein_expression_of_brain_weight_associated_genes.txt";
my $out_file7 = "$out_dir/data_for_making_barplot_of_protein_expression_of_HAR_associated_genes.txt";
my $out_file8 = "$out_dir/data_for_making_barplot_of_protein_expression_of_2_or_more_explicit_genes.txt";
my $out_file9 = "$out_dir/data_for_making_barplot_of_protein_expression_of_3_or_more_explicit_genes.txt";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!; open(OUT4,'>',$out_file4) or die $!; open(OUT5,'>',$out_file5) or die $!; open(OUT6,'>',$out_file6) or die $!; open(OUT7,'>',$out_file7) or die $!; open(OUT8,'>',$out_file8) or die $!; open(OUT9,'>',$out_file9) or die $!;
print OUT1 "Tissue\tNumber\n"; print OUT2 "Tissue\tNumber\n"; print OUT3 "Tissue\tNumber\n"; print OUT4 "Tissue\tNumber\n"; print OUT5 "Tissue\tNumber\n"; print OUT6 "Tissue\tNumber\n"; print OUT7 "Tissue\tNumber\n"; print OUT8 "Tissue\tNumber\n"; print OUT9 "Tissue\tNumber\n";

# STORE THE SET OF PROTEINS EXPRESSED PER TISSUE
my %protein_expression = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_name = $line[1]; my $tissue = $line[2]; my $cell_type = $line[3]; my $level = $line[4]; my $reliability = $line[5]; # note: these column numbers are only accurate if looking at v23 data
	  $cell_type =~ s/\'//g;
	  next if ($reliability eq 'Uncertain'); # CHECKPOINT: retain only those proteins with status of Enhanced, Supported, or Approved
	  if (($level eq 'Low') or ($level eq 'Medium') or ($level eq 'High'))
		{ $protein_expression{$tissue}{$gene_id}++; }
	}
close(IN) or die $!;

# STORE THE SET OF BRAIN-ASSOCIATED GENES FOR EACH PHENOTYPE
my %brain_genes = (); my %category_per_gene = ();
open(IN,$in_file2) or die $!;
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
	  print "$category: $total_key_genes genes\n";
	  my @tissues = ();
	  while((my $tissue,my $irrel)=each(%protein_expression))
		{ push(@tissues,$tissue); }
	  my @sorted_tissues = sort {$a cmp $b} @tissues;
	  foreach my $tissue (@sorted_tissues)
		{ my $num = 0;
		  while((my $gene_id,my $irrel)=each(%{$brain_genes{$category}}))
			{ if (exists($protein_expression{$tissue}{$gene_id}))
				{ $num++; }
			}
		  if 	($x == 1) { print OUT1 "$tissue\t$num\n"; }
		  elsif ($x == 2) { print OUT2 "$tissue\t$num\n"; }
		  elsif ($x == 3) { print OUT3 "$tissue\t$num\n"; }
		  elsif ($x == 4) { print OUT4 "$tissue\t$num\n"; }
		  elsif ($x == 5) { print OUT5 "$tissue\t$num\n"; }
		  elsif ($x == 6) { print OUT6 "$tissue\t$num\n"; }
		  elsif ($x == 7) { print OUT7 "$tissue\t$num\n"; }
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
	  my @tissues = ();
	  while((my $tissue,my $irrel)=each(%protein_expression))
		{ push(@tissues,$tissue); }
	  my @sorted_tissues = sort {$a cmp $b} @tissues;
	  foreach my $tissue (@sorted_tissues)
		{ my $num = 0;
		  while((my $gene_id,my $irrel)=each(%genes))
			{ if (exists($protein_expression{$tissue}{$gene_id}))
				{ $num++; }
			}
		  if 	($x == 0) { print OUT8 "$tissue\t$num\n"; }
		  elsif ($x == 1) { print OUT9 "$tissue\t$num\n"; }
		}
	}
close(OUT8) or die $!; close(OUT9) or die $!;
exit 1;