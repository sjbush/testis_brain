=head

AFTER USAGE, RUN R:

# the following code produces Supplementary Figure 2: a box plot of the expression level of the differentially expressed brain-associated genes per testicular cell cluster of our own single-cell expression atlas

library(tidyverse)
library(ggpubr)
library(grid)
theme_set(theme_bw())

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_boxplot_of_avg_germline_exp_levels_of_brain_associated_genes.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
my_comparisons <- list( c("undiff SPG", "diff SPG"), c("undiff SPG", "spermatocyte"), c("undiff SPG", "early spermatid 1"), c("undiff SPG", "early spermatid 2"), c("undiff SPG", "late spermatid 1"), c("undiff SPG", "late spermatid 2") ) # c("undiff SPG", "myoid/Leydig"), c("undiff SPG", "Sertoli"), c("undiff SPG", "endothelia")
figS2a<-ggplot(df, aes(x=as.factor(Cell.type),y=Avg..expression)) + geom_boxplot(notch=TRUE,notchwidth=0.1,varwidth=TRUE,outlier.size=1) + scale_y_log10() + xlab('Cell cluster') + ylab('Avg. expression across all cells in the cluster') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test', fontsize = 0.2) + ggtitle('A. Inclusive set of brain-associated genes (n=7193)')

df<-read.table('C:/Users/User/Desktop/testis_brain/data_for_making_boxplot_of_avg_germline_exp_levels_of_brain_associated_genes.conservative.txt',header=T,sep='\t')
df$Cell.type <- factor(df$Cell.type, levels = c("undiff SPG","diff SPG","spermatocyte","early spermatid 1","early spermatid 2","late spermatid 1","late spermatid 2","myoid/Leydig","Sertoli","endothelia"))
my_comparisons <- list( c("undiff SPG", "diff SPG"), c("undiff SPG", "spermatocyte"), c("undiff SPG", "early spermatid 1"), c("undiff SPG", "early spermatid 2"), c("undiff SPG", "late spermatid 1"), c("undiff SPG", "late spermatid 2") ) # c("undiff SPG", "myoid/Leydig"), c("undiff SPG", "Sertoli"), c("undiff SPG", "endothelia")
figS2b<-ggplot(df, aes(x=as.factor(Cell.type),y=Avg..expression)) + geom_boxplot(notch=TRUE,notchwidth=0.1,varwidth=TRUE,outlier.size=1) + scale_y_log10() + xlab('Cell cluster') + ylab('Avg. expression across all cells in the cluster') + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(legend.position = "none") + stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test', fontsize = 0.2) + ggtitle('B. Exclusive set of brain-associated genes (n=5212)')

png(file="C:/Users/User/Desktop/testis_brain/Documents/Supplementary Figure 2.png",width=16,height=10,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=1,ncol=2)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(figS2a, vp = define_region(row=1,col=1))
print(figS2b, vp = define_region(row=1,col=2))
dev.off()

=cut

use strict;
use warnings;

# PARAMETERS
my $conservative = 'no'; # 'yes'; # conservative criteria are to restrict analysis only to those genes associated with a defined brain phenotype: macro/megalencephaly, autism, schizophrenia, epilepsy, brain weight

# REQUIREMENTS
my $in_file1 = ''; # from 2.summarise_germline_protein_expression_of_brain_genes.pl
if ($conservative eq 'no')
	{ $in_file1 = 'C:/Users/User/Desktop/testis_brain/germline_expression_at_protein_level_of_brain_associated_genes.txt'; } # from 2.summarise_germline_protein_expression_of_brain_genes.pl
elsif ($conservative eq 'yes')
	{ $in_file1 = 'C:/Users/User/Desktop/testis_brain/germline_expression_at_protein_level_of_brain_associated_genes.conservative.txt'; } # from 2.summarise_germline_protein_expression_of_brain_genes.pl
my $in_file2 = 'C:/Users/User/Desktop/human_adult_SSCs/results/human_adult.atlas.txt';
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_file = '';
if ($conservative eq 'no')
	{ $out_file = "data_for_making_boxplot_of_avg_germline_exp_levels_of_brain_associated_genes.txt"; }
elsif ($conservative eq 'yes')
	{ $out_file = "data_for_making_boxplot_of_avg_germline_exp_levels_of_brain_associated_genes.conservative.txt"; }
open(OUT,'>',$out_file) or die $!;
print OUT "Gene ID\tGene name\tCell type\tAvg. expression\n";

# STORE THE SET OF BRAIN-ASSOCIATED GENES
my %key_genes = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if (!(defined($line[4])));
	  my $gene_id = $line[1];
	  $key_genes{$gene_id}++;
	}
close(IN) or die $!;
my $total_key_genes = scalar keys %key_genes;
print "there are $total_key_genes key genes in total\n";

# PARSE THE WHOLE TESTIS SINGLE CELL ATLAS TO DETERMINE THE EXPRESSION PROFILES OF BRAIN-ASSOCIATED GENES
my %tpms = (); my %gene_names = ();
open(IN,$in_file2) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0]; my $ens_id = $line[1];
	  $gene_names{$ens_id} = $gene_name;
	  next if (!(exists($key_genes{$ens_id}))); # CHECKPOINT: restrict analysis only to brain-developmental genes
	  my $exp_sertoli 	 = $line[38];
	  my $exp_leydig 	 = $line[39];
	  my $exp_endothelia = $line[40];
	  my $exp_undiff_spg = $line[41];
	  my $exp_diff_spg   = $line[42];
	  my $exp_scyte 	 = $line[43];
	  my $exp_early_std1 = $line[44];
	  my $exp_early_std2 = $line[45];
	  my $exp_late_std1  = $line[46];
	  my $exp_late_std2  = $line[47];
	  $tpms{$ens_id}{'undiff SPG'} 		  = $exp_undiff_spg;
	  $tpms{$ens_id}{'diff SPG'} 		  = $exp_diff_spg;
	  $tpms{$ens_id}{'spermatocyte'} 	  = $exp_undiff_spg;
	  $tpms{$ens_id}{'early spermatid 1'} = $exp_early_std1;
	  $tpms{$ens_id}{'early spermatid 2'} = $exp_early_std2;
	  $tpms{$ens_id}{'late spermatid 1'}  = $exp_late_std1;
	  $tpms{$ens_id}{'late spermatid 2'}  = $exp_late_std2;
	  $tpms{$ens_id}{'endothelia'} 		  = $exp_endothelia;
	  $tpms{$ens_id}{'myoid/Leydig'} 	  = $exp_leydig;
	  $tpms{$ens_id}{'Sertoli'} 		  = $exp_sertoli;
	}
close(IN) or die $!;

my @gene_ids = ();
while((my $gene_id,my $irrel)=each(%tpms))
	{ push(@gene_ids,$gene_id); }
my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
foreach my $gene_id (@sorted_gene_ids)
	{ my $gene_name = 
	  my @cell_types = ();
	  while((my $cell_type,my $irrel)=each(%{$tpms{$gene_id}}))
		{ push(@cell_types,$cell_type); }
	  my @sorted_cell_types = sort {$a cmp $b} @cell_types;
	  foreach my $cell_type (@sorted_cell_types)
		{ my $tpm = $tpms{$gene_id}{$cell_type};
		  print OUT "$gene_id\t$gene_names{$gene_id}\t$cell_type\t$tpm\n";
		}
	}

close(OUT) or die $!;
exit 1;