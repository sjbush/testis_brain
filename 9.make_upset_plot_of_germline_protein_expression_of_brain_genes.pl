use strict;
use warnings;

# PARAMETERS
my $conservative = 'yes'; # 'no'; # conservative criteria are to restrict analysis only to those genes associated with a defined brain phenotype: macro/megalencephaly, autism, schizophrenia, epilepsy, brain weight

# REQUIREMENTS
my $in_file = '';
if ($conservative eq 'no')
	{ $in_file = 'germline_expression_at_protein_level_of_brain_associated_genes.txt'; } # from 2.summarise_germline_protein_expression_of_brain_genes.pl
elsif ($conservative eq 'yes')
	{ $in_file = 'germline_expression_at_protein_level_of_brain_associated_genes.conservative.txt'; } # from 2.summarise_germline_protein_expression_of_brain_genes.pl
if (!(-e($in_file))) { print "ERROR: cannot find $in_file\n"; exit 1; }

# OUTPUT
my $out_file = '';
if ($conservative eq 'no')
	{ $out_file = 'dataframe_for_making_upset_plot_of_germline_expression_of_brain_associated_genes.txt'; } # NOTE: for use with the R package ComplexUpset, following the instructions https://krassowski.github.io/complex-upset/articles/Examples_R.html
elsif ($conservative eq 'yes')
	{ $out_file = 'dataframe_for_making_upset_plot_of_germline_expression_of_brain_associated_genes.conservative.txt'; }
open(OUT,'>',$out_file) or die $!;

my %associations_per_gene = ();
open(IN,$in_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if (!(defined($line[4])));
	  my $gene_name = $line[0]; my $gene_id = $line[1]; my $specific_staining_performed = $line[4]; my $reliability = $line[5];
	  next if (($reliability eq 'uncertain') or ($reliability eq 'no data')); # CHECKPOINT: restrict analysis only to genes where we can have confidence in the staining pattern
	  next unless ($specific_staining_performed eq 'yes'); # CHECKPOINT: restrict analysis only to genes where antibody staining for specific germ cells was performed
	  my $ducts = $line[8]; my $leydig = $line[9]; my $sertoli = $line[10]; my $myoid = $line[11]; my $spermatogonia = $line[12]; my $preleptotene = $line[13]; my $pachytene = $line[14]; my $round_spermatid = $line[15]; my $elongated_spermatid = $line[16];
#	  $associations_per_gene{$gene_id}{'cells in seminiferous ducts'}++	 unless (($ducts 			   eq 'not detected') or ($ducts 			   eq ''));
	  $associations_per_gene{$gene_id}{'Leydig cells'}++ 				 unless (($leydig 			   eq 'not detected') or ($leydig 			   eq ''));
	  $associations_per_gene{$gene_id}{'Sertoli cells'}++ 				 unless (($sertoli 			   eq 'not detected') or ($sertoli 			   eq ''));
	  $associations_per_gene{$gene_id}{'peritubular myoid cells'}++ 	 unless (($myoid 			   eq 'not detected') or ($myoid 			   eq ''));
	  $associations_per_gene{$gene_id}{'spermatogonia'}++ 	    	     unless (($spermatogonia 	   eq 'not detected') or ($spermatogonia 	   eq ''));
	  $associations_per_gene{$gene_id}{'preleptotene spermatocytes'}++   unless (($preleptotene 	   eq 'not detected') or ($preleptotene 	   eq ''));
	  $associations_per_gene{$gene_id}{'pachytene spermatocytes'}++ 	 unless (($pachytene           eq 'not detected') or ($pachytene           eq ''));
	  $associations_per_gene{$gene_id}{'round or early spermatids'}++    unless (($round_spermatid 	   eq 'not detected') or ($round_spermatid 	   eq ''));
	  $associations_per_gene{$gene_id}{'elongated or late spermatids'}++ unless (($elongated_spermatid eq 'not detected') or ($elongated_spermatid eq ''));
	}
close(IN) or die $!;

my $num_genes = scalar keys %associations_per_gene;
print "number of brain-associated genes processed for upset plot: $num_genes\n"; # note that this is NOT the same as the number of genes where $specific_staining_performed eq 'yes' and $reliability ne 'uncertain' (n = 1502) because not all of these genes will actually be detectable in at least one cell type. As it happens, neither LRP1B nor PSMB10 are found, so we have n = 1500 genes here.

# OUTPUT A DATA FRAME FOR MAKING AN UPSET PLOT USING THE R PACKAGE ComplexUpset, FOLLOWING THE VIGNETTE AT https://krassowski.github.io/complex-upset/articles/Examples_R.html
my @associations = ("Sertoli cells","Leydig cells","peritubular myoid cells","spermatogonia","preleptotene spermatocytes","pachytene spermatocytes","round or early spermatids","elongated or late spermatids");
my $associations_line = join("\t",@associations);
print OUT "Gene\t$associations_line\n";	
my @gene_ids = ();
while((my $gene_id,my $irrel)=each(%associations_per_gene))
	{ push(@gene_ids,$gene_id); }
my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
foreach my $gene_id (@sorted_gene_ids)
	{ my $out_line = '';
	  foreach my $association (@associations)
		{ my $true_or_false = 'FALSE'; # the default
		  if (exists($associations_per_gene{$gene_id}{$association}))
			{ $true_or_false = 'TRUE'; }
		  $out_line .= "$true_or_false\t";
		}
	  $out_line =~ s/\t$//;
	  print OUT "$gene_id\t$out_line\n";
	}
close(OUT) or die $!;

print "COMPLETE. Did we use a conservative gene list? $conservative\n";

exit 1;