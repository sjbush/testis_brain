use strict;
use warnings;

# REQUIREMENTS
my $in_file1 = 'HPA/normal_tissue.tsv'; # from https://v23.proteinatlas.org/about/download
my $in_file2 = 'Ens112.gene_annotations.txt'; # from Ensembl BioMart (GRCh38.p14): Gene stable ID, Gene name, Gene description, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand, Gene type, Phenotype description, Source name
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_file1 = 'num_of_germline_proteins_also_found_in_brain.txt';
my $out_file2 = 'num_of_proteins_only_found_in_germline_and_brain.txt';
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;
print OUT1 "No. of genes where protein expression was assessed in the male germline, i.e. in one of five germ cell types\tNo. of genes whose protein was detected in a male germ cell\tOther tissue\tNo. of cell types for which protein expression was assessed\tCell types for which protein expression was assessed\tNo. of genes with a protein detected in both a germ cell and a cell of another tissue\t% of genes with a protein detected in both a germ cell and a cell of another tissue\n";
print OUT2 "Gene ID\tGene name\tBrain cells in which this protein is expressed\tMale germ cells in which this protein is expressed\n";

# STORE PROTEIN-LEVEL EXPRESSION DATA PER CELL TYPE
my %levels = (); my %gene_ids_seen = (); my %tissues_seen = (); my %cell_types_seen = (); my %cell_types_seen_per_tissue = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_name = $line[1]; my $tissue = $line[2]; my $cell_type = $line[3]; my $level = $line[4]; my $reliability = $line[5];
	  $cell_type =~ s/\'//g;
	  next if ($reliability eq 'Uncertain'); # CHECKPOINT: retain only those proteins with status of Enhanced, Supported, or Approved
	  next if ($tissue eq 'N/A'); # CHECKPOINT: exclude entries which do not have an assigned tissue
	  $gene_ids_seen{$gene_id}++;
	  $levels{$tissue}{$cell_type}{$gene_id} = $level;
	  $tissues_seen{$tissue}++; $cell_types_seen{$cell_type}++; $cell_types_seen_per_tissue{$tissue}{"$tissue/$cell_type"}++;
	}
close(IN) or die $!;

my $num_genes = scalar keys %gene_ids_seen;
my $num_tissues = scalar keys %tissues_seen;
my $num_cell_types = scalar keys %cell_types_seen;

print "$in_file1 contains data on the (reliable) protein products of $num_genes in $num_cell_types cell types from $num_tissues tissues (not every gene is assessed in each cell type)\n";

# STORE GENE ANNOTATION DATA SO WE CAN RESTRICT EACH GENE LIST ONLY TO KNOWN PROTEIN-CODING GENES
my %approved_gene_ids = ();
open(IN,$in_file2) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0];
	  next if (!(defined($line[1]))); # CHECKPOINT: discard genes which have not been annotated with an HGNC name
	  my $gene_name = $line[1]; my $gene_type = $line[7];
	  next if ($gene_type ne 'protein_coding'); # CHECKPOINT: discard genes which are not protein-coding
	  $approved_gene_ids{$gene_id} = $gene_name;
	}
close(IN) or die $!;

# IS THIS GENE FOUND IN THE BRAIN?
my $in_brain = 0;
my @brain_tissues = ("caudate","cerebellum","cerebral cortex","choroid plexus","dorsal raphe","hippocampus","substantia nigra");
my %brain_tissues = map {$_ => 1} @brain_tissues;
my %genes_in_brain = ();
while((my $brain_tissue,my $irrel)=each(%brain_tissues))
	{ while((my $cell_type,my $irrel)=each(%{$levels{$brain_tissue}}))
		{ while((my $gene_id,my $irrel)=each(%{$levels{$brain_tissue}{$cell_type}}))
			{ my $level = $levels{$brain_tissue}{$cell_type}{$gene_id};
			  if (($level eq 'Low') or ($level eq 'Medium') or ($level eq 'High'))
				{ $genes_in_brain{$gene_id}{level} = $level;
				  $genes_in_brain{$gene_id}{location}{"$brain_tissue/$cell_type"}++;
				}
			}
		}
	}

# IS THIS GENE FOUND IN THE GERMLINE AND *NOT* IN TESTICULAR SOMATIC CELLS? NOTE THAT ONLY A LIMITED SUBSET OF GENES HAVE PROTEIN EXPRESSION DIRECTLY QUANTIFIED FOR GERM CELL TYPES, AS OPPOSED TO THE BROADER CATEGORY OF "CELLS IN SEMINIFEROUS DUCTS".
my %genes_in_germline = (); my %genes_where_germline_directly_assessed = ();
my %genes_in_spg = (); my %genes_in_preleptotene = (); my %genes_in_pachytene = (); my %genes_in_early_spermatid = (); my %genes_in_late_spermatid = ();
my %genes_in_testis_somatic = ();
while((my $cell_type,my $irrel)=each(%{$levels{'testis'}}))
	{ if (($cell_type eq 'Leydig cells') or ($cell_type eq 'peritubular cells') or ($cell_type eq 'sertoli cells'))
		{ while((my $gene_id,my $irrel)=each(%{$levels{'testis'}{$cell_type}}))
			{ my $level = $levels{'testis'}{$cell_type}{$gene_id};
			  if (($level eq 'Low') or ($level eq 'Medium') or ($level eq 'High'))
				{ $genes_in_testis_somatic{$gene_id}++;
				}
			}
		}
	  if (($cell_type eq 'spermatogonia cells') or ($cell_type eq 'preleptotene spermatocytes') or ($cell_type eq 'pachytene spermatocytes') or ($cell_type eq 'round or early spermatids') or ($cell_type eq 'elongated or late spermatids'))
		{ while((my $gene_id,my $irrel)=each(%{$levels{'testis'}{$cell_type}}))
			{ $genes_where_germline_directly_assessed{$gene_id}++;
			  my $level = $levels{'testis'}{$cell_type}{$gene_id};
			  if (($level eq 'Low') or ($level eq 'Medium') or ($level eq 'High'))
				{ $genes_in_germline{$gene_id}{level} = $level;
				  $genes_in_germline{$gene_id}{location}{$cell_type}++;
				  if    ($cell_type eq 'spermatogonia cells')	       { $genes_in_spg{$gene_id} 			 = $level; }
				  elsif ($cell_type eq 'preleptotene spermatocytes')   { $genes_in_preleptotene{$gene_id} 	 = $level; }
				  elsif ($cell_type eq 'pachytene spermatocytes')	   { $genes_in_pachytene{$gene_id} 		 = $level; }
				  elsif ($cell_type eq 'round or early spermatids')	   { $genes_in_early_spermatid{$gene_id} = $level; }
				  elsif ($cell_type eq 'elongated or late spermatids') { $genes_in_late_spermatid{$gene_id}  = $level; }
				}
			}
		}
	}
	
# IS THIS GENE FOUND IN ANY OTHER TISSUE BESIDES BRAIN?
my %genes_in_other = (); my %genes_in_other_excl_separate_brain = ();
while((my $tissue,my $irrel)=each(%levels))
	{ next if ($tissue eq 'testis'); # CHECKPOINT: skip testis (we have parsed this already)
	  while((my $cell_type,my $irrel)=each(%{$levels{$tissue}}))
		{ while((my $gene_id,my $irrel)=each(%{$levels{$tissue}{$cell_type}}))
			{ my $level = $levels{$tissue}{$cell_type}{$gene_id};
			  if (($level eq 'Low') or ($level eq 'Medium') or ($level eq 'High'))
				{ $genes_in_other{$gene_id}{$tissue} = $level;
				  
				  # for the file OUT2, we are going to print those genes which are only found in a male germ cell and a brain cell, and no other (see code below)
				  # however, to do this properly, we must not include in the count of 'other' the SEPARATE brain tissues
				  if (!(exists($brain_tissues{$tissue})))
					{ $genes_in_other_excl_separate_brain{$gene_id}{$tissue} = $level;
					}
				}
			}
		}
	}

# SUMMARISE RESULTS: OF ALL THE PROTEINS WHOSE EXPRESSION WAS ASSESSED IN THE GERMLINE, HOW MANY DO YOU ALSO FIND IN THE BRAIN - AND ELSEWHERE?
# We should note that protein expression is only assessed in these germ cells in the first place if there is high expression in the testis (otherwise, the HPA resorts to the cruder analysis of testing for expression in "cells in seminiferous ducts").
my $num_genes_where_germline_directly_assessed = scalar keys %genes_where_germline_directly_assessed;
print "no. of proteins where germline expression directly assessed = $num_genes_where_germline_directly_assessed\n";
my $num_genes_in_germline = scalar keys %genes_in_germline;
print "no. of proteins expressed in germline = $num_genes_in_germline\n";

my $num_genes_in_brain_and_germline = 0; my $num_genes_in_brain_and_spg = 0; my $num_genes_in_brain_and_preleptotene = 0; my $num_genes_in_brain_and_pachytene = 0; my $num_genes_in_brain_and_early_spermatid = 0; my $num_genes_in_brain_and_late_spermatid = 0;
my %num_genes_in_germline_and_other = (); my %num_genes_only_in_germline_and_brain = ();
while((my $gene_id,my $irrel)=each(%genes_in_germline))
	{ if (exists($genes_in_brain{$gene_id}))
		{ $num_genes_in_brain_and_germline++;
		  if (exists($genes_in_spg{$gene_id})) 			   { $num_genes_in_brain_and_spg++; 			}
		  if (exists($genes_in_preleptotene{$gene_id}))    { $num_genes_in_brain_and_preleptotene++; 	}
		  if (exists($genes_in_pachytene{$gene_id})) 	   { $num_genes_in_brain_and_pachytene++; 		}
		  if (exists($genes_in_early_spermatid{$gene_id})) { $num_genes_in_brain_and_early_spermatid++; }
		  if (exists($genes_in_late_spermatid{$gene_id}))  { $num_genes_in_brain_and_late_spermatid++;  }
		}
	  if (exists($genes_in_other{$gene_id}))
		{ while((my $other_tissue,my $irrel)=each(%{$genes_in_other{$gene_id}}))
			{ $num_genes_in_germline_and_other{$other_tissue}++;
			}
		}
	  
	  # is this gene only expressed in a brain cell and a male germ cell and no other?
	  if ( (exists($genes_in_brain{$gene_id})) and (!(exists($genes_in_other_excl_separate_brain{$gene_id}))) and (!(exists($genes_in_testis_somatic{$gene_id}))) )
		{ my $gene_name = $approved_gene_ids{$gene_id};
		  my @locations = ();
		  while((my $cell_type,my $irrel)=each(%{$genes_in_brain{$gene_id}{location}}))
			{ push(@locations,$cell_type); }
		  my @sorted_locations = sort {$a cmp $b} @locations;
		  my $brain_locations = join(", ",@sorted_locations);
		  @locations = ();
		  while((my $cell_type,my $irrel)=each(%{$genes_in_germline{$gene_id}{location}}))
			{ push(@locations,$cell_type); }
		  @sorted_locations = sort {$a cmp $b} @locations;
		  my $testis_locations = join(", ",@sorted_locations);
		  print OUT2 "$gene_id\t$gene_name\t$brain_locations\t$testis_locations\n";
		  $num_genes_only_in_germline_and_brain{$gene_id}++;
		}
	}
	
# PRINT THE PROPORTION OF MALE GERMLINE-EXPRESSED PROTEINS FOUND IN ANOTHER TISSUE
my @other_tissue = ();
while((my $other_tissue,my $irrel)=each(%num_genes_in_germline_and_other))
	{ push(@other_tissue,$other_tissue); }
my @sorted_other_tissue = sort {$a cmp $b} @other_tissue;
foreach my $other_tissue (@sorted_other_tissue)
	{ my $num_cell_types_seen_per_tissue = scalar keys %{$cell_types_seen_per_tissue{$other_tissue}};
	  my @cell_types_seen_per_tissue = ();
	  while((my $cell_type,my $irrel)=each(%{$cell_types_seen_per_tissue{$other_tissue}}))
		{ push(@cell_types_seen_per_tissue,$cell_type); }
	  my @sorted_cell_types_seen_per_tissue = sort {$a cmp $b} @cell_types_seen_per_tissue;
	  my $cell_types_seen_per_tissue = join(", ",@sorted_cell_types_seen_per_tissue);
	  
	  my $num_genes_in_germline_and_other = $num_genes_in_germline_and_other{$other_tissue};
	  my $pct_genes_in_germline_also_in_other = sprintf("%.2f",(($num_genes_in_germline_and_other/$num_genes_in_germline)*100));
	  print OUT1 "$num_genes_where_germline_directly_assessed\t$num_genes_in_germline\t$other_tissue\t$num_cell_types_seen_per_tissue\t$cell_types_seen_per_tissue\t$num_genes_in_germline_and_other\t$pct_genes_in_germline_also_in_other\n";
	}

# PRINT SUMMARY STATISTICS OF THE NUMBER OF MALE GERM CELL-EXPRESSED PROTEINS FOUND IN A BRAIN CELL (BROKEN DOWN BY GERM CELL TYPE)
my $pct_genes_in_germline_also_in_brain = sprintf("%.2f",(($num_genes_in_brain_and_germline/$num_genes_in_germline)*100));
my $num_genes_only_in_germline_and_brain = scalar keys %num_genes_only_in_germline_and_brain;
my $num_cell_types_in_brain = 0;
my @cell_types_seen_in_brain = ();
while((my $brain_tissue,my $irrel)=each(%brain_tissues))
	{ while((my $cell_type,my $irrel)=each(%{$cell_types_seen_per_tissue{$brain_tissue}}))
		{ push(@cell_types_seen_in_brain,$cell_type);
		  $num_cell_types_in_brain++;
		}
	}
my @sorted_cell_types_seen_in_brain = sort {$a cmp $b} @cell_types_seen_in_brain;
my $cell_types_seen_in_brain = join(", ",@sorted_cell_types_seen_in_brain);
print OUT1 "$num_genes_where_germline_directly_assessed\t$num_genes_in_germline\tbrain\t$num_cell_types_in_brain\t$cell_types_seen_in_brain\t$num_genes_in_brain_and_germline\t$pct_genes_in_germline_also_in_brain\n";
print OUT1 "\n";
print OUT1 "Summary:\n";
print OUT1 "$num_genes_in_brain_and_germline proteins ($pct_genes_in_germline_also_in_brain%) which are expressed in at least one of five male germ cell types are also found in the brain (that is, in at least one cell type in at least one of seven brain tissue types: caudate, cerebellum, cerebral cortex, choroid plexus, dorsal raphe, hippocampus, substantia nigra)\n";
my $pct_genes_in_germline_also_in_spg = sprintf("%.2f",(($num_genes_in_brain_and_spg/$num_genes_in_germline)*100));
print OUT1 "$num_genes_in_brain_and_spg proteins ($pct_genes_in_germline_also_in_spg%) are found in both a brain cell and spermatogonia\n";
my $pct_genes_in_germline_also_in_preleptotene = sprintf("%.2f",(($num_genes_in_brain_and_preleptotene/$num_genes_in_germline)*100));
print OUT1 "$num_genes_in_brain_and_preleptotene proteins ($pct_genes_in_germline_also_in_preleptotene%) are found in both a brain cell and preleptotene spermatocytes\n";
my $pct_genes_in_germline_also_in_pachytene = sprintf("%.2f",(($num_genes_in_brain_and_pachytene/$num_genes_in_germline)*100));
print OUT1 "$num_genes_in_brain_and_pachytene proteins ($pct_genes_in_germline_also_in_pachytene%) are found in both a brain cell and pachytene spermatocytes\n";
my $pct_genes_in_germline_also_in_early_spermatid = sprintf("%.2f",(($num_genes_in_brain_and_early_spermatid/$num_genes_in_germline)*100));
print OUT1 "$num_genes_in_brain_and_early_spermatid proteins ($pct_genes_in_germline_also_in_early_spermatid%) are found in both a brain cell and early/elongated spermatids\n";
my $pct_genes_in_germline_also_in_late_spermatid = sprintf("%.2f",(($num_genes_in_brain_and_late_spermatid/$num_genes_in_germline)*100));
print OUT1 "$num_genes_in_brain_and_late_spermatid proteins ($pct_genes_in_germline_also_in_late_spermatid%) are found in both a brain cell and late/round spermatids\n";
print OUT1 "\n";
print OUT1 "$num_genes_only_in_germline_and_brain proteins are ONLY found BOTH in a male germ cell and brain cell, and no other cell in any other tissue\n";

close(OUT1) or die $!; close(OUT2) or die $!;
exit 1;