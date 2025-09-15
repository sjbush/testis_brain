# PURPOSE: summarise germline expression, at the protein-level, of genes associated with human-specific single nucleotide changes.

use strict;
use warnings;

# REQUIREMENTS
my $in_file1 = 'HPA/rna_tissue_hpa.tsv'; 	  # from https://v23.proteinatlas.org/about/download (transcript-level expression per tissue: gives nTPM in the testis)
my $in_file2 = 'HPA/normal_tissue.tsv'; 	  # from https://v23.proteinatlas.org/about/download (protein-level expression per cell type)
my $in_file3 = 'Ens112.gene_annotations.txt'; # from Ensembl BioMart (GRCh38.p14): Gene stable ID, Gene name, Gene description, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand, Gene type, Phenotype description, Source name
my $in_file4 = 'genes_with_human_specific_changes.txt'; # from Table S1 ("HHMCs": human-lineage high-frequency missense changes) of Kulwhilm 2019: a catalog of single nucleotide changes distinguishing modern humans from archaic hominins
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }
if (!(-e($in_file3))) { print "ERROR: cannot find $in_file3\n"; exit 1; }
if (!(-e($in_file4))) { print "ERROR: cannot find $in_file3\n"; exit 1; }

# OUTPUT
my $out_file = 'germline_expression_at_protein_level_of_genes_with_human_specific_changes.txt';
open(OUT,'>',$out_file) or die $!;
print OUT "Gene name\tEnsembl gene ID\tDescription\t";
print OUT "RNA expression in testis (nTPM)\t";
print OUT "Was antibody staining performed for specific germ cell types? (yes/no)\t";
print OUT "Reliability classification for protein expression data\t";
print OUT "Protein expression in the germline? (yes/probably/no/no data)\t";
print OUT "Protein expression in somatic testicular cells? (yes/no/no data)\t";
print OUT "Protein expression in cells in seminiferous ducts?\tProtein expression in Leydig cells?\t";
print OUT "Protein expression in Sertoli cells?\tProtein expression in peritubular cells?\tProtein expression in spermatogonial cells?\tProtein expression in preleptotene spermatocytes?\tProtein expression in pachytene spermatocytes?\tProtein expression in early or rounded spermatids?\tProtein expression in late or elongated spermatids?\t";
print OUT "If protein is only expressed in one testicular cell type (somatic or germline), what is it? (given only when staining was performed for specific germ cell types)\t";
print OUT "If protein is only expressed in one germ cell type, what is it? (given only when staining was performed for specific germ cell types)\n";

# STORE TRANSCRIPT-LEVEL EXPRESSION PER TISSUE ("nTPM")
my %tpms = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_name = $line[1]; my $tissue = $line[2]; my $ntpm = $line[5];
	  $tpms{$gene_id}{$tissue} = $ntpm;
	}
close(IN) or die $!;

# STORE PROTEIN-LEVEL EXPRESSION PER CELL TYPE
my %levels = (); my %reliabilities = ();
open(IN,$in_file2) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_name = $line[1]; my $tissue = $line[2]; my $cell_type = $line[3]; my $level = lc($line[4]); my $reliability = lc($line[5]);
	  $cell_type =~ s/\'//g;
	  next if ($reliability eq 'uncertain'); # CHECKPOINT: retain only those proteins with status of Enhanced, Supported, or Approved
	  next if ($tissue eq 'N/A'); # CHECKPOINT: exclude entries which do not have an assigned tissue
	  if ($cell_type eq 'sertoli cells') { $cell_type = 'Sertoli cells'; }
	  $levels{$gene_id}{$tissue}{$cell_type} = $level;
	  $reliabilities{$gene_id} = $reliability;
	}
close(IN) or die $!;

# STORE GENE DESCRIPTIONS AND RESTRICT EACH GENE LIST ONLY TO KNOWN PROTEIN-CODING GENES ON THE MAJOR AUTOSOMES + XY (I.E. NO CHROMOSOMAL PATCHES)
my %annotations = ();
my @acceptable_chrs = (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT/);
my %acceptable_chrs = map {$_ => 1} @acceptable_chrs;
open(IN,$in_file3) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_name = '';
	  $gene_id =~ s/\s//g;
	  if (!(defined($line[1]))) { $gene_name = $gene_id; } else { $gene_name = $line[1]; }
	  if ($gene_name eq '') { $gene_name = $gene_id; }
	  my $gene_desc = $line[2]; my $chr = $line[3]; my $gene_type = $line[7];
	  next if (!(exists($acceptable_chrs{$chr})));
	  $annotations{$gene_id}{gene_name} = $gene_name;
	  $annotations{$gene_id}{gene_desc} = $gene_desc;
	  $annotations{$gene_id}{gene_type} = $gene_type;
	}
close(IN) or die $!;

# STORE A LIST OF GENES CONTAINING HUMAN-LINEAGE HIGH-FREQUENCY MISSENSE CHANGES
my %genes = ();
open(IN,$in_file4) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^gene\=(ENSG.*?)\;.+$/)
		{ my $gene_id = $1;
		  if (!(exists($annotations{$gene_id}{gene_name})))
			{ print "ERROR: unable to find $gene_id in the latest Ensembl annotation\n"; }
		  else
			{ $genes{$gene_id}++; }
		}
	}
close(IN) or die $!;

my $total_genes = scalar keys %genes;
print "there are $total_genes genes which are associated with human-specific single nucleotide changes\n";

my %out_lines = (); my %gene_ids_seen = (); my %gene_ids_seen_excl_uncertain = ();
my $no_with_no_reliable_data = 0;
my $no_stained_for_germ_cells = 0; my $no_stained_for_ducts = 0; my $no_stained_for_ducts_or_germ_cells = 0; my $no_expressed_in_soma = 0; my $no_prob_or_def_in_germline = 0; my $no_prob_in_germline = 0; my $no_expressed_in_germline = 0; my $no_expressed_only_in_germline = 0; my $no_expressed_in_all_germ_cells = 0;
my $no_expressed_in_spg = 0; my $no_expressed_in_scyte = 0; my $no_expressed_in_stid = 0;
while((my $gene_id,my $irrel)=each(%genes))
	{ my $gene_name = $annotations{$gene_id}{gene_name};
	  my $gene_type = $annotations{$gene_id}{gene_type};
	  my $gene_desc = $annotations{$gene_id}{gene_desc};
	  if (exists($gene_ids_seen{$gene_id})) { print "ERROR: $gene_id occurs more than once\n"; exit 1; }
	  my $testis_ntpm = 'not available';
	  if (exists($tpms{$gene_id}{'testis'})) { $testis_ntpm = $tpms{$gene_id}{'testis'}; }
	  my $reliability = 'no data';
	  if (exists($reliabilities{$gene_id})) { $reliability = $reliabilities{$gene_id}; }
	  
	  next if ($reliability eq 'no data'); # CHECKPOINT: exclude proteins with no data
	  
	  my $was_specific_germ_cell_staining_performed = 'no';
	  my $in_germline = 'no data'; my $in_somatic = 'no data';
	  my $exp_duct = ''; my $exp_Leydig = ''; my $exp_Sertoli = ''; my $exp_peritubular = ''; my $exp_spg = ''; my $exp_prelep = ''; my $exp_pachy = ''; my $exp_std_early = ''; my $exp_std_late = '';
	  my $unique_cell_type = ''; my $unique_germ_cell_type = '';
	  
	  if ( (exists($levels{$gene_id}{'testis'}{'cells in seminiferous ducts'})) or (exists($levels{$gene_id}{'testis'}{'spermatogonia cells'})) ) # either non-specific OR specific staining for germ cells
		{ $no_stained_for_ducts_or_germ_cells++ unless ($reliability eq 'uncertain');
		  if (exists($levels{$gene_id}{'testis'}{'cells in seminiferous ducts'}))
			{ my $exp_duct = $levels{$gene_id}{'testis'}{'cells in seminiferous ducts'};
			  if (($exp_duct eq 'low') or ($exp_duct eq 'medium') or ($exp_duct eq 'high'))
				{ $no_prob_or_def_in_germline++ unless ($reliability eq 'uncertain'); }
			}
		  elsif ( (exists($levels{$gene_id}{'testis'}{'spermatogonia cells'})) and (exists($levels{$gene_id}{'testis'}{'preleptotene spermatocytes'})) and (exists($levels{$gene_id}{'testis'}{'pachytene spermatocytes'})) and (exists($levels{$gene_id}{'testis'}{'round or early spermatids'})) and (exists($levels{$gene_id}{'testis'}{'elongated or late spermatids'})) )
			{ my $exp_spg       = $levels{$gene_id}{'testis'}{'spermatogonia cells'};
			  my $exp_prelep    = $levels{$gene_id}{'testis'}{'preleptotene spermatocytes'};
			  my $exp_pachy     = $levels{$gene_id}{'testis'}{'pachytene spermatocytes'};
			  my $exp_std_early = $levels{$gene_id}{'testis'}{'round or early spermatids'};
			  my $exp_std_late  = $levels{$gene_id}{'testis'}{'elongated or late spermatids'};
			  if ( (($exp_spg 	    eq 'low') or ($exp_spg 	     eq 'medium') or ($exp_spg 	     eq 'high')) or 
				   (($exp_prelep 	eq 'low') or ($exp_prelep 	 eq 'medium') or ($exp_prelep 	 eq 'high')) or 
				   (($exp_pachy     eq 'low') or ($exp_pachy     eq 'medium') or ($exp_pachy     eq 'high')) or 
				   (($exp_std_early eq 'low') or ($exp_std_early eq 'medium') or ($exp_std_early eq 'high')) or 
				   (($exp_std_late  eq 'low') or ($exp_std_late  eq 'medium') or ($exp_std_late  eq 'high')) )
				{ $no_prob_or_def_in_germline++ unless ($reliability eq 'uncertain'); }
			}
		}
	  
	  if ( (exists($levels{$gene_id}{'testis'}{'cells in seminiferous ducts'})) and (exists($levels{$gene_id}{'testis'}{'Leydig cells'})) ) # non-specific staining for germ cells
		{ $no_stained_for_ducts++ unless ($reliability eq 'uncertain');
		  $exp_duct   = $levels{$gene_id}{'testis'}{'cells in seminiferous ducts'};
		  $exp_Leydig = $levels{$gene_id}{'testis'}{'Leydig cells'};
		  
		  if (($exp_duct   eq 'low') or ($exp_duct   eq 'medium') or ($exp_duct   eq 'high'))
			{ $in_germline = 'probably'; $no_prob_in_germline++ unless ($reliability eq 'uncertain'); }
		  elsif ($exp_duct eq 'not detected')
			{ $in_germline = 'no'; }
		  
		  if (($exp_Leydig eq 'low') or ($exp_Leydig eq 'medium') or ($exp_Leydig eq 'high'))
			{ $in_somatic  = 'yes'; $no_expressed_in_soma++ unless ($reliability eq 'uncertain'); }
		  elsif ($exp_Leydig eq 'not detected')
			{ $in_somatic  = 'no'; }
		}
	  elsif ( (exists($levels{$gene_id}{'testis'}{'Leydig cells'})) and (exists($levels{$gene_id}{'testis'}{'Sertoli cells'})) and (exists($levels{$gene_id}{'testis'}{'peritubular cells'})) and (exists($levels{$gene_id}{'testis'}{'spermatogonia cells'})) and (exists($levels{$gene_id}{'testis'}{'preleptotene spermatocytes'})) and (exists($levels{$gene_id}{'testis'}{'pachytene spermatocytes'})) and (exists($levels{$gene_id}{'testis'}{'round or early spermatids'})) and (exists($levels{$gene_id}{'testis'}{'elongated or late spermatids'})) ) # specific staining for germ cells
		{ $no_stained_for_germ_cells++ unless ($reliability eq 'uncertain');
		  $exp_Leydig 	   = $levels{$gene_id}{'testis'}{'Leydig cells'};
		  $exp_Sertoli 	   = $levels{$gene_id}{'testis'}{'Sertoli cells'};
		  $exp_peritubular = $levels{$gene_id}{'testis'}{'peritubular cells'};
		  $exp_spg 		   = $levels{$gene_id}{'testis'}{'spermatogonia cells'};
		  $exp_prelep 	   = $levels{$gene_id}{'testis'}{'preleptotene spermatocytes'};
		  $exp_pachy 	   = $levels{$gene_id}{'testis'}{'pachytene spermatocytes'};
		  $exp_std_early   = $levels{$gene_id}{'testis'}{'round or early spermatids'};
		  $exp_std_late    = $levels{$gene_id}{'testis'}{'elongated or late spermatids'};
		  $was_specific_germ_cell_staining_performed = 'yes';
		  if ( (($exp_spg 	    eq 'low') or ($exp_spg 	     eq 'medium') or ($exp_spg 	     eq 'high')) or 
		       (($exp_prelep 	eq 'low') or ($exp_prelep 	 eq 'medium') or ($exp_prelep 	 eq 'high')) or 
		       (($exp_pachy     eq 'low') or ($exp_pachy     eq 'medium') or ($exp_pachy     eq 'high')) or 
			   (($exp_std_early eq 'low') or ($exp_std_early eq 'medium') or ($exp_std_early eq 'high')) or 
			   (($exp_std_late  eq 'low') or ($exp_std_late  eq 'medium') or ($exp_std_late  eq 'high')) )
			{ $in_germline = 'yes'; $no_expressed_in_germline++ unless ($reliability eq 'uncertain');
			  if ($exp_spg ne 'not detected')
				{ $no_expressed_in_spg++   unless ($reliability eq 'uncertain'); }
			  if (($exp_prelep ne 'not detected') or ($exp_pachy ne 'not detected'))
				{ $no_expressed_in_scyte++ unless ($reliability eq 'uncertain'); }
			  if (($exp_std_early ne 'not detected') or ($exp_std_late ne 'not detected'))
				{ $no_expressed_in_stid++  unless ($reliability eq 'uncertain'); }
			  if (($exp_spg ne 'not detected') and ($exp_prelep ne 'not detected') and ($exp_pachy ne 'not detected') and ($exp_std_early ne 'not detected') and ($exp_std_late ne 'not detected'))
				{ $no_expressed_in_all_germ_cells++ unless ($reliability eq 'uncertain'); }
			}
		  elsif ( (($exp_spg 	   eq 'not detected') and ($exp_spg 	  eq 'not detected') and ($exp_spg 	     eq 'not detected')) and 
				  (($exp_prelep    eq 'not detected') and ($exp_prelep 	  eq 'not detected') and ($exp_prelep 	 eq 'not detected')) and 
		          (($exp_pachy     eq 'not detected') and ($exp_pachy     eq 'not detected') and ($exp_pachy     eq 'not detected')) and
			      (($exp_std_early eq 'not detected') and ($exp_std_early eq 'not detected') and ($exp_std_early eq 'not detected')) and
			      (($exp_std_late  eq 'not detected') and ($exp_std_late  eq 'not detected') and ($exp_std_late  eq 'not detected')) )
			{ $in_germline = 'no'; }
		  if ( (($exp_Leydig 	  eq 'low') or ($exp_Leydig 	 eq 'medium') or ($exp_Leydig 	   eq 'high')) or 
		       (($exp_Sertoli 	  eq 'low') or ($exp_Sertoli 	 eq 'medium') or ($exp_Sertoli 	   eq 'high')) or 
		       (($exp_peritubular eq 'low') or ($exp_peritubular eq 'medium') or ($exp_peritubular eq 'high')) )
			{ $in_somatic = 'yes'; $no_expressed_in_soma++ unless ($reliability eq 'uncertain'); }
		  elsif ( (($exp_Leydig 	  eq 'not detected') and ($exp_Leydig 	   eq 'not detected') and ($exp_Leydig 	    eq 'not detected')) and 
				  (($exp_Sertoli 	  eq 'not detected') and ($exp_Sertoli 	   eq 'not detected') and ($exp_Sertoli 	eq 'not detected')) and 
				  (($exp_peritubular  eq 'not detected') and ($exp_peritubular eq 'not detected') and ($exp_peritubular eq 'not detected')) )
			{ $in_somatic = 'no'; }
		  if (($in_germline eq 'yes') and ($in_somatic eq 'no'))
			{ $no_expressed_only_in_germline++ unless ($reliability eq 'uncertain'); }
		  my %cell_types_gene_is_expressed_in = (); my %germ_cell_types_gene_is_expressed_in = ();
		  if ($reliability ne 'uncertain')
			{ if ($exp_Leydig 	   ne 'not detected') { $cell_types_gene_is_expressed_in{'Leydig cells'}++; 				}
			  if ($exp_Sertoli 	   ne 'not detected') { $cell_types_gene_is_expressed_in{'Sertoli cells'}++; 				}
			  if ($exp_peritubular ne 'not detected') { $cell_types_gene_is_expressed_in{'peritubular cells'}++; 			}
			  if ($exp_spg 		   ne 'not detected') { $cell_types_gene_is_expressed_in{'spermatogonia cells'}++; 			$germ_cell_types_gene_is_expressed_in{'spermatogonia cells'}++; 		 }
			  if ($exp_prelep 	   ne 'not detected') { $cell_types_gene_is_expressed_in{'preleptotene spermatocytes'}++; 	$germ_cell_types_gene_is_expressed_in{'preleptotene spermatocytes'}++; 	 }
			  if ($exp_pachy 	   ne 'not detected') { $cell_types_gene_is_expressed_in{'pachytene spermatocytes'}++; 		$germ_cell_types_gene_is_expressed_in{'pachytene spermatocytes'}++; 	 }
			  if ($exp_std_early   ne 'not detected') { $cell_types_gene_is_expressed_in{'round or early spermatids'}++; 	$germ_cell_types_gene_is_expressed_in{'round or early spermatids'}++;    }
			  if ($exp_std_late    ne 'not detected') { $cell_types_gene_is_expressed_in{'elongated or late spermatids'}++; $germ_cell_types_gene_is_expressed_in{'elongated or late spermatids'}++; }
			  my $no_cell_types_gene_is_expressed_in = scalar keys %cell_types_gene_is_expressed_in;
			  if ($no_cell_types_gene_is_expressed_in == 1)
				{ while((my $cell_type,my $irrel)=each(%cell_types_gene_is_expressed_in))
					{ $unique_cell_type = $cell_type; }
				}
			  my $no_germ_cell_types_gene_is_expressed_in = scalar keys %germ_cell_types_gene_is_expressed_in;
			  if ($no_germ_cell_types_gene_is_expressed_in == 1)
				{ while((my $cell_type,my $irrel)=each(%germ_cell_types_gene_is_expressed_in))
					{ $unique_germ_cell_type = $cell_type; }
				}
			}
		}
	  
	  if ($reliability ne 'uncertain')
		{ $gene_ids_seen_excl_uncertain{$gene_id}++; }
	  elsif ($reliability eq 'uncertain')
		{ $no_with_no_reliable_data++; } # this will be incremented only if there is uncertain data
	  
	  my $out_line = "$gene_name\t$gene_id\t$gene_desc\t$testis_ntpm\t$was_specific_germ_cell_staining_performed\t$reliability\t$in_germline\t$in_somatic\t$exp_duct\t$exp_Leydig\t$exp_Sertoli\t$exp_peritubular\t$exp_spg\t$exp_prelep\t$exp_pachy\t$exp_std_early\t$exp_std_late\t$unique_cell_type\t$unique_germ_cell_type";
	  $out_lines{$gene_name}{$gene_id} = $out_line;
	  $gene_ids_seen{$gene_id}++;
	}
	
my $no_of_genes_with_uncertain = scalar keys %gene_ids_seen;
my $no_of_genes_excl_uncertain = scalar keys %gene_ids_seen_excl_uncertain;

# PRINT OUTPUT
my @gene_names = ();
while((my $gene_name,my $irrel)=each(%out_lines))
	{ push(@gene_names,$gene_name); }
my @sorted_gene_names = sort {$a cmp $b} @gene_names;
foreach my $gene_name (@sorted_gene_names)
	{ my @gene_ids = ();
	  while((my $gene_id,my $irrel)=each(%{$out_lines{$gene_name}}))
		{ push(@gene_ids,$gene_id); }
	  my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
	  my $num_gene_ids = @sorted_gene_ids;
	  foreach my $gene_id (@sorted_gene_ids)
		{ if ($num_gene_ids > 1) { print "$gene_name is associated with ID $gene_id\n"; }
		  print OUT "$out_lines{$gene_name}{$gene_id}\n";
		}
	}
	
my $pc_with_no_reliable_data 	   = sprintf("%.2f",(($no_with_no_reliable_data/$no_of_genes_with_uncertain)*100));
my $pc_stained_for_germ_cells 	   = sprintf("%.2f",(($no_stained_for_germ_cells/$no_of_genes_excl_uncertain)*100));
my $pc_stained_for_ducts 	  	   = sprintf("%.2f",(($no_stained_for_ducts/$no_of_genes_excl_uncertain)*100));
my $pc_expressed_in_soma 	  	   = sprintf("%.2f",(($no_expressed_in_soma/$no_of_genes_excl_uncertain)*100));
my $pc_prob_or_def_in_germline	   = sprintf("%.2f",(($no_prob_or_def_in_germline/$no_of_genes_excl_uncertain)*100));
my $pc_prob_in_germline 	  	   = sprintf("%.2f",(($no_prob_in_germline/$no_stained_for_ducts)*100));
my $pc_expressed_in_germline  	   = sprintf("%.2f",(($no_expressed_in_germline/$no_stained_for_germ_cells)*100));
my $pc_expressed_only_in_germline  = sprintf("%.2f",(($no_expressed_only_in_germline/$no_stained_for_germ_cells)*100));
my $pc_expressed_in_all_germ_cells = sprintf("%.2f",(($no_expressed_in_all_germ_cells/$no_stained_for_germ_cells)*100));
my $pc_expressed_in_spg 		   = sprintf("%.2f",(($no_expressed_in_spg/$no_stained_for_germ_cells)*100));
my $pc_expressed_in_scyte 		   = sprintf("%.2f",(($no_expressed_in_scyte/$no_stained_for_germ_cells)*100));
my $pc_expressed_in_stid 		   = sprintf("%.2f",(($no_expressed_in_stid/$no_stained_for_germ_cells)*100));

print OUT "\n";
print OUT "no. of genes seen:\t$no_of_genes_with_uncertain\n";
print OUT "no. of genes with whose protein expression data is usable, i.e. not classified as 'uncertain':\t$no_with_no_reliable_data ($pc_with_no_reliable_data%)\n";
print OUT "no. of genes seen, excluding those with no usable data:\t$no_of_genes_excl_uncertain\n";
print OUT "no. of genes where antibody staining was performed for at least one of five germ cell types (because the gene has testis-enriched transcript expression):\t$no_stained_for_germ_cells ($pc_stained_for_germ_cells%)\n";
print OUT "no. of genes where antibody staining was only performed for the vaguer category of 'seminiferous ducts':\t$no_stained_for_ducts ($pc_stained_for_ducts%)\n";
print OUT "no. of genes with detectable somatic expression (applies regardless of whether we stain for ducts or germ cells as we look at Leydig cells in both):\t$no_expressed_in_soma ($pc_expressed_in_soma%)\n";
print OUT "no. of genes with probable OR detectable germline expression (counting all genes irrespective of whether staining was performed only for ducts or for specific germ cells):\t$no_prob_or_def_in_germline ($pc_prob_or_def_in_germline%)\n";
print OUT "no. of genes with probable germline expression (only applicable when we stain for ducts):\t$no_prob_in_germline ($pc_prob_in_germline%)\n";
print OUT "no. of genes with detectable germline expression (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_germline ($pc_expressed_in_germline%)\n";
print OUT "no. of genes with detectable germline expression that are only expressed in the germline, not the soma (only applicable when we stain explicitly for germ cells):\t$no_expressed_only_in_germline ($pc_expressed_only_in_germline%)\n";
print OUT "no. of genes with detectable germline expression that are expressed in all germ cell types (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_all_germ_cells ($pc_expressed_in_all_germ_cells%)\n";
print OUT "no. of genes with detectable germline expression that are expressed in spermatogonia (i.e. pre-meiosis) (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_spg ($pc_expressed_in_spg%)\n";
print OUT "no. of genes with detectable germline expression that are expressed in spermatocytes (i.e. meiosis) (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_scyte ($pc_expressed_in_scyte%)\n";
print OUT "no. of genes with detectable germline expression that are expressed in spermatids (i.e. post-meiosis) (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_stid ($pc_expressed_in_stid%)\n";

close(OUT) or die $!;
exit 1;