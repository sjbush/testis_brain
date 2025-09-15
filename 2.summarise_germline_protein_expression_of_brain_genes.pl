# PURPOSE: summarise germline expression, at the protein-level, of genes associated with brain growth, development, and evolution.

use strict;
use warnings;

# PARAMETERS
my $conservative = 'no'; # 'yes'; # conservative criteria are to restrict analysis only to those genes associated with a defined brain phenotype: macro/megalencephaly, autism, schizophrenia, epilepsy, brain weight

# REQUIREMENTS
my $in_file1 = 'HPA/rna_tissue_hpa.tsv'; 	  # from https://v23.proteinatlas.org/about/download # transcript-level expression per tissue: gives nTPM in the testis
my $in_file2 = 'HPA/normal_tissue.tsv'; 	  # from https://v23.proteinatlas.org/about/download # protein-level expression per cell type
my $in_file3 = 'Ens112.gene_annotations.txt'; # from Ensembl BioMart (GRCh38.p14): Gene stable ID, Gene name, Gene description, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand, Gene type, Phenotype description, Source name
my $in_file4 = '';
if ($conservative eq 'no')
	{ $in_file4 = 'brain_associated_genes.txt'; } # from 1.create_list_of_brain_associated_genes.pl
elsif ($conservative eq 'yes')
	{ $in_file4 = 'brain_associated_genes.conservative.txt'; } # from 1.create_list_of_brain_associated_genes.pl
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }
if (!(-e($in_file3))) { print "ERROR: cannot find $in_file3\n"; exit 1; }
if (!(-e($in_file4))) { print "ERROR: cannot find $in_file4\n"; exit 1; }

# OUTPUT
my $out_file = '';
if ($conservative eq 'no')
	{ $out_file = 'germline_expression_at_protein_level_of_brain_associated_genes.txt'; }
elsif ($conservative eq 'yes')
	{ $out_file = 'germline_expression_at_protein_level_of_brain_associated_genes.conservative.txt'; }
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
print OUT "If protein is only expressed in one germ cell type, what is it? (given only when staining was performed for specific germ cell types)\t";
print OUT "Reference(s)\t";
print OUT "Associated with macrocephaly or megalencephaly (in Bastos 2022 or DeCasien 2022)?\t";
print OUT "Associated with autism (in Qiu 2022, Rylaarsdam 2019, Sanders 2015, Satterstrom 2020, or SFARI 2025)?\t";
print OUT "Associated with schizophrenia (in Ripke 2014, Singh 2022 or Owen 2023)?\t";
print OUT "Associated with epilepsy (in Perucca 2020, Thakran 2020 or Rastin 2023)?\t";
print OUT "Associated with intelligence or educational attainment (in Lee 2018 or Savage 2018)?\t"; # if $conservative eq 'yes', this column will always be empty
print OUT "Associated with brain weight (in Boddy 2017 or Seidlitz 2023)?\t";
print OUT "Associated with a human-accelerated region (in Wei 2019)?\t"; # if $conservative eq 'yes', this column will always be empty
print OUT "No. of SPECIFIC associations\n";

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
#	  next if ($reliability eq 'uncertain'); # CHECKPOINT: retain only those proteins with status of Enhanced, Supported, or Approved
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

# ITERATE THROUGH A MANUALLY-CURATED LIST OF BRAIN-ASSOCIATED GENES TO DETERMINE THE LEVEL OF PROTEIN EXPRESSION IN THE GERMLINE
my %out_lines = (); my %gene_ids_seen = (); my %gene_ids_seen_excl_uncertain = ();
my %genes_with_2plus_brain = (); my %genes_with_3plus_brain = (); my %genes_with_4plus_brain = ();
my %genes_with_2plus_brain_specific_staining = (); my %genes_with_3plus_brain_specific_staining = (); my %genes_with_4plus_brain_specific_staining = ();
my $no_with_no_data_or_no_certain_data = 0; my $no_stained_for_germ_cells = 0; my $no_stained_for_ducts = 0; my $no_stained_for_ducts_or_germ_cells = 0; my $no_expressed_in_soma = 0; my $no_prob_or_def_in_germline = 0; my $no_prob_in_germline = 0; my $no_expressed_in_germline = 0; my $no_expressed_only_in_germline = 0; my $no_expressed_in_all_germ_cells = 0;
my $no_expressed_in_spg 		   = 0; my $no_expressed_in_scyte 		   = 0; my $no_expressed_in_stid 		   = 0;
my $no_only_expressed_in_spg 	   = 0; my $no_only_expressed_in_scyte 	   = 0; my $no_only_expressed_in_stid 	   = 0;
my $no_brain_genes_2_in_germline   = 0; my $no_brain_genes_3_in_germline   = 0; my $no_brain_genes_4_in_germline   = 0;
my $no_brain_genes_2_in_spg 	   = 0; my $no_brain_genes_3_in_spg 	   = 0; my $no_brain_genes_4_in_spg 	   = 0;
my $no_brain_genes_2_in_scyte 	   = 0; my $no_brain_genes_3_in_scyte 	   = 0; my $no_brain_genes_4_in_scyte 	   = 0;
my $no_brain_genes_2_in_stid 	   = 0; my $no_brain_genes_3_in_stid 	   = 0; my $no_brain_genes_4_in_stid 	   = 0;
my $no_brain_genes_2_only_in_spg   = 0; my $no_brain_genes_3_only_in_spg   = 0; my $no_brain_genes_4_only_in_spg   = 0;
my $no_brain_genes_2_only_in_scyte = 0; my $no_brain_genes_3_only_in_scyte = 0; my $no_brain_genes_4_only_in_scyte = 0;
my $no_brain_genes_2_only_in_stid  = 0; my $no_brain_genes_3_only_in_stid  = 0; my $no_brain_genes_4_only_in_stid  = 0;
open(IN,$in_file4) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  next if ($line eq '');
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; $gene_id =~ s/\s//g; my $refs = $line[1];
	  my $associated_macro   = $line[2];
	  my $associated_autism  = $line[3];
	  my $associated_schizo  = $line[4];
	  my $associated_epilep  = $line[5];
	  my $associated_iq_ea   = $line[6];
	  my $associated_weight  = $line[7];
	  my $associated_har  	 = $line[8];
	  my $no_of_associations = $line[9];
	  if (!(exists($annotations{$gene_id}))) { print "ERROR: unable to verify $gene_id: it is not recorded in the latest GRCh38 annotation\n"; }
	  next if (!(exists($annotations{$gene_id}))); # CHECKPOINT: skip genes which are not present in the latest GRCh38 annotation
	  my $gene_type = $annotations{$gene_id}{gene_type};
	  if ($gene_type ne 'protein_coding') { print "NOTE: discarding $gene_id because it is not protein-coding; it's $gene_type\n"; }
	  next if ($gene_type ne 'protein_coding'); # CHECKPOINT: discard genes which are not protein-coding
	  if (exists($gene_ids_seen{$gene_id})) { print "ERROR: $gene_id occurs more than once\n"; exit 1; }
	  my $gene_name = $annotations{$gene_id}{gene_name};
	  my $gene_desc = $annotations{$gene_id}{gene_desc};
	  my $testis_ntpm = 'not available';
	  if (exists($tpms{$gene_id}{'testis'})) { $testis_ntpm = $tpms{$gene_id}{'testis'}; }
	  my $reliability = 'no data';
	  if (exists($reliabilities{$gene_id})) { $reliability = $reliabilities{$gene_id}; }
	  
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
			{ $in_germline  = 'yes'; $no_expressed_in_germline++ unless ($reliability eq 'uncertain');
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
			{ $in_germline  = 'no'; }
		  if ( (($exp_Leydig 	  eq 'low') or ($exp_Leydig 	 eq 'medium') or ($exp_Leydig 	   eq 'high')) or 
		       (($exp_Sertoli 	  eq 'low') or ($exp_Sertoli 	 eq 'medium') or ($exp_Sertoli 	   eq 'high')) or 
		       (($exp_peritubular eq 'low') or ($exp_peritubular eq 'medium') or ($exp_peritubular eq 'high')) )
			{ $in_somatic  = 'yes'; $no_expressed_in_soma++ unless ($reliability eq 'uncertain'); }
		  elsif ( (($exp_Leydig 	  eq 'not detected') and ($exp_Leydig 	   eq 'not detected') and ($exp_Leydig 	    eq 'not detected')) and 
				  (($exp_Sertoli 	  eq 'not detected') and ($exp_Sertoli 	   eq 'not detected') and ($exp_Sertoli 	eq 'not detected')) and 
				  (($exp_peritubular  eq 'not detected') and ($exp_peritubular eq 'not detected') and ($exp_peritubular eq 'not detected')) )
			{ $in_somatic  = 'no'; }
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
#	  else
#		{ $no_with_no_data_or_no_certain_data++; } # this will be incremented if there is literally no data (but NOT if there is data, but it's uncertain)

	  if (($reliability eq 'uncertain') or ($reliability eq 'no data'))
		{ $no_with_no_data_or_no_certain_data++; } # this will be incremented only if there is uncertain data
	  
	  if (($reliability ne 'uncertain') and ($reliability ne 'no data')) # DERIVE SUMMARY STATISTICS ONLY FOR PROTEINS WITH RELIABLE EXPRESSION DATA
		{ $gene_ids_seen_excl_uncertain{$gene_id}++;
		  
		  # how many genes have proteins only expressed either before, during or after meiosis?
		  if ($unique_cell_type eq 'spermatogonia cells')
			{ $no_only_expressed_in_spg++; }
		  elsif (($unique_cell_type eq 'preleptotene spermatocytes') or ($unique_cell_type eq 'pachytene spermatocytes'))
			{ $no_only_expressed_in_scyte++; }
		  elsif (($unique_cell_type eq 'round or early spermatids') or ($unique_cell_type eq 'elongated or late spermatids'))
			{ $no_only_expressed_in_stid++; }
		  
		  # for genes with multiple specific associations with the brain, count how many have proteins detected in the germline, as well as spermatogonia/spermatocytes/spermatids in particular
		  if (($no_of_associations >= 2) and ($in_germline ne 'no data'))
			{ $genes_with_2plus_brain{$gene_id}++;
			  if (($in_germline eq 'yes') or ($in_germline eq 'probably'))
				{ $no_brain_genes_2_in_germline++; }
			  if ($was_specific_germ_cell_staining_performed eq 'yes')
				{ $genes_with_2plus_brain_specific_staining{$gene_id}++;
				  
				  # how many genes are found in SPG/SCYTE/STID?
				  if ($exp_spg ne 'not detected')
					{ $no_brain_genes_2_in_spg++; }
				  if (($exp_prelep ne 'not detected') or ($exp_pachy ne 'not detected'))
					{ $no_brain_genes_2_in_scyte++; }
				  if (($exp_std_early ne 'not detected') or ($exp_std_late ne 'not detected'))
					{ $no_brain_genes_2_in_stid++; }
				  
				  # how many genes are ONLY found in SPG/SCYTE/STID?
				  if (($exp_spg ne 'not detected') and ($unique_cell_type eq 'spermatogonia cells'))
					{ $no_brain_genes_2_only_in_spg++; }
				  if ( (($exp_prelep ne 'not detected') or ($exp_pachy ne 'not detected')) and (($unique_cell_type eq 'preleptotene spermatocytes') or ($unique_cell_type eq 'pachytene spermatocytes')) )
					{ $no_brain_genes_2_only_in_scyte++; }
				  if ( (($exp_std_early ne 'not detected') or ($exp_std_late ne 'not detected')) and (($unique_cell_type eq 'round or early spermatids') or ($unique_cell_type eq 'elongated or late spermatids')) )
					{ $no_brain_genes_2_only_in_stid++; }
				}
			}
		  if (($no_of_associations >= 3) and ($in_germline ne 'no data'))
			{ $genes_with_3plus_brain{$gene_id}++;
			  if (($in_germline eq 'yes') or ($in_germline eq 'probably'))
				{ $no_brain_genes_3_in_germline++; }
			  if ($was_specific_germ_cell_staining_performed eq 'yes')
				{ $genes_with_3plus_brain_specific_staining{$gene_id}++;
				  
				  # how many genes are found in SPG/SCYTE/STID?
				  if ($exp_spg ne 'not detected')
					{ $no_brain_genes_3_in_spg++; }
				  if (($exp_prelep ne 'not detected') or ($exp_pachy ne 'not detected'))
					{ $no_brain_genes_3_in_scyte++; }
				  if (($exp_std_early ne 'not detected') or ($exp_std_late ne 'not detected'))
					{ $no_brain_genes_3_in_stid++; }
				  
				  # how many genes are ONLY found in SPG/SCYTE/STID?
				  if (($exp_spg ne 'not detected') and ($unique_cell_type eq 'spermatogonia cells'))
					{ $no_brain_genes_3_only_in_spg++; }
				  if ( (($exp_prelep ne 'not detected') or ($exp_pachy ne 'not detected')) and (($unique_cell_type eq 'preleptotene spermatocytes') or ($unique_cell_type eq 'pachytene spermatocytes')) )
					{ $no_brain_genes_3_only_in_scyte++; }
				  if ( (($exp_std_early ne 'not detected') or ($exp_std_late ne 'not detected')) and (($unique_cell_type eq 'round or early spermatids') or ($unique_cell_type eq 'elongated or late spermatids')) )
					{ $no_brain_genes_3_only_in_stid++; }
				}
			}
		  if (($no_of_associations >= 4) and ($in_germline ne 'no data'))
			{ $genes_with_4plus_brain{$gene_id}++;
			  if (($in_germline eq 'yes') or ($in_germline eq 'probably'))
				{ $no_brain_genes_4_in_germline++; }
			  if ($was_specific_germ_cell_staining_performed eq 'yes')
				{ $genes_with_4plus_brain_specific_staining{$gene_id}++;
				  
				  # how many genes are found in SPG/SCYTE/STID?
				  if ($exp_spg ne 'not detected')
					{ $no_brain_genes_4_in_spg++; }
				  if (($exp_prelep ne 'not detected') or ($exp_pachy ne 'not detected'))
					{ $no_brain_genes_4_in_scyte++; }
				  if (($exp_std_early ne 'not detected') or ($exp_std_late ne 'not detected'))
					{ $no_brain_genes_4_in_stid++; }
				  
				  # how many genes are ONLY found in SPG/SCYTE/STID?
				  if (($exp_spg ne 'not detected') and ($unique_cell_type eq 'spermatogonia cells'))
					{ $no_brain_genes_4_only_in_spg++; }
				  if ( (($exp_prelep ne 'not detected') or ($exp_pachy ne 'not detected')) and (($unique_cell_type eq 'preleptotene spermatocytes') or ($unique_cell_type eq 'pachytene spermatocytes')) )
					{ $no_brain_genes_4_only_in_scyte++; }
				  if ( (($exp_std_early ne 'not detected') or ($exp_std_late ne 'not detected')) and (($unique_cell_type eq 'round or early spermatids') or ($unique_cell_type eq 'elongated or late spermatids')) )
					{ $no_brain_genes_4_only_in_stid++; }
				}
			}
		}

	  my $out_line = "$gene_name\t$gene_id\t$gene_desc\t$testis_ntpm\t$was_specific_germ_cell_staining_performed\t$reliability\t$in_germline\t$in_somatic\t$exp_duct\t$exp_Leydig\t$exp_Sertoli\t$exp_peritubular\t$exp_spg\t$exp_prelep\t$exp_pachy\t$exp_std_early\t$exp_std_late\t$unique_cell_type\t$unique_germ_cell_type\t$refs\t$associated_macro\t$associated_autism\t$associated_schizo\t$associated_epilep\t$associated_iq_ea\t$associated_weight\t$associated_har\t$no_of_associations";
	  $out_lines{$gene_name}{$gene_id} = $out_line;
	  $gene_ids_seen{$gene_id}++;
	}
close(IN) or die $!;

my $no_of_genes_with_uncertain = scalar keys %gene_ids_seen;
my $no_of_genes_excl_uncertain = scalar keys %gene_ids_seen_excl_uncertain;

# SANITY TEST: CONFIRM THAT WE HAVE PRINTED OUTPUT FOR EVERY GENE THAT WE HAVE SEEN
my $num_printed = 0;
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
		  $num_printed++;
		}
	}
if ($no_of_genes_with_uncertain != $num_printed) { print "ERROR: seen $no_of_genes_with_uncertain Ensembl gene IDs but printed $num_printed\n"; exit 1; }

my $pc_with_no_data 		  	   = sprintf("%.2f",(($no_with_no_data_or_no_certain_data/$no_of_genes_with_uncertain)*100));
my $no_of_genes_seen_excl_no_data  = $no_of_genes_excl_uncertain; # $no_of_genes_with_uncertain-$no_with_no_data_or_no_certain_data;
my $pc_stained_for_germ_cells 	   = sprintf("%.2f",(($no_stained_for_germ_cells/$no_of_genes_seen_excl_no_data)*100));
my $pc_stained_for_ducts 	  	   = sprintf("%.2f",(($no_stained_for_ducts/$no_of_genes_seen_excl_no_data)*100));
my $pc_expressed_in_soma 	  	   = sprintf("%.2f",(($no_expressed_in_soma/$no_of_genes_seen_excl_no_data)*100));
my $pc_prob_or_def_in_germline	   = sprintf("%.2f",(($no_prob_or_def_in_germline/$no_of_genes_seen_excl_no_data)*100));
my $pc_prob_in_germline 	  	   = sprintf("%.2f",(($no_prob_in_germline/$no_stained_for_ducts)*100));
my $pc_expressed_in_germline  	   = sprintf("%.2f",(($no_expressed_in_germline/$no_stained_for_germ_cells)*100));
my $pc_expressed_only_in_germline  = sprintf("%.2f",(($no_expressed_only_in_germline/$no_stained_for_germ_cells)*100));
my $pc_expressed_in_all_germ_cells = sprintf("%.2f",(($no_expressed_in_all_germ_cells/$no_stained_for_germ_cells)*100));
my $pc_expressed_in_spg 		   = sprintf("%.2f",(($no_expressed_in_spg/$no_stained_for_germ_cells)*100));
my $pc_expressed_in_scyte 		   = sprintf("%.2f",(($no_expressed_in_scyte/$no_stained_for_germ_cells)*100));
my $pc_expressed_in_stid 		   = sprintf("%.2f",(($no_expressed_in_stid/$no_stained_for_germ_cells)*100));
my $pc_only_expressed_in_spg 	   = sprintf("%.2f",(($no_only_expressed_in_spg/$no_stained_for_germ_cells)*100));
my $pc_only_expressed_in_scyte 	   = sprintf("%.2f",(($no_only_expressed_in_scyte/$no_stained_for_germ_cells)*100));
my $pc_only_expressed_in_stid 	   = sprintf("%.2f",(($no_only_expressed_in_stid/$no_stained_for_germ_cells)*100));
my $no_of_genes_with_2plus_brain   = scalar keys %genes_with_2plus_brain;
my $no_of_genes_with_3plus_brain   = scalar keys %genes_with_3plus_brain;
my $no_of_genes_with_4plus_brain   = scalar keys %genes_with_4plus_brain;
my $pc_brain_genes_2_in_germline   = 0; if ($no_of_genes_with_2plus_brain > 0) { $pc_brain_genes_2_in_germline = sprintf("%.2f",(($no_brain_genes_2_in_germline/$no_of_genes_with_2plus_brain)*100)); }
my $pc_brain_genes_3_in_germline   = 0; if ($no_of_genes_with_3plus_brain > 0) { $pc_brain_genes_3_in_germline = sprintf("%.2f",(($no_brain_genes_3_in_germline/$no_of_genes_with_3plus_brain)*100)); }
my $pc_brain_genes_4_in_germline   = 0; if ($no_of_genes_with_4plus_brain > 0) { $pc_brain_genes_4_in_germline = sprintf("%.2f",(($no_brain_genes_4_in_germline/$no_of_genes_with_4plus_brain)*100)); }
my $no_of_genes_with_2plus_brain_specific_staining = scalar keys %genes_with_2plus_brain_specific_staining;
my $no_of_genes_with_3plus_brain_specific_staining = scalar keys %genes_with_3plus_brain_specific_staining;
my $no_of_genes_with_4plus_brain_specific_staining = scalar keys %genes_with_4plus_brain_specific_staining;
my $pc_brain_genes_2_in_spg   	   = 0; if ($no_of_genes_with_2plus_brain_specific_staining > 0) { $pc_brain_genes_2_in_spg 	   = sprintf("%.2f",(($no_brain_genes_2_in_spg/$no_of_genes_with_2plus_brain_specific_staining)*100)); 		   }
my $pc_brain_genes_3_in_spg   	   = 0; if ($no_of_genes_with_3plus_brain_specific_staining > 0) { $pc_brain_genes_3_in_spg 	   = sprintf("%.2f",(($no_brain_genes_3_in_spg/$no_of_genes_with_3plus_brain_specific_staining)*100)); 		   }
my $pc_brain_genes_4_in_spg   	   = 0; if ($no_of_genes_with_4plus_brain_specific_staining > 0) { $pc_brain_genes_4_in_spg 	   = sprintf("%.2f",(($no_brain_genes_4_in_spg/$no_of_genes_with_4plus_brain_specific_staining)*100)); 		   }
my $pc_brain_genes_2_in_scyte      = 0; if ($no_of_genes_with_2plus_brain_specific_staining > 0) { $pc_brain_genes_2_in_scyte 	   = sprintf("%.2f",(($no_brain_genes_2_in_scyte/$no_of_genes_with_2plus_brain_specific_staining)*100)); 	   }
my $pc_brain_genes_3_in_scyte      = 0; if ($no_of_genes_with_3plus_brain_specific_staining > 0) { $pc_brain_genes_3_in_scyte 	   = sprintf("%.2f",(($no_brain_genes_3_in_scyte/$no_of_genes_with_3plus_brain_specific_staining)*100)); 	   }
my $pc_brain_genes_4_in_scyte      = 0; if ($no_of_genes_with_4plus_brain_specific_staining > 0) { $pc_brain_genes_4_in_scyte 	   = sprintf("%.2f",(($no_brain_genes_4_in_scyte/$no_of_genes_with_4plus_brain_specific_staining)*100)); 	   }
my $pc_brain_genes_2_in_stid       = 0; if ($no_of_genes_with_2plus_brain_specific_staining > 0) { $pc_brain_genes_2_in_stid       = sprintf("%.2f",(($no_brain_genes_2_in_stid/$no_of_genes_with_2plus_brain_specific_staining)*100)); 	   }
my $pc_brain_genes_3_in_stid       = 0; if ($no_of_genes_with_3plus_brain_specific_staining > 0) { $pc_brain_genes_3_in_stid       = sprintf("%.2f",(($no_brain_genes_3_in_stid/$no_of_genes_with_3plus_brain_specific_staining)*100)); 	   }
my $pc_brain_genes_4_in_stid       = 0; if ($no_of_genes_with_4plus_brain_specific_staining > 0) { $pc_brain_genes_4_in_stid       = sprintf("%.2f",(($no_brain_genes_4_in_stid/$no_of_genes_with_4plus_brain_specific_staining)*100)); 	   }
my $pc_brain_genes_2_only_in_spg   = 0; if ($no_of_genes_with_2plus_brain_specific_staining > 0) { $pc_brain_genes_2_only_in_spg   = sprintf("%.2f",(($no_brain_genes_2_only_in_spg/$no_of_genes_with_2plus_brain_specific_staining)*100));    }
my $pc_brain_genes_3_only_in_spg   = 0; if ($no_of_genes_with_3plus_brain_specific_staining > 0) { $pc_brain_genes_3_only_in_spg   = sprintf("%.2f",(($no_brain_genes_3_only_in_spg/$no_of_genes_with_3plus_brain_specific_staining)*100));    }
my $pc_brain_genes_4_only_in_spg   = 0; if ($no_of_genes_with_4plus_brain_specific_staining > 0) { $pc_brain_genes_4_only_in_spg   = sprintf("%.2f",(($no_brain_genes_4_only_in_spg/$no_of_genes_with_4plus_brain_specific_staining)*100));    }
my $pc_brain_genes_2_only_in_scyte = 0; if ($no_of_genes_with_2plus_brain_specific_staining > 0) { $pc_brain_genes_2_only_in_scyte = sprintf("%.2f",(($no_brain_genes_2_only_in_scyte/$no_of_genes_with_2plus_brain_specific_staining)*100));  }
my $pc_brain_genes_3_only_in_scyte = 0; if ($no_of_genes_with_3plus_brain_specific_staining > 0) { $pc_brain_genes_3_only_in_scyte = sprintf("%.2f",(($no_brain_genes_3_only_in_scyte/$no_of_genes_with_3plus_brain_specific_staining)*100));  }
my $pc_brain_genes_4_only_in_scyte = 0; if ($no_of_genes_with_4plus_brain_specific_staining > 0) { $pc_brain_genes_4_only_in_scyte = sprintf("%.2f",(($no_brain_genes_4_only_in_scyte/$no_of_genes_with_4plus_brain_specific_staining)*100));  }
my $pc_brain_genes_2_only_in_stid  = 0; if ($no_of_genes_with_2plus_brain_specific_staining > 0) { $pc_brain_genes_2_only_in_stid  = sprintf("%.2f",(($no_brain_genes_2_only_in_stid/$no_of_genes_with_2plus_brain_specific_staining)*100));   }
my $pc_brain_genes_3_only_in_stid  = 0; if ($no_of_genes_with_3plus_brain_specific_staining > 0) { $pc_brain_genes_3_only_in_stid  = sprintf("%.2f",(($no_brain_genes_3_only_in_stid/$no_of_genes_with_3plus_brain_specific_staining)*100));   }
my $pc_brain_genes_4_only_in_stid  = 0; if ($no_of_genes_with_4plus_brain_specific_staining > 0) { $pc_brain_genes_4_only_in_stid  = sprintf("%.2f",(($no_brain_genes_4_only_in_stid/$no_of_genes_with_4plus_brain_specific_staining)*100));   }
print "no. of genes seen: $no_of_genes_with_uncertain\n";
print "no. of genes with no protein expression data or no protein expression data classified other than 'uncertain': $no_with_no_data_or_no_certain_data ($pc_with_no_data%)\n";
print "no. of genes seen, excluding those with no data or data classified as 'uncertain': $no_of_genes_seen_excl_no_data\n";
print "no. of genes where antibody staining was performed for at least one of five germ cell types (because the gene has testis-enriched transcript expression): $no_stained_for_germ_cells ($pc_stained_for_germ_cells%)\n";
print "no. of genes where antibody staining was only performed for the vaguer category of 'seminiferous ducts': $no_stained_for_ducts ($pc_stained_for_ducts%)\n";
print "no. of genes with detectable somatic expression (applies regardless of whether we stain for ducts or germ cells as we look at Leydig cells in both): $no_expressed_in_soma ($pc_expressed_in_soma%)\n";
print "no. of genes with probable OR detectable germline expression (counting all genes irrespective of whether staining was performed only for ducts or for specific germ cells): $no_prob_or_def_in_germline ($pc_prob_or_def_in_germline%)\n";
print "no. of genes with probable germline expression (only applicable when we stain for ducts): $no_prob_in_germline ($pc_prob_in_germline%)\n";
print "no. of genes with detectable germline expression (only applicable when we stain explicitly for germ cells): $no_expressed_in_germline ($pc_expressed_in_germline%)\n";
print "no. of genes with detectable germline expression that are only expressed in the germline, not the soma (only applicable when we stain explicitly for germ cells): $no_expressed_only_in_germline ($pc_expressed_only_in_germline%)\n";
print "no. of genes with detectable germline expression that are expressed in all germ cell types (only applicable when we stain explicitly for germ cells): $no_expressed_in_all_germ_cells ($pc_expressed_in_all_germ_cells%)\n";
print "\n";
print "no. of genes with detectable germline expression that are expressed in spermatogonia (i.e. pre-meiosis) (only applicable when we stain explicitly for germ cells): $no_expressed_in_spg ($pc_expressed_in_spg%)\n";
print "no. of genes with detectable germline expression that are expressed in spermatocytes (i.e. meiosis) (only applicable when we stain explicitly for germ cells): $no_expressed_in_scyte ($pc_expressed_in_scyte%)\n";
print "no. of genes with detectable germline expression that are expressed in spermatids (i.e. post-meiosis) (only applicable when we stain explicitly for germ cells): $no_expressed_in_stid ($pc_expressed_in_stid%)\n";
print "\n";
print "no. of genes with detectable germline expression that are ONLY expressed in spermatogonia (i.e. pre-meiosis) (only applicable when we stain explicitly for germ cells): $no_only_expressed_in_spg ($pc_only_expressed_in_spg%)\n";
print "no. of genes with detectable germline expression that are ONLY expressed in spermatocytes (i.e. meiosis) (only applicable when we stain explicitly for germ cells): $no_only_expressed_in_scyte ($pc_only_expressed_in_scyte%)\n";
print "no. of genes with detectable germline expression that are ONLY expressed in spermatids (i.e. post-meiosis) (only applicable when we stain explicitly for germ cells): $no_only_expressed_in_stid ($pc_only_expressed_in_stid%)\n";
print "\n";
print "no. of genes with protein expression data and 2 or more specific associations (with the brain or a human-accelerated region): $no_of_genes_with_2plus_brain\n";
print "no. of genes with protein expression data and 2 or more specific associations (with the brain or a human-accelerated region) that have probable or confirmed germline expression: $no_brain_genes_2_in_germline ($pc_brain_genes_2_in_germline%)\n";
print "no. of genes with protein expression data and 3 or more specific associations (with the brain or a human-accelerated region): $no_of_genes_with_3plus_brain\n";
print "no. of genes with protein expression data and 3 or more specific associations (with the brain or a human-accelerated region) that have probable or confirmed germline expression: $no_brain_genes_3_in_germline ($pc_brain_genes_3_in_germline%)\n";
print "no. of genes with protein expression data and 4 or more specific associations (with the brain or a human-accelerated region): $no_of_genes_with_4plus_brain\n";
print "no. of genes with protein expression data and 4 or more specific associations (with the brain or a human-accelerated region) that have probable or confirmed germline expression: $no_brain_genes_4_in_germline ($pc_brain_genes_4_in_germline%)\n";
print "\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region): $no_of_genes_with_2plus_brain_specific_staining\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatogonia: $no_brain_genes_2_in_spg ($pc_brain_genes_2_in_spg%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatocytes: $no_brain_genes_2_in_scyte ($pc_brain_genes_2_in_scyte%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatids: $no_brain_genes_2_in_stid ($pc_brain_genes_2_in_stid%)\n";
print "\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region): $no_of_genes_with_3plus_brain_specific_staining\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatogonia: $no_brain_genes_3_in_spg ($pc_brain_genes_3_in_spg%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatocytes: $no_brain_genes_3_in_scyte ($pc_brain_genes_3_in_scyte%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatids: $no_brain_genes_3_in_stid ($pc_brain_genes_3_in_stid%)\n";
print "\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region): $no_of_genes_with_4plus_brain_specific_staining\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatogonia: $no_brain_genes_4_in_spg ($pc_brain_genes_4_in_spg%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatocytes: $no_brain_genes_4_in_scyte ($pc_brain_genes_4_in_scyte%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatids: $no_brain_genes_4_in_stid ($pc_brain_genes_4_in_stid%)\n";
print "\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatogonia: $no_brain_genes_2_only_in_spg ($pc_brain_genes_2_only_in_spg%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatogonia: $no_brain_genes_3_only_in_spg ($pc_brain_genes_3_only_in_spg%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatogonia: $no_brain_genes_4_only_in_spg ($pc_brain_genes_4_only_in_spg%)\n";
print "\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatocytes: $no_brain_genes_2_only_in_scyte ($pc_brain_genes_2_only_in_scyte%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatocytes: $no_brain_genes_3_only_in_scyte ($pc_brain_genes_3_only_in_scyte%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatocytes: $no_brain_genes_4_only_in_scyte ($pc_brain_genes_4_only_in_scyte%)\n";
print "\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatids: $no_brain_genes_2_only_in_stid ($pc_brain_genes_2_only_in_stid%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatids: $no_brain_genes_3_only_in_stid ($pc_brain_genes_3_only_in_stid%)\n";
print "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatids: $no_brain_genes_4_only_in_stid ($pc_brain_genes_4_only_in_stid%)\n";
print OUT "\n";
print OUT "no. of genes seen:\t$no_of_genes_with_uncertain\n";
print OUT "no. of genes with no protein expression data or no protein expression data classified other than 'uncertain':\t$no_with_no_data_or_no_certain_data ($pc_with_no_data%)\n";
print OUT "no. of genes seen, excluding those with no data or data classified as 'uncertain':\t$no_of_genes_seen_excl_no_data\n";
print OUT "no. of genes where antibody staining was performed for at least one of five germ cell types (because the gene has testis-enriched transcript expression):\t$no_stained_for_germ_cells ($pc_stained_for_germ_cells%)\n";
print OUT "no. of genes where antibody staining was only performed for the vaguer category of 'seminiferous ducts':\t$no_stained_for_ducts ($pc_stained_for_ducts%)\n";
print OUT "no. of genes with detectable somatic expression (applies regardless of whether we stain for ducts or germ cells as we look at Leydig cells in both):\t$no_expressed_in_soma ($pc_expressed_in_soma%)\n";
print OUT "no. of genes with probable OR detectable germline expression (counting all genes irrespective of whether staining was performed only for ducts or for specific germ cells):\t$no_prob_or_def_in_germline ($pc_prob_or_def_in_germline%)\n";
print OUT "no. of genes with probable germline expression (only applicable when we stain for ducts):\t$no_prob_in_germline ($pc_prob_in_germline%)\n";
print OUT "no. of genes with detectable germline expression (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_germline ($pc_expressed_in_germline%)\n";
print OUT "no. of genes with detectable germline expression that are only expressed in the germline, not the soma (only applicable when we stain explicitly for germ cells):\t$no_expressed_only_in_germline ($pc_expressed_only_in_germline%)\n";
print OUT "no. of genes with detectable germline expression that are expressed in all germ cell types (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_all_germ_cells ($pc_expressed_in_all_germ_cells%)\n";
print OUT "\n";
print OUT "no. of genes with detectable germline expression that are expressed in spermatogonia (i.e. pre-meiosis) (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_spg ($pc_expressed_in_spg%)\n";
print OUT "no. of genes with detectable germline expression that are expressed in spermatocytes (i.e. meiosis) (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_scyte ($pc_expressed_in_scyte%)\n";
print OUT "no. of genes with detectable germline expression that are expressed in spermatids (i.e. post-meiosis) (only applicable when we stain explicitly for germ cells):\t$no_expressed_in_stid ($pc_expressed_in_stid%)\n";
print OUT "\n";
print OUT "no. of genes with detectable germline expression that are ONLY expressed in spermatogonia (i.e. pre-meiosis) (only applicable when we stain explicitly for germ cells):\t$no_only_expressed_in_spg ($pc_only_expressed_in_spg%)\n";
print OUT "no. of genes with detectable germline expression that are ONLY expressed in spermatocytes (i.e. meiosis) (only applicable when we stain explicitly for germ cells):\t$no_only_expressed_in_scyte ($pc_only_expressed_in_scyte%)\n";
print OUT "no. of genes with detectable germline expression that are ONLY expressed in spermatids (i.e. post-meiosis) (only applicable when we stain explicitly for germ cells):\t$no_only_expressed_in_stid ($pc_only_expressed_in_stid%)\n";
print OUT "\n";
print OUT "no. of genes with protein expression data and 2 or more specific associations (with the brain or a human-accelerated region):\t$no_of_genes_with_2plus_brain\n";
print OUT "no. of genes with protein expression data and 2 or more specific associations (with the brain or a human-accelerated region) that have probable or confirmed germline expression:\t$no_brain_genes_2_in_germline ($pc_brain_genes_2_in_germline%)\n";
print OUT "no. of genes with protein expression data and 3 or more specific associations (with the brain or a human-accelerated region):\t$no_of_genes_with_3plus_brain\n";
print OUT "no. of genes with protein expression data and 3 or more specific associations (with the brain or a human-accelerated region) that have probable or confirmed germline expression:\t$no_brain_genes_3_in_germline ($pc_brain_genes_3_in_germline%)\n";
print OUT "no. of genes with protein expression data and 4 or more specific associations (with the brain or a human-accelerated region):\t$no_of_genes_with_4plus_brain\n";
print OUT "no. of genes with protein expression data and 4 or more specific associations (with the brain or a human-accelerated region) that have probable or confirmed germline expression:\t$no_brain_genes_4_in_germline ($pc_brain_genes_4_in_germline%)\n";
print OUT "\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region):\t$no_of_genes_with_2plus_brain_specific_staining\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatogonia:\t$no_brain_genes_2_in_spg ($pc_brain_genes_2_in_spg%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatocytes:\t$no_brain_genes_2_in_scyte ($pc_brain_genes_2_in_scyte%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatids:\t$no_brain_genes_2_in_stid ($pc_brain_genes_2_in_stid%)\n";
print OUT "\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) and protein expression data:\t$no_of_genes_with_3plus_brain_specific_staining\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatogonia:\t$no_brain_genes_3_in_spg ($pc_brain_genes_3_in_spg%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatocytes:\t$no_brain_genes_3_in_scyte ($pc_brain_genes_3_in_scyte%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatids:\t$no_brain_genes_3_in_stid ($pc_brain_genes_3_in_stid%)\n";
print OUT "\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region):\t$no_of_genes_with_4plus_brain_specific_staining\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatogonia:\t$no_brain_genes_4_in_spg ($pc_brain_genes_4_in_spg%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatocytes:\t$no_brain_genes_4_in_scyte ($pc_brain_genes_4_in_scyte%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) have protein-level expression in spermatids:\t$no_brain_genes_4_in_stid ($pc_brain_genes_4_in_stid%)\n";
print OUT "\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatogonia:\t$no_brain_genes_2_only_in_spg ($pc_brain_genes_2_only_in_spg%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatogonia:\t$no_brain_genes_3_only_in_spg ($pc_brain_genes_3_only_in_spg%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatogonia:\t$no_brain_genes_4_only_in_spg ($pc_brain_genes_4_only_in_spg%)\n";
print OUT "\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatocytes:\t$no_brain_genes_2_only_in_scyte ($pc_brain_genes_2_only_in_scyte%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatocytes:\t$no_brain_genes_3_only_in_scyte ($pc_brain_genes_3_only_in_scyte%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatocytes:\t$no_brain_genes_4_only_in_scyte ($pc_brain_genes_4_only_in_scyte%)\n";
print OUT "\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 2 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatids:\t$no_brain_genes_2_only_in_stid ($pc_brain_genes_2_only_in_stid%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 3 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatids:\t$no_brain_genes_3_only_in_stid ($pc_brain_genes_3_only_in_stid%)\n";
print OUT "no. of genes with protein expression data, explicit staining for germ cells, and 4 or more specific associations (with the brain or a human-accelerated region) that have protein-level expression ONLY in spermatids:\t$no_brain_genes_4_only_in_stid ($pc_brain_genes_4_only_in_stid%)\n";
close(OUT) or die $!;

print "COMPLETE. Did we use a conservative gene list? $conservative\n";

exit 1;