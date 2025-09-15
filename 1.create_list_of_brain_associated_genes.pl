use strict;
use warnings;

# REQUIREMENTS
my $in_dir   = 'gene_lists_used_to_make_master_list_of_brain_genes'; # manually created; see notes to Supplementary Table 1 for sources
my $in_file1 = 'Ens112.gene_annotations.txt'; # from Ensembl BioMart (GRCh38.p14): Gene stable ID, Gene name, Gene description, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand, Gene type, Phenotype description, Source name
my $in_file2 = 'Ens112.gene_synonyms.txt'; # from Ensembl BioMart (GRCh38.p14): Gene name, Gene synonym
if (!(-d($in_dir)))   { print "ERROR: cannot find $in_dir\n";   exit 1; }
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# PARAMETERS
my $conservative = 'yes'; # 'no'; # conservative criteria are to restrict analysis only to those genes associated with a defined brain phenotype: macro/megalencephaly, autism, schizophrenia, epilepsy, brain weight
my @conservative = ("Bastos 2022","DeCasien 2022","Dhindsa 2025 (autism)","Dhindsa 2025 (epilepsy)","Qiu 2022","Rylaarsdam 2019","Sanders 2015","Satterstrom 2020","Ripke 2014","Singh 2022","SFARI 2025","Owen 2023","Perucca 2020","Thakran 2020","Rastin 2023","Boddy 2017","Seidlitz 2023");
my %conservative = map {$_ => 1} @conservative;

# OUTPUT
my $out_file1 = ''; my $out_file2 = '';
if ($conservative eq 'no')
	{ $out_file1 = 'brain_associated_genes.txt';
	  $out_file2 = 'dataframe_for_making_upset_plot_of_brain_associated_genes.txt'; # NOTE: for use with the R package ComplexUpset, following the instructions https://krassowski.github.io/complex-upset/articles/Examples_R.html
	}
elsif ($conservative eq 'yes')
	{ $out_file1 = 'brain_associated_genes.conservative.txt';
	  $out_file2 = 'dataframe_for_making_upset_plot_of_brain_associated_genes.conservative.txt';
	}
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;
print OUT1 "Gene ID\tReference(s)\t";
print OUT1 "Associated with macrocephaly or megalencephaly (in Bastos 2022 or DeCasien 2022)?\t";
print OUT1 "Associated with autism (in Dhindsa 2025, Qiu 2022, Rylaarsdam 2019, Sanders 2015, Satterstrom 2020, or SFARI 2025)?\t";
print OUT1 "Associated with schizophrenia (in Ripke 2014, Singh 2022 or Owen 2023)?\t";
print OUT1 "Associated with epilepsy (in Dhindsa 2025, Perucca 2020, Thakran 2020 or Rastin 2023)?\t";
print OUT1 "Associated with intelligence or educational attainment (in Lee 2018 or Savage 2018)?\t";
print OUT1 "Associated with brain weight, either positively or negatively (in Boddy 2017 or Seidlitz 2023)?\t";
print OUT1 "Associated with a human-accelerated region (in Wei 2019)?\t";
print OUT1 "No. of SPECIFIC associations (that is, with any of the seven aforementioned categories)\n";

# STORE GENE ANNOTATION DATA SO WE CAN (A) RESTRICT EACH GENE LIST ONLY TO KNOWN PROTEIN-CODING GENES, AND (B) ASSIGN EACH GENE ITS ENSEMBL ID, IF ONE HAS NOT OTHERWISE BEEN MADE AVAILABLE
# WE RESTRICT EACH GENE LIST ONLY TO KNOWN PROTEIN-CODING GENES ON THE MAJOR AUTOSOMES + XY (I.E. NO CHROMOSOMAL PATCHES)
my @acceptable_chrs = (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT/);
my %acceptable_chrs = map {$_ => 1} @acceptable_chrs;
my %approved_gene_ids = (); my %approved_gene_names = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0];
	  next if (!(defined($line[1]))); # CHECKPOINT: discard genes which have not been annotated with an HGNC name
	  my $gene_name = $line[1]; my $chr = $line[3]; my $gene_type = $line[7];
	  next if (!(exists($acceptable_chrs{$chr}))); # CHECKPOINT: skip genes which are not on an autosome or sex chromosome (i.e. discard all chromosome patches)
	  next if ($gene_type ne 'protein_coding'); # CHECKPOINT: discard genes which are not protein-coding
	  $approved_gene_ids{$gene_id} = $gene_name; $approved_gene_names{$gene_name}{$gene_id}++;
	}
close(IN) or die $!;

# STORE GENE SYNONYMS SO THAT WHEN PARSING THE LISTS OF BRAIN-ASSOCIATED GENES (SOME OF WHICH ARE QUITE OLD), WE CAN TEST TO SEE WHETHER THE NAMES REPORTED IN EACH LIST ARE ACTUALLY SYNONYMS FOR A MORE RECENT (UPDATED) GENE NAME
my %synonyms = ();
open(IN,$in_file2) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if ( (!(defined($line[0]))) or (!(defined($line[1]))) ); # CHECKPOINT: discard genes which have not been annotated with an HGNC name
	  my $gene_name = $line[0]; my $gene_synonym = $line[1];
	  $synonyms{$gene_synonym}{$gene_name}++;
	}
close(IN) or die $!;

# CREATE A NON-REDUNDANT LIST OF BRAIN-ASSOCIATED GENES
my %gene_ids = ();
opendir(DIR,$in_dir) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
foreach my $file (@files)
	{ next if (($file eq '.') or ($file eq '..'));
	  my $list_source;
	  if ($file =~ /^(.*?)\.txt$/) { $list_source = $1; } else { print "ERROR: cannot parse the name of $file\n"; exit 1; }
	  next if ( ($conservative eq 'yes') and (!(exists($conservative{$list_source}))) ); # CHECKPOINT: if we are restricting the dataset to a conservative set of brain-associated genes, exclude all other gene lists
	  
	  # for the Boddy 2017 file, we wish to only keep those gene IDs whose names appear in the list TWICE
	  # this is because the file aggregates data from both Tables S7 and S8 of their paper and what we need are the genes which appear in both. These tables list genes with a coevolutionary relationship with brain mass and encephalisation quotient, respectively.
	  my %gene_count_boddy2017 = ();
	  if ($list_source eq 'Boddy 2017')
		{ open(IN,"$in_dir/$file") or die $!;
		  while(<IN>)
			{ my $val = $_; chomp($val);
			  $val =~ s/^\s+//; $val =~ s/\s+$//;
			  if ($val =~ s/^(.*?)\s+$//g) { $val = $1; }
			  if ($val =~ /^ENSG/) # if the gene already has an Ensembl ID, we don't need to do anything else other than confirm it's on the approved list
				{ if (exists($approved_gene_ids{$val}))
					{ $gene_count_boddy2017{$val}++; }
				}
			  else
				{ if (exists($approved_gene_names{$val})) # else if it's a gene name, we first need to convert it to a gene ID
					{ my $num_possible_gene_ids = scalar keys %{$approved_gene_names{$val}};
					  if ($num_possible_gene_ids == 1) # we can only convert a gene name to a gene ID if there is a one-to-one correspondence between the two. This would exclude some possibilities IF we had not already removed those genes only included on chromosome patches, e.g. the Wei 2019 list contains ZC3H3 which could be either ENSG00000014164 (on chr 8) or ENSG00000282684 (on HSCHR8_3_CTG7). Other ambiguities persist regardless of chromosome patches: CD99 is either ENSG00000002586 (on chr X) or ENSG00000292348 (on chr Y).
						{ while((my $gene_id,my $irrel)=each(%{$approved_gene_names{$val}}))
							{ if (exists($approved_gene_ids{$gene_id}))
								{ $gene_count_boddy2017{$gene_id}++; }
							}
						}
					  else
						{ print "discarding $val because it has $num_possible_gene_ids gene IDs\n"; }
					}
				}
			}
		  close(IN) or die $!;	
		}
	  
	  open(IN,"$in_dir/$file") or die $!;
	  while(<IN>)
		{ my $val = $_; chomp($val);
		  $val =~ s/^\s+//; $val =~ s/\s+$//;
		  if ($val =~ s/^(.*?)\s+$//g) { $val = $1; }
		  my $gene_id = '';
		  if ($val =~ /^ENSG/) # if the gene already has an Ensembl ID, we don't need to do anything else other than confirm it's on the approved list
			{ if (exists($approved_gene_ids{$val}))
				{ $gene_id = $val; }
			}
		  else
			{ if (exists($approved_gene_names{$val})) # else if it's a gene name, we first need to convert it to a gene ID
				{ my $num_possible_gene_ids = scalar keys %{$approved_gene_names{$val}};
				  if ($num_possible_gene_ids == 1) # we can only convert a gene name to a gene ID if there is a one-to-one correspondence between the two. This would exclude some possibilities IF we had not already removed those genes only included on chromosome patches, e.g. the Wei 2019 list contains ZC3H3 which could be either ENSG00000014164 (on chr 8) or ENSG00000282684 (on HSCHR8_3_CTG7). Other ambiguities persist regardless of chromosome patches: CD99 is either ENSG00000002586 (on chr X) or ENSG00000292348 (on chr Y).
					{ while((my $id,my $irrel)=each(%{$approved_gene_names{$val}}))
						{ if (exists($approved_gene_ids{$id}))
							{ $gene_id = $id; }
						}
					}
				  else
					{ print "discarding $val because it has $num_possible_gene_ids gene IDs\n"; }
				}
			}
		  if ($gene_id eq '') # if we cannot assign a gene ID, it may be because the gene name used in the list is out-of-date, and instead of a synonym of an updated one - so let's update that where appropriate
			{ print "$file: discarding $val because we cannot obtain an approved Ensembl gene ID for it\n";
			  if (exists($synonyms{$val}))
				{ my $num_poss_synonyms = scalar keys %{$synonyms{$val}};
				  while((my $poss_name,my $irrel)=each(%{$synonyms{$val}}))
					{ print "... but note that $val is a synonym for $poss_name\n";
					  if (exists($approved_gene_names{$poss_name})) # else if it's a gene name, we first need to convert it to a gene ID
						{ my $num_possible_gene_ids = scalar keys %{$approved_gene_names{$poss_name}};
						  if (($num_poss_synonyms == 1) and ($num_possible_gene_ids == 1)) # we can only convert a gene name to a gene ID if there is a one-to-one correspondence between the two. This would exclude some possibilities IF we had not already removed those genes only included on chromosome patches, e.g. the Wei 2019 list contains ZC3H3 which could be either ENSG00000014164 (on chr 8) or ENSG00000282684 (on HSCHR8_3_CTG7). Other ambiguities persist regardless of chromosome patches: CD99 is either ENSG00000002586 (on chr X) or ENSG00000292348 (on chr Y).
							{ while((my $id,my $irrel)=each(%{$approved_gene_names{$poss_name}}))
								{ print "... ... $poss_name updated to $id ($approved_gene_ids{$id}))\n";
								  if (exists($approved_gene_ids{$id}))
									{ $gene_id = $id; }
								}
							}
						  else
							{ print "discarding $val because it has $num_possible_gene_ids gene IDs\n"; }
						}
					}
				}
			}
		  next if ($gene_id eq ''); # CHECKPOINT: skip if we do not have an approved gene ID
		  next if (!(exists($approved_gene_ids{$gene_id}))); # CHECKPOINT: skip if we do not have an approved gene ID (secondary check; does the same as the above)
		  next if ( ($list_source eq 'Boddy 2017') and (!(exists($gene_count_boddy2017{$gene_id}))) ); # CHECKPOINT: skip if we are parsing the Boddy 2017 gene list and do not recognise the gene name
		  next if ( ($list_source eq 'Boddy 2017') and ($gene_count_boddy2017{$gene_id} != 2) ); # CHECKPOINT: skip if we are parsing the Boddy 2017 gene list and the gene ID did not appear twice (a necessity, for reasons given above)
		  $gene_ids{$gene_id}{$list_source}++;
		}
	  close(IN) or die $!;
	}

my @gene_ids = ();
while((my $gene_id,my $irrel)=each(%gene_ids))
	{ push(@gene_ids,$gene_id); }
my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
my %associations_per_gene = ();
foreach my $gene_id (@sorted_gene_ids)
	{ my @list_sources = (); my %sources = ();
	  while((my $list_source,my $irrel)=each(%{$gene_ids{$gene_id}}))
		{ push(@list_sources,$list_source);
		  $sources{$list_source}++;
		}
	  my @sorted_list_sources = sort {$a cmp $b} @list_sources;
	  my $list_sources = join(", ",@sorted_list_sources);
	  my $associated_macro = ''; my $associated_autism = ''; my $associated_schizo = ''; my $associated_epilepsy = ''; my $associated_IQ_or_EA = ''; my $associated_weight = ''; my $associated_har = '';
	  my $no_of_associations = 0;
	  if ( (exists($sources{'Bastos 2022'})) or (exists($sources{'DeCasien 2022'})) )
		{ $associated_macro = 'yes';    $no_of_associations++; $associations_per_gene{$gene_id}{'macro/megalencephaly'}++; 	    }
	  if ( (exists($sources{'Dhindsa 2025 (autism)'})) or (exists($sources{'Qiu 2022'})) or (exists($sources{'Rylaarsdam 2019'})) or (exists($sources{'Sanders 2015'})) or (exists($sources{'Satterstrom 2020'})) or (exists($sources{'SFARI 2025'})) )
		{ $associated_autism = 'yes';   $no_of_associations++; $associations_per_gene{$gene_id}{'autism'}++;            }
	  if ( (exists($sources{'Ripke 2014'})) or (exists($sources{'Singh 2022'})) or (exists($sources{'Owen 2023'})) )
		{ $associated_schizo = 'yes';   $no_of_associations++; $associations_per_gene{$gene_id}{'schizophrenia'}++;     }
	  if ( (exists($sources{'Dhindsa 2025 (epilepsy)'})) or (exists($sources{'Perucca 2020'})) or (exists($sources{'Thakran 2020'})) or (exists($sources{'Rastin 2023'})) )
		{ $associated_epilepsy = 'yes'; $no_of_associations++; $associations_per_gene{$gene_id}{'epilepsy'}++;          }
	  if ( (exists($sources{'Lee 2018'})) or (exists($sources{'Savage 2018'})) )
		{ $associated_IQ_or_EA = 'yes'; $no_of_associations++; $associations_per_gene{$gene_id}{'IQ or attainment'}++;  }
	  if ( (exists($sources{'Boddy 2017'})) or (exists($sources{'Seidlitz 2023'})) )
		{ $associated_weight = 'yes';   $no_of_associations++; $associations_per_gene{$gene_id}{'brain weight'}++;      }
	  if (exists($sources{'Wei 2019'}))
		{ $associated_har = 'yes';      $no_of_associations++; $associations_per_gene{$gene_id}{'human accelerated'}++; }
	  
	  if ($no_of_associations == 0)
		{ $associations_per_gene{$gene_id}{'other'}++;
		  if ($conservative eq 'yes')
		   { print "ERROR: we are making a conservative list of brain-associated genes and therefore do not expect to see any gene (like $gene_id) which has no clear functional association\n"; exit 1; }
		}
	  
	  print OUT1 "$gene_id\t$list_sources\t$associated_macro\t$associated_autism\t$associated_schizo\t$associated_epilepsy\t$associated_IQ_or_EA\t$associated_weight\t$associated_har\t$no_of_associations\n";
	}
close(OUT1) or die $!;

my $num_genes = scalar keys %associations_per_gene;
print "number of brain-associated genes processed for upset plot: $num_genes\n";

# OUTPUT A DATA FRAME FOR MAKING AN UPSET PLOT USING THE R PACKAGE ComplexUpset, FOLLOWING THE VIGNETTE AT https://krassowski.github.io/complex-upset/articles/Examples_R.html
my @associations = ();
if ($conservative eq 'no')
	{ @associations = ("macro/megalencephaly","brain weight","autism","schizophrenia","epilepsy","IQ or attainment","human accelerated","other"); }
elsif ($conservative eq 'yes')
	{ @associations = ("macro/megalencephaly","brain weight","autism","schizophrenia","epilepsy"); }
my $associations_line = join("\t",@associations);
print OUT3 "Gene\t$associations_line\n";	
@gene_ids = ();
while((my $gene_id,my $irrel)=each(%associations_per_gene))
	{ push(@gene_ids,$gene_id); }
@sorted_gene_ids = sort {$a cmp $b} @gene_ids;
foreach my $gene_id (@sorted_gene_ids)
	{ my $out_line = '';
	  foreach my $association (@associations)
		{ my $true_or_false = 'FALSE'; # the default
		  if (exists($associations_per_gene{$gene_id}{$association}))
			{ $true_or_false = 'TRUE'; }
		  $out_line .= "$true_or_false\t";
		}
	  $out_line =~ s/\t$//;
	  print OUT2 "$gene_id\t$out_line\n";
	}
close(OUT2) or die $!;

print "COMPLETE. Did we generate a conservative gene list? $conservative\n";

exit 1;