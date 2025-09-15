use strict;
use warnings;

# PARAMETERS
my $conservative = 'no'; # 'yes'; # conservative criteria are to restrict analysis only to those genes associated with a defined brain phenotype: macro/megalencephaly, autism, schizophrenia, epilepsy, brain weight
my $total_num_of_cells_in_atlas = 60427; # of which 405 are Sertoli, 2964 are endothelia, and 7419 are myoid/Leydig = 10,788 of which are somatic cells

# REQUIREMENTS
my $in_file1 = ''; # from 8.summarise_germline_protein_expression_of_brain_genes.pl
if ($conservative eq 'no')
	{ $in_file1 = 'C:/Users/User/Desktop/testis_brain/germline_expression_at_protein_level_of_brain_associated_genes.txt'; } # from 2.summarise_germline_protein_expression_of_brain_genes.pl
elsif ($conservative eq 'yes')
	{ $in_file1 = 'C:/Users/User/Desktop/testis_brain/germline_expression_at_protein_level_of_brain_associated_genes.conservative.txt'; } # from 2.summarise_germline_protein_expression_of_brain_genes.pl
my $in_file2 = 'results/human_adult.atlas.txt'; # produced by https://github.com/sjbush/spg_atlas/12a.create_summary_table_of_whole_testes_atlas.pl
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_file1 = ''; my $out_file2 = ''; my $out_file3 = '';
if ($conservative eq 'no')
	{ $out_file1 = 'num_of_brain_developmental_genes_detectable_in_whole_testis_atlas.summary.txt';
	  $out_file2 = 'num_of_brain_developmental_genes_detectable_in_whole_testis_atlas.full.txt';
	  $out_file3 = 'data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.txt';
	}
elsif ($conservative eq 'yes')
	{ $out_file1 = 'num_of_brain_developmental_genes_detectable_in_whole_testis_atlas.summary.conservative.txt';
	  $out_file2 = 'num_of_brain_developmental_genes_detectable_in_whole_testis_atlas.full.conservative.txt';
	  $out_file3 = 'data_for_making_barplot_of_germline_transcription_of_brain_associated_genes.conservative.txt';
	}
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!;
print OUT1 "Gene name\tEnsembl gene ID\t% of somatic cells in which this gene is detected\t% of germline cells in which this gene is detected\tRatio of germline to somatic cells in which this gene is detected\tCell cluster(s) in which this gene is differentially expressed\n";
print OUT3 "Gene category\tNumber\n";

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
my $de_in_somatic = 0; my $de_in_germline = 0; my $de_in_myoid = 0; my $de_in_sertoli = 0; my $de_in_endothelia = 0; my $de_in_undiff_SPG = 0; my $de_in_diff_SPG = 0; my $de_in_scyte = 0; my $de_in_early_stid1 = 0; my $de_in_early_stid2 = 0; my $de_in_late_stid1 = 0; my $de_in_late_stid2 = 0;
my $no_expressed_in_gt1pc = 0; my $de_in_spg = 0; my $de_in_stid = 0; my $only_de_in_spg = 0; my $only_de_in_scyte = 0; my $only_de_in_stid = 0;
my %key_genes_seen = ();
open(IN,$in_file2) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($. == 1) { print OUT2 "$line\n"; }
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0]; my $ens_id = $line[1];
	  next if (!(exists($key_genes{$ens_id}))); # CHECKPOINT: restrict analysis only to brain-developmental genes
	  my $num_of_non_germline_cells = $line[16]+$line[17]+$line[18];
	  my $pct_of_non_germline_cells = sprintf("%.2f",(($num_of_non_germline_cells/$total_num_of_cells_in_atlas)*100));
	  my $num_of_germline_cells = $line[19]+$line[20]+$line[21]+$line[22]+$line[23]+$line[24]+$line[25];
	  my $pct_of_germline_cells = sprintf("%.2f",(($num_of_germline_cells/$total_num_of_cells_in_atlas)*100));
	  my $ratio_of_germline_to_non_germline_cells = 'NA';
	  if ($num_of_non_germline_cells > 0) { $ratio_of_germline_to_non_germline_cells = sprintf("%.2f",($num_of_germline_cells/$num_of_non_germline_cells)); }
	  my $clusters = $line[49]; my $proporp_in_this_cluster = $line[50]; my $proporp_in_all_other_clusters = $line[51]; my $abs_diff_in_proporp = $line[52]; my $max_abs_diff = $line[53]; my $avg_log2fc = $line[54]; my $adj_p = $line[56];
	  my $is_biomarker_of = $line[57];
	  my %de_in = ();
	  if ($clusters =~ /Cluster 0 \(late spermatid 1\)/)  { $de_in_late_stid1++;  $de_in_germline++; $de_in{'spermatid'}++; 	}
	  if ($clusters =~ /Cluster 1 \(early spermatid 1\)/) { $de_in_early_stid1++; $de_in_germline++; $de_in{'spermatid'}++; 	}
	  if ($clusters =~ /Cluster 2 \(early spermatid 2\)/) { $de_in_early_stid2++; $de_in_germline++; $de_in{'spermatid'}++; 	}
	  if ($clusters =~ /Cluster 3 \(myoid\/Leydig\)/) 	  { $de_in_myoid++; 	  $de_in_somatic++;  $de_in{'myoid/Leydig'}++; 	}
	  if ($clusters =~ /Cluster 4 \(spermatocyte\)/) 	  { $de_in_scyte++; 	  $de_in_germline++; $de_in{'spermatocyte'}++;	}
	  if ($clusters =~ /Cluster 5 \(late spermatid 2\)/)  { $de_in_late_stid2++;  $de_in_germline++; $de_in{'spermatid'}++; 	}
	  if ($clusters =~ /Cluster 6 \(SSC\)/) 			  { $de_in_undiff_SPG++;  $de_in_germline++; $de_in{'spermatogonia'}++; }
	  if ($clusters =~ /Cluster 7 \(endothelia\)/) 		  { $de_in_endothelia++;  $de_in_somatic++;  $de_in{'endothelia'}++;	}
	  if ($clusters =~ /Cluster 8 \(spermatogonia\)/)	  { $de_in_diff_SPG++;    $de_in_germline++; $de_in{'spermatogonia'}++; }
	  if ($clusters =~ /Cluster 9 \(Sertoli\)/) 		  { $de_in_sertoli++;     $de_in_somatic++;  $de_in{'Sertoli'}++; 		}
	  if ($is_biomarker_of =~ /SSC/) 		   { $is_biomarker_of = 'undiff SPG'; 			  }
	  if ($is_biomarker_of =~ /endothelia/)    { $is_biomarker_of = 'endothelia/macrophage';  }
	  if ($is_biomarker_of =~ /spermatogonia/) { $is_biomarker_of = 'diff SPG/early meiosis'; }
	  $clusters =~ s/Cluster 0 \(late spermatid 1\)/late spermatid 1/;
	  $clusters =~ s/Cluster 1 \(early spermatid 1\)/early spermatid 1/;
	  $clusters =~ s/Cluster 2 \(early spermatid 2\)/early spermatid 2/;
	  $clusters =~ s/Cluster 3 \(myoid\/Leydig\)/myoid\/Leydig/;
	  $clusters =~ s/Cluster 4 \(spermatocyte\)/spermatocyte/;
	  $clusters =~ s/Cluster 5 \(late spermatid 2\)/late spermatid 2/;
	  $clusters =~ s/Cluster 6 \(SSC\)/undiff SPG/;
	  $clusters =~ s/Cluster 7 \(endothelia\)/endothelia\/macrophage/;
	  $clusters =~ s/Cluster 8 \(spermatogonia\)/diff SPG\/early meiosis/;
	  $clusters =~ s/Cluster 9 \(Sertoli\)/Sertoli/;
	  $line[49] = $clusters; $line[57] = $is_biomarker_of;
	  $line = join("\t",@line);
	  print OUT1 "$gene_name\t$ens_id\t$pct_of_non_germline_cells\t$pct_of_germline_cells\t$ratio_of_germline_to_non_germline_cells\t$clusters\n";
	  print OUT2 "$line\n";
	  $key_genes_seen{$ens_id}++;
	  if ($pct_of_germline_cells >= 1) { $no_expressed_in_gt1pc++; }
	  my $num_de_in = scalar keys %de_in;
	  if 	(exists($de_in{'spermatogonia'})) { $de_in_spg++;  } # IMPORTANT: 'spermatogonia' and 'spermatid' both comprise multiple clusters (e.g. 'spermatocyte' comprises 'undiff SPG' and 'diff SPG') which is why we increment these counters here. We don't need to do this for spermatocytes as there is only one spermatocyte cluster to begin with
	  elsif (exists($de_in{'spermatid'}))	  { $de_in_stid++; }
	  if ($num_de_in == 1)
		{ if 	(exists($de_in{'spermatogonia'})) { $only_de_in_spg++;   }
		  elsif (exists($de_in{'spermatocyte'}))  { $only_de_in_scyte++; }
		  elsif (exists($de_in{'spermatid'}))	  { $only_de_in_stid++;  }
		}
	}
close(IN) or die $!;

my $pc_expressed_in_gt1pc = sprintf("%.2f",(($no_expressed_in_gt1pc/$total_key_genes)*100));
my $pc_de_in_spg     	  = sprintf("%.2f",(($de_in_spg/$total_key_genes)*100));
my $pc_de_in_scyte   	  = sprintf("%.2f",(($de_in_scyte/$total_key_genes)*100));
my $pc_de_in_stid    	  = sprintf("%.2f",(($de_in_stid/$total_key_genes)*100));
my $pc_only_de_in_spg     = sprintf("%.2f",(($only_de_in_spg/$total_key_genes)*100));
my $pc_only_de_in_scyte   = sprintf("%.2f",(($only_de_in_scyte/$total_key_genes)*100));
my $pc_only_de_in_stid    = sprintf("%.2f",(($only_de_in_stid/$total_key_genes)*100));

print OUT1 "\n";
print OUT1 "No. of genes expressed in >= 1% of germ cells\t$no_expressed_in_gt1pc ($pc_expressed_in_gt1pc%)\n";
print OUT1 "No. of genes differentially expressed in spermatogonia\t$de_in_spg ($pc_de_in_spg%)\n";
print OUT1 "No. of genes differentially expressed in spermatocytes\t$de_in_scyte ($pc_de_in_scyte%)\n";
print OUT1 "No. of genes differentially expressed in spermatids\t$de_in_stid ($pc_de_in_stid%)\n";
print OUT1 "No. of genes only differentially expressed in spermatogonia\t$only_de_in_spg ($pc_only_de_in_spg%)\n";
print OUT1 "No. of genes only differentially expressed in spermatocytes\t$only_de_in_scyte ($pc_only_de_in_scyte%)\n";
print OUT1 "No. of genes only differentially expressed in spermatids\t$only_de_in_stid ($pc_only_de_in_stid%)\n";

my $pct_in_somatic     = sprintf("%.2f",(($de_in_somatic/$total_key_genes)*100));
my $pct_in_myoid 	   = sprintf("%.2f",(($de_in_myoid/$total_key_genes)*100));
my $pct_in_sertoli     = sprintf("%.2f",(($de_in_sertoli/$total_key_genes)*100));
my $pct_in_endothelia  = sprintf("%.2f",(($de_in_endothelia/$total_key_genes)*100));
my $pct_in_germline    = sprintf("%.2f",(($de_in_germline/$total_key_genes)*100));
my $pct_in_undiff_SPG  = sprintf("%.2f",(($de_in_undiff_SPG/$total_key_genes)*100));
my $pct_in_diff_SPG    = sprintf("%.2f",(($de_in_diff_SPG/$total_key_genes)*100));
my $pct_in_scyte 	   = sprintf("%.2f",(($de_in_scyte/$total_key_genes)*100));
my $pct_in_early_stid1 = sprintf("%.2f",(($de_in_early_stid1/$total_key_genes)*100));
my $pct_in_early_stid2 = sprintf("%.2f",(($de_in_early_stid2/$total_key_genes)*100));
my $pct_in_late_stid1  = sprintf("%.2f",(($de_in_late_stid1/$total_key_genes)*100));
my $pct_in_late_stid2  = sprintf("%.2f",(($de_in_late_stid2/$total_key_genes)*100));

print "number of genes differentially expressed in a somatic cell cluster = $de_in_somatic ($pct_in_somatic%)\n";
print "number of genes differentially expressed in myoid/Leydig cells = $de_in_myoid ($pct_in_myoid%)\n";
print "number of genes differentially expressed in Sertoli cells = $de_in_sertoli ($pct_in_sertoli%)\n";
print "number of genes differentially expressed in endothelia/macrophages = $de_in_endothelia ($pct_in_endothelia%)\n";
print "\n";
print "number of genes differentially expressed in a germ cell cluster = $de_in_germline ($pct_in_germline%)\n";
print "number of genes differentially expressed in undiff SPG = $de_in_undiff_SPG ($pct_in_undiff_SPG%)\n";
print "number of genes differentially expressed in diff SPG = $de_in_diff_SPG ($pct_in_diff_SPG%)\n";
print "number of genes differentially expressed in spermatocyte = $de_in_scyte ($pct_in_scyte%)\n";
print "number of genes differentially expressed in early spermatid 1 = $de_in_early_stid1 ($pct_in_early_stid1%)\n";
print "number of genes differentially expressed in early spermatid 2 = $de_in_early_stid2 ($pct_in_early_stid2%)\n";
print "number of genes differentially expressed in late spermatid 1 = $de_in_late_stid1 ($pct_in_late_stid1%)\n";
print "number of genes differentially expressed in late spermatid 2 = $de_in_late_stid2 ($pct_in_late_stid2%)\n";

# which genes on the checklist aren't seen in $in_file2?
while((my $gene_name,my $irrel)=each(%key_genes))
	{ if (!(exists($key_genes_seen{$gene_name})))
		{ print "missing $gene_name\n"; }
	}

# SUMMARISE THE KEY FIGURES IN ORDER TO MAKE A BARPLOT SUMMARISING THE GERMLINE TRANSCRIPTION OF BRAIN-ASSOCIATED GENES, I.E. NUMBER OF GENES DIFFERENTIALLY EXPRESSED IN A GIVEN CLUSTER
print OUT3 "undiff SPG\t$de_in_undiff_SPG\n";
print OUT3 "diff SPG\t$de_in_diff_SPG\n";
print OUT3 "spermatocyte\t$de_in_scyte\n";
print OUT3 "early spermatid 1\t$de_in_early_stid1\n";
print OUT3 "early spermatid 2\t$de_in_early_stid2\n";
print OUT3 "late spermatid 1\t$de_in_late_stid1\n";
print OUT3 "late spermatid 2\t$de_in_late_stid2\n";
print OUT3 "endothelia\t$de_in_endothelia\n";
print OUT3 "myoid/Leydig\t$de_in_myoid\n";
print OUT3 "Sertoli\t$de_in_sertoli\n";

close(OUT1) or die $!; close(OUT2) or die $!; close(OUT3) or die $!;

print "COMPLETE. Did we use a conservative gene list? $conservative\n";

exit 1;