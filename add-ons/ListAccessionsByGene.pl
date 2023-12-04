#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

sub ScriptDir { return './' if ($0 !~ /\//); $0 =~ /^(.*\/)[^\/]+$/; return $1; }
use lib ScriptDir();
use bureaucracy;


if (@ARGV != 1) {
    die "\n  USAGE: ./ListAccesssionsByGene.pl [Mirage-Results]\n\n";
}

my %SpeciesGeneToAccession;

my $species_msas_dir_name = ConfirmDirectory(ConfirmDirectory($ARGV[0]).'Species-MSAs');
my $AllSpeciesDir = OpenDirectory($species_msas_dir_name);

my $longest_species_name_len = 0;
my $longest_gene_name_len = 0;

while (my $species = readdir($AllSpeciesDir)) {

    next if ($species =~ /^\./);
    $species =~ s/\/$//;

    my $species_ali_dir_name = $species_msas_dir_name.$species.'/alignments/';
    next if (!(-d $species_ali_dir_name));

    if (length($species) > $longest_species_name_len) {
	$longest_species_name_len = length($species);
    }

    my $SpeciesAliDir = OpenDirectory($species_ali_dir_name);

    while (my $fname = readdir($SpeciesAliDir)) {

	next if ($fname !~ /^(\S+)\.afa$/);
	my $gene = lc($1);

	if (length($gene) > $longest_gene_name_len) {
	    $longest_gene_name_len = length($gene);
	}

	my $has_accessions = 0;
	my %Accessions;

	$fname = $species_ali_dir_name.$fname;

	my $AccessionGrep = OpenSystemCommand("grep -i 'accession:' \"$fname\"");
	while (my $line = <$AccessionGrep>) {

	    $line = uc($line);

	    next if ($line !~ /ACCESSION:(\S+)/);
	    my $accession = $1;

	    $accession =~ s/\-\d+$//;
	    $Accessions{$accession} = 1;

	    $has_accessions = 1;

	}
	close($AccessionGrep);

	next if (!$has_accessions);

	my $accessions_list_str = '';
	foreach my $accession (keys %Accessions) {
	    $accessions_list_str = $accessions_list_str.'/'.$accession;
	}
	$accessions_list_str =~ s/^\///;

	$SpeciesGeneToAccession{$species.'|'.$gene} = $accessions_list_str;
	
    }
    closedir($SpeciesAliDir);
    
}
closedir($AllSpeciesDir);

foreach my $species_gene_pair (sort keys %SpeciesGeneToAccession) {

    $species_gene_pair =~ /^([^\|]+)\|(\S+)$/;
    my $species = $1;
    my $gene = $2;

    my $accession = $SpeciesGeneToAccession{$species_gene_pair};

    while (length($species) < $longest_species_name_len) {
	$species = $species.' ';
    }

    while (length($gene) < $longest_gene_name_len) {
	$gene = $gene.' ';
    }

    print "$species $gene $accession\n";
    
}

1;













# EOF
