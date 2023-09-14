#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

sub ScriptDir { return './' if ($0 !~ /\//); $0 =~ /^(.*\/)[^\/]+$/; return $1; }
use lib ScriptDir();
use bureaucracy;


if (@ARGV != 2) {
    die "\n  USAGE: ./GetLengthAndStoppiness.pl [Extracted-DCEs] [Species-Guide]\n\n";
}


my $sfetch = '.'.(ScriptDir()).'../hsi/sfetch';
die "\n  Failed to locate sfetch (looking for '$sfetch')\n\n" if (!(-e $sfetch));


my %IsStopCodon;
$IsStopCodon{'TAG'} = 'amber';
$IsStopCodon{'TAA'} = 'ochre';
$IsStopCodon{'TGA'} = 'opal';


my $in_dir_name = ConfirmDirectory($ARGV[0]);

my $SpeciesGuide = OpenInputFile($ARGV[1]);
my %SpeciesToGenome;
while (my $line = <$SpeciesGuide>) {
    next if ($line !~ /^(\S+)\s+(\S+)\s*\S*\s*$/);
    $SpeciesToGenome{$1} = $2;
}
close($SpeciesGuide);

my @Species;
my $InDir = OpenDirectory($in_dir_name);
while (my $fname = readdir($InDir)) {

    next if ($fname !~ /^(\S+)\.csv$/);
    my $species = $1;

    next if (!$SpeciesToGenome{$species});

    push(@Species,$species);
    
}
closedir($InDir);


if (scalar(@Species) == 0) {
    die "\n  ERROR:  No recognized species found in directory '$in_dir_name'\n\n";
}


my $out_dir_name = CreateDirectory('Length-And-Stoppiness');


foreach my $species (@Species) {
    
    my $species_dir_name = ConfirmDirectory($in_dir_name.$species);
    my $genome = $SpeciesToGenome{$species};

    my $InCSV = OpenInputFile($in_dir_name.$species.'.csv');
    <$InCSV>; # Burn the header
    my %IndexGroupIsNearCTerm;
    while (my $line = <$InCSV>) {

	next if ($line !~ /Insufficient Nucls \(\S+\)\, \S+\s*$/);

	$line =~ /^(\d+)\, \S+\, (\S+)\,/;
	my $index = $1;
	my $group_list_str = $2;

	foreach my $group (split(/\//,$group_list_str)) {
	    $IndexGroupIsNearCTerm{$index.'|'.$group} = 1;
	}

    }
    close($InCSV);

    my $OutCSV = OpenOutputFile($out_dir_name.$species.'.csv');
    print $OutCSV "Index, Gene, Group, Nucl.s in Dual-Coding Region, Num. Remaining Nucl.s, Post-Coding Triple, Is a Stop Codon\n";

    my $SpeciesDir = OpenDirectory($species_dir_name);
    while (my $fname = readdir($SpeciesDir)) {

	next if ($fname !~ /^(\d+)\.(\S+)\.out$/);
	my $index = $1;
	my $gene  = $2;

	$fname = $species_dir_name.$fname;

	my $GeneFile = OpenInputFile($fname);

	my $chr = <$GeneFile>;
	while ($chr !~ /Chromosome\s+\:\s+(\S+)/) {
	    $chr = <$GeneFile>;
	}
	$chr = $1;

	my $revcomp = 0;
	if ($chr =~ /\[revcomp\]/) {
	    $chr =~ s/\[revcomp\]//;
	    $revcomp = 1;
	}
	
	while (my $line = <$GeneFile>) {

	    next if ($line !~ /\> Group (\d+)/);
	    my $group = $1;
	    
	    while ($line !~ /\|\s+([^\|]+)\s+\|\s+([^\|]+)\s*$/) {
		$line = <$GeneFile>;
	    }

	    my $dual_coding_coords_str = $1;
	    my $c_terminal_coords_str = $2;

	    $dual_coding_coords_str =~ s/\s//g;
	    $c_terminal_coords_str  =~ s/\s//g;

	    my $dual_coding_len = 0;
	    foreach my $dual_coding_run (split(/\//,$dual_coding_coords_str)) {

		$dual_coding_run =~ /^(\d+)\.\.(\d+)$/;
		my $exon_len = $2 - $1;

		$exon_len *= -1 if ($revcomp);
		$exon_len += 1;

		$dual_coding_len += $exon_len;
		
	    }

	    print $OutCSV "$index,$gene,$group,$dual_coding_len,";

	    # If we aren't approaching the end of the protein, jump ship!
	    if ($c_terminal_coords_str ne 'C-TERM' && !$IndexGroupIsNearCTerm{$index.'|'.$group}) {
		print $OutCSV "-,-,-\n"; 
		next;
	    }

	    # How far to the end of the protein (nucleotides)?
	    # What are the next 3 nucls after the end?
	    # Are the a stop codon?

	    my $dist_to_end   = 0;
	    my $exit_codon    = '';
	    my $is_stop_codon = 0;

	    my $final_coding_coord;
	    if ($c_terminal_coords_str ne 'C-TERM') {

		foreach my $c_terminal_run (split(/\//,$c_terminal_coords_str)) {
		    
		    $c_terminal_run =~ /^(\d+)\.\.(\d+)$/;
		    my $exon_len = $2 - $1;
		    
		    $exon_len *= -1 if ($revcomp);
		    $exon_len += 1;
		    
		    $dist_to_end += $exon_len;
		    
		}

		$c_terminal_coords_str =~ /(\d+)$/;
		$final_coding_coord = $1;
		
	    } else {

		$dual_coding_coords_str =~ /(\d+)$/;
		my $final_coding_coord = $1;
		
	    }

	    my $sfetch_range;
	    if ($revcomp) {
		$sfetch_range = ($final_coding_coord-1).'..'.($final_coding_coord-3);
	    } else {
		$sfetch_range = ($final_coding_coord+1).'..'.($final_coding_coord+3);
	    }
	    
	    my $SFetch = OpenSystemCommand($sfetch.' -range '.$sfetch_range.' '.$genome.' '.$chr);
	    $exit_codon = <$SFetch>;
	    $exit_codon = <$SFetch>;
	    close($SFetch);
	    
	    $exit_codon =~ s/\n|\r//g;
	    $exit_codon = uc($exit_codon);
	    
	    if ($IsStopCodon{$exit_codon}) {
		$is_stop_codon = 1;
	    }

	    print $OutCSV "$dist_to_end,$exit_codon,$is_stop_codon\n";
	    
	}
	close($GeneFile);

    }
    closedir($SpeciesDir);
    
}



1;













# EOF
