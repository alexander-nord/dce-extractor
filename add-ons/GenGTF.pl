#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

sub ScriptDir { return './' if ($0 !~ /\//); $0 =~ /^(.*\/)[^\/]+$/; return $1; }
use lib ScriptDir();
use bureaucracy;


if (@ARGV != 1) {
    die "\n  USAGE: ./GenGTFs.pl [Extracted-DCEs]\n\n";
}


my $base_dir_name = ConfirmDirectory($ARGV[0]);
my $BaseDir = OpenDirectory($base_dir_name);

while (my $csv_name = readdir($BaseDir)) {

    next if ($csv_name !~ /^(\S+)\.csv$/);

    my $species = $1;
    my $species_dir_name = $base_dir_name.$species.'/';
    
    next if (!(-d $species_dir_name));

    my $gtf_name = $base_dir_name.$species.'.gtf';
    RunSystemCommand("rm \"$gtf_name\"") if (-e $gtf_name);

    my $GTF = OpenOutputFile($gtf_name);

    my $SpeciesDir = OpenDirectory($species_dir_name);
    while (my $dce_file_name = readdir($SpeciesDir)) {

	next if ($dce_file_name !~ /^\d+\.\S+\.out$/);

	$dce_file_name = $species_dir_name.$dce_file_name;

	ExtractDCEsToGTF($dce_file_name,$GTF);

    }
    closedir($SpeciesDir);
    close($GTF);
    
}
closedir($BaseDir);



1;











######################################################################
#
#  Function:  ExtractDCEsToGTF
#
sub ExtractDCEsToGTF
{

    my $dce_file_name = shift;
    my $GTF = shift;

    my $DCEFile = OpenInputFile($dce_file_name);

    my $species_line    = <$DCEFile>;
    my $gene_fam_line   = <$DCEFile>;
    my $dce_index_line  = <$DCEFile>;
    my $chromosome_line = <$DCEFile>;

    $gene_fam_line =~ /Gene Family \: (\S+)/;
    my $gene_fam = $1;

    $dce_index_line =~ /DCE Index   \: (\d+)/;
    my $dce_index = $1;

    $chromosome_line =~ /Chromosome  \: (\S+)/;
    my $chr = $1;

    my $strand = '+';
    my $revcomp = 0;
    if ($chr =~ /\[revcomp\]/) {
	$strand = '-';
	$revcomp = 1;
	$chr =~ s/\[revcomp\]//;
    }

    
    my %RangesToGroups;
    
    my $line = <$DCEFile>;
    while (!eof($DCEFile) && $line !~ /^\-\-\-\-\-\-\-\-\-/) {

	if ($line !~ /^\> Group (\d+)\:/) {
	    $line = <$DCEFile>;
	    next;
	}

	my $group_id = $1;

	while ($line !~ /Genome Ranges/) {
	    $line = <$DCEFile>;
	}

	$line =~ /\|([^\|]+)\|/;

	my $ranges = $1;
	$ranges =~ s/\s//g;

	if ($RangesToGroups{$ranges}) {
	    $RangesToGroups{$ranges} = $RangesToGroups{$ranges}.'/'.$group_id;
	} else {
	    $RangesToGroups{$ranges} = $group_id;
	}

    }
    close($DCEFile);

    foreach my $range_list_str (keys %RangesToGroups) {

	my $group_str = $RangesToGroups{$range_list_str};

	my $exon_num = 0;
	foreach my $range (split(/\//,$range_list_str)) {

	    $exon_num++;
	    
	    $range =~ /^(\d+)\.\.(\d+)$/;
	    my $start = $1;
	    my $end = $2;

	    if ($revcomp) {
		my $swap = $end;
		$end = $start;
		$start = $swap;
	    }

	    print $GTF "$chr";
	    print $GTF "\tExtractedDCEs\texon";
	    print $GTF "\t$start\t$end\t.\t$strand\t.";
	    print $GTF "\tdce_index \"$dce_index\";";
	    print $GTF "\tgene_name \"$gene_fam\";";
	    print $GTF "\tgroup \"$group_str\";";
	    print $GTF "\texon_number \"$exon_num\";";
	    print $GTF "\n";
	    
	}
	
    }
    
}









# EOF
