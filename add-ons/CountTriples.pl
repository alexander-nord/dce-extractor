#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) { die "\n  USAGE:  ./CountTriples.pl [Extracted-DCEs/[species].IREs-removed.csv]\n\n"; }


my $in_file_name = $ARGV[0];
die "\n  ERROR:  Failed to locate input CSV file '$in_file_name'\n\n"
	if (!(-e $in_file_name));
open(my $InFile,'<',$in_file_name)
	|| die "\n  ERROR:  Failed to open input CSV file '$in_file_name'\n\n";

<$InFile>; # Header

my %IdToFrames;
my %IdToGene;
while (my $line = <$InFile>)
{
	if ($line =~ /^(\d+),\s*([^,]+),[^,]+,\s*(\d),/)
	{
		my $id    = $1;
		my $gene  = $2;
		my $frame = $3;
		$IdToGene {$id} = $gene;
		if ($IdToFrames{$id}) { $IdToFrames{$id} = $IdToFrames{$id}.'/'.$frame; }
		else                  { $IdToFrames{$id} =                      $frame; }
	}
}

close($InFile);


my $num_triples = 0;

foreach my $id (sort {$a <=> $b} keys %IdToFrames)
{
	my %Frames;
	foreach my $frame (split(/\//,$IdToFrames{$id}))
	{
		$Frames{$frame} = 1;
	}
	if ($Frames{'1'} && $Frames{'2'} && $Frames{'3'})
	{
		my $gene = $IdToGene{$id};
		print " $id / $gene\n";
		$num_triples++;
	}
}

print "\n  Total triples: $num_triples\n\n";

1;
