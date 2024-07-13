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
while (my $line = <$InFile>)
{
	if ($line =~ /^(\d+),[^,]+,[^,]+,\s*(\d),/)
	{
		my $id    = $1;
		my $frame = $2;
		if ($IdToFrames{$id}) { $IdToFrames{$id} = $IdToFrames{$id}.'/'.$frame; }
		else                  { $IdToFrames{$id} =                      $frame; }
	}
}

close($InFile);


foreach my $id (sort {$a <=> $b} keys %IdToFrames)
{
	my %Frames;
	foreach my $frame (split(/\//,$IdToFrames{$id}))
	{
		$Frames{$frame} = 1;
	}
	if ($Frames{'1'} && $Frames{'2'} && $Frames{'3'})
	{
		print " $id\n";
	}
}

1;
