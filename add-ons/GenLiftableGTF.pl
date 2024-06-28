#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub WriteLiftableEntry;
sub MergeOverlappingRanges;
sub ScanForOverlaps;
sub RangesOverlap;


if (@ARGV != 1) { die "\n  USAGE:  ./GenLiftableGTF.pl [Extracted-DCEs]\n\n"; }


my $top_dir_name = $ARGV[0];
die "\n  ERROR:  Failed to locate directory '$top_dir_name'\n\n"
	if (!(-d $top_dir_name));
$top_dir_name = $top_dir_name.'/' if ($top_dir_name !~ /\/$/);


opendir(my $TopDir, $top_dir_name)
	|| die "\n  ERROR:  Failed to open directory '$top_dir_name'\n\n";

my @Species;
while (my $file_name = readdir($TopDir))
{
	next if ($file_name !~ /^([^\.]+)\.csv$/);
	push(@Species,$1) if (-d $top_dir_name.$1);
}

closedir($TopDir);


die "\n  ERROR:  No valid species detected in directory '$top_dir_name'\n\n"
	if (scalar(@Species) == 0);


foreach my $species (sort @Species)
{
	my $species_dir_name = $top_dir_name.$species.'/';
	opendir(my $SpeciesDir, $species_dir_name)
		|| die "\n  ERROR:  Failed to open species directory '$species_dir_name'\n\n";

	my @DCEFiles;
	my $max_id = 0;

	while (my $file_name = readdir($SpeciesDir))
	{
		
		next if ($file_name !~ /^(\d+)\.[^\.]+\.out$/);
		my $dce_id = $1;

		$DCEFiles[$dce_id] = $species_dir_name.$file_name;
		$max_id = $dce_id if ($dce_id > $max_id);

	}

	closedir($SpeciesDir);

	next if ($max_id == 0);


	my $out_file_name = $top_dir_name.$species.'.liftable.gtf';
	open(my $OutFile,'>',$out_file_name)
		|| die "\n  ERROR:  Failed to open output file '$out_file_name'\n\n";

	for (my $dce_id = 0; $dce_id <= $max_id; $dce_id++)
	{
		next if (!$DCEFiles[$dce_id]);
		WriteLiftableEntry($DCEFiles[$dce_id],$OutFile);
	}

	close($OutFile);

}


1;








#################################################################################
#
#  Subroutine:  WriteLiftableEntry
#
sub WriteLiftableEntry
{

	my $in_file_name = shift;
	my $OutFile = shift;

	$in_file_name =~ /([^\/]+)$/;
	my $in_file_base_name = $1;

	$in_file_base_name =~ /^(\d+)\.([^\.]+)\.out$/;
	my $dce_id = $1;
	my $gene = $2;

	open(my $InFile,'<',$in_file_name)
		|| die "\n  ERROR:  Failed to open input file '$in_file_name'\n\n";

	my $chr = <$InFile>;
	$chr = <$InFile> while ($chr !~ /^Chromosome\s+:\s+(\S+)/);
	$chr = $1;

	my $revcomp = 0;
	my $strand = '+';
	if ($chr =~ /\[revcomp\]/)
	{
		$chr =~ s/\[revcomp\]//g;
		$revcomp = 1;
		$strand = '-';
	}


	my @AllRanges;
	while (my $line = <$InFile>)
	{
		if ($line =~ /Genome Ranges:/)
		{

			$line =~ /\| (\d+\.[^\|]*\.\d+) \|/;
			my $range_list_str = $1;

			foreach my $group_exon_range (split(/ \/ /,$range_list_str))
			{
				if ($revcomp)
				{
					$group_exon_range =~ /^(\d+)\.\.(\d+)$/;
					push(@AllRanges,$2.'..'.$1);
				}
				else 
				{
					push(@AllRanges,$group_exon_range);					
				}
			}

		}
	}

	close($InFile);


	my $final_ranges_ref = MergeOverlappingRanges(\@AllRanges);
	my @FinalRanges = @{$final_ranges_ref};


	foreach my $range (@FinalRanges)
	{
		$range =~ /^(\d+)\.\.(\d+)$/;
		my $start = $1;
		my $end   = $2;
		print $OutFile "$chr\tGenLiftableGTF\tCodingRegion\t$start\t$end\t.\t$strand\tgene_name \"$gene\";\tdce_id \"$dce_id\";\n";
	}


}







#################################################################################
#
#  Subroutine:  MergeOverlappingRanges
#
sub MergeOverlappingRanges
{
	my $range_list_ref = shift;
	my @RangeList = @{$range_list_ref};

	my @FinalList;

	my $active_ranges = scalar(@RangeList);
	while ($active_ranges)
	{

		my $range_id = 0;
		$range_id++ while (!$RangeList[$range_id]);

		$RangeList[$range_id] =~ /^(\d+)\.\.(\d+)$/;
		my $range_start = $1;
		my $range_end   = $2;

		$RangeList[$range_id] = 0;

		my $overlaps_found;
		($range_list_ref,$overlaps_found,$range_start,$range_end) 
			= ScanForOverlaps(\@RangeList,$range_start,$range_end);
		@RangeList = @{$range_list_ref};

		while ($overlaps_found)
		{
			($range_list_ref,$overlaps_found,$range_start,$range_end) 
				= ScanForOverlaps(\@RangeList,$range_start,$range_end);
			@RangeList = @{$range_list_ref};
		}


		push(@FinalList,$range_start.'..'.$range_end);


		$active_ranges = 0;
		for (my $i=0; $i<scalar(@RangeList); $i++)
		{
			if ($RangeList[$i])
			{
				$active_ranges = 1;
				last;
			}
		}

	}


	return \@FinalList;

}







#################################################################################
#
#  Subroutine:  ScanForOverlaps
#
sub ScanForOverlaps
{
	my $range_list_ref = shift;
	my $range_start = shift;
	my $range_end = shift;

	my @RangeList = @{$range_list_ref};

	my $overlaps_found = 0;

	for (my $i=0; $i<scalar(@RangeList); $i++)
	{
		next if (!$RangeList[$i]);

		$RangeList[$i] =~ /^(\d+)\.\.(\d+)$/;
		my $next_start = $1;
		my $next_end   = $2;

		if (RangesOverlap($range_start,$next_start,$range_end,$next_end))
		{
			$range_start = $next_start if ($next_start < $range_start);
			$range_end   = $next_end   if ($next_end   > $range_end  );
			$overlaps_found++;

			$RangeList[$i] = 0;
		}
	}

	return (\@RangeList,$overlaps_found,$range_start,$range_end);

}






#################################################################################
#
#  Subroutine:  RangesOverlap
#
sub RangesOverlap
{
	my $start1 = shift;
	my $start2 = shift;
	my $end1 = shift;
	my $end2 = shift;

	return 1 if ($start1 <= $start2 && $end1 >= $start2);
	return 1 if ($start1 <= $end2   && $end1 >= $end2  );
	
	return 1 if ($start2 <= $start1 && $end2 >= $start1);
	return 1 if ($start2 <= $end1   && $end2 >= $end1  );

	return 1 if ($start1 >= $start2 && $end1 <= $end2  );
	return 1 if ($start2 >= $start1 && $end2 <= $end1  );

	return 0;
}


