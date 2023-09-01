#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub ScriptDir { return './' if ($0 !~ /\//); $0 =~ /^(.*\/)[^\/]+$/; return $1; }
use lib ScriptDir();
use bureaucracy;


my $sfetch = ScriptDir().'hsi/sfetch';

if (!(-e $sfetch)) {
    die "\n\n  >>>  Run 'lazy_maker.sh' in 'hsi/' (please and thank you!)  <<<\n\n\n";
}


sub HasDCEs;
sub ExtractDCEs;
sub RemoveEmptyColumns;
sub RecordWindowsToCSV;
sub VisDualCodingRegion;
sub CheckForMidExonIntrons;
sub GenStrongWindowFastas;
sub GenMidExonIntronFastas;

# DEBUGGING OUTPUT
sub DumpSpeciesGeneData;



if (@ARGV != 2) {
    die "\n  USAGE: ./DCE-Extractor.pl [Mirage-Results] [Species-Guide]\n\n";
}


my $mirage_results_dirname = ConfirmDirectory($ARGV[0]);

my $SpeciesGuide = OpenInputFile($ARGV[1]);
my %SpeciesToGenomes;

my $num_mapped_species = 0;

while (my $line = <$SpeciesGuide>) {
    if ($line =~ /^\s*(\S+)\s+(\S+)/) {

	my $species = lc($1);
	my $genome  = $2;

	if (-e $genome) {
	    $SpeciesToGenomes{$species} = $genome;
	    $num_mapped_species++;
	}
	
    }
}

close($SpeciesGuide);


if ($num_mapped_species == 0) {
    die "\n  ERROR:  Failed to locate the genomes listed in species guide file '$ARGV[1]'\n\n";
}


my $out_dirname = CreateDirectory('Extracted-DCEs');

my %GappyMappingSeqs; # Global hash for sequences that have gappy mappings (from Spaln)

my $num_dces = 0;

my $all_species_dirname = ConfirmDirectory($mirage_results_dirname.'Species-MSAs');

my $mei_fname = $out_dirname.'Mid-Exon-Intron-Warnings.out';
my $MEIFile = OpenOutputFile($mei_fname);

my $AllSpeciesDir = OpenDirectory($all_species_dirname);
while (my $species = readdir($AllSpeciesDir)) {

    
    $species =~ s/\/$//;
    
    next if (!$SpeciesToGenomes{$species});

    
    # For searching, we'll take our great big (nasty) outputs and extract
    # the critical (16-nucl|16-nucl) windows on the left and right sides
    # and write them to a CSV file
    my $species_csv_name = $out_dirname.$species.'.csv';
    my $SpeciesCSV = OpenOutputFile($species_csv_name);
    print $SpeciesCSV "Index, Gene, Group(s), Frame(s), Left Window, Right Window, Trans Ali Pct ID\n";

    
    # Hold onto the starting value of num_dces to determine if we recorded anything...
    my $species_dce_start = $num_dces;


    # Do some data prep
    my $genome = $SpeciesToGenomes{$species};

    my $species_in_dirname  = ConfirmDirectory($all_species_dirname.$species);
    my $species_ali_dirname = ConfirmDirectory($species_in_dirname.'alignments');
    my $species_map_dirname = ConfirmDirectory($species_in_dirname.'mappings');
    
    my $species_out_dirname = CreateDirectory($out_dirname.$species);


    my $SpeciesAliDir = OpenDirectory($species_ali_dirname);
    while (my $gene = readdir($SpeciesAliDir)) {

	next if ($gene !~ /^(\S+)\.afa$/);
	$gene = $1;

	my $gene_ali_fname = $species_ali_dirname.$gene.'.afa';
	my $gene_map_fname = $species_map_dirname.$gene.'.out';
	
	# Yee haw?
	next if (!HasDCEs($gene_ali_fname));

	# YEE HAW!!!
	my $gene_start_dce = $num_dces;
	$num_dces = ExtractDCEs($species_out_dirname,$gene_ali_fname,
				$gene_map_fname,$genome,$num_dces);

	# Record this gene's DCE(s) to our CSV
	while ($gene_start_dce < $num_dces) {

	    $gene_start_dce++;
	    my $dce_fname = $species_out_dirname.$gene_start_dce.'.'.$gene.'.out';

	    VisDualCodingRegion($dce_fname);
	    CheckForMidExonIntrons($MEIFile,$dce_fname); # This also records frame num
	    RecordWindowsToCSV($SpeciesCSV,$dce_fname);

	}
	
    }
    closedir($SpeciesAliDir);


    # If we didn't have any DCEs added for this species, clear its output
    # directory and move along
    if ($num_dces - $species_dce_start == 0) {
	RunSystemCommand("rm -rf \"$species_out_dirname\"");
	RunSystemCommand("rm \"$species_csv_name\"");
	next;
    }


}
closedir($AllSpeciesDir);

close($MEIFile);
if (!(-s $mei_fname)) { RunSystemCommand("rm $mei_fname"); }


# If we had any gappy mappings (indicating low-quality), report them
if (scalar(keys %GappyMappingSeqs)) {

    # We'll want to order them by their DCE indices
    my %GappySeqsByIndex;
    my $max_gappy_index = 0;
    foreach my $gms_key (keys %GappyMappingSeqs) {

	my @GMSData = split(/\&/,$gms_key);

	my $seq_name  = $GMSData[0];
	my $dce_index = $GMSData[1];

	if ($GappySeqsByIndex{$dce_index}) {
	    $GappySeqsByIndex{$dce_index} = $GappySeqsByIndex{$dce_index}.'&'.$seq_name;
	} else {
	    $GappySeqsByIndex{$dce_index} = $seq_name;
	}

    }

    my $GappyOutFile = OpenOutputFile($out_dirname.'Gappy-Alignment-Warnings.out');
    foreach my $dce_index (sort { $a <=> $b } keys %GappySeqsByIndex) {

	my $seq_name = $GappySeqsByIndex{$dce_index};
	$seq_name =~ s/\&.+$//;
	$seq_name =~ /^([^\|]+)\|([^\|]+)\|/;

	my $species = lc($1);
	my $gene_fam_list = lc($2);

	$gene_fam_list =~ /^([^\/]+)/;
	my $gene_fam = $1;
	
	print $GappyOutFile "$species: $dce_index ($gene_fam)\n";
	    
    }
    close($GappyOutFile);
    
}


# Now, let's make some FASTA files!
GenStrongWindowFastas($out_dirname);
GenMidExonIntronFastas($out_dirname);

1;





#############################################################################
#
#  Function:  HasDCEs
#
sub HasDCEs
{
    my $fname = shift;

    my $InFile = OpenInputFile($fname);

    my $has_dces = 0;
    while (my $line = <$InFile>) {

	next if ($line !~ /^\>/);

	if ($line =~ / ARFs?\:(\d+)/) {
	    $has_dces = 1;
	    last;
	}
	
    }
    close($InFile);

    return $has_dces;
    
}







#############################################################################
#
#  Function:  ExtractDCEs
#
sub ExtractDCEs
{
    my $out_dirname = shift;
    my $ali_fname   = shift;
    my $map_fname   = shift;
    my $genome      = shift;
    my $dce_index   = shift;

    $out_dirname =~ /\/([^\/]+)\/$/;
    my $species = $1;

    $ali_fname =~ /\/([^\/]+)\.afa$/;
    my $gene = $1;

    my $MapFile = OpenInputFile($map_fname);

    my %MappedSeqs;
    my $num_seqs = 0;
    
    my $canon_chr = <$MapFile>;
    $canon_chr =~ /Canonical Chromosome\:\s+(\S+)/;
    $canon_chr = $1;

    while (my $line = <$MapFile>) {

	next if ($line !~ /Sequence ID\: (\S+)/);
	my $seq_name = $1;

	$line = <$MapFile>;

	$line = <$MapFile>;
	$line =~ /Chromosome \: (\S+)/;

	my $chr = $1;

	if ($chr eq $canon_chr) {
	    $num_seqs++;
	    $MappedSeqs{$seq_name} = $num_seqs;
	}

    }
    close($MapFile);
    
    my $revcomp = 0;
    if ($canon_chr =~ /\[revcomp\]/) {
	$revcomp = 1;
	$canon_chr =~ s/\[revcomp\]//;
    }

    my $AliFile = OpenInputFile($ali_fname);

    my @MSA; # NOTE: MSA[0] is empty!
    my @SeqNames;
    my @ARFLabels;
    my $msa_len;
    my $seq_num;
    my $is_mapped = 0;

    while (my $line = <$AliFile>) {

	$line =~ s/\n|\r//g;
	next if (!$line);

	if ($line =~ /^\>(\S+)/) {

	    my $seq_name = $1;

	    if ($MappedSeqs{$seq_name}) {

		$seq_num = $MappedSeqs{$seq_name};

		my $arf_label = 'none';
		if ($line =~ / ARFs?\:(\S+)/) {
		    $arf_label = $1;
		}

		$SeqNames[$seq_num]  = $seq_name;
		$ARFLabels[$seq_num] = $arf_label;
	    
		$msa_len = 0;

		$is_mapped = 1;
		
	    } else {

		$is_mapped = 0;

	    }
	    
	} elsif ($is_mapped) {

	    foreach my $char (split(//,uc($line))) {
		$MSA[$seq_num][$msa_len] = $char;
		$msa_len++;
	    }
	    
	}
	
    }
    close($AliFile);

    my $msa_ref;
    my $num_exons;
    ($msa_ref,$msa_len,$num_exons) = RemoveEmptyColumns(\@MSA,$num_seqs,$msa_len);

    @MSA = @{$msa_ref};


    # At this point, all empty (all-gap) columns have been removed,
    # and each position of MSA[0][0..msa_len-1] is either:
    #
    #    (a.) a splice-site indicator ('/')
    #    (b.) the numerical index of the exon (from 1 to num_exons)
    #
    # Yay!

    
    $MapFile = OpenInputFile($map_fname);

    my @MapMSA;
    my @MapRangesMSA; # Not super data-efficient, but WHAT*EVER*

    while (my $line = <$MapFile>) {

	next if ($line !~ /Sequence ID\: (\S+)/);
	my $seq_name = $1;

	next if (!$MappedSeqs{$seq_name}); # Non-canonical exon usage check
	
	$seq_num = $MappedSeqs{$seq_name};

	$line = <$MapFile>; # Map method
	$line = <$MapFile>; # Chromosome

	$line = <$MapFile>;
	$line =~ /Num Exons\s+\: (\d+)/;
	my $num_seq_exons = $1;

	my $msa_pos = 0;
	
	for (my $exon_num = 0; $exon_num < $num_seq_exons; $exon_num++) {

	    $line = <$MapFile>;

	    $line =~ /\:(\d+\.\.\d+)\s*$/;
	    my $exon_map_range = $1;

	    $line = <$MapFile>;
	    $line =~ s/\n|\r//g;
	    my @ExonMapCoords = split(/\,/,$line);

	    my $exon_map_coord_pos = 0;
	    while ($exon_map_coord_pos < scalar(@ExonMapCoords)) {

		if ($MSA[$seq_num][$msa_pos] =~ /[A-Z]/) {

		    $MapMSA[$seq_num][$msa_pos] = $ExonMapCoords[$exon_map_coord_pos];
		    $MapRangesMSA[$seq_num][$msa_pos] = $exon_map_range;

		    $exon_map_coord_pos++;
		    
		} else {

		    $MapMSA[$seq_num][$msa_pos]       = $MSA[$seq_num][$msa_pos];
		    $MapRangesMSA[$seq_num][$msa_pos] = $MSA[$seq_num][$msa_pos];

		}

		$msa_pos++;
		
	    }
	    
	}
	
    }

    close($MapFile);


    # The very last bit of prep work we're going to do is make an MSA indicating
    # whether an amino acid is recognized as Alternative or Standard
    my @AminoStatusMSA;
    for (my $i=1; $i<=$num_seqs; $i++) {

	my $msa_pos = 0;

	my $arf_label = $ARFLabels[$i];

	if ($arf_label ne 'none') {

	    my @ARFRanges = split(/\,/,$arf_label);

	    my $amino_pos = 0;
	    foreach my $arf_range (@ARFRanges) {

		$arf_range =~ /(\d+)\.\.(\d+)/;
		my $start_amino = $1-1;
		my $end_amino   = $2; # This let's us just do strict '<'
		
		while ($amino_pos < $start_amino) {

		    if ($MSA[$i][$msa_pos] =~ /[A-Z]/) {
			$AminoStatusMSA[$i][$msa_pos] = 'S';
			$amino_pos++;
		    } else {
			$AminoStatusMSA[$i][$msa_pos] = $MSA[$i][$msa_pos];
		    }
		    
		    $msa_pos++;
		    
		}
		
		while ($amino_pos < $end_amino) {

		    if ($MSA[$i][$msa_pos] =~ /[A-Z]/) {
			$AminoStatusMSA[$i][$msa_pos] = 'A';
			$amino_pos++;
		    } else {
			$AminoStatusMSA[$i][$msa_pos] = $MSA[$i][$msa_pos];
		    }

		    $msa_pos++;
		    
		}
	    
	    }

	}
	
	while ($msa_pos < $msa_len) {

	    if ($MSA[$i][$msa_pos] =~ /[A-Z]/) {
		$AminoStatusMSA[$i][$msa_pos] = 'S';
	    } else {
		$AminoStatusMSA[$i][$msa_pos] = $MSA[$i][$msa_pos];
	    }

	    $msa_pos++;

	}

    }


    # Now we can get down to the real business!
    #
    # For the sake of clarity, here are the data we've compiled so far:
    #
    #    + MSA            : The amino acid matrix, with MSA[0][0..msa_len-1]
    #                         indicating which exon number a give column corresponds
    #                         to (in the range from 1..num_exons).
    #
    #    + SeqNames       : The names of the sequences, without comment data
    #
    #    + ARFLabels      : What ARF coordinates are associated with this sequence?
    #
    #    + MapMSA         : What center-nucleotide mapping coordinates are associated
    #                         with each amino acid?
    #
    #    + MapRangesMSA   : What range of precise nucleotide mapping coordinates are
    #                         associated with the exon that each amino acid belongs to?
    #
    #    + AminoStatusMSA : Is each amino acid 'S'tandard or 'A'lternative?
    #


    # DATA DUMP CALL (DEBUGGING)
    #
    #DumpSpeciesGeneData($out_dirname.$gene.'.out',$canon_chr,$revcomp,\@MSA,\@SeqNames,
    #$num_seqs,$msa_len,\@MapMSA,\@MapRangesMSA,\@AminoStatusMSA);


    # In case multiple sequences have the same 'A'-labelled content (and are
    # effectively identical), we'll try not to double-count.
    my %ObservedAltContent;
    
    # We'll go through our MSA sequence-by-sequence looking for
    # 'A's in our AminoStatusMSA
    for (my $seq_id=1; $seq_id<$num_seqs; $seq_id++) {

	my $msa_pos=1;
	while ($msa_pos < $msa_len-1) {

	    if ($AminoStatusMSA[$seq_id][$msa_pos] ne 'A') {
		$msa_pos++;
		next;
	    }

	    my $alt_start_col = $msa_pos;
	    my $alt_end_col = $msa_pos;

	    my $num_alt_aminos = 1;

	    $msa_pos++;
	    while ($msa_pos < $msa_len-1 && $AminoStatusMSA[$seq_id][$msa_pos] ne 'S') {

		if ($AminoStatusMSA[$seq_id][$msa_pos] eq 'A') {
		    $alt_end_col = $msa_pos;
		    $num_alt_aminos++;
		}
		$msa_pos++;
		
	    }

	    next if ($num_alt_aminos < 7);

	    my $alt_start_exon = $MSA[0][$alt_start_col];
	    my $alt_end_exon   = $MSA[0][$alt_end_col];

	    my $alt_msa_range_start = $alt_start_col;
	    while ($MSA[0][$alt_msa_range_start-1] ne '/') {
		$alt_msa_range_start--;
	    }

	    my $alt_msa_range_end = $alt_end_col;
	    while ($MSA[0][$alt_msa_range_end+1] ne '/') {
		$alt_msa_range_end++;
	    }

	    # We'll make a list of all sequences that have some representation within
	    # the range of interest, as well as a list of MSA coordinates where they
	    # have aminos and what those aminos are.
	    my @ReppedSeqIDs;          # The sequence ID wrt the full MSA
	    my $focal_rep_id;          # The representative ID of the ARF-labelled seq
	    my @ReppedCorePosListStrs; # The (A/S)RF of interest's columns in the MSA
	    my @ReppedSeqStrs;         # The (A/S)RF of interest's sequence
	    for (my $repped_seq_id = 1; $repped_seq_id <= $num_seqs; $repped_seq_id++) {

		my $repped_seq_pos_list_str = '';
		my $repped_seq_str = '';
		my $num_aminos_in_range = 0;

		for (my $repped_seq_pos = $alt_msa_range_start;
		     $repped_seq_pos <= $alt_msa_range_end;
		     $repped_seq_pos++)
		{

		    if ($MSA[$repped_seq_id][$repped_seq_pos] =~ /[A-Z]/) {

			$repped_seq_pos_list_str
			    = $repped_seq_pos_list_str.','.$repped_seq_pos;

			$repped_seq_str
			    = $repped_seq_str.$MSA[$repped_seq_id][$repped_seq_pos];

			$num_aminos_in_range++;
			
		    } elsif ($MSA[$repped_seq_id][$repped_seq_pos] eq '/') {

			$repped_seq_str = $repped_seq_str.'/';
			
		    }
		    
		}

		next if ($num_aminos_in_range < 7);


		$repped_seq_pos_list_str =~ s/^\,//;

		
		if ($repped_seq_id == $seq_id) {
		    $focal_rep_id = scalar(@ReppedSeqIDs);
		}
		
		push(@ReppedSeqIDs,$repped_seq_id);
		push(@ReppedSeqStrs,$repped_seq_str);
		push(@ReppedCorePosListStrs,$repped_seq_pos_list_str);

		
	    }


	    # It's conceivable that we won't have any competing frames that had
	    # six aminos (unlikely, but conceivable), so if we don't have multiple
	    # sequences to consider here we'll jump ship.
	    my $num_repped = scalar(@ReppedSeqIDs);
	    next if ($num_repped < 2);

	    
	    # Next, for each sequence that had representation in the zone with our
	    # alternative frame of interest, we'll want to get a list of the MSA
	    # positons where they have their next seven aminos
	    my @ReppedLeftPosListStrs;  # The upstream content (up to 7 aminos)
	    my @ReppedRightPosListStrs; # The downstream content (up to 7 aminos)
	    for (my $rep_id=0; $rep_id<$num_repped; $rep_id++) {

		my $repped_seq_id = $ReppedSeqIDs[$rep_id];

		# LEFT SIDE
		$ReppedCorePosListStrs[$rep_id] =~ /^(\d+)\,/;
		my $repped_seq_pos = $1 - 1;

		my $remaining_positions  = 7;
		my $repped_left_list_str = '';
		my $repped_left_seq_str  = '';
		
		while ($repped_seq_pos && $remaining_positions) {

		    if ($MSA[$repped_seq_id][$repped_seq_pos] =~ /[A-Z]/) {

			$repped_left_list_str
			    = $repped_seq_pos.','.$repped_left_list_str;

			$repped_left_seq_str
			    = $MSA[$repped_seq_id][$repped_seq_pos].$repped_left_seq_str;
			
			$remaining_positions--;
			
		    } elsif ($MSA[$repped_seq_id][$repped_seq_pos] eq '/') {

			$repped_left_seq_str = '/'.$repped_left_seq_str;
			
		    }

		    $repped_seq_pos--;
		    
		}

		# If the alt frame was N-terminal, we might not have anything here...
		if ($repped_left_list_str) { $repped_left_list_str =~ s/\,$//; }

		$ReppedLeftPosListStrs[$rep_id] = $repped_left_list_str;


		# RIGHT SIDE
		$ReppedCorePosListStrs[$rep_id] =~ /\,(\d+)$/;
		$repped_seq_pos = $1 + 1;
		
		$remaining_positions  = 7;
		my $repped_right_list_str = '';
		my $repped_right_seq_str  = '';
		
		while ($repped_seq_pos < $msa_len && $remaining_positions) {

		    if ($MSA[$repped_seq_id][$repped_seq_pos] =~ /[A-Z]/) {

			$repped_right_list_str
			    = $repped_right_list_str.','.$repped_seq_pos;

			$repped_right_seq_str
			    = $repped_right_seq_str.$MSA[$repped_seq_id][$repped_seq_pos];
			
			$remaining_positions--;
			
		    } elsif ($MSA[$repped_seq_id][$repped_seq_pos] eq '/') {

			$repped_right_seq_str = $repped_right_seq_str.'/';
			
		    }

		    $repped_seq_pos++;
		    
		}

		# If the alt frame was C-terminal, we might not have anything here...
		if ($repped_right_list_str) { $repped_right_list_str =~ s/^\,//; }

		$ReppedRightPosListStrs[$rep_id] = $repped_right_list_str;


		# Before we move onto the next represented sequence, let's update the
		# full amino sequence info!
		while ($repped_left_seq_str =~ /\/\//) {
		    $repped_left_seq_str =~ s/\/\//\//;
		}
		if ($repped_left_seq_str eq '/') {
		    $repped_left_seq_str = '';
		}

		while ($repped_right_seq_str =~ /\/\//) {
		    $repped_right_seq_str =~ s/\/\//\//;
		}
		if ($repped_right_seq_str eq '/') {
		    $repped_right_seq_str = '';
		}

		$ReppedSeqStrs[$rep_id]
		    = $repped_left_seq_str.'|'.$ReppedSeqStrs[$rep_id].'|'.$repped_right_seq_str;

		

	    }

	    
	    # OH, BOI, WE R GETTIN' CLOSE!

	    # The next step is to make a list of the genome ranges for the implicated
	    # exons for each sequence.
	    
	    my @ReppedGenomeRanges;
	    for (my $rep_id=0; $rep_id<$num_repped; $rep_id++) {
		
		my $repped_seq_id = $ReppedSeqIDs[$rep_id];

		# LEFT
		my $left_ranges_str = '';
		my $current_range   = '';
		foreach my $pos (split(/\,/,$ReppedLeftPosListStrs[$rep_id])) {

		    if ($MapRangesMSA[$repped_seq_id][$pos] ne $current_range) {
			$current_range   = $MapRangesMSA[$repped_seq_id][$pos];
			$left_ranges_str = $left_ranges_str.'/'.$current_range;
		    }
		    
		}

		if (!$left_ranges_str) { $left_ranges_str = 'N-TERM'; }
		else                   { $left_ranges_str =~ s/^\///; }

		
		# OVERLAP WITH ARF
		my $mid_ranges_str = '';
		$current_range     = '';
		foreach my $pos (split(/\,/,$ReppedCorePosListStrs[$rep_id])) {

		    if ($MapRangesMSA[$repped_seq_id][$pos] ne $current_range) {
			$current_range  = $MapRangesMSA[$repped_seq_id][$pos];
			$mid_ranges_str = $mid_ranges_str.'/'.$current_range;
		    }
		    
		}

		$mid_ranges_str =~ s/^\///;

		
		# RIGHT
		my $right_ranges_str = '';
		$current_range    = '';
		foreach my $pos (split(/\,/,$ReppedRightPosListStrs[$rep_id])) {

		    if ($MapRangesMSA[$repped_seq_id][$pos] ne $current_range) {
			$current_range    = $MapRangesMSA[$repped_seq_id][$pos];
			$right_ranges_str = $right_ranges_str.'/'.$current_range;
		    }
		    
		}

		if (!$right_ranges_str) { $right_ranges_str = 'C-TERM'; }
		else                    { $right_ranges_str =~ s/^\///; }

		
		# Put it all together
		my $genome_ranges_str
		    = $left_ranges_str.'|'.$mid_ranges_str.'|'.$right_ranges_str;
		push(@ReppedGenomeRanges,$genome_ranges_str);

	    }


	    # Quick check -- have we seen this 'alt' content before?
	    my $double_count_check = $ReppedGenomeRanges[$focal_rep_id].'&'.$ReppedSeqStrs[$focal_rep_id];
	    next if ($ObservedAltContent{$double_count_check});
	    
	    # It's a new one!
	    $ObservedAltContent{$double_count_check} = 1;
	    

	    # Since we have the amino acid and genome range data for each of the
	    # repped sequences, we're in a position to say which of these are
	    # using the exact same content (with respect to the DCE of interest),
	    # so let's link up anything that would be redundant.
	    #
	    # UPDATE: An obvious issue is that focal_rep_id==0 makes all sequences
	    #         report as 'identical' -- for a quick fix I'm going to add a
	    #         second flag 'IsUniqueSeq' that's either 1 or 0
	    #
	    my @IsGroupLeader;
	    my @IsIdenticalTo;

	    # We'll do a preliminary check for sequences that are identical to the
	    # sequence of interest
	    for (my $rep_id=0; $rep_id<$num_repped; $rep_id++) {

		$IsGroupLeader[$rep_id] = 1;
		$IsIdenticalTo[$rep_id] = 0;

		next if ($rep_id == $focal_rep_id);
		
		if ($ReppedGenomeRanges[$rep_id] eq $ReppedGenomeRanges[$focal_rep_id]
		    && $ReppedSeqStrs[$rep_id] eq $ReppedSeqStrs[$focal_rep_id]) {

		    $IsGroupLeader[$rep_id] = 0;
		    $IsIdenticalTo[$rep_id] = $focal_rep_id;
		    
		}
		
	    }

	    # Next, we'll consider other identities...
	    for (my $rep_id=0; $rep_id<$num_repped; $rep_id++) {

		if ($ReppedSeqIDs[$rep_id] == $seq_id || $IsIdenticalTo[$rep_id]) {
		    next;
		}
		
		for (my $comp_rep_id=0; $comp_rep_id < $rep_id; $comp_rep_id++) {

		    if ($ReppedGenomeRanges[$rep_id] eq $ReppedGenomeRanges[$comp_rep_id]
			&& $ReppedSeqStrs[$rep_id] eq $ReppedSeqStrs[$comp_rep_id]) {

			$IsGroupLeader[$rep_id] = 0;
			$IsIdenticalTo[$rep_id] = $comp_rep_id;
			last;

		    }
		    
		}
		
	    }

	    
	    
	    # Awesome! Now it's time to actually pull all that stinky data together!
	    my @OutAminoSeqs;
	    my @OutNuclSeqs;
	    my @OutTransSeqs;
	    for (my $rep_id=0; $rep_id<$num_repped; $rep_id++) {

		next if ($IsIdenticalTo[$rep_id]);

		my $repped_seq_id  = $ReppedSeqIDs[$rep_id];
		my $repped_seq_str = $ReppedSeqStrs[$rep_id];
		$repped_seq_str =~ s/\/|\|//g;

		my $repped_seq_pos_list_str   = $ReppedCorePosListStrs[$rep_id];
		my $repped_left_pos_list_str  = $ReppedLeftPosListStrs[$rep_id];
		my $repped_right_pos_list_str = $ReppedRightPosListStrs[$rep_id];

		my @CorePositions  = split(/\,/,$ReppedCorePosListStrs[$rep_id]);
		my @LeftPositions  = split(/\,/,$ReppedLeftPosListStrs[$rep_id]);
		my @RightPositions = split(/\,/,$ReppedRightPosListStrs[$rep_id]);

		my $repped_genome_ranges_str = $ReppedGenomeRanges[$rep_id];

		$ReppedGenomeRanges[$rep_id] =~ /^([^\|]+)\|([^\|]+)\|([^\|]+)$/;
		my  $left_ranges_str = $1;
		my  $core_ranges_str = $2;
		my $right_ranges_str = $3;

		my $full_nucl_str = '';

		
		#
		# UPSTREAM!
		#
		if ($left_ranges_str ne 'N-TERM') {

		    # If we aren't N-terminal, we'll want to be sure that we have the
		    # true first nucleotide for our first codon
		    $left_ranges_str =~ /^(\d+)\./;
		    my $first_range_start = $1;
		    
		    my $first_codon_center = $MapMSA[$repped_seq_id][$LeftPositions[0]];
		    
		    if ($first_codon_center == $first_range_start) {

			# Backtrack!
			my $pos = $LeftPositions[0]-1;
			while ($MSA[$repped_seq_id][$pos] !~ /[A-Z]/) {
			    $pos--;
			}

			$MapRangesMSA[$repped_seq_id][$pos] =~ /\.\.(\d+)$/;
			my $first_nucl_range = $1.'..'.$1;

			$left_ranges_str = $first_nucl_range.'/'.$left_ranges_str;
			
		    } else {

			my $first_codon_start = $first_codon_center;		    

			if ($revcomp) { $first_codon_start++; }
			else          { $first_codon_start--; }

			$left_ranges_str =~ s/^\d+\./$first_codon_start\./;
			
		    }

		    # Now we can party!
		    foreach my $range (split(/\//,$left_ranges_str)) {

			my $sfetch_cmd = $sfetch.' -range '.$range.' '.$genome.' '.$canon_chr;
			my $Fetched = OpenSystemCommand($sfetch_cmd);

			my $line = <$Fetched>; # Eat the header
			while ($line = <$Fetched>) {

			    $line =~ s/\n|\r//g;
			    next if (!$line);

			    $full_nucl_str = $full_nucl_str.lc($line);

			}

			close($Fetched);

			$full_nucl_str = $full_nucl_str.'/';
			
		    }

		    $full_nucl_str =~ s/\/$//;
		    
		}

		
		$full_nucl_str = $full_nucl_str.'|';

		
		#
		# THE CORE!
		#
		foreach my $range (split(/\//,$core_ranges_str)) {
		    
		    my $sfetch_cmd = $sfetch.' -range '.$range.' '.$genome.' '.$canon_chr;
		    my $Fetched = OpenSystemCommand($sfetch_cmd);
		    
		    my $line = <$Fetched>; # Eat the header
		    while ($line = <$Fetched>) {
			
			$line =~ s/\n|\r//g;
			next if (!$line);
			
			$full_nucl_str = $full_nucl_str.uc($line);
			
		    }
		    
		    close($Fetched);
		    
		    $full_nucl_str = $full_nucl_str.'/';
		    
		}
		
		
		$full_nucl_str =~ s/\/$/\|/;

		
		#
		# DOWNSTREAM!
		#
		if ($right_ranges_str ne 'C-TERM') {

		    # If we aren't C-terminal, we'll want to be sure that we have the
		    # true last nucleotide for our first codon
		    $right_ranges_str =~ /\.(\d+)$/;
		    my $last_range_end = $1;

		    my $final_pos = $RightPositions[scalar(@RightPositions)-1];
		    my $last_codon_center = $MapMSA[$repped_seq_id][$final_pos];
		    
		    if ($last_codon_center == $last_range_end) {

			# Forwardtrack!
			my $pos = $final_pos+1;
			while ($MSA[$repped_seq_id][$pos] !~ /[A-Z]/) {
			    $pos++;
			}

			$MapRangesMSA[$repped_seq_id][$pos] =~ /^(\d+)\.\./;
			my $last_nucl_range = $1.'..'.$1;

			$right_ranges_str = $right_ranges_str.'/'.$last_nucl_range;
			
		    } else {

			my $last_codon_end = $last_codon_center;		    
			
			if ($revcomp) { $last_codon_end--; }
			else          { $last_codon_end++; }

			$right_ranges_str =~ s/\.\d+$/\.$last_codon_end/;
			
		    }

		    # Now we can party!
		    foreach my $range (split(/\//,$right_ranges_str)) {

			my $sfetch_cmd = $sfetch.' -range '.$range.' '.$genome.' '.$canon_chr;
			my $Fetched = OpenSystemCommand($sfetch_cmd);

			my $line = <$Fetched>; # Eat the header
			while ($line = <$Fetched>) {

			    $line =~ s/\n|\r//g;
			    next if (!$line);

			    $full_nucl_str = $full_nucl_str.lc($line);

			}

			close($Fetched);

			$full_nucl_str = $full_nucl_str.'/';
			
		    }

		    $full_nucl_str =~ s/\/$//;
		    
		}


		# Make sure any updates to our nucleotide ranges are preserved
		$ReppedGenomeRanges[$rep_id] = $left_ranges_str.'|'.$core_ranges_str.'|'.$right_ranges_str;


		# As a final bit of fun, let's translate the nucleotides
		# (this will be helpful for debugging and quality assurance)
		my @ReppedSeq = split(//,$repped_seq_str);
		my $amino_num = 0;
		my @Amino;
		my @Trans;
		my $codon = '';
		my $last_nucl_pos;
		my @Nucls = split(//,$full_nucl_str);
		for (my $pos = 0; $pos < scalar(@Nucls); $pos++) {

		    my $nucl = uc($Nucls[$pos]);

		    if ($nucl =~ /\/|\|/) {
			push(@Amino,$nucl);
			push(@Trans,$nucl);
			next;
		    }

		    $codon = $codon.$nucl;
		    if (length($codon) == 3) {

			$Trans[$last_nucl_pos] = TranslateCodon($codon);
			$codon = '';
			
			# Gappy mappings from Spaln can create issues,
			# so we'll want to be able to handle cases where
			# we've exhausted our 'ReppedSeq'
			if ($amino_num < scalar(@ReppedSeq)) {
			
			    $Amino[$last_nucl_pos] = $ReppedSeq[$amino_num++];
			    
			    # Use lowercase to signify a mismatch
			    if ($Amino[$last_nucl_pos] ne $Trans[$last_nucl_pos]) {
				$Amino[$last_nucl_pos] = lc($Amino[$last_nucl_pos]);
			    }

			} else {

			    # We'll provide a warning about this sequence being gappy
			    my $gms_key = $SeqNames[$repped_seq_id].'&'.($dce_index+1);
			    $GappyMappingSeqs{$gms_key} = 1;
			    
			}
			
		    }

		    $last_nucl_pos = scalar(@Amino);
		    
		    push(@Amino,' ');
		    push(@Trans,' ');
			
		}

		my $final_amino_str = '';
		my $final_trans_str = '';
		for (my $pos = 0; $pos < scalar(@Nucls); $pos++) {
		    $final_amino_str = $final_amino_str.$Amino[$pos];
		    $final_trans_str = $final_trans_str.$Trans[$pos];
		}

		$OutAminoSeqs[$rep_id] = $final_amino_str;
		$OutNuclSeqs[$rep_id]  = $full_nucl_str;
		$OutTransSeqs[$rep_id] = $final_trans_str;
		

	    }


	    # AAAAAAND, now we do Output!!!
	    $dce_index++;

	    my $out_filename = $out_dirname.$dce_index.'.'.$gene.'.out';
	    my $OutFile = OpenOutputFile($out_filename);

	    print $OutFile "Species     : $species\n";
	    print $OutFile "Gene Family : $gene\n";
	    print $OutFile "DCE Index   : $dce_index\n";
	    print $OutFile "Chromosome  : $canon_chr";
	    print $OutFile "[revcomp]" if ($revcomp);
	    print $OutFile "\n\n";
		
	    my $group_number = 0;


	    for (my $rep_id=0; $rep_id<$num_repped; $rep_id++) {

		next if (!$IsGroupLeader[$rep_id]);
		
		$group_number++;
		my $group_title = "> Group $group_number:";

		print $OutFile "\n";
		print $OutFile "$group_title $SeqNames[$ReppedSeqIDs[$rep_id]]\n";

		for (my $check_id=$rep_id+1; $check_id<$num_repped; $check_id++) {

		    if (!$IsGroupLeader[$check_id] && $IsIdenticalTo[$check_id] == $rep_id) {
			for (my $i=0; $i<length($group_title); $i++) {
			    print $OutFile " ";
			}
			print $OutFile " $SeqNames[$ReppedSeqIDs[$check_id]]\n";
		    }
		    
		}

		print $OutFile "\n";

		
		# What genomic coordinates are implicated with this group?
		my $genome_ranges = $ReppedGenomeRanges[$rep_id];
		$genome_ranges =~ s/\// \/ /g;
		$genome_ranges =~ s/\|/ \| /g;		
		print $OutFile "  Genome Ranges: $genome_ranges\n";

		print $OutFile "\n";


		# How would a human enjoy this data?
		my @Amino = split(//,$OutAminoSeqs[$rep_id]);
		my @Nucls = split(//,$OutNuclSeqs[$rep_id]);
		my @Trans = split(//,$OutTransSeqs[$rep_id]);
		my $out_len = scalar(@Amino);

		my $line_start_len = 15;
		my $line_start     = ' ';
		while (length($line_start) < $line_start_len) {
		    $line_start = $line_start.' ';
		}

		my $line_len = 60;
		for (my $pos=0; $pos<$out_len; $pos+=$line_len) {

		    print $OutFile "\n";
		    
		    # Aminos
		    print $OutFile "$line_start";
		    print $OutFile "Protein Seq.  :   ";

		    for (my $i=$pos; $i<Min($out_len,$pos+$line_len); $i++) {
			print $OutFile "$Amino[$i]";
		    }
		    print $OutFile "\n";

		    
		    # DNA
		    print $OutFile "$line_start";
		    print $OutFile "Coding Nucl.s :   ";

		    for (my $i=$pos; $i<Min($out_len,$pos+$line_len); $i++) {
			print $OutFile "$Nucls[$i]";
		    }
		    print $OutFile "\n";

		    
		    # Translation
		    print $OutFile "$line_start";
		    print $OutFile "Translation   :   ";


		    for (my $i=$pos; $i<Min($out_len,$pos+$line_len); $i++) {
			print $OutFile "$Trans[$i]";
		    }
		    print $OutFile "\n";


		    # Line down!
		    print $OutFile "\n";
		    
		}
		
		print $OutFile "\n";
		
	    }
	    
	    close($OutFile);
	    
	    # THAT'S ANOTHER DCE DOWN!
	    
	}
	
    }
    
    
    # Species/gene DCEs examined!
    return $dce_index;

}







#############################################################################
#
#  Function:  RemoveEmptyColumns
#
sub RemoveEmptyColumns
{
    my $orig_msa_ref = shift;
    my $num_seqs     = shift;
    my $orig_msa_len = shift;

    my @OrigMSA = @{$orig_msa_ref};

    my @MSA;
    my $msa_len = 0;

    my $num_exons = 0;

    for (my $j=0; $j<$orig_msa_len; $j++) {

	# Splice column
	if ($OrigMSA[1][$j] eq '/') {

	    $num_exons++;

	    for (my $i=0; $i<=$num_seqs; $i++) {
		$MSA[$i][$msa_len] = '/';
	    }
	    $msa_len++;

	    next;
	    
	}


	# Not a splice column (but maybe empty?)
	my $empty_col = 1;

	$MSA[0][$msa_len] = $num_exons;
	
	for (my $i=1; $i<=$num_seqs; $i++) {

	    $MSA[$i][$msa_len] = $OrigMSA[$i][$j];

	    if ($OrigMSA[$i][$j] =~ /[A-Z]/) {
		$empty_col = 0;
	    }

	}

	$msa_len++ if ($empty_col == 0);
	
    }

    return (\@MSA,$msa_len,$num_exons);
}







#############################################################################
#
#  Function: RecordWindowsToCSV
#
sub RecordWindowsToCSV
{
    my $CSV = shift;
    my $in_filename = shift;

    $in_filename =~ /\/(\S+)\/(\d+)\.([^\/]+)\.out$/;
    my $species = $1;
    my $dce_index = $2;
    my $gene = $3;

    my $InFile = OpenInputFile($in_filename);

    my %WindowToGroup;
    my @WindowsInOrder;

    my $line = <$InFile>;
    while ($line !~ /^\s+Percents identity/) {

	if ($line !~ /\> Group (\d+)/) {
	    $line = <$InFile>;
	    next;
	}
	my $group = $1;

	# Advance to the start of this group's alignment visualization
	$line = <$InFile> while ($line !~ /Protein Seq\./);

	my $group_nucl_str = '';
	while ($line =~ /Protein Seq\./) {

	    $line = <$InFile>; # Coding nucleotides

	    $line =~ /Coding Nucl\.s\s+\:\s+(\S+)/;
	    $group_nucl_str = $group_nucl_str.$1;

	    $line = <$InFile>; # Translation
	    $line = <$InFile>; # Blank Line 1
	    $line = <$InFile>; # Blank Line 2
	    $line = <$InFile>; # Blank Line 3  *or*  Protein Seq.
	    
	}

	my @Nucls = split(//,$group_nucl_str);

	my $splice_in_pos = 0;
	while ($Nucls[$splice_in_pos] ne '|') {
	    $splice_in_pos++;
	}

	my $splice_out_pos = $splice_in_pos+1;
	while ($Nucls[$splice_out_pos] ne '|') {
	    $splice_out_pos++;
	}

	
	# Build up the left window
	my $L_L_str = ''; # nucleotides left  of the upstream splice junction
	my $L_R_str = ''; # nucleotides right of the upstream splice junction

	my $pos = $splice_in_pos-1;
	while ($pos >= 0 && length($L_L_str) < 16) {

	    if ($Nucls[$pos] =~ /[A-Za-z]/) {
		$L_L_str = $Nucls[$pos].$L_L_str;
	    }

	    $pos--;
	    
	}

	$pos = $splice_in_pos+1;
	while ($Nucls[$pos] ne '|' && length($L_R_str) < 16) {

	    if ($Nucls[$pos] =~ /[A-Za-z]/) {
		$L_R_str = $L_R_str.$Nucls[$pos];
	    }

	    $pos++;
	    
	}

	my $L_str = $L_L_str.$L_R_str;

	if (length($L_str) != 32) {
	    if (!$L_L_str) { $L_str = "N-Terminal";                  }
	    else           { $L_str = "Insufficient Nucls ($L_str)"; }
	}


	# Build up the right window
	my $R_L_str = ''; # nucleotides left  of the downstream splice junction
	my $R_R_str = ''; # nucleotides right of the downstream splice junction

	$pos = $splice_out_pos-1;
	while ($Nucls[$pos] ne '|' && length($R_L_str) < 16) {

	    if ($Nucls[$pos] =~ /[A-Za-z]/) {
		$R_L_str = $Nucls[$pos].$R_L_str;
	    }

	    $pos--;
	    
	}

	$pos = $splice_out_pos+1;
	while ($pos < scalar(@Nucls) && length($R_R_str) < 16) {

	    if ($Nucls[$pos] =~ /[A-Za-z]/) {
		$R_R_str = $R_R_str.$Nucls[$pos];
	    }

	    $pos++;
	    
	}

	my $R_str = $R_L_str.$R_R_str;

	if (length($R_str) != 32) {
	    if (!$R_R_str) { $R_str = "C-Terminal";                  }
	    else           { $R_str = "Insufficient Nucls ($R_str)"; }
	}


	# Put the left and the right together, and you've got yourself some windows!
	my $window = $L_str.'&'.$R_str;

	if ($WindowToGroup{$window}) {
	    $WindowToGroup{$window} = $WindowToGroup{$window}.'/'.$group;
	} else {
	    $WindowToGroup{$window} = $group;
	}
	
	push(@WindowsInOrder,$window);
	
    }

    
    # Before printing to the CSV, we'll also want to associate each window with
    # its alignment percent identity
    my %GroupToPctID;
    my %GroupToFrame;
    while ($line !~ /^\s+Overlaid alignment/) {

	if ($line !~ /Groups?\s+(\S+)\s+\:\=\s+(\S+)\s+\[frames?\:([^\]]+)\]/) {
	    $line = <$InFile>;
	    next;
	}
	
	my $group  = $1;
	my $pct_id = $2;
	my $frame  = $3;

	$group =~ s/\,/\//g;
	$frame =~ s/\,/\//g;

	$GroupToPctID{$group} = $pct_id;
	$GroupToFrame{$group} = $frame;

	# We'll also want all individual groups to have their percents ID
	# recorded, in case they don't group w.r.t. 32-nucleotide windows
	foreach my $subgroup (split(/\//,$group)) {
	    $GroupToPctID{$subgroup} = $pct_id;
	    $GroupToFrame{$subgroup} = $frame;
	}

	$line = <$InFile>;

    }


    # Now we just need to write our windows out to the CSV!
    foreach my $window (@WindowsInOrder) {

	next if (!$WindowToGroup{$window});

	my $group = $WindowToGroup{$window};

	# Do we have funniness with the boundary def.s?
	my $pct_id;
	my $frame;
	if (!$GroupToPctID{$group}) {

	    $group =~ /^(\d+)\//;
	    my $lead_fam = $1;

	    $pct_id = $GroupToPctID{$lead_fam};
	    $frame  = $GroupToFrame{$lead_fam};

	} else {

	    $pct_id = $GroupToPctID{$group};
	    $frame  = $GroupToFrame{$group};

	}


	$window =~ /^([^\&]+)\&([^\&]+)$/;
	my $left  = $1;
	my $right = $2;

	print $CSV "$dce_index, $gene, $group, $frame, $left, $right, $pct_id\n";

	$WindowToGroup{$window} = 0;
	
    }
    
    close($InFile);
    
}








#############################################################################
#
#  Function:  VisDualCodingRegion
#
sub VisDualCodingRegion
{

    my $fname = shift;
    
    my $InFile = OpenInputFile($fname);

    my $line = <$InFile>; # Species
    $line = <$InFile>;    # Gene
    $line = <$InFile>;    # DCE Index
    $line = <$InFile>;    # Chromosome
    
    my $strand = 1;
    if ($line =~ /\[revcomp\]/) {
	$strand = -1;
    }

    my $num_groups = 0;

    # It's a real sick thing what we're doing here
    my %CoordToGroupAmino;
    my %CoordToNucl;
    my %GroupToIndex;
    my @GroupNames;
    
    $line = <$InFile>; # Buffer line 1
    $line = <$InFile>; # Buffer line 2

    $line = <$InFile>; # Group start line
    while ($line =~ /\> Group (\d+)\:/) {

	my $group = $1;
	$num_groups++;

	$GroupToIndex{$group} = $num_groups;
	$GroupNames[$num_groups] = $group;

	while ($line !~ /Genome Ranges\:/) {
	    $line = <$InFile>;
	}

	# Grab the key coord. ranges
	$line =~ /\|([^\|]+)\|/;
	my $coord_ranges_str = $1;
	$coord_ranges_str =~ s/\s//g;

	my @CoordRanges = split(/\//,$coord_ranges_str);

	$line = <$InFile>; # Dead line 1
	$line = <$InFile>; # Dead line 2

	my $prot_str = '';
	my $nucl_str = '';

	$line = <$InFile>; # Protein Seq.
	while ($line =~ /Protein Seq\./) {

	    $line =~ s/\n|\r//g;
	    $line =~ /\:   (.*)$/;

	    $prot_str = $prot_str.$1;

	    $line = <$InFile>; # Coding Nucl.s
	    $line =~ s/\n|\r//g;
	    $line =~ /\:   (.*)$/;

	    $nucl_str = $nucl_str.$1;

	    $line = <$InFile>; # Translation
	    $line = <$InFile>; # Dead line 1
	    $line = <$InFile>; # Dead line 2

	    last if (eof($InFile));

	    $line = <$InFile>; # Protein Seq.?
	    
	}


	$prot_str =~ s/\///g;
	$prot_str =~ /\|([^\|]+)\|/;
	$prot_str = $1;

	$nucl_str =~ s/\///g;
	$nucl_str =~ /\|([^\|]+)\|/;
	$nucl_str = $1;

	my @Aminos = split(//,$prot_str);
	my @Nucls  = split(//,$nucl_str);

	my $pos = 0;
	
	foreach my $nucl_range (@CoordRanges) {

	    $nucl_range =~ /^(\d+)\.\.(\d+)$/;
	    my $range_start = $1;
	    my $range_end   = $2;
	    
	    while ($range_start != $range_end + $strand) {

		my $amino = $Aminos[$pos];
		my $nucl  = $Nucls[$pos];
		$pos++;

		$CoordToNucl{$range_start} = $nucl;

		if ($amino =~ /[A-Za-z]/) {

		    if ($CoordToGroupAmino{$range_start}) {
			$CoordToGroupAmino{$range_start} =
			    $CoordToGroupAmino{$range_start}.'&'.$group.':'.$amino;
		    } else {
			$CoordToGroupAmino{$range_start} = $group.':'.$amino;
		    }

		}

		$range_start += $strand;
		
	    }
	    
	}

	if (!eof($InFile)) {
	    $line = <$InFile>;
	    next;
	}
	
    }
    close($InFile);


    my @CoordList;
    if ($strand == 1) { @CoordList = sort { $a <=> $b } keys %CoordToNucl; }
    else              { @CoordList = sort { $b <=> $a } keys %CoordToNucl; }


    my @Ali;
    my $ali_len = 0;
    for (my $coord_id = 0; $coord_id < scalar(@CoordList); $coord_id++) {
	
	my $coord = $CoordList[$coord_id];

	# Intron catch
	if ($coord_id && abs($coord-$CoordList[$coord_id-1]) > 1) {
	    for (my $i=0; $i<=$num_groups; $i++) {
		$Ali[$i][$ali_len] = '/';
	    }
	    $ali_len++;
	}

	# Now, back to work!
	$Ali[0][$ali_len] = $CoordToNucl{$coord};
	
	for (my $i=1; $i<=$num_groups; $i++) {
	    $Ali[$i][$ali_len] = ' ';
	}

	# If this coordinate isn't the center of any codon, move along
	if (!$CoordToGroupAmino{$coord}) {
	    $ali_len++;
	    next;
	}

	foreach my $group_amino_pair (split(/\&/,$CoordToGroupAmino{$coord})) {

	    $group_amino_pair =~ /^([^\:]+)\:(\S)$/;
	    my $group = $1;
	    my $amino = $2;
	    
	    my $group_id = $GroupToIndex{$group};

	    $Ali[$group_id][$ali_len] = $amino;
	    
	}

	# That's another column down!
	$ali_len++;
	
    }


    # Let's translate out each frame! FUN!
    my @TransAli;
    for (my $frame = 0; $frame < 3; $frame++) {

	# Prime the start of the translations
	for (my $pos=0; $pos<$frame; $pos++) {
	    $TransAli[$frame][$pos] = ' ';
	}

	# Now to the party!
	my $codon = '';
	my $last_nucl_pos;
	for (my $pos = $frame; $pos < $ali_len; $pos++) {

	    $TransAli[$frame][$pos] = ' ';
	    
	    if ($Ali[0][$pos] =~ /[A-Za-z]/) {
		
		$codon = $codon.$Ali[0][$pos];

		if (length($codon) == 3) {
		    
		    $TransAli[$frame][$last_nucl_pos] = TranslateCodon($codon);
		    
		    $codon = '';
		    
		}

		$last_nucl_pos = $pos;
		
	    }
	    
	}
	
    }


    # Before we output, we'll see how each of our mapped amino sequences
    # matches with the translated nucleotides...
    my @GroupMatches;
    my @GroupMismatches;
    my @GroupPctsID;
    for (my $group_id=1; $group_id<=$num_groups; $group_id++) {

	$GroupMatches[$group_id]    = 0;
	$GroupMismatches[$group_id] = 0;

	for (my $col=0; $col<$ali_len; $col++) {

	    if (!$Ali[$group_id][$col]) {
		$Ali[$group_id][$col] = ' ';
		next;
	    }

	    if ($Ali[$group_id][$col] =~ /[A-Za-z]/) {

		my $trans_amino;
		if    ($TransAli[0][$col] =~ /[A-Za-z]/) { $trans_amino = uc($TransAli[0][$col]); }
		elsif ($TransAli[1][$col] =~ /[A-Za-z]/) { $trans_amino = uc($TransAli[1][$col]); }
		elsif ($TransAli[2][$col] =~ /[A-Za-z]/) { $trans_amino = uc($TransAli[2][$col]); }

		# Uppercase if we have a match (or nothing to compare),
		# lowercase for a mismatch
		if (!$trans_amino || $trans_amino eq uc($Ali[$group_id][$col])) {

		    $Ali[$group_id][$col] = uc($Ali[$group_id][$col]); #    match

		    $GroupMatches[$group_id]++;

		} else {

		    $Ali[$group_id][$col] = lc($Ali[$group_id][$col]); # mismatch

		    $GroupMismatches[$group_id]++;

		}
		
	    }

	}

	# DEBUGGING
	if ($GroupMatches[$group_id] + $GroupMismatches[$group_id] == 0) {
	    print "\n  APPARENT MAPPING ERROR: $fname\n\n";
	    return;
	}

	$GroupPctsID[$group_id] = int(1000.0 * $GroupMatches[$group_id] / ($GroupMatches[$group_id] + $GroupMismatches[$group_id])) / 10.0;
	
	$GroupPctsID[$group_id] = $GroupPctsID[$group_id].'.0' if ($GroupPctsID[$group_id] !~ /\.\d$/);
	$GroupPctsID[$group_id] = $GroupPctsID[$group_id].'%';
	while (length($GroupPctsID[$group_id]) < length('100.0%')) {
	    $GroupPctsID[$group_id] = ' '.$GroupPctsID[$group_id];
	}
	
    }


    # The very last thing we'll do is consolidate any groups that are
    # identical within this region.
    # Of course, it would make more sense to have done this earlier, but
    # I'm LAZY
    my @FinalGroups;
    my @FinalGroupNames;
    my $final_num_groups = 0;

    for (my $group_id=1; $group_id<=$num_groups; $group_id++) {
	$FinalGroups[$group_id] = 0;
    }

    for (my $group_id=1; $group_id<=$num_groups; $group_id++) {

	next if ($FinalGroups[$group_id]);

	$final_num_groups++;
	$FinalGroups[$group_id] = $final_num_groups;

	$FinalGroupNames[$final_num_groups] = $GroupNames[$group_id];

	for (my $check_id=$group_id+1; $check_id<=$num_groups; $check_id++) {

	    my $is_identical = 1;
	    for (my $pos=0; $pos<$ali_len; $pos++) {
		if ($Ali[$group_id][$pos] ne $Ali[$check_id][$pos]) {
		    $is_identical = 0;
		    last;
		}
	    }

	    next if (!$is_identical);
	    
	    $FinalGroups[$check_id] = $final_num_groups;
	    $FinalGroupNames[$final_num_groups] =
		$FinalGroupNames[$final_num_groups].','.$check_id;
	    
	}
	
    }


    my @FinalGroupPctsID;

    # Finally time to prep for output!  We'll want to figure out
    # the length of the longest group name for left-side buffering
    my $longest_name_len = 0;
    for (my $group_id=1; $group_id<=$final_num_groups; $group_id++) {

	my $group_name = $FinalGroupNames[$group_id];

	if ($group_name =~ /\,/) {
	    $group_name = 'Groups '.$group_name.'    ';
	} else {
	    $group_name = 'Group  '.$group_name.'    ';
	}

	$FinalGroupNames[$group_id] = $group_name;

	if (length($group_name) > $longest_name_len) {
	    $longest_name_len = length($group_name);
	}

	$group_name =~ /Groups?\s+(\d+)/;
	$FinalGroupPctsID[$group_id] = $GroupPctsID[$1];

    }

    for (my $group_id=1; $group_id<=$final_num_groups; $group_id++) {
	while (length($FinalGroupNames[$group_id]) < $longest_name_len) {
	    $FinalGroupNames[$group_id] = $FinalGroupNames[$group_id].' ';
	}
    }


    # COULD IT BE?! WE'RE OPENING OUR OUTPUT FILE TO PRINT TO?!?!
    open(my $OutFile,'>>',$fname)
	|| die "\n  ERROR:  File append open failed (VDCR:'$fname')\n\n";


    print $OutFile "\n\n-------------------------------------------------------";
    print $OutFile "-------------------------------------------------------\n\n";

    print $OutFile "\n   Percents identity for alignments of translated dual-coding nucleotides to group aminos\n\n";
    for (my $meta_group_id=1; $meta_group_id<=$final_num_groups; $meta_group_id++) {
	print $OutFile "   $FinalGroupNames[$meta_group_id]:= $FinalGroupPctsID[$meta_group_id]\n";
    }
    
    print $OutFile "\n\n\n   Overlaid alignment of the dual-coding regions for each group\n";

    my $line_len = 60;
    for (my $line_start = 0; $line_start < $ali_len; $line_start += $line_len) {

	my $line_end = Min($line_start+$line_len,$ali_len);

	
	print $OutFile "\n";
	

	# Print each of the 'final' groups
	for (my $meta_group_id = $final_num_groups; $meta_group_id > 0; $meta_group_id--) {

	    print $OutFile "   $FinalGroupNames[$meta_group_id]";

	    # Find a representative member
	    my $group_id = $num_groups;
	    while ($FinalGroups[$group_id] != $meta_group_id) {
		$group_id--;
	    }
	    
	    for (my $pos = $line_start; $pos < $line_end; $pos++) {
		print $OutFile "$Ali[$group_id][$pos]";
	    }
	    print $OutFile "\n";
	    
	}

	
	# Print the nucleotide sequence
	my $nucl_name = 'Nucl.s';
	while (length($nucl_name) < $longest_name_len) {
	    $nucl_name = $nucl_name.' ';
	}

	print $OutFile "   $nucl_name";

	for (my $pos = $line_start; $pos < $line_end; $pos++) {
	    print $OutFile "$Ali[0][$pos]";
	}
	print $OutFile "\n";
	print $OutFile "\n"; # Extra break before translated frames

	
	# And, finally, print the translations!
	for (my $frame = 0; $frame < 3; $frame++) {

	    my $fmtd_name = 'Frame '.($frame+1);
	    while (length($fmtd_name) < $longest_name_len) {
		$fmtd_name = $fmtd_name.' ';
	    }
	    
	    print $OutFile "   $fmtd_name";
	    
	    for (my $pos = $line_start; $pos < $line_end; $pos++) {
		print $OutFile "$TransAli[$frame][$pos]";
	    }
	    print $OutFile "\n";
	    
	}


	# Add a lil' extra for prettiness
	print $OutFile "\n\n";

	
    }

    print $OutFile "\n";
    close($OutFile);
    
}








#############################################################################
#
#  Function:  CheckForMidExonIntrons
#
#  About:  This subroutine checks whether the dual coding region is created
#          by intron insertion w.r.t. genomic sequence that is exonic in one
#          of the other sequences (indicating that the reported left/right
#          window may need adjustment)
#
sub CheckForMidExonIntrons
{
    my $MidExonIntronFile = shift;
    my $dce_filename = shift;
    
    $dce_filename =~ /\/([^\/]+)\/(\d+)\.([^\/]+)\.out$/;
    my $species = $1;
    my $dce_index = $2;
    my $gene = $3;

    my $InFile = OpenInputFile($dce_filename);

    while (my $line = <$InFile>) { 
	last if ($line =~ /Percents identity/); 
    }


    my %GroupNamesToFormatted;
    my $longest_name_len = 0;
    while (my $line = <$InFile>) {

	last if ($line =~ /Overlaid alignment/);

	if ($line =~ /Groups?\s+(\S+)/) {

	    my $group_name = $1;
	    if (length($group_name) > $longest_name_len) {
		$longest_name_len = length($group_name);
	    }

	    $GroupNamesToFormatted{$group_name} = $group_name;

	}

    }

    # We may need to do some special spacing to make sure we start
    # grabbing characters at the right place to determine frame
    #
    foreach my $group_name (keys %GroupNamesToFormatted) {

	my $formatted_name = $GroupNamesToFormatted{$group_name};
	while (length($formatted_name) < $longest_name_len) {
	    $formatted_name = $formatted_name.' ';
	}

	$GroupNamesToFormatted{$group_name} = $formatted_name;

    }

    my @GroupMSA;
    my %GroupNamesToIDs;
    my $num_groups = 0;
    while (my $line = <$InFile>) {

	$line =~ s/\n|\r//g;

	next if ($line !~ /Groups?\s+(\S+)/);
	my $group_name = $GroupNamesToFormatted{$1};

	$line =~ /Groups?\s+$group_name    (.+)$/;
	my $ali_chars = $1;

	my $group_id;
	if (!$GroupNamesToIDs{$group_name}) {
	    $GroupNamesToIDs{$group_name} = ++$num_groups;
	    $group_id = $num_groups;
	} else {
	    $group_id = $GroupNamesToIDs{$group_name};
	}

	my $group_msa_len;
	if ($GroupMSA[$group_id]) {
	    $group_msa_len = scalar(@{$GroupMSA[$group_id]});
	} else {
	    $group_msa_len = 0;
	}

	foreach my $char (split(//,$ali_chars)) {
	    if ($char ne '/') {
		$GroupMSA[$group_id][$group_msa_len++] = $char;
	    }
	}
	
    }

    my $multi_frame_observed = 0;
    my @FramesUsed;
    for (my $group_id=1; $group_id<=$num_groups; $group_id++) {

	$FramesUsed[$group_id][0] = 0;
	$FramesUsed[$group_id][1] = 0;
	$FramesUsed[$group_id][2] = 0;
	
	for (my $pos=0; $pos<scalar(@{$GroupMSA[$group_id]}); $pos++) {
	    if ($GroupMSA[$group_id][$pos] =~ /[A-Za-z]/) {
		$FramesUsed[$group_id][$pos % 3] = 1;
	    }
	}

	if ($FramesUsed[$group_id][0] +
	    $FramesUsed[$group_id][1] +
	    $FramesUsed[$group_id][2] > 1) {
	    $multi_frame_observed = 1;
	}
	
    }
    
    close($InFile);

    print $MidExonIntronFile "$species: $dce_index ($gene)\n"
	if ($multi_frame_observed);


    # The last thing that we'll do is add the frames used to the dce info file
    my $tmp_filename = $dce_filename;
    $tmp_filename =~ s/\.out$/\.tmp/;

    $InFile = OpenInputFile($dce_filename);
    my $TmpOutFile = OpenOutputFile($tmp_filename);
    while (my $line = <$InFile>) {

	$line =~ s/\n|\r//g;

	if ($line =~ /Groups?\s+(\S+\s*)    \:\=\s+\S+\%/) {

	    my $group_id = $GroupNamesToIDs{$1};

	    my $frames_used = '';
	    for (my $frame=0; $frame<3; $frame++) {
		if ($FramesUsed[$group_id][$frame]) {
		    $frames_used = $frames_used.','.($frame+1);
		}
	    }
	    $frames_used =~ s/^\,//;

	    $line = $line.'  [frame';
	    $line = $line.'s' if ($frames_used =~ /\,/);
	    $line = $line.':'.$frames_used.']';

	}

	print $TmpOutFile "$line\n";
	
    }
    close($InFile);
    close($TmpOutFile);

    RunSystemCommand("mv \"$tmp_filename\" \"$dce_filename\"");
    
}








#############################################################################
#
#  Function: GenStrongWindowFastas
#
#  About: This function compiles a FASTA file for each species, containing
#         the left and/or right windows nucleotide windows.  It only does this
#         for sequences from DCE indices that are not implicated in either
#         'Gappy-Alignment-Warnings' or 'Mid-Exon-Intron-Warnings'
#
sub GenStrongWindowFastas
{
    my $dirname = shift;

    my %WeakDCEIndices;

    # DCE indices where there are low-quality mappings
    if (-e($dirname.'Gappy-Alignment-Warnings.out')) {

	my $InFile = OpenInputFile($dirname.'Gappy-Alignment-Warnings.out');

	while (my $line = <$InFile>) {
	    if ($line =~ /\:\s+(\d+)\s+/) {
		$WeakDCEIndices{$1} = 1;
	    }
	}
	close($InFile);
	
    }

    # DCE indices where the alt. frame has a funny look...
    if (-e ($dirname.'Mid-Exon-Intron-Warnings.out')) {

	my $InFile = OpenInputFile($dirname.'Mid-Exon-Intron-Warnings.out');

	while (my $line = <$InFile>) {
	    if ($line =~ /\:\s+(\d+)\s+/) {
		$WeakDCEIndices{$1} = 1;
	    }
	}
	close($InFile);
	
    }

    # Now let's pull the right/left windows!  We'll add an extra check for
    # alignments where there's mapping identity less than 95%, since that
    # would imply the mapping isn't strong enough to trust for RNA-Seq verification
    my $ResultsDir = OpenDirectory($dirname);
    while (my $in_filename = readdir($ResultsDir)) {

	next if ($in_filename !~ /^(\S+)\.csv$/);
	my $species = $1;

	$in_filename     = $dirname.$in_filename;
	my $out_filename = $dirname.$species.'.fa';

	my $InFile  = OpenInputFile($in_filename);
	my $OutFile = OpenOutputFile($out_filename);

	my $current_dce_index = 0;
	my $num_current_index_seqs = 0;
	my @IndexLines;
	while (my $line = <$InFile>) {

	    next if ($line !~ /^(\d+)\,/);
	    my $next_dce_index = $1;

	    if ($next_dce_index != $current_dce_index) {

		# Wooo! let's write this collection of seqs!
		if ($num_current_index_seqs > 1 && !$WeakDCEIndices{$current_dce_index}) {
		    
		    for (my $i=0; $i<$num_current_index_seqs; $i++) {

			my $index_line = $IndexLines[$i];
			$index_line =~ /^[^\,]+\, ([^\,]+)\, ([^\,]+)\, ([^\,]+)\, ([^\,]+)\, ([^\,]+)\,/;
			my $gene  = $1;
			my $group = $2;
			my $frame = $3;
			my $left  = $4;
			my $right = $5;

			if ($left !~ /Terminal|Insufficient/) {

			    my $seq_name = '>'.$current_dce_index.'_'.$gene.'_'.$group;
			    $seq_name = $seq_name.'_frame'.$frame.'_left';
			    
			    print $OutFile "$seq_name\n$left\n";
			    
			}
			
			if ($right !~ /Terminal|Insufficient/) {
			    
			    my $seq_name = '>'.$current_dce_index.'_'.$gene.'_'.$group;
			    $seq_name = $seq_name.'_frame'.$frame.'_right';

			    print $OutFile "$seq_name\n$right\n";
			    
			}
			
		    }
		}

		$current_dce_index = $next_dce_index;
		$num_current_index_seqs = 0;
		
	    }
	    
	    $line =~ /\, ([^\,]+)\, ([^\,]+)\, ([^\,]+)\%\s*$/;
	    my $left   = $1;
	    my $right  = $2;
	    my $pct_id = $3;
	    
	    if ($pct_id > 95.0 && ($left !~ /Terminal|Insufficient/ || $right !~ /Terminal|Insufficient/)) {
		$IndexLines[$num_current_index_seqs++] = $line;
	    }
	    
	}
	
    }
    closedir($ResultsDir);
    
    
}










#############################################################################
#
#  Function: GenMidExonIntronFastas
#
sub GenMidExonIntronFastas
{
    my $result_dirname = shift;

    my $mei_fname = $result_dirname.'Mid-Exon-Intron-Warnings.out';
    return if (!(-e $mei_fname));

    my $out_dirname = CreateDirectory($result_dirname.'Mid-Exon-Intron-FASTAs');
    my $num_outputs = 0;

    my $MEIFile = OpenInputFile($mei_fname);
    while (my $mei_line = <$MEIFile>) {

	next if ($mei_line !~ /^(\S+)\:\s+(\d+)\s+\((\S+)\)/);

	my $species = $1;
	my $dce_index = $2;
	my $gene = $3;

	my $dce_fname = $result_dirname.$species.'/'.$dce_index.'.'.$gene.'.out';

	my $DCEFile = OpenInputFile($dce_fname);

	my %GroupToAli;
	while (my $dce_line = <$DCEFile>) {

	    last if ($dce_line =~ /Percents identity/);

	    next if ($dce_line !~ /\> Group (\d+)\:/);
	    my $group_id = $1;

	    while ($dce_line !~ /Genome Ranges/) {
		$dce_line = <$DCEFile>;
	    }
	    $dce_line = <$DCEFile>; # blank
	    $dce_line = <$DCEFile>; # blank
	    
	    $dce_line = <$DCEFile>; # * Protein Seq.
	    $dce_line =~ s/\n|\r//g;

	    my $protein_ali_str = '';
	    my $nucl_ali_str    = '';

	    while ($dce_line =~ /Protein Seq\.  \:   (.*)$/) {

		$protein_ali_str = $protein_ali_str.$1;

		$dce_line = <$DCEFile>; # * Coding Nucl.s
		$dce_line =~ s/\n|\r//g;

		$dce_line =~ /Coding Nucl\.s \:   (.*)$/;
		
		$nucl_ali_str = $nucl_ali_str.$1;

		$dce_line = <$DCEFile>; # Translation
		$dce_line = <$DCEFile>; # blank
		$dce_line = <$DCEFile>; # blank

		$dce_line = <$DCEFile>; # * Protein Seq.?
		$dce_line =~ s/\n|\r//g;
		
	    }

	    $protein_ali_str =~ s/^[^\|]*\|//;
	    $nucl_ali_str    =~ s/^[^\|]*\|//;

	    $protein_ali_str =~ s/\|.*$//;
	    $nucl_ali_str    =~ s/\|.*$//;

	    $protein_ali_str =~ s/\///g;
	    $nucl_ali_str    =~ s/\///g;

	    $GroupToAli{$group_id} = $protein_ali_str.'&'.$nucl_ali_str;
	    
	}

	
	# Now that we have the alignments for each group, we can see which groups
	# end up using multiple reading frames (i.e., have an intron in the middle
	# of another sequence's exon).
	my %GroupIsMultiFrame;
	my %GroupLeaderToFull;
	my $longest_group_name_len = 0;
	while (my $dce_line = <$DCEFile>) {

	    last if ($dce_line =~ /Overlaid alignment/);

	    if ($dce_line =~ /Groups?\s+(\S+)/) {

		my $group_name_len = length($1);
		if ($group_name_len > $longest_group_name_len) {
		    $longest_group_name_len = $group_name_len;
		}
		
		# Capture the lead group for any multi-frames
		if ($dce_line =~ /Groups?\s+(\S+)\s+\:\=\s+\S+\s+\[frames/) {

		    my $full_group = $1;

		    $full_group =~ /^(\d+)/;
		    my $group_lead = $1;

		    $full_group =~ s/\,/\//g;

		    $GroupIsMultiFrame{$group_lead} = 1;
		    $GroupLeaderToFull{$group_lead} = $full_group;

		}

	    }
	    
	}

	
	# Next up, we'll build up an MSA representing the group alignment
	my @GroupLeadNums;
	my @GroupMSAStrs;
	my $num_groups = -1;
	my $full_ali_nucls = '';
	while (my $dce_line = <$DCEFile>) {

	    next if ($dce_line !~ /Groups?/);

	    $num_groups = 0;
	    while ($dce_line =~ /Groups?\s+(\S+)/) {

		my $full_group_name = $1;

		$full_group_name =~ /^(\d+)/;
		my $lead_group_num = $1;

		while (length($full_group_name) < $longest_group_name_len) {
		    $full_group_name = $full_group_name.' ';
		}

		$dce_line =~ s/\n|\r//g;
		$dce_line =~ /Groups?\s+$full_group_name    (.*)$/;

		my $ali_str = $1;

		$GroupLeadNums[$num_groups] = $lead_group_num;

		if ($GroupMSAStrs[$num_groups]) {
		    $GroupMSAStrs[$num_groups] = $GroupMSAStrs[$num_groups].$ali_str;
		} else {
		    $GroupMSAStrs[$num_groups] = $ali_str;
		}

		$num_groups++;

		$dce_line = <$DCEFile>; # Either 'Group' or 'Nucl.s'

	    }

	    $dce_line =~ /Nucl\.s\s+(\S+)/;
	    $full_ali_nucls = $full_ali_nucls.$1;

	    $dce_line = <$DCEFile>; # blank
	    $dce_line = <$DCEFile>; # Frame 1
	    $dce_line = <$DCEFile>; # Frame 2
	    $dce_line = <$DCEFile>; # Frame 3
	    $dce_line = <$DCEFile>; # blank
	    $dce_line = <$DCEFile>; # blank
	    $dce_line = <$DCEFile>; # blank

	}

	# That's it for you, DCE File!
	close($DCEFile);

	# Quick sanity check...
	if ($num_groups == -1) {
	    print "\n  WARNING: No groups found in alignment for $species $gene\n\n";
	    next;
	}
	
	my @FullAliNucls = split(//,$full_ali_nucls);
	
	# For each multi-frame group, find the index of the last amino acid before
	# a long gap (indicating the inserted intron)
	for (my $i=0; $i<$num_groups; $i++) {

	    my $group_id = $GroupLeadNums[$i];
	    next if (!$GroupIsMultiFrame{$group_id});

	    my @AliChars = split(//,$GroupMSAStrs[$i]);

	    my $pre_gap_amino_index = 0;
	    my $num_adjacent_gaps   = 0;

	    # Find the index in the alignment corresponding to the start of the
	    # big gap (i.e., the intron)
	    my $char_index = 0;
	    my $pre_gap_amino_coord;
	    while ($char_index < scalar(@AliChars) && $num_adjacent_gaps < 4) {

		if ($AliChars[$char_index] =~ /[A-Za-z]/) {
		    $num_adjacent_gaps = 0;
		    $pre_gap_amino_index++;
		    $pre_gap_amino_coord = $char_index;
		} else {
		    $num_adjacent_gaps++;
		}

		$char_index++;
		
	    }

	    # Weird if this happens, but better catch...
	    next if ($char_index == scalar(@AliChars));

	    $pre_gap_amino_index--;

	    $GroupToAli{$group_id} =~ /^([^\&]+)\&([^\&]+)$/;
	    my $group_prot_ali_str = $1;
	    my $group_nucl_ali_str = $2;

	    my @GroupAminos = split(//,$group_prot_ali_str);
	    my @GroupNucls  = split(//,$group_nucl_ali_str);

	    my $amino_index = 0;
	    my $group_pos = 0;
	    while ($amino_index < $pre_gap_amino_index) {

		if ($GroupAminos[$group_pos] =~ /[A-Za-z]/) {
		    $amino_index++;
		}

		$group_pos++;
		
	    }

	    # There should be a splice marker close ahead...
	    while ($GroupAminos[$group_pos] !~ /\//) {
		$group_pos++;
	    }

	    # Great! Now let's just grab the 16 nucl.s in each direction!
	    my $group_nucls_L = '';
	    for (my $pos = $group_pos-1; $pos >= 0; $pos--) {
		if ($GroupNucls[$pos] =~ /[A-Za-z]/) {
		    $group_nucls_L = $GroupNucls[$pos].$group_nucls_L;
		    last if (length($group_nucls_L) == 16);
		}
	    }

	    my $group_nucls_R = '';
	    for (my $pos = $group_pos+1; $pos < scalar(@GroupNucls); $pos++) {
		if ($GroupNucls[$pos] =~ /[A-Za-z]/) {
		    $group_nucls_R = $group_nucls_R.$GroupNucls[$pos];
		    last if (length($group_nucls_R) == 16);
		}
	    }


	    # If we don't have enough nucleotides to be searchable, we'll move
	    # on with our lives!
	    my $group_nucls = $group_nucls_L.$group_nucls_R;
	    next if (length($group_nucls) < 32);

	    # OOOOH LALA!  We've got something cookin'!  Let's really quickly
	    # build up the sequences that we'll want to report for this group
	    my $group_name = '>'.$dce_index.'_'.$gene.'_'.$GroupLeaderToFull{$group_id}.'_IntronInsertion';

	    
	    # The last bit of work we'll need to do is determine nucleotide
	    # windows for sequences that don't have the "bonus" intron

	    #
	    # LEFT SIDE OF "INTRON"
	    #

	    $char_index = $pre_gap_amino_coord + 1; # Last nucl in common codon

	    my $comp_nucls_LL = '';
	    while (length($comp_nucls_LL) < 16 && $char_index >= 0) {

		if ($FullAliNucls[$char_index] =~ /[A-Za-z]/) {
		    $comp_nucls_LL = $FullAliNucls[$char_index].$comp_nucls_LL;
		}
		$char_index--;
		
	    }
	    
	    
	    $char_index = $pre_gap_amino_coord + 2; # First nucl in 'intron'

	    my $comp_nucls_LR = '';
	    while (length($comp_nucls_LR) < 16 && $char_index < scalar(@FullAliNucls)) {

		if ($FullAliNucls[$char_index] =~ /[A-Za-z]/) {
		    $comp_nucls_LR = $comp_nucls_LR.$FullAliNucls[$char_index];
		}
		$char_index++;
		
	    }

	    my $comp_nucls_L = $comp_nucls_LL.$comp_nucls_LR;
	    

	    #
	    # RIGHT SIDE OF "INTRON"
	    #

	    # Walk up to the end of the group of interest's gappy run
	    my $post_gap_amino_coord = $pre_gap_amino_coord+2;

	    while ($post_gap_amino_coord < scalar(@AliChars) &&
		   $AliChars[$post_gap_amino_coord] !~ /[A-Za-z]/) {
		$post_gap_amino_coord++;
	    }
	    

	    $char_index = $post_gap_amino_coord - 2; # Last nucl in 'intron'
	    
	    my $comp_nucls_RL = '';
	    while (length($comp_nucls_RL) < 16 && $char_index >= 0) {

		if ($FullAliNucls[$char_index] =~ /[A-Za-z]/) {
		    $comp_nucls_RL = $FullAliNucls[$char_index].$comp_nucls_RL;
		}
		$char_index--;
		
	    }
	    
	    
	    $char_index = $post_gap_amino_coord - 1; # First nucl in follow-up exon

	    my $comp_nucls_RR = '';
	    while (length($comp_nucls_RR) < 16 && $char_index < scalar(@FullAliNucls)) {

		if ($FullAliNucls[$char_index] =~ /[A-Za-z]/) {
		    $comp_nucls_RR = $comp_nucls_RR.$FullAliNucls[$char_index];
		}
		$char_index++;
		
	    }

	    my $comp_nucls_R = $comp_nucls_RL.$comp_nucls_RR;


	    # GREAT! Name those suckas!
	    my $comp_nucls_L_name = $group_name;
	    $comp_nucls_L_name =~ s/IntronInsertion/LeftSpliceSite/;
	    
	    my $comp_nucls_R_name = $group_name;
	    $comp_nucls_R_name =~ s/IntronInsertion/RightSpliceSite/;


	    # Time to sing some beautiful songs!!!
	    open(my $OutFile,'>>',$out_dirname.$species.'.mid-exon-introns.fa')
		|| die "\n  ERROR:  Failed to open MEI output file ($species)!\n\n";

	    print $OutFile "$group_name\n$group_nucls\n";

	    if (length($comp_nucls_L) == 32) {
		print $OutFile "$comp_nucls_L_name\n$comp_nucls_L\n";
	    }

	    if (length($comp_nucls_R) == 32) {
		print $OutFile "$comp_nucls_R_name\n$comp_nucls_R\n";
	    }

	    close($OutFile);

	    # Mark it!
	    $num_outputs++;

	}

    }
    close($MEIFile);

    if ($num_outputs == 0) { RunSystemCommand("rm -rf \"$out_dirname\""); }
    
}












#############################################################################
#
#  Function: DumpSpeciesGeneData
#
sub DumpSpeciesGeneData
{
    my $out_filename = shift;

    my $canon_chr = shift;
    my $revcomp = shift;
    
    my $msa_ref = shift;
    my $seq_names_ref = shift;
    my $num_seqs = shift;
    my $msa_len = shift;
    
    my $map_msa_ref = shift;
    my $map_ranges_msa_ref = shift;
    my $amino_status_msa_ref = shift;

    my @MSA = @{$msa_ref};
    my @SeqNames = @{$seq_names_ref};
    my @MapMSA = @{$map_msa_ref};
    my @MapRangesMSA = @{$map_ranges_msa_ref};
    my @AminoStatusMSA = @{$amino_status_msa_ref};

    my $num_exons = $MSA[0][$msa_len-2];

    my $OutFile = OpenOutputFile($out_filename);

    $out_filename =~ /\/([^\/]+)\.out$/;
    my $gene = $1;

    $canon_chr = $canon_chr.'[revcomp]' if ($revcomp);

    print $OutFile "Gene Family Name : $gene\n";
    print $OutFile "Canon Chromosome : $canon_chr\n";

    for (my $i=1; $i<$num_seqs; $i++) {

	my $fmt_seq_num = $i;
	while (length("$fmt_seq_num") < length("$num_seqs")) {
	    $fmt_seq_num = ' '.$fmt_seq_num;
	}
	
	print $OutFile "\n---------------------------------------------------------\n";
	print $OutFile "\n>$SeqNames[$i]\n\n";

	my $fmt_exon_num;
	
	for (my $j=0; $j<$msa_len-1; $j++) {

	    if ($MSA[$i][$j] eq '/') {

		if ($fmt_exon_num) {
		    
		    print $OutFile "  //\n";
		    print $OutFile "  // Intron!\n";
		    print $OutFile "  //\n";
		    
		}
		
		$fmt_exon_num = $MSA[0][$j+1];
		while (length("$fmt_exon_num") < length("$num_exons")) {
		    $fmt_exon_num = ' '.$fmt_exon_num;
		}

		next;

	    }

	    my $fmt_col = $j;
	    while (length("$fmt_col") < length("$msa_len")) {
		$fmt_col = ' '.$fmt_col;
	    }

	    print $OutFile "  seq $fmt_seq_num, exon $fmt_exon_num, col $fmt_col: ";

	    if ($MSA[$i][$j] eq '-') {
		print $OutFile "   -       -\n";
		next;
	    }

	    if ($AminoStatusMSA[$i][$j] eq 'S') { print $OutFile "  std   :  "; }
	    else                                { print $OutFile "  ALT   :  "; }

	    print $OutFile "$MSA[$i][$j]  :=  $MapMSA[$i][$j]  (exon:$MapRangesMSA[$i][$j])\n";
	    
	}

	print $OutFile "\n";
	
    }
    
    close($OutFile);

}







# EOF
