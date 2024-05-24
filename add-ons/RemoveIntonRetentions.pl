#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

sub ScriptDir { return './' if ($0 !~ /\//); $0 =~ /^(.*\/)[^\/]+$/; return $1; }
use lib ScriptDir();
use bureaucracy;


if (@ARGV != 1) {
    die "\n  USAGE: ./RemoveIntronRetentionsFromCSV.pl [Extracted-DCEs/]\n\n";
}


my $in_dir_name = ConfirmDirectory($ARGV[0]);

my $ire_file_name = $in_dir_name.'Intron-Retention-Events.out';
if (!(-e $ire_file_name)) 
{
    print "\n  No intron retention events reported (couldn't find file '$ire_file_name')\n\n";
    exit(0);
}


my $IREFile = OpenInputFile($ire_file_name);
my %Species;
while (my $line = <$IREFile>)
{
    if ($line =~ /^(\S+):/)
    {
        $Species{$1} = 1;
    }
}
close($IREFile);


foreach my $species (keys %Species)
{
    
    $IREFile = OpenInputFile($ire_file_name);
    my %SpeciesIREs;
    while (my $line = <$IREFile>)
    {
        if ($line =~ /$species: (\d+) \((\S+)\)/)
        {
            $SpeciesIREs{$1.'/'.$2} = 1;
        }
    }
    close($IREFile);


    my $species_file_name = $in_dir_name.$species.'.csv';
    my $species_no_ires_file_name = $species_file_name;
    $species_no_ires_file_name =~ s/\.csv$/\.IREs-removed\.csv/;


    my $SpeciesInFile = OpenInputFile($species_file_name);
    my $SpeciesNoIREs = OpenOutputFile($species_no_ires_file_name);

    my $header = <$SpeciesInFile>;
    print $SpeciesNoIREs "$header";

    while (my $line = <$SpeciesInFile>)
    {
        if ($line =~ /^(\d+), (\S+),/)
        {
            if (!$SpeciesIREs{$1.'/'.$2})
            {
                print $SpeciesNoIREs "$line";
            }
        }
    }
    close($SpeciesInFile);
    close($SpeciesNoIREs);


    $species_file_name = $in_dir_name.$species.'.fa';
    $species_no_ires_file_name = $species_file_name;
    $species_no_ires_file_name =~ s/\.fa$/\.IREs-removed\.fa/;


    $SpeciesInFile = OpenInputFile($species_file_name);
    $SpeciesNoIREs = OpenOutputFile($species_no_ires_file_name);
    while (my $line = <$SpeciesInFile>)
    {
        if ($line =~ /^>(\d+)_([^_]+)_/)
        {
            if ($SpeciesIREs{$1.'/'.$2})
            {
                <$SpeciesInFile>;
            }
            else
            {
                print $SpeciesNoIREs "$line";
                $line = <$SpeciesInFile>;
                print $SpeciesNoIREs "$line";
            }
        }
    }
    close($SpeciesInFile);
    close($SpeciesNoIREs);

}


1;













# EOF
