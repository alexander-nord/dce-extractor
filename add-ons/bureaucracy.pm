#
#  bureaucracy.pm - Alex Nord, 2022
#
use warnings;
use strict;
use POSIX;
use Time::HiRes;


# GROUP 1: Unavoidables
sub Min;
sub Max;


# GROUP 2: System-related management subroutines
sub RunSystemCommand;
sub OpenSystemCommand;


# GROUP 3: Directory-related management subroutines
sub NameDirectory;
sub OpenDirectory;
sub CreateDirectory;
sub DirectoryExists;     # 1 or 0 (doesn't kill)
sub ConfirmDirectory;    # Returns directory name with terminal '/' (kills if doesn't exist)
sub RemoveDirectory;     # Warning: 'rm -rf' (drops a NUKE!)


# GROUP 4: File-related management subroutines
sub OpenInputFile;
sub NameOutputFile;
sub OpenOutputFile;      # Optional 2nd variable flags overwrite (forces requrested name)
sub AppendToOutputFile;
sub FileExists;          # 1 or 0 (doesn't kill)
sub ConfirmFile;         # Returns filename (kills if doesn't exist)
sub RemoveFile;


# GROUP 5: Timing-related subroutines
sub StartTimer;
sub GetElapsedSeconds;
sub GetYearMonthDayNum;
sub GetYearMonthDayStr;


# GROUP 6: Process parallelization subroutines
sub SpawnProcesses;
sub RegroupProcesses;






#####################
#                   #
#    SUBROUTINES    #
#                   #
#####################






#############################
#                           #
#   GROUP 1: Unavoidables   #
#                           #
#############################






########################################################################
#
#   Function:   Min
#
sub Min
{
	my $val1 = shift;
	my $val2 = shift;
	return $val1 if ($val1 < $val2);
	return $val2;
}







########################################################################
#
#   Function:   Max
#
sub Max
{
	my $val1 = shift;
	my $val2 = shift;
	return $val1 if ($val1 > $val2);
	return $val2;
}







#############################
#                           #
#   GROUP 2: System Calls   #
#                           #
#############################







########################################################################
#
#   Function:   RunSystemCommand
#
sub RunSystemCommand
{
	my $command = shift;
	if (system($command)) {
		die "\n  ERROR:  System command '$command' failed\n\n";
	}
}







########################################################################
#
#   Function:   OpenSystemCommand
#
sub OpenSystemCommand
{
	my $command = shift;

	$command = $command.' |' if ($command !~ /\|\s*$/);
	open(my $CommandOutput,$command)
		|| die "\n  ERROR:  Failed to open output from system command '$command'\n\n";
	
	return $CommandOutput;

}







############################
#                          #
#   GROUP 3: Directories   #
#                          #
############################







########################################################################
#
#   Function:   NameDirectory
#
sub NameDirectory
{
	my $dir_name = shift;

	$dir_name =~ s/([^\\]) /$1\\ /g;
	$dir_name =~ s/\/$//;

	my $check_name = $dir_name;
	my $attempt_num = 1;

	while (DirectoryExists($check_name)) {		

		$attempt_num++;
		$check_name = $dir_name.'-'.$attempt_num;

	}

	return $check_name.'/';

}







########################################################################
#
#   Function:   OpenDirectory
#
sub OpenDirectory
{
	my $dir_name = shift;
	ConfirmDirectory($dir_name);	
	opendir(my $Dir,$dir_name) || 
		die "\n  ERROR:  Failed to open directory '$dir_name'\n\n";
	return $Dir;
}






########################################################################
#
#   Function:   CreateDirectory
#
sub CreateDirectory
{
	my $dir_name = shift;
	my $force_name = shift;

	if ($force_name && DirectoryExists($dir_name)) {
		RemoveDirectory($dir_name);
	}

	$dir_name = NameDirectory($dir_name);
	RunSystemCommand("mkdir \"$dir_name\"");

	return ConfirmDirectory($dir_name);

}







########################################################################
#
#   Function:   DirectoryExists
#
sub DirectoryExists
{
	my $dir_name = shift;
	return 1 if (-d $dir_name);
	return 0;
}







########################################################################
#
#   Function:   ConfirmDirectory
#
sub ConfirmDirectory
{
	my $dir_name = shift;

	$dir_name =~ s/([^\\]) /$1\\ /g;
	$dir_name = $dir_name.'/' if ($dir_name !~ /\/$/);

	if (!DirectoryExists($dir_name)) {
		die "\n  ERROR:  Failed to locate directory '$dir_name'\n\n";
	}

	return $dir_name;
}







########################################################################
#
#   Function:   RemoveDirectory
#
sub RemoveDirectory
{
	my $dir_name = shift;

	if (!DirectoryExists($dir_name)) {
		print "  Warning:  Attempt to delete directory non-existent '$dir_name' cancelled\n";
		return;
	}

	$dir_name = ConfirmDirectory($dir_name);

	RunSystemCommand("rm -rf \"$dir_name\"");

}






######################
#                    #
#   GROUP 4: Files   #
#                    #
######################







########################################################################
#
#   Function:   OpenInputFile
#
sub OpenInputFile
{
	my $file_name = shift;
	$file_name = ConfirmFile($file_name);
	
	open(my $InFile,'<',$file_name)
		|| die "  ERROR:  Failed to open input file '$file_name'\n";
	return $InFile;
}







########################################################################
#
#   Function:   NameOutputFile
#
sub NameOutputFile
{
	my $file_name = shift;

	my $file_base = $file_name;
	my $file_ext = '';

	my $attempt_num = 1;
	my $check_name  = $file_name;

	if ($file_name =~ /^(.+)\.([^\.]+)$/) {

		$file_base = $1;
		$file_ext  = $2;

		while (FileExists($check_name)) {
			$attempt_num++;
			$check_name = $file_base.'-'.$attempt_num.'.'.$file_ext;
		}

	} else {

		while (FileExists($check_name)) {
			$attempt_num++;
			$check_name = $file_name.'-'.$attempt_num;
		}

	}

	return $check_name;

}







########################################################################
#
#   Function:   OpenOutputFile
#
sub OpenOutputFile
{
	my $file_name = shift;
	my $force_name = shift;

	if ($force_name && FileExists($file_name)) {
		RemoveFile($file_name);
	}

	$file_name = NameOutputFile($file_name);

	open(my $OutFile,'>',$file_name)
		|| die "  ERROR:  Failed to open file '$file_name' for writing\n";

	return $OutFile;

}







########################################################################
#
#   Function:   AppendToOutputFile
#
sub AppendToOutputFile
{
	my $file_name = shift;
	open(my $OutFile,'>>',$file_name)
		|| die "  ERROR:  Failed to re-open file '$file_name' for writing\n";
	return $OutFile;
}






########################################################################
#
#   Function:   FileExists
#
sub FileExists
{
	my $filename = shift;
	return 1 if (-e $filename);
	return 0;
}







########################################################################
#
#   Function:   ConfirmFile
#
sub ConfirmFile
{
	my $file_name = shift;
	$file_name =~ s/([^\\]) /$1\\ /g;
	if (!FileExists($file_name)) {
		die "  ERROR:  File '$file_name' does not exist\n";
	}
	return $file_name;
}







########################################################################
#
#   Function:   RemoveFile
#
sub RemoveFile
{
	my $file_name = shift;

	if (!FileExists($file_name)) {
		print "  Warning:  Attempt to delete non-existent file '$file_name' cancelled\n";
		return;
	}

	RunSystemCommand("rm \"$file_name\"");
}








#######################
#                     #
#   GROUP 5: Timing   #
#                     #
#######################









########################################################################
#
#   Function:   StartTimer
#
sub StartTimer
{
    return [Time::HiRes::gettimeofday()];
}







########################################################################
#
#   Function:   GetElapsedSeconds
#
sub GetElapsedSeconds
{
	my $timer = shift;
	return Time::HiRes::tv_interval($timer);
}







########################################################################
#
#   Function:   GetYearMonthDayNum
#
sub GetYearMonthDayNum
{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	return (1900+$year,$mon,$mday);	
}






########################################################################
#
#   Function:   GetYearMonthDayStr
#
sub GetYearMonthDayStr
{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

	$year += 1900;

	my @Months = qw( January February March April May June July August September October November December );
	$mon = $Months[$mon];

	if ($mday == 1 || $mday == 21 || $mday == 31) {
		$mday = $mday.'st';
	} elsif ($mday == 2 || $mday == 22) {
		$mday = $mday.'nd';
	} elsif ($mday == 3 || $mday == 23) {
		$mday = $mday.'rd';
	} else {
		$mday = $mday.'th';
	}

	return $year.', '.$mon.' '.$mday;

}








########################################
#                                      #
#   GROUP 6: Process Parallelization   #
#                                      #
########################################









########################################################################
#
#   Function:   SpawnProcesses
#
sub SpawnProcesses
{
	my $requested_threads = shift;

	my $thread_id = 0;
	my $num_threads = 1;
	while ($num_threads < $requested_threads) {

		if (my $pid = fork) {
            
            if (not defined $pid) { 
            	die "  ERROR:  Fork failed (only $num_threads of $requested_threads spawned successfully)\n";
            }

            $num_threads++;
        
        } else {
        
            $thread_id = $num_threads;
            last;
        
        }

	}

	return $thread_id;
}





########################################################################
#
#   Function:   RegroupProcesses
#
sub RegroupProcesses
{
	my $thread_id = shift;
	
	if ($thread_id) {
		exit(0);
	}

	while (wait() != -1) {}

}





1; # EOF





