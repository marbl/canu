#!/usr/local/bin/perl

use strict;
use Math::BigFloat;
use TIGR::Foundation;
use TIGR::ConfigFile;
use IO::File;
use File::Basename;
use File::Copy;
use File::Find;
use File::Path;
use POSIX qw(uname);
use POSIX qw(:signal_h :errno_h :sys_wait_h);
use sigtrap 'handler', \&SIGHANDLER, 'normal-signals', 'error-signals',
  'old-interface-signals';

my $HELPTEXT = qq~

Retrieve assembly results from the aserver marshalling base.

   aget [options] <request_id> <resdir>

      <request_id>  The request id for the assembly job.

      <resdir>  The directory in user space where the
                assembly results are copied.

      options:

       -mode <mode_num> The mode for copying the files(default mode=2)

       -aconf <aserver_config> Specify a config file for the aserver
                               (by default use the config file in
                               /usr/local/packages/aserver/etc/aserver.conf)

       -L <list_file> Provide a file containing the suffixes of the files
                      and directories to be copied.

       -du Get only the disk space usage of all the files and directories to
           be copied.

   aget copies the assembly output files from the aserver marshalling base
   to the user space. It accepts three different modes(0,1,2) that determine
   the number of files copied. Mode 0 copies all the files(skip no file) and
   mode 1 mode 2 copy a small subset of the files.
~;

# General info
my $VERSION = " Version 1.85 ";

my $cancel_track = 0;
my $resdir       = undef;

my $MAX_MODE = 0;    #copy all the files
my $MED_MODE = 1;    #copy the files with suffixes in %KEEP_FILES1
my $MIN_MODE = 2;    #copy the files with suffixes in %KEEP_FILES2


#the suffix of files to be kept when clean level is 1
my %KEEP_FILES1 = (
	'log'        => 0,
	'out'        => 0,
	'info'       => 0,
	'asm'        => 0,
	'contig'     => 0,
	'fasta'      => 0,
	'feat'       => 0,
	'seqs'       => 0,
	'frg'        => 0,
	'gz'         => 0,
	'qc'         => 0,
	'tasm'       => 0,
	'cgi'        => 0,
	'clone'      => 0,
	'insert'     => 0,
	'frgStore'   => 0,
	'gkpStore'   => 0,
	'ofg'        => 0,
	'ovlStore'   => 0,
	'autoEditor' => 0,
	'upload'     => 0,
	'autoFinish' => 0,
	'scaff'      => 0,
	'cga.0'      => 0,
	'bcp'        => 0,
	'sq.context' => 0,
	'sq.contigs' => 0,
	'9-terminator' => 0,
);

#the suffix of files to be kept when clean level is 2
my %KEEP_FILES2 = (
	'log'        => 0,
	'out'        => 0,
	'info'       => 0,
	'asm'        => 0,
	'contig'     => 0,
	'fasta'      => 0,
	'feat'       => 0,
	'seqs'       => 0,
	'frg'        => 0,
	'gz'         => 0,
	'qc'         => 0,
	'tasm'       => 0,
	'clone'      => 0,
	'insert'     => 0,
	'autoEditor' => 0,
	'upload'     => 0,
	'autoFinish' => 0,
	'scaff'      => 0,
	'cga.0'      => 0,
	'bcp'        => 0,
	'sq.context' => 0,
	'sq.contigs' => 0,
	'9-terminator' => 0,
);

my @MY_DEPENDS = ( "TIGR::Foundation", "TIGR::ConfigFile" );

my $debug               = 0;
my $aserver_config_file = "/usr/local/aserver/etc/aserver.conf";
my $tf_object           = new TIGR::Foundation;
my $basedir             = $tf_object->getProgramInfo("env_path");
$tf_object->setHelpInfo($HELPTEXT);
$tf_object->setVersionInfo($VERSION);
$tf_object->addDependInfo(@MY_DEPENDS);
my $aget_log = "$basedir/aget.log";
$tf_object->setLogFile($aget_log);
my $newDir = 0;
my @createdDirs;
my @copiedFiles;
my $file_size=0;

# this function takes in a string representing a relative or absolute
# path name. It returns the absolute path name for the string.
sub relPaths($) {
	my $path = shift;
	return File::Spec->rel2abs($path);
}

#read the user provided list of files to be kept in the assembly
sub readList($) {
	my $file_list = shift;
	my %copy_hash = ();
	my $rl        = new IO::File $file_list
	  or $tf_object->bail("Cannot open $file_list: $!");
	my $line = undef;

	while ( defined( $line = <$rl> ) ) {
		chomp $line;
		$line =~ s/\s*$//;
		$copy_hash{$line} = 0;
	}
	return \%copy_hash;
}

#Name:   SIGHANDLER
#Input:  none
#Output: none
#Usage:  this functions acts as the signal handler for the external signals
#        that aget might receive. It indicates that an external signal has
#        been received by changing the status field in the Request table to
#        "C" (cancel).

sub SIGHANDLER {

	$tf_object->logLocal( "aget received a cancel request", 1 );
	print("aget received a cancel request\n");
	$cancel_track = 1;
	exit(2);
}

MAIN:
{
	my $copy_mode  = $MIN_MODE;
	my %copy_files = ();
	my $file_list;
	my $disk_space = undef;
	my $result     =
	  $tf_object->TIGR_GetOptions( 'mode=i', \$copy_mode, 'aconf=s',
		\$aserver_config_file, 'L=s', \$file_list, 'du', \$disk_space );

	if ( !defined($result) ) {
		$tf_object->bail("The options could not be read");
	}

	# Update debug level for logging if it is not specified by the user
	if ( !defined $tf_object->getDebugLevel() ) {
		$tf_object->setDebugLevel($debug);
	}

	my $request_id = $ARGV[0];

	$tf_object->bail(
		"Please provide a valid request_id for the " . "assembly job" )
	  if ( !defined $request_id );

	$tf_object->logLocal( "The request id is $request_id", 2 );

	$resdir = $ARGV[1];

	if ( !defined $disk_space ) {
		$tf_object->bail( "Please provide a name for a directory in "
			  . "user space for the assembly output files" )
		  if ( !defined $resdir );
	}

	$tf_object->logLocal( "The output directory in user space is $resdir", 2 );

	my $acf = undef;
	if ( defined $aserver_config_file ) {
		$acf = new TIGR::ConfigFile($aserver_config_file);
	}

	if ( !defined $acf ) {
		$tf_object->bail("Could not initialize config file object");
	}

	my $marshalling_base = $acf->getOption("marshalling_base");
	my $outdir           = "$marshalling_base/$request_id";

	if ( !( -d $outdir ) ) {
		$tf_object->bail( "The output directory $outdir corresponding to "
			  . "request_id $request_id does not exist" );
	}

	if ( !defined $disk_space ) {
		my $temp_resdir = $resdir;
		$resdir = relPaths($temp_resdir);

		if ( !-e $resdir ) {
			$tf_object->logLocal(
				"The directory $resdir does not "
				  . "exist in user space. Making dir...",
				2
			);
			if ( !mkdir($resdir) ) {
				$tf_object->bail( "Failed to create results directory "
					  . "\'$resdir\' ($!)" );
			}
			$newDir = 1;
		}

		if ( !isWritableDir($resdir) ) {
			$tf_object->bail( "The output directory \'$resdir\' does "
				  . "not have write permissions" );
		}

		$tf_object->logLocal( "Copying the intermediate files", 3 );
	}

	if ( defined $file_list ) {
		my $copy_hash = readList($file_list);
		%copy_files = %$copy_hash;
		$copy_mode = 3;  #Set to other than MED, MAX, MIN mode
	}

	if ( $copy_mode == $MED_MODE ) {
		%copy_files = %KEEP_FILES1;
	}
	elsif ( $copy_mode == $MIN_MODE ) {
		%copy_files = %KEEP_FILES2;
	}
	my $error_check = undef;

	if ( !chdir $outdir ) {
		$tf_object->bail("Failed to cd to $outdir");
	}

	my @files = `ls -1 $outdir`;
	my $dir_size_;
	my $sizeOfFile_;
	my $outdir_mbs=0;
	foreach my $file (@files) {
		chomp $file;

		my @name_parts = split( /\./, $file );
		my $suffix     = $name_parts[ scalar(@name_parts) - 1 ];
		my $suffix1    = $file;
		if ( $file =~ /\./ig ) { $suffix1 = $'; }

		if (   ( defined( $copy_files{$suffix} ) )
			|| ( defined( $copy_files{$suffix1} ) )
			|| ( $copy_mode == $MAX_MODE ) )
		{    #copy the files

			if ($disk_space) {
			    $outdir_mbs += dirSize ($outdir, $file);
				find(\&fileSize,"$file");
			}
			else {
				find( \&filecopy, "$file" );
				if ( $error_check == 1 ) {
					removeFailedCopyRemnants( $resdir, $tf_object, $newDir,
						\@createdDirs, \@copiedFiles );
				}
			}
		}
	}
	print("File $outdir takes $outdir_mbs MB\n");

	if ($disk_space) {
		my $mbs=Math::BigFloat->bceil(($file_size/(1024*1024)));
		print("Total Mega bytes taken = $mbs MB\n");
	}


	sub fileSize {
		my $filename = $File::Find::name;
		my $sizeOfFile=0;
		my $file_;
		my $dir_size=0;
		if( -l "$outdir/$filename" ){
			$sizeOfFile = (lstat("$outdir/$filename"))[7];
			$file_size += $sizeOfFile;
		}
		elsif ( -d "$outdir/$filename" ) {
			opendir(BIN, "$outdir/$filename") or die "Can't open $outdir/$filename: $!";
			while( defined ($file_ = readdir BIN) ) {
			    #$dir_size += dirSize("$outdir/$filename",$file_);
				if( -l "$outdir/$filename/$file_" ){
					$sizeOfFile = (lstat("$outdir/$filename/$file_"))[7];
					$dir_size += $sizeOfFile;
				}
				else{
					$sizeOfFile = (-s "$outdir/$filename/$file_");
					$dir_size += $sizeOfFile;
				}

			}
			closedir(BIN);
			my $f = ($dir_size/(1024*1024));
			my $dir_mbs=Math::BigFloat->new("$f");
			$dir_mbs = $dir_mbs->bround(2);
			print("File $outdir/$filename takes $dir_mbs MB\n");
		}
		else{
			$sizeOfFile = (-s $_);
			$file_size += $sizeOfFile;
		}

	}

	sub filecopy {
		my $filename = $File::Find::name;

		if ( -d "$outdir/$filename" ) {
			if ( ( !-e "$resdir/$filename" ) && ( !mkdir "$resdir/$filename" ) )
			{
				$tf_object->logError("Could not create $resdir/$filename");
				$error_check = 1;
				return;
			}
			push @createdDirs, "$resdir/$filename";
		}
		else {

			# test to see if it is a symlink and use symlink() instead of copy(): Added by Prabhu - 01/17/07
			# this is the skip-no-files mode(0) where symlinks are copied as symlinks, the rest as real files
			if ( ( $copy_mode == $MAX_MODE ) && ( -l "$outdir/$filename" ) ) {
				my $link_read = readlink "$outdir/$filename";
				if ( index( $link_read, "/" ) == 0 ) {    ## absolute symlink
					my $abs2rel          = File::Spec->abs2rel($link_read);
					my $sym_link_success = symlink $abs2rel,"$resdir/$filename";
					testSymLinkSuccess( $sym_link_success, $resdir, $outdir,$filename );
				}
				else {
					my $sym_link_pass = symlink $link_read, "$resdir/$filename";
					testSymLinkSuccess( $sym_link_pass, $resdir, $outdir,$filename );
				}
				return;
			}
			#Copy and preserve file attributes (-p)
			if ( $tf_object->runCommand( "cp -p $outdir/$filename $resdir/$filename" ) ) {
				$tf_object->logError(
					    "Could not copy the file $outdir/$filename "
					  . "to $resdir/$filename" );
				print(  "Could not copy the file $outdir/$filename "
					  . "to $resdir/$filename\n" );
				$error_check = 1;
				return;
			}
			push @copiedFiles, "$resdir/$filename";
		}
	}

}

END {
	if ( ( $cancel_track == 1 ) && ( defined $resdir ) ) {
		if ( -e $resdir ) {
			removeFailedCopyRemnants( $resdir, $tf_object, $newDir,\@createdDirs,
\@copiedFiles );
		}
	}
}

sub removeFailedCopyRemnants {
	my $resdir      = shift;
	my $tf_object   = shift;
	my $newDir      = shift;
	my $createdDirs = shift;
	my $copiedFiles = shift;

	if ($newDir) {
		$tf_object->logLocal( "Removing $resdir", 2 );
		my $rm_cmd = "rm -r $resdir";
		my $bad    = $tf_object->runCommand($rm_cmd);
		$tf_object->bail( "Could not remove $resdir after " . "copy failed" )
		  if ($bad);
	}
	else {
		for my $tmp_file (@$copiedFiles) {
			my $rm_cmd = "rm $tmp_file";
			my $bad    = $tf_object->runCommand($rm_cmd);
			$tf_object->bail(
				"Could not remove $tmp_file after " . "copy failed" )
			  if ($bad);

		}
		for my $tmp_dir (@$createdDirs) {

			my $rm_cmd = "rm -r $tmp_dir";
			my $bad    = $tf_object->runCommand($rm_cmd);
			$tf_object->bail(
				"Could not remove $tmp_dir after " . "copy failed" )
			  if ($bad);
		}

	}
}

sub testSymLinkSuccess {
	my $sym_link_success = shift;
	my $resdir           = shift;
	my $outdir           = shift;
	my $filename         = shift;
	my $error_check;

	if ($sym_link_success) {
		push @copiedFiles, "$resdir/$filename";
	}
	else {
		$tf_object->logError(
			    "Could not copy the symlink file $outdir/$filename "
			  . "to $resdir/$filename" );
		print(  "Could not copy the symlink file $outdir/$filename "
			  . "to $resdir/$filename\n" );
		$error_check = 1;

	}
}





sub dirSize {
	my $outdir      = shift;
	my $file_ = shift;
	my $dir_size_ ;
	my $sizeOfFile_;

	if ( -l "$outdir/$file_" ) {
		$sizeOfFile_ = ( lstat("$outdir/$file_") )[7];
		$dir_size_ += $sizeOfFile_;
	}
	else {
		$sizeOfFile_ = ( -s "$outdir/$file_" );
		$dir_size_ += $sizeOfFile_;
	}
	my $f = ( $dir_size_ / ( 1024 * 1024 ) );
	my $dir_mbs = Math::BigFloat->new("$f");
	$dir_mbs = $dir_mbs->bround(2);

	return $dir_mbs;
}
