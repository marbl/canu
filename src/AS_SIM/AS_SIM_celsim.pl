#!/usr/local/bin/perl
#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#
eval 'exec /usr/opt/PERL5004/bin/perl -S $0 ${1+"$@"}'
    if $running_under_some_shell;
			# this emulates #! processing on NIH machines.
			# (remove #! line above if indigestible)

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;
			# process any FOO=bar switches

use Math::BigInt;

use FindBin;
my $bin = "$FindBin::Bin";

$[ = 1;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator
$| = 1;                 # set autoflush

$notdone = 1;
srand(join($,,@Fields));
$keepDNA = 0;
# Don't nuke the dna file
$fileComment = '(none)';
$sawfileseed = 0;

# Have we seen an (optional) .seed section?
$sawcmdlineseed = 0;
$dna_length = -1;

# Have we seen a command line seed
$seed = int(1000 * rand(1));
$notations = 0;
$comment = 1;
$fastaOutput = 0;
$uniform = 0;
$base = 1;
$prompt = 0;
$protoflag = '';

# the initial acnum can be overridden at the command line
# $acnum = Math::BigInt->new("10 000 000 000 000 000");
# this is the UID of the first drosophila read listed in the first input file
$acnum = 17000001585819;

# this is the time stamp of the batch message in the first drosophila file
# it may not be historically accurate...
$etm = 915170460;

$ubac_type = 1;
$fbac_type = 2;
$lbac_type = 3;

# the folowing arrays are used in add_bactig
$acgt{0} = 'a';
$acgt{1} = 'c';
$acgt{2} = 'g';
$acgt{3} = 't';

$base_index{a} = 0;
$base_index{c} = 1;
$base_index{g} = 2;
$base_index{t} = 3;

# $dummy is justa device to get around CVS mechanisms.  If there is just a $ then
# Perl complains.  If you put \$ and check-in CVS removes the \.  So $dummy fools
# them both.
@revFields = split(' ', "\$Revision: 1.5 $dummy", 9999);
$revision = $revFields[2];
@dateFields = split(' ', "\$Date: 2006-10-16 03:59:30 $dummy", 9999);
$date = $dateFields[2];
print $STDERR . 'date = ' . $date . ' revision = ' . $revision;
while ($base < ($#ARGV+1) - 1) {
    if ($ARGV[$base] =~ '-F') {
	$comment = 0;
    }
    elsif ($ARGV[$base] =~ '-s' && $base + 1 < ($#ARGV+1) - 1) {
	$cmdlineseed = $ARGV[$base + 1];
	$sawcmdlineseed = 1;
	$base += 1;
    }
    elsif ($ARGV[$base] =~ '-a') {
	print "* Setting initial acnum to $ARGV[$base + 1]\n";
	$acnum = $ARGV[$base + 1];
	$base += 1;
    }
    elsif ($ARGV[$base] =~ '-d') {
	print "* Saving DNA File\n";
	$keepDNA = 1;
    }
    elsif ($ARGV[$base] =~ '-p') {
	print "* Will Prompt before clobbering files\n";
	$prompt = 1;
    }
    elsif ($ARGV[$base] =~ '-u') {
	print "* Uniform Distance Distribution (Normal is default)\n";
	$uniform = 1;
    }
    elsif ($ARGV[$base] =~ '-f') {
	print "* FASTA output (Assembler ProtoI/O is default)\n";
	$fastaOutput = 1;
	$protoflag = ' ';
    }
    elsif ($ARGV[$base] =~ '-n') {
	print "* Normal Distance Distribution (default)\n";
	$uniform = 0;
    }
    elsif ($ARGV[$base] =~ '-b') { # && $base + 1 < ($#ARGV+1) - 1) {
	# $cmdlinebatchsize = $ARGV[$base + 1];
	$sawcmdlinebatchflag = 1;
	# $base += 1;
    }
    else {
	print " *** Illegal argument '%s'", $ARGV[$base];
	$notdone = 0;
	$ExitValue = (1);
	goto line;  # MJF was: last line
    }
    $base += 1;
}
if ($base != ($#ARGV+1) - 1) {
    die ' *** Usage: celsim [-F] [-s #] [-d] [-p] [-u] [-f] [-n] [-b <batch size>] file\n';
}

$pfile = $ARGV[1] = $ARGV[($#ARGV+1) - 1];
$#ARGV = 1;  # MJF ($#ARGV+1) = 2;
$numFields = (@fields = split(/\./, $pfile, 9999));
if ($numFields == 1) {
    # No Extension
    $tryfile = $pfile . '.sim';
    $prefix = $pfile;
}
else {
    # print "file extension is " fields[numFields];
    #      if("spec" != fields[numFields] &&
    # "sim" != fields[numFields]){
    $n;#	print "File extension " fields[numFields] " is invalid -- try .sim\n";
    #	exit(1);
    #      }else{

    $prefix = '';
    for ($i = 1; $i < $numFields; $i++) {
	$prefix = $prefix . $fields[$i];
    }
}

$frgfileCount = 1;
$frgfile = sprintf "%s_A%05.d.frg", $prefix, $frgfileCount;
$outfile = $prefix . '.cms';
$notfile = $prefix . '.not.' . $seed;
$adlfile = $prefix . '.adl';
$comfile = $prefix . '.cms.' . $seed;
$qltfile = $prefix . '.qlt.' . $seed;
#CMM    dnafile = ".dna." seed;
$dnafile = $prefix . '.dna';
$bactig_file = $prefix . '.bacs';

$phase = -1;
$npoly = 0;
$ptotl = 0;
$nlibs = 0;

line: while (<>) {

    chop;	# strip record separator
    @Fields = split(' ', $_, 9999);

    if (/^\.seed/ && $notdone) {
	if ($phase != -1) {
	    print ' *** Seed specification must come first.';
	    &cleanup(1);
	}
	# Do not change phase...we are still before .dna
	$sawfileseed = 1;

	# read in the seed right here
	if ($phase == -1 && $sawfileseed == 1) {
	    $_ = <>;
	    chop;
	    $fileseed = $_;
	    $sawfileseed = 2;
	    # Read only one seed in this section
	    # print "Seed = " seed "\n";
	}
	goto line;
    }

    if (/^\.comment/ && $notdone) {
	$fileComment = '';
	for ($i = 2; $i <= $#Fields; $i++) {
	    $fileComment = $fileComment . ' ' . $Fields[$i];
	}
	print 'FileComment' . $fileComment;
	goto line;
    }

    if (/^\.dna/ && $notdone) {
	if ($phase != -1) {
	    print ' *** Dna specification must come first.';
	    &cleanup(1);
	}
	$phase = 0;
	# Cmd line uber alles
	if ($sawcmdlineseed == 1) {
	    $seed = $cmdlineseed;
	    print 'Command line seed won: seed = ' . $seed;
	    if ($cmdlineseed == 0) {
		srand( time() ^ ($$ + ($$ << 15)));
	    }
	}
	elsif ($sawfileseed == 2) {
	    #file .seed section over rand
	    $seed = $fileseed;
	    print 'File seed won: seed = ' . $seed;
	}
	goto line;
    }

    if (/^\.poly([\t ]+[0-9.]+)+$/ && $notdone) {
	if ($phase == -1) {
	    print ' *** Missing dna specification.';
	    &cleanup(1);
	}
	if ($phase == 2) {
	    print ' *** Poly specification after sample spec.';
	    &cleanup(1);
	}
	if ($phase == 0) {
	    &make_dna();
	}
	else {
	    &make_polys($pcopy);
	}
	$phase = 1;
	for ($i = 2; $i <= $#Fields; $i++) {
	    $pweight{$npoly + ($i - 1)} = $Fields[$i];
	}
	$pcopy = $#Fields - 1;
	$lines = 0;
	goto line;
    }

    if (/^\.sample$/ && $notdone) {
	$bacends = 0;
	if ($phase == -1) {
	    print ' *** Missing dna and polymorphism specifications.';
	    &cleanup(1);
	}
	if ($phase == 0) {
	    &make_dna();
	    $npoly = 1;
	    $ptotl = 1;
	    $pweight{1} = 1;
	    if ($notations > 0) {
		delete $opened{$notfile} && close($notfile);
	    }
	    &Pick('>>', $outfile) &&
		(print $fh '# NO POLYMORPHISMS APPLIED, USED BASE SEQUENCE');
	    &Pick('>>', $outfile) &&
		(print $fh '#');
	    $cmd = 'cp ' . $dnafile . ' ' . $prefix . '.poly.1';
	    # print 'Executing: ' . $cmd;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 1: " . $cmd;
	    system($cmd);
	    &read_sample_section();
	    &make_sample();
	}
	elsif ($phase == 1) {
	    &make_polys($pcopy);
	}
	else {
	    print '* Invoking make_sample from .sample';
	    print "\$fragspec: $fragspec";
	    &read_sample_section();
	    &make_sample();
	}
	$phase = 2;
	$items = 0;
	goto line;
    }

    if (/^\.bacends$/ && $notdone) {
	$bacends = 1;
	if ($phase == -1) {
	    print ' *** Missing dna and polymorphism specifications.';
	    &cleanup(1);
	}
	if ($phase == 0) {
	    &make_dna();
	    $npoly = 1;
	    $ptotl = 1;
	    $pweight{1} = 1;
	    if ($notations > 0) {
		delete $opened{$notfile} && close($notfile);
	    }
	    &Pick('>>', $outfile) &&
		(print $fh '# NO POLYMORPHISMS APPLIED, USED BASE SEQUENCE');
	    &Pick('>>', $outfile) &&
		(print $fh '#');
	    $cmd = 'cp ' . $dnafile . ' ' . $prefix . '.poly.1';
	    # print 'Executing: ' . $cmd;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 2: " . $cmd;
	    system($cmd);
	    &read_sample_section();
	    &make_sample();
	}
	elsif ($phase == 1) {
	    &make_polys($pcopy);
	}
	else {
	    # print "* Invoking make_sample from .bacends";
	    &read_sample_section();
	    &make_sample();
	}
	$phase = 2;
	$items = 0;
	goto line;
    }

    if (/^\.lbac$/ && $notdone) {
	print "in .lbac, \$phase = ", $phase;
	$lbac = 1;
	if ($phase == -1) {
	    print ' *** Missing dna and polymorphism specifications.';
	    &cleanup(1);
	}
	if ($phase == 0) {
	    &make_dna();
	    $npoly = 1;
	    $ptotl = 1;
	    $pweight{1} = 1;
	    if ($notations > 0) {
		delete $opened{$notfile} && close($notfile);
	    }
	    &Pick('>>', $outfile) &&
		(print $fh '# NO POLYMORPHISMS APPLIED, USED BASE SEQUENCE');
	    &Pick('>>', $outfile) &&
		(print $fh '#');
	    $cmd = 'cp ' . $dnafile . ' ' . $prefix . '.poly.1';
	    # print 'Executing: ' . $cmd;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 2: " . $cmd;
	    system($cmd);
	    &read_lbac_section();
	    $shred = 0;
	    &make_uflbacs($lbac_type, $shred, $min_num_lbacs, $max_num_lbacs, $min_lbac_length, $max_lbac_length,
			  $sequencing_error, $insert_error_rate, $delete_error_rate);
	}
	elsif ($phase == 1) {
	    &make_polys($pcopy);
	}
	else {
	    # print "* Invoking make_sample from .lbac";
	    &read_lbac_section();
	    $shred = 0;
	    &make_uflbacs($lbac_type, $shred, $min_num_lbacs, $max_num_lbacs, $min_lbac_length, $max_lbac_length,
			  $sequencing_error, $insert_error_rate, $delete_error_rate);
	}
	$phase = 2;
	$items = 0;
	goto line;
    }

    if (/^\.fbac$/ && $notdone) {
	if ($phase == -1) {
	    print ' *** Missing dna and polymorphism specifications.';
	    &cleanup(1);
	}
	if ($phase == 0) {
	    &make_dna();
	    $npoly = 1;
	    $ptotl = 1;
	    $pweight{1} = 1;
	    if ($notations > 0) {
		delete $opened{$notfile} && close($notfile);
	    }
	    &Pick('>>', $outfile) &&
		(print $fh '# NO POLYMORPHISMS APPLIED, USED BASE SEQUENCE');
	    &Pick('>>', $outfile) &&
		(print $fh '#');
	    $cmd = 'cp ' . $dnafile . ' ' . $prefix . '.poly.1';
	    # print 'Executing: ' . $cmd;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 2: " . $cmd;
	    system($cmd);
	    &read_fbac_section();
	    &make_uflbacs($fbac_type, $shred, $min_num_fbacs, $max_num_fbacs, $min_fbac_length, $max_fbac_length,
			  $sequencing_error, $insert_error_rate, $delete_error_rate);
	}
	elsif ($phase == 1) {
	    &make_polys($pcopy);
	}
	else {
	    # print "* Invoking make_sample from .fbac";
	    &read_fbac_section();
	    &make_uflbacs($fbac_type, $shred, $min_num_fbacs, $max_num_fbacs, $min_fbac_length, $max_fbac_length,
			  $sequencing_error, $insert_error_rate, $delete_error_rate);
	}
	$phase = 2;
	$items = 0;
	goto line;
    }

    if (/^\.ubac$/ && $notdone) {
	if ($phase == -1) {
	    print ' *** Missing dna and polymorphism specifications.';
	    &cleanup(1);
	}
	if ($phase == 0) {
	    &make_dna();
	    $npoly = 1;
	    $ptotl = 1;
	    $pweight{1} = 1;
	    if ($notations > 0) {
		delete $opened{$notfile} && close($notfile);
	    }
	    &Pick('>>', $outfile) &&
		(print $fh '# NO POLYMORPHISMS APPLIED, USED BASE SEQUENCE');
	    &Pick('>>', $outfile) &&
		(print $fh '#');
	    $cmd = 'cp ' . $dnafile . ' ' . $prefix . '.poly.1';
	    # print 'Executing: ' . $cmd;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 2: " . $cmd;
	    system($cmd);
	    &read_ubac_section();
	    &make_uflbacs($ubac_type, $shred, $min_num_ubacs, $max_num_ubacs, $min_ubac_length, $max_ubac_length,
			  $sequencing_error, $insert_error_rate, $delete_error_rate);
	}
	elsif ($phase == 1) {
	    &make_polys($pcopy);
	}
	else {
	    # print "* Invoking make_sample from .ubac";
	    &read_ubac_section();
	    print "after reading ubac section, \$acnum = $acnum";
	    &make_uflbacs($ubac_type, $shred, $min_num_ubacs, $max_num_ubacs, $min_ubac_length, $max_ubac_length,
			  $sequencing_error, $insert_error_rate, $delete_error_rate);
	}
	$phase = 2;
	$items = 0;
	goto line;
    }

    if (/^\.length$/) 
    {
	$dna_length = <>;
	print "read \$dna_length from file: $dna_length";
    }

    &readbody();
}

if ($notdone) {
    if ($phase != 2) {
	print ' *** Generated no data! (phase = ' . $phase . ')';
	print ' *** Generated no data! (phase = ' . $Fields[$phase] . ')';
	&cleanup(1);
    }
    # &make_sample();
    if ($npoly == 1 && $notations > 0 && !$fastaOutput) {
	delete $opened{$frgfile} && close($frgfile);
	# close the file, so we can annotate it

	# do the header file first
	{
	    $frgfile = sprintf "%s_A%05.d.frg", $prefix, 1;
	    $perlcmd = "perl $bin/AS_SIM_labelfrag.pl " . $outfile . ' '
		. $frgfile . ' ' . $notfile;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 3: " . $perlcmd;
	    if (system($perlcmd) != 0) {
		&Pick('>>', $outfile) &&
		    (print $fh "*** Error in fragment labeling\n");
		&cleanup(1);
		exit (1);
	    }
	}

	# now do everybody past the header
	for ($frgfileNum = 2; $frgfileNum < $frgfileCount; $frgfileNum++)
	{
	    $frgfile = sprintf "%s_B%05.d.frg", $prefix, $frgfileNum;
	    $perlcmd = "perl $bin/AS_SIM_labelfrag.pl " . $outfile . ' '
		. $frgfile . ' ' . $notfile;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 3: " . $perlcmd;
	    if (system($perlcmd) != 0) {
		&Pick('>>', $outfile) &&
		    (print $fh "*** Error in fragment labeling\n");
		&cleanup(1);
		exit (1);
	    }
	}
    }
    elsif ($notations > 0) {
	if (!$fastaOutput) {
	    print "*** Poly > 1 -- No fragment labelling applied\n";
	}
    }
    
    # if (!$fastaOutput && $sawcmdlinebatchflag) {
	# print "\n %%%%%%%%%%%%%%%%%%%%% command 4: " . 'tobatches -b ' . $cmdlinebatchsize . ' ' . $frgfile;
	# system('tobatches -b ' . $cmdlinebatchsize . ' ' . $frgfile);
    # }

    &cleanup(1);
}   
exit $ExitValue;    


sub readbody
{
    # print "in readbody 1, \$_: ", $_, "| notdone: ", $notdone;

    if (/^([^.]|\.[0-9])/ && $notdone) {
	# print "in readbody 2, \$_: ", $_, "notdone: ", $notdone;

	# Process a section body .......
	if ($phase == 0) {
	    if ($Fields[1] eq '@') {
		$notations++;
		# track number of notations
		&Pick('>', $notfile) &&
		    (print $fh $Fields[2]);
		$Fields[1] = '';
		# Strip the @ sign, so frag doesn't see it
	    }
	    if ($Fields[1] eq 'Z') {
		$dna_length = $Fields[3];
	    }
	    &Pick('>', $comfile) &&
		(print $fh join($,,@Fields));
	}
	elsif ($phase == 1) {
	    $lines += 1;
	    if ($phase == 1 && $lines == 1) {
		&Pick('>', $comfile) &&
		    (print $fh $dnafile);
	    }
	    &Pick('>', $comfile) &&
		(print $fh join($,,@Fields));
	}
	elsif ($phase == 2) {
	    if ($items == 0 && $#Fields >= 1) {
		$j = 2;
		$fnum = $Fields[1];
		$items++;
		#	  print STDERR "items " items " fnum = " fnum ;
	    }
	    else {
		$j = 1;
	    }
	    for ($i = $j; $i <= $#Fields; $i++) {
		#	  print STDERR "i = " i " NF = " NF " items = " items " text = " $i;
		if ($items < 4 || $items > 7) {
		    $fragspec = $fragspec . ' ' . $Fields[$i];
		}
		else {
		    $massageSpec = $massageSpec . ' ' . $Fields[$i];
		}
		$items++;
	    }
	}
    }

    # specifications
    # massageSpec
    #   1  min prefix unclear range
    #   2  max prefix unclear range
    #   3  min suffix unclear range
    #   4  max suffix unclear range

    # fragspec gets split as follows:
    #   1  min frag len
    #   2  max frag len
    #   3  F/R odds
    #   4  Error ramp Min
    #   5  Error ramp Max
    #   6  Ins odds - if there is an error, the odds it's an insertion
    #   7  Del odds - if there is an error, the odds it's a deletion
    #   8  Single Odds - the odds of a frag not having a mate pair
    #   9  min insert length
    #   10 max insert length
    #   11 False mate rate
    #   Note: 1.0 - "Ins odds" - "Del odds" is implicitly taken as the odds of a substitution
}

# read_sample_section is used for both samples and ebacs
sub read_sample_section {

    my $icnt;

    $#fragspec = 0;
    $#massageSpec = 0;
    
    $icnt = 0;
    while ($icnt < 5)
    {
	$_ = <>;
	chop;	# strip record separator
	@Fields = split(' ', $_, 9999);
	# print @Fields;
	if ($Fields[0] ne '#' && substr($Fields[0], 1, 1) ne '#')  # skip comments
	{
	    if ($icnt == 0) {
		$fnum = $Fields[1];
	    }
	    elsif ($icnt == 2) {
		@massageSpec = @Fields;
	    }
	    else {
		push (@fragspec, @Fields);
	    }
	    $icnt++;
	}
    }
    $items = 13;  # this is needed later in make_sample to trigger things
}

sub read_lbac_section {
    my $icnt;

    $#fragspec = 0;
    $#massageSpec = 0;
    
    $icnt = 0;
    while ($icnt < 6)
    {
	$_ = <>;
	chop;	# strip record separator
	@Fields = split(' ', $_, 9999);
	# print @Fields;
	if ($Fields[0] ne '#' && substr($Fields[0], 1, 1) ne '#')  # skip comments
	{
	    if ($icnt == 0) {
		$lbac_coverage = $Fields[1];
	    }
	    elsif ($icnt == 2) {
		@massageSpec = @Fields;
	    }
	    elsif ($icnt == 5) {
		($min_num_lbacs, $max_num_lbacs, $min_lbac_length, $max_lbac_length) = @Fields;
	    }
	    else {
		push (@fragspec, @Fields);
	    }
	    $icnt++;
	}
    }
    # these perl variables are not relevant to lbacs, since frag is called on the lbacs
    $insert_error_rate = 0;
    $delete_error_rate = 0;
    $sequencing_error = 0;

    $items = 13;  # this is needed later in make_sample to trigger things
}

sub read_fbac_section {

    my $icnt;

    $#fragspec = 0;
    $#massageSpec = 0;

    $icnt = 0;
    while ($icnt < 3)
    {
	$_ = <>;
	chop;	# strip record separator
	@Fields = split(' ', $_, 9999);
	print @Fields;

	if ($Fields[0] ne '#' && substr($Fields[0], 1, 1) ne '#')  # skip comments
	{
	    if ($icnt == 0) {
		$shred = $Fields[0];
	    }
	    elsif ($icnt == 1) {
		($min_num_fbacs, $max_num_fbacs, $min_fbac_length, $max_fbac_length,
		 $sequencing_error, $insert_error_rate, $delete_error_rate) = @Fields;
	    }
	    else {
		@massageSpec = @Fields;
	    }
	    $icnt++;
	}
    }
    $items = 13;  # this is needed later in make_fbacs to trigger things
}

sub read_ubac_section {

    my $icnt;

    # $#massageSpec = 0;
    # @massageSpec = join ' ', 0, 0, 0, 0;

    $icnt = 0;
    while ($icnt < 3)
    {
	$_ = <>;
	chop;	# strip record separator
	@Fields = split(' ', $_, 9999);
	print @Fields;

	if ($Fields[0] ne '#' && substr($Fields[0], 1, 1) ne '#')  # skip comments
	{
	    if ($icnt == 0) {
		($shred, $ubac_coverage) = @Fields;
	    }
	    elsif ($icnt == 1) {
		($min_num_ubacs, $max_num_ubacs, $min_ubac_length, $max_ubac_length,
		 $sequencing_error, $insert_error_rate, $delete_error_rate) = @Fields;
	    }
	    else {
		@massageSpec = @Fields;
	    }
	    $icnt++;
	}
    }

    if ($dna_length == -1)
    {
	print "\$dna_length not defined (did you forget to say Z = dna_length?)";
	exit;
    }	
    if ($min_ubac_length > $dna_length) {
	print "$min_ubac_length must be greater than $dna_length";
	exit;
    }
    if ($max_ubac_length > $dna_length) {
	$max_ubac_length = $dna_length;
    }
}

sub make_dna {
    delete $opened{$comfile} && close($comfile);
    print ' +++ Building DNA sequence';
    print "command: $bin/AS_SIM_frag -s " . $seed . $protoflag . ' ' . $comfile . ' >' . $dnafile;
    print "\n %%%%%%%%%%%%%%%%%%%%% command 5: " . "$bin/AS_SIM_frag -s " . $seed . $protoflag . ' ' . $comfile  . ' >' . $dnafile;
    if (system("$bin/AS_SIM_frag -s " . $seed . $protoflag . ' ' . $comfile . ' >' . $dnafile) != 0) {
	print ' *** Dna spec error';
	&cleanup(1);
    }
    if ($uniform) {
	$uni = ' (uniform) ';
    }
    else {
	$uni = ' (gaussian) ';
    }
    if (!$fastaOutput) {
	&Pick('>', $adlfile) &&
	    (print $fh 'Celsim ' . $revision . $uni . ' ' . $date);
	&Pick('>>', $adlfile) &&
	    (print $fh $fileComment);
	close($adlfile);  # force the write of the previous two lines
	$perlcmd = "perl $bin/AS_SIM_extractLength.pl  < " . $dnafile
	  . ' >> ' . $adlfile;
	# print STDERR perlcmd;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 6: " . $perlcmd;
	if (system($perlcmd) != 0) {
	    print ' *** Error in extractLength.pl';
	    &cleanup(1);
	}
	# print "*** Celsim Input was:\n" >> adlfile;
	$perlcmd = "echo \" ********Celsim Input ******** \">> " . $adlfile;

	#      print STDERR perlcmd;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 7: " . $perlcmd;
	if (system($perlcmd) != 0) {
	    print $STDERR . ' *** Error outputing comment';
	    &cleanup(1);
	}
	$perlcmd = 'cat ' . $pfile . ' >>' . $adlfile;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 8: " . $perlcmd;
	if (system($perlcmd) != 0) {
	    print $STDERR . ' *** Error in cat input to adlfile';
	    &cleanup(1);
	}
	$perlcmd = "echo \" ********End Celsim Input ******** \">> " .
	  $adlfile;

	print "\n %%%%%%%%%%%%%%%%%%%%% command 9: " . $perlcmd;
	if (system($perlcmd) != 0) {
	    print $STDERR . ' *** Error outputting comment';
	    &cleanup(1);
	}
    }

    &write_param_file();

    if ($comment) {
	if (&overwrite($outfile)) {
	    return;
	}
	&Pick('>>', $outfile) &&
	    (print $fh '# Celsim V' . $revision . ' ' . $fileComment . ' ');
# solves a buffering problems
	select $fh; $| = 1; select(STDOUT);   
	&Pick('>>', $outfile) &&
	    (print $fh '# DNA SEQUENCE GENERATION');
	print "\n %%%%%%%%%%%%%%%%%%%%% command 10: " . "awk '/^#/ && (\$0 !~ /^# DNA Sequence:/)' " . $dnafile . ' >>' . $outfile;
	system("awk '/^#/ && (\$0 !~ /^# DNA Sequence:/)' " . $dnafile . ' >>' . $outfile);
    }
}

sub write_param_file
{
    my $param_file = $prefix . '.param';
    
    open PF, ">$param_file" or die "Can't open file to write: $param_file\n";
    print PF 'gatekeeper.project_name = ' . $prefix;
    print PF 'cgb.estimated_genome_length = ' . $dna_length;
    close PF;
}

sub make_polys {
    local($pcopy) = @_;
    delete $opened{$comfile} && close($comfile);
    print "in make_polys, pcopy = ", $pcopy;
    for ($i = 1; $i <= $pcopy; $i++) {
	$npoly += 1;
	$ptotl += $pweight{$npoly};
	print ' +++ Making polymorphism ' . $npoly;
	$command = "$bin/AS_SIM_poly -s " . ($seed + $npoly) . ' ' .

	  $comfile . ' > ' . $prefix . '.poly.' . $npoly;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 11: " . $command;
	$v = system($command);
	if ($v != 0) {
	    print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
	    &cleanup(1);
	}
	if ($comment) {
	    if ($npoly != 1) {
		&Pick('>>', $outfile) &&
		    (print $fh '#');
	    }
	    &Pick('>>', $outfile) &&
		(print $fh '# POLYMORPHISM GENERATION: WEIGHT = ' .
		  $pweight{$npoly});
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 12: " . "awk '/^#/' " . $prefx . '.poly.' . $npoly . ' >>' . $outfile;
	    system("awk '/^#/' " . $prefx . '.poly.' . $npoly . ' >>' . $outfile);
	}
    }
}

sub make_sample {
    $nlibs += 1;
    if ($nlibs == 1) {
	if (&overwrite($frgfile)) {
	    return;
	}
    }

    # we always output an ADT record on the first go round
    if ($nlibs == 1 && !$fastaOutput) 
    {
	$command = "$bin/AS_SIM_outputADT -b " . $acnum++ . " < " . $adlfile . '  > ' . $frgfile;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 13: " . $command;
	$v = system($command);
	if ($v != 0) 
	{
	    print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
	    &cleanup(1);
	}

	# if we are doing things by sample start a new frag file, the first one is special, a header as it were
	if ($sawcmdlinebatchflag)
	{
	    $frgfile = sprintf "%s_B%05.d.frg", $prefix, ++$frgfileCount;
	}
    }

    # now we only output an ADT record if we are not going to cat all the frgfiles together at end
    if ($sawcmdlinebatchflag && !$fastaOutput) {
	system( "echo celsim > $adlfile" . '.' . $frgfileCount );
	system( "echo $dna_length >> $adlfile" . '.' . $frgfileCount );
	$command = "$bin/AS_SIM_outputADT -b " . $acnum++ . " < " . $adlfile . '.' . $frgfileCount . '  > ' . $frgfile;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 13: " . $command;
	$v = system($command);
	if ($v != 0) {
	    print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
	    &cleanup(1);
	}
    }

    {
	@splitFragSpec = @fragspec;
	@splitMassageSpec = @massageSpec;

	$med = int(($splitFragSpec[10] + $splitFragSpec[9]) / 2);
	$dta = int(($splitFragSpec[10] - $splitFragSpec[9]) / 2);

	if ($items > 12 && !$fastaOutput) {
	    $command = "$bin/AS_SIM_outputDST " . $acnum . ' ' . $med . ' ' . $dta . ' >>' . $frgfile;

	    # print "fragspec " fragspec " \nend of fragpsec\n";
	    # print command;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 14: " . $command;
	    $v = system($command);
	    if ($v != 0) {
		print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
		&cleanup(1);
	    }
	}
	$dist_id = $acnum;
	$acnum += 1;
    }

    $cwght = 0;
    $lind = 0;
    for ($i = 1; $i <= $npoly; $i++) {
	&Pick('>', $comfile) &&
	    (print $fh '< ' . $prefix . '.poly.' . $i);
	$cwght += $pweight{$i};
	$nind = int($fnum * $cwght / $ptotl);
	$nfrag{$i} = $nind - $lind;
	&Pick('>', $comfile) &&
	    (print $fh $nfrag{$i});
	$lind = $nind;
	&Pick('>', $comfile) &&
	    (print $fh @fragspec);
	delete $opened{$comfile} && close($comfile);
	$c = $seed + $nlibs * $npoly + $i;
	# frag could have a command line argument with a label
	# for the particular polymorphism we're generating from.
	# This would facilitate labeling the fragments with their poly
	# origin.
	if ($uniform) {
	    $uni = ' -U ';
	}
	else {
	    $uni = ' -N ';
	}
	if ($bacends) {
	    $bac = ' -e ';
	}
	else {
	    $bac = ' -r ';
	    # print "### uni = " uni;
	    ;
	}

	print "-----------------------------------------------------------";
	print @massageSpec;
	print "-----------------------------------------------------------";

	if (!$fastaOutput) {
	    $command = "$bin/AS_SIM_frag -s " . $c . $protoFlag . $uni .
	      ' -F -N ' . $comfile . " | $bin/AS_SIM_massage -q " . $qltfile . ' ';
	    $command = $command . $uni . $bac . $dist_id . ' ' . $acnum . ' ' .
	      "@massageSpec" . ' >>' . $frgfile;

	    #       print STDERR "command: " command;
	    #	print STDERR massageSpec;
	    #	print STDERR fragspec;
	}
	else {
	    $command = "$bin/AS_SIM_frag -s " . $c . $uni . ' -F -N ' .
	      $comfile . ' >> ' . $frgfile;
	}
	print ' +++ Generating library ' . $nlibs . ' from poly ' . $i;
	print $STDERR . 'command: ' . $command;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 15: " . $command;
	$v = system($command);
	$acnum += $nfrag{$i};
	if ($v != 0) {
	    print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
	    &cleanup(1);
	}
    }
    if ($comment) {
	&Pick('>>', $outfile) &&
	    (printf $fh "#\n");
	&Pick('>>', $outfile) &&
	    (printf $fh "# FRAGMENT LIBRARY %d\n", $nlibs);
	&Pick('>>', $outfile) &&
	    (printf $fh "#\n");
	&Pick('>>', $outfile) &&
	    (printf $fh "# Fragments   Seed   Poly Source 1\n");
	&Pick('>>', $outfile) &&
	    (printf $fh "# ---------   ----   -----------\n");
	for ($i = 1; $i <= $npoly; $i++) {
	    &Pick('>>', $outfile) &&
		(printf $fh "#  %8d  %5d   %s.poly.%d.%d\n", 
	    $nfrag{$i}, $seed + $nlibs * $npoly + $i, $prefix, $i, $seed);
	}
	&Pick('>>', $outfile) &&
	    (printf $fh "# ---------\n#  %8d\n", $fnum);
	&Pick('>>', $outfile) &&
	    (printf $fh "#\n");
	&Pick('>>', $outfile) &&
	    (printf $fh "#    Length Range = [%d,%d], F/R odds = %.2f/%.2f\n",

	      
	$splitFragSpec[1], $splitFragSpec[2], $splitFragSpec[3],

	  1 - $splitFragSpec[3]);
	&Pick('>>', $outfile) &&
	    (printf $fh "#\n");
	&Pick('>>', $outfile) &&
	    (printf $fh "# Clear Range Characteristics:\n");
	&Pick('>>', $outfile) &&
	    (printf $fh "#    Prefix Unclear [%d,%d]\n", $splitMassageSpec[1],

	      $splitMassageSpec[2]);
	&Pick('>>', $outfile) &&
	    (printf $fh "#    Suffix Unclear [%d,%d]\n", $splitMassageSpec[3],

	      $splitMassageSpec[4]);
	&Pick('>>', $outfile) &&
	    (printf $fh "# Edit Characteristics:\n");
	&Pick('>>', $outfile) &&
	    (printf $fh '#    Error Ramp = %.2f->%.2f, ', $splitFragSpec[4],

	      $splitFragSpec[5]);
	&Pick('>>', $outfile) &&
	    (printf $fh "Ins/Del/Sub Odds = %.2f/%.2f/%.2f\n", 
	$splitFragSpec[6], $splitFragSpec[7],

	  1 - ($splitFragSpec[6] + $splitFragSpec[7]));
	&Pick('>>', $outfile) &&
	    (printf $fh "#\n");
	if ($items > 12) {
	    &Pick('>>', $outfile) &&
		(printf $fh "# Dual-End Inserts:\n");
	    &Pick('>>', $outfile) &&
		(printf $fh "#    Single Odds = %.2f\n", $splitFragSpec[8]);
	    &Pick('>>', $outfile) &&
		(printf $fh "#    Insert Range = [%d,%d]\n",
		  $splitFragSpec[9], $splitFragSpec[10]);
	    &Pick('>>', $outfile) &&
		(printf $fh "#    Pairing Error Rate = %.2f\n",
		  $splitFragSpec[11]);
	    &Pick('>>', $outfile) &&
		(printf $fh "#\n");
	}
    }
    $fragspec = '';
    $massageSpec = '';

    # if we are doing things by sample start a new frag file
    if ($sawcmdlinebatchflag)
    {
	$frgfile = sprintf "%s_B%05.d.frg", $prefix, ++$frgfileCount;
    }
}

sub make_uflbacs {
    my $bac_type = shift @_;
    my $shred = shift @_;
    my $min_num_bacs = shift @_;
    my $max_num_bacs = shift @_;
    my $min_bac_length = shift @_;
    my $max_bac_length = shift @_;
    my $sequencing_error = shift @_;
    my $insert_error_rate = shift @_;
    my $delete_error_rate = shift @_;

    $nlibs += 1;
    if ($nlibs == 1) {
	if (&overwrite($frgfile)) {
	    return;
	}
    }

    # we always output an ADT record on the first go round
    if ($nlibs == 1 && !$fastaOutput) 
    {
	$command = "$bin/AS_SIM_outputADT -b " . $acnum++ . " < " . $adlfile . '  > ' . $frgfile;
	print $STDERR . 'adlfile = ' . $adlfile . ' frgfile= ' . $frgfile . ' command = ' . $command;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 132: " . $command;
	$v = system($command);
	if ($v != 0) 
	{
	    print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
	    &cleanup(1);
	}

	# if we are doing things by sample start a new frag file, the first one is special, a header as it were
	if ($sawcmdlinebatchflag)
	{
	    $frgfile = sprintf "%s_B%05.d.frg", $prefix, ++$frgfileCount;
	}
    }

    # now we only output an ADT record if we are not going to cat all the frgfiles together at end
    if ($sawcmdlinebatchflag && !$fastaOutput) {
	system( "echo celsim > $adlfile" . '.' . $frgfileCount );
	system( "echo $dna_length >> $adlfile" . '.' . $frgfileCount );
	$command = "$bin/AS_SIM_outputADT -b " . $acnum++ . " < " . $adlfile . '.' . $frgfileCount . '  > ' . $frgfile;
	print "\n %%%%%%%%%%%%%%%%%%%%% command 13: " . $command;
	$v = system($command);
	if ($v != 0) {
	    print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
	    &cleanup(1);
	}
    }

    $num_bacs = $min_num_bacs + int rand ($max_num_bacs - $min_num_bacs);

    # determine length of bacs
    $total_bac_bases = 0;
    for ($bac_cnt = 1; $bac_cnt <= $num_bacs; $bac_cnt++) 
    {
	$total_bac_length{$bac_cnt} = $min_bac_length + int rand ($max_bac_length - $min_bac_length);
	$total_bac_bases += $total_bac_length{$bac_cnt};
    }

    # now determine size of gaps between them
    $total_gap_bases = $dna_length - $total_bac_bases;

    # for each gap generate a random number, sum them, normalize, and assign gap size based on percent of sum
    $gap_sum = 0;

    # if we have more bac bases than genome bases the first and last gaps cannot be negative
    if ($total_gap_bases < 0)
    {
	$gap_percent{0} = 0;
	$gap_percent{$num_bacs} = 0;
    }
    else
    {
	$gap_percent{0} = rand();
	$gap_sum += $gap_percent{0};
	$gap_percent{$num_bacs} = rand();
	$gap_sum += $gap_percent{$num_bacs};
    }

    # gap 0 is to left of first bac
    for ($gap_cnt = 1; $gap_cnt < $num_bacs; $gap_cnt++) 
    {
	$gap_percent{$gap_cnt} = rand();
	# $gap_percent{$gap_cnt} = 1.0;  #for tiling, any non-zero number will do
	$gap_sum += $gap_percent{$gap_cnt};
    }

    # determine gap sizes
    $check_gap_bases = 0;
    for ($gap_cnt = 0; $gap_cnt <= $num_bacs; $gap_cnt++) 
    {
	$gap_size{$gap_cnt} = int (($gap_percent{$gap_cnt} / $gap_sum) * $total_gap_bases);
	$check_gap_bases += $gap_size{$gap_cnt};
    }
    #print "********* \$check_gap_bases: $check_gap_bases";
    #print "********* \$num_bacs: $num_bacs";
    #print "********* \$gap_size{0}: $gap_size{0}";

    $bac_start{1} = $gap_size{0};
    for ($bac_cnt = 2; $bac_cnt <= $num_bacs; $bac_cnt++) 
    {
	$bac_start{$bac_cnt} = $bac_start{$bac_cnt - 1} + $total_bac_length{$bac_cnt - 1} + $gap_size{$bac_cnt - 1};
	print "\$bac_start{" . ($bac_cnt - 1) . "}: " . $bac_start{$bac_cnt - 1}, 
	      "\t bac_end: " . ($bac_start{$bac_cnt - 1} + $total_bac_length{$bac_cnt - 1}),
	      "\t\$gap_size{" . ($bac_cnt - 1) . "}: " . $gap_size{$bac_cnt - 1};
    }
    print "\$bac_start{" . $num_bacs . "}: " . $bac_start{$num_bacs}, 
          "\t bac_end: " . ($bac_start{$bac_cnt - 1} + $total_bac_length{$bac_cnt - 1}),
	  "\t\$gap_size{" . $num_bacs . "}: " . $gap_size{$num_bacs};

    for ($bac_cnt = 1; $bac_cnt <= $num_bacs; $bac_cnt++) 
    {
	$bac_length = 0;

	# this is the file we write frags for this bac to
	$temp_bac_file = "fastaBacs";
	system("rm $temp_bac_file");

	# determine number of bactigs
	if ($bac_type == $ubac_type)
	{
	    # total_bac_length is not what we put into the DST record
	    # that quantity is the sum of the bactig sizes
	    # this is the size of the insert over which we are going to create bactigs

	    $num_bactigs = build_bactigs($ubac_coverage, $total_bac_length{$bac_cnt});
	}
	else
	{
	    $num_bactigs = 1;
	    $bactig_length{1} = $total_bac_length{$bac_cnt};
	}

	# get all the acnums needed for the DST and BAC records
	$bac_dist_id = $acnum++;
	$locale = $acnum++;  
	$seq_id = $acnum++;

	# determine size of bac not including gaps
	for ($bactig_cnt = 1; $bactig_cnt <= $num_bactigs; $bactig_cnt++) 
	{
	    $bactig_acnum{$bactig_cnt} = $acnum++;

	    $bac_length += $bactig_length{$bactig_cnt};
	}

	# put out a distance record for this bac
	&emitDSTRecord($frgfile, $bac_length, int ($bac_length / 3), $bac_dist_id);

	# now we have all the info needed for a BAC record
        # we need $first_frag_acnum for massage
	if ($shred == 0 && $bac_type == $ubac_type) {
	    $first_frag_acnum = $bactig_acnum{1};  # accession numbers for unshredded UBACs are the bactig ID
	}
	elsif ($shred == 0 && $bac_type == $fbac_type) {
	    $first_frag_acnum = $seq_id;  # accession numbers for fullBacs are the sequence id
	}
	else {
	    $first_frag_acnum = $acnum;  # frag has yet to be created but will have this accession number
	}

	# now determine placement of bac
	# $bac_start = int rand ($dna_length - $bac_length);
	# $bactig_start_current = $bac_start;

	# create and output the bactigs
	for ($bactig_cnt = 1; $bactig_cnt <= $num_bactigs; $bactig_cnt++) 
	{
	    #print "\$bactig $bactig_cnt starts at $bactig_start{$bactig_cnt} and ends at " . 
		#($bactig_start{$bactig_cnt} + $bactig_length{$bactig_cnt});

	    if ($shred == 1) 
	    {
		$fragsize = 500;
		$coverage = 2;
	    }
	    else  # not shredding
	    {
		$fragsize = $bactig_length{$bactig_cnt};
		$coverage = 1;
	    }

	    # now create the bactig, shred if necessary, and write results to a temp file ($temp_bac_file)

	    print "in make_uflbacs, before adding bactig, \$bactig_length{$bactig_cnt}: $bactig_length{$bactig_cnt}";
	    print "in make_uflbacs, before adding bactig, \$bac_start{$bac_cnt}: $bac_start{$bac_cnt}";
	    print "in make_uflbacs, before adding bactig, \$bactig_start{$bactig_cnt}: $bactig_start{$bactig_cnt}";

	    &add_bactig($locale, $bactig_acnum{$bactig_cnt}, $bac_start{$bac_cnt} + $bactig_start{$bactig_cnt}, 
			$bactig_length{$bactig_cnt},
			$fragsize, $coverage, $sequencing_error, $insert_error_rate, $delete_error_rate, $temp_bac_file);
	    $bactig_length{$bactig_cnt} = $bactig_length_after_idas;

	    #print "in make_uflbacs, after adding bactig, \$acnum = $acnum";
	    #print "in make_uflbacs, after adding bactig, \$first_frag_acnum = $first_frag_acnum";
	}

	# open the temporary file that will hold the bactig information
	open BF, ">$bactig_file" or die "Can't open file to write: $bactig_file\n";
	print BF $num_bactigs;

	for ($bactig_cnt = 1; $bactig_cnt <= $num_bactigs; $bactig_cnt++) 
	{
	    # output an accession number and a length
	    print BF $bactig_acnum{$bactig_cnt}, $bactig_length{$bactig_cnt};
	}
	close BF;

	if ($uniform) {
	    $uni = ' -U ';
	}
	else {
	    $uni = ' -N ';
	}
	
	if ($bac_type == $ubac_type || $bac_type == $fbac_type)
	{
	    if ($bac_type == $ubac_type)
	    {
		if ($shred == 0) {
		    $bac = '-b'; # bactigs
		}
		else {
		    $bac = '-u'; # shredded ubacs
		}
	    }
	    elsif ($bac_type == $fbac_type)
	    {
		if ($shred == 0) {
		    $bac = '-c'; # fullbacs
		}
		else {
		    $bac = '-f'; # shredded fbacs
		}
	    }
	    print 'fastaOutput: ' . $fastaOutput;
	    print ' +++ Generating library ' . $nlibs . ' from poly ' . $i;
	    
	    # massage and output the fragments
	    $command = "$bin/AS_SIM_massage " . 
		$uni . $bac . " $bactig_file" . ' ' . $bac_dist_id . ' ' . $locale . ' ' . $seq_id . ' ' . 
		    ' ' . $first_frag_acnum . ' ' .
			"@massageSpec" . ' >> ' . $frgfile . ' < ' . $temp_bac_file;
	    print $STDERR . 'command: ' . $command;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 21: " . $command;
	    $v = system($command);
	    if ($v != 0) {
		print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
		&cleanup(1);
	    }
	}
	elsif ($bac_type == $lbac_type)
	{
	    $bac = '-l ';

	    # calculate number of frags we need to cover lbac to specified depth
	    $fnum = $lbac_coverage * int ($bac_length / 500.0);

	    # for lbacs we need to call frag
	    # first we write the file for frag to operate on
	    $frag_specs = 'fastaBacs.spec';
	    open FRAG_IN, ">$frag_specs" or die "Can't open file to read: $frag_specs";
	    print FRAG_IN "< $temp_bac_file";
	    print FRAG_IN "$fnum";
	    print FRAG_IN @fragspec;
	    close FRAG_IN;

	    # frag the "bactig" that is the lbac
	    $lbac_frags = 'lbac_frags_temp';
	    $command = "$bin/AS_SIM_frag -s " . $c . $protoFlag . $uni . ' -F -N ' . $frag_specs . "> " . $lbac_frags;

	    print $STDERR . 'command: ' . $command;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 25: " . $command;
	    $v = system($command);
	    if ($v != 0) {
		print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
		&cleanup(1);
	    }

	    # $bac_start is where the lbac starts
	    # need to use the quality file created by frag in $temp_bac_file . "qlt"
	    $temp_quality_file = $temp_bac_file . ".qlt";
	    $command = "$bin/AS_SIM_massage -q " . $temp_quality_file . ' ' .
	               $uni . $bac . $dist_id . ' ' . $locale . ' ' . $bac_start{$bac_cnt} . ' ' .$acnum . ' ' .
	               "@massageSpec" . ' >> ' . $frgfile  . ' < ' . $lbac_frags;
	    print $STDERR . 'command: ' . $command;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 26: " . $command;
	    $v = system($command);
	    if ($v != 0) {
		print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
		&cleanup(1);
	    }
	    # frag has split the "bactig" into $fnum fragments, need to adjust $acnum
	    $acnum += $fnum;
	}
	else
	{
	    print "unknown bac type";
	    exit;
	}
    }
    $fragspec = '';
    $massageSpec = '';

    # if we are doing things by sample start a new frag file
    if ($sawcmdlinebatchflag)
    {
	$frgfile = sprintf "%s_B%05.d.frg", $prefix, ++$frgfileCount;
    }
}

sub cleanup {
    local($ecode) = @_;

    exit(ecode); # TEMPORARY, for debug
    if ($phase >= 0) {
	if ($keepDNA == 0) {
	    system('rm -f ' . $dnafile);
	}
	if ($notations > 0) {
	    system('rm ' . $notfile);
	}
    }
    if ($phase >= 1) {
	for ($i = 1; $i <= $npoly; $i++) {
	    system('rm ' . $prefix . '.poly.' . $i);
	}
    }
    system('rm ' . $comfile);
    system('rm ' . $qltfile);
    if (!$fastaOutput) {
	system('rm ' . $adlfile);
    }

    if (!$sawcmdlinebatchflag)
    {
	# first remove the header
	{
	    $frgfile = sprintf "%s_B%05.d.frg", $prefix, $frgfileNum;
	    $perlcmd = "rm " . $frgfile;
	    system($perlcmd);
	}

	# now remove everybody else
	for ($frgfileNum = 2; $frgfileNum < $frgfileCount; $frgfileNum++)
	{
	    $frgfile = sprintf "%s_B%05.d.frg", $prefix, $frgfileNum;
	    $perlcmd = "rm " . $frgfile;
	    system($perlcmd);
	}
    }
    $notdone = 0;
    $ExitValue = ($ecode);
    exit;  # MJF was: last line
}

sub overwrite {
    local($name) = @_;
    if (system('test -e ' . $name) == 0) {
	if ($prompt) {
	    printf '  Overwrite %s [y/n]? ', $name;
	    $answer = &Getline3('/dev/stdin');
	    if ($answer =~ "^y\$") {
		system('rm ' . $name);
	    }
	    else {
		return (1);
	    }
	}
	else {
	    system('rm ' . $name);
	}
    }
    (0);
}

sub Getline3 {
    &Pick('',@_);
    local($_);
    if ($getline_ok = (($_ = <$fh>) ne '')) {
	chop;	# strip record separator
    }
    $_;
}

sub Pick {
    local($mode,$name,$pipe) = @_;
    $fh = $name;
    open($name,$mode.$name.$pipe) unless $opened{$name}++;
}

sub min {
    ($a, $b) = @_;
    (($a) < ($b)) ? return ($a) : return ($b);
}

sub max {
    ($a, $b) = @_;
    (($a) > ($b)) ? return ($a) : return ($b);
}

# _______________________________________________________________________________
# Below is an altered version of what used to be "addfbac" written by Daniel 
# Huson.  

sub add_bactig
{
    my $locale = shift @_;
    my $bactig_acnum = shift @_;
    my $startpos = shift @_;
    my $bactig_length_in = shift @_;
    my $fragsize = shift @_;
    my $coverage = shift @_;
    my $seq_error = shift @_;
    my $insert_error = shift @_;
    my $delete_error = shift @_;
    my $outfile = shift @_;  # always the fragment file, so far

    # because of the use of a2p, the array are 1-based in the main routine (see $[ above)
    # want the 0-based here
    $[ = 0;

    # script switches:
    $verbose = 0;

    if ($verbose == 1)
    {
	print "------add_bactig------";
	print "locale           = $locale";
	print "startpos         = $startpos";
	print "bactig length    = $bactig_length_in";
	print "fragsize         = $fragsize";
	print "coverage         = $coverage";
	print "outfile          = $outfile";
    }

    $quality = chr(48+30);  # quality value for fbac fragments, here: 30

    if($verbose == 1)
    {
	print "accession = $acnum";
    }
    
    # We need to compute the increment that will yield the given coverage for
    # the given fragsize
    
    $increment = int $fragsize / $coverage;
    if($verbose == 1)
    {
	print "increment = $increment";
    }
    
    if (!$read{$dnafile})
    {
	# Read the dna:
	$seq = "";
	$seqlen = 0;    

	$read{$dnafile} = 1;

	open FP, "$dnafile" or die "Can't open file to read: $dnafile";

	while(<FP>)
	{
	    chomp;
	    if ($_ eq "> 0") # we leave this state after finding ">"
	    {
		last;
	    }
	}
	
	$old_input_record_separator = $/; # current input record separator
	$/ = "";
	
	while (<FP>)
	{
	    $seq = $_;    # suck everything remaining into $seq
	}
	close FP;
	
	$/ = $old_input_record_separator;
	
	$seq =~ s/\s//g;  # get rid of newlines

	$seq =~ tr /A-Z/a-z/; # canonicalize to lower case

	$seqlen = length $seq;
    
	if($verbose)
	{
	    print "Sequence  = [0, ", $seqlen-1, " ], length = $seqlen";
	}
    }
    # we cut out the fbac from the dna:
    
    $fbac = substr($seq, $startpos, $bactig_length_in);
    $fbaclen = length $fbac;
    
    if ((rand) < 0.5)
    {
	$reversed = 0;
	$direction = '';
    }
    else
    {
	$reversed = 1;
	$direction = "REVERSED";
	$fbac = reverse $fbac;
	$fbac =~ tr /aA/%/;
	$fbac =~ tr /tT/a/;
	$fbac =~ tr /%/t/;
	$fbac =~ tr /cC/%/;
	$fbac =~ tr /gG/c/;
	$fbac =~ tr /%/g/;
    }
    
    if($verbose)
    {
	print "bactig    = [$startpos,", $startpos + $fbaclen - 1, "], length = $fbaclen";
    }

    # now introduce the sequencing error
    for ($icnt = 0; $icnt < $fbaclen; $icnt++)
    {
	if (rand(1) < $seq_error)
	{
	    $error_type = rand(1);
	    if ($error_type < $insert_error) # insertion
	    {
		$new_base = $acgt{int(rand(4))};
		$fbac = substr($fbac, 0, $icnt + 1) . $new_base . substr($fbac, $icnt + 1, $fbaclen - $icnt - 1);
		$fbaclen++;
		$icnt++;
	    }
	    elsif ($error_type < $insert_error + $delete_error)  # deletion
	    {
		$fbac = substr($fbac, 0, $icnt) . substr($fbac, $icnt + 1, $fbaclen - $icnt - 1);
		$fbaclen--;
	    }
	    else #substitution
	    {
		$current_base = substr($fbac, $icnt, 1);
		$rand_index = ($base_index{$current_base} + 1 + int(rand(3))) % 4;  #choose randomly amongst the other bases
		substr($fbac, $icnt, 1) = $acgt{$rand_index};
	    }
	}
    }
    
    # open the outfile:
    open FP, ">>$outfile" or die "Can't open file to write: $outfile";

    # shred the bac and write the output:
    $done = 0;
    $count = 0;
    $pos = 0;
    $linelength = 50;

    # insertions and deletions can change $fbaclen
    if ($shred == 0)
    {
	$fragsize = $fbaclen;
    }

    while($done == 0)
    {
	$count++;
	if ($pos + $fragsize >= $fbaclen)	# last fragment, make it fit and after
	                                        # printing it, we're done!
	{
	    $pos = max(0, $fbaclen - $fragsize);  # sometimes $fbaclen is smaller than $fragsize
	    $done = 1;
	}
	$frag = substr($fbac, $pos, $fragsize);
	$fraglen = length $frag;
	$fasta = 1;
	
	if (!$fasta)
	{
	    print FP "{FRG";
	    print FP "act:A";
	    print FP "acc:",$acnum++,"";
	    print FP "typ:F";
	    print FP "loc:$locale";
	    print FP "pos:$pos,",$pos + $fraglen - 1,"";
	    print FP "src:";
	    print FP "f";
	    print FP "[",$startpos+$pos,",",$startpos+$pos+$fraglen-1,"].";
	    print FP "etm:",$etm++,"";
	    print FP "seq:";
	}
	else
	{
	    if ($pos < 0)
	    {
		print "uh oh, \$pos: $pos";
	    }

	    if ($reversed == 0)
	    {
		$start = $pos;
		$end = $pos + $fraglen;
		print FP "> " . $acnum++ . ' Interval: [' . ($startpos + $start) . ", " . ($startpos + $end) . "] $direction" .
		    " bactig_id:" . $bactig_acnum . " locpos[" . $pos . "," . ($pos + $fraglen) . "]" .
		    " genomepos[" . $startpos . "," . ($startpos + $bactig_length_in) . "]";
	    }
	    else # ($reversed == 1)
	    {
		$start = $fbaclen - $pos;
		$end = $fbaclen - $pos - $fraglen;
		print FP "> " . $acnum++ . ' Interval: [' . ($startpos + $end) . ", " . ($startpos + $start) . "] $direction" .
		    " bactig_id:" . $bactig_acnum . " locpos[" . $pos . "," . ($pos + $fraglen) . "]" .
		    " genomepos[" . $startpos . "," . ($startpos + $bactig_length_in) . "]";
	    }
	}

	for ($i = 0; $i < $fraglen; $i += $linelength)
	{
	    print FP substr($frag, $i, $linelength);
	}

	if (!$fasta)
	{
	    print FP ".";
	    print FP "qlt:";
	    for($i = 0; $i < $fraglen; $i++)
	    {
		if( ( ($i) % $linelength) == 0)
		{
		    print FP "";
		}
		print FP $quality;
	    }
	    print FP ".";
	    print FP "clr:0,",$fraglen-1,"";
	    print FP "}";
	}
	$pos += $increment;
    }
    if($verbose)
    {
	print "Wrote $count fragments to file: $outfile";
	print "-------end-------";
    }
    close FP;

    # arrays are 1-based in main routine
    $[ = 1;

    print "\$fbaclen: $fbaclen";
    $bactig_length_after_idas = $fbaclen;
}

sub emitDSTRecord
{
    ($file, $length, $stddev, $dist_id) = @_;

    print STDERR "============================== emitting DST (acc: $dist_id) record to $file";

    if (0)
    {
	open FF, ">>$file" or die "Can't open file to write: FF\n";
	print FF "{DST";
	print FF "act:A";
	print FF "acc:$dist_id";
	print FF "mea:$length";
	print FF "std:$stddev";
	print FF "}";
	close FF;
    }
    else
    {
	$command = "$bin/AS_SIM_outputDST $dist_id $length $stddev >> $file";
	    print $STDERR . 'command: ' . $command;
	    print "\n %%%%%%%%%%%%%%%%%%%%% command 17: " . $command;
	    $v = system($command);
	    if ($v != 0) {
		print $STDERR . ' Error Invoking ' . $command . ' ...exiting';
		&cleanup(1);
	    }
    }
}


# this routine uses simulated fragment reads to build bactigs
# used only for building bactig lengths and placement for ubacs
sub build_bactigs
{
    my $coverage = shift @_;
    my $bac_length = shift @_;
    my $icnt;
    my ($frag_begin, $frag_start, $frag_end, $frag_terminus);
    my ($bactig_end, $bactig_count);

    my $num_fragments = int ($coverage * $bac_length / 500);

    #print "in build_bactigs, \$coverage: $coverage\n";
    #print "in build_bactigs, \$bac_length: $bac_length\n";
    #print "in build_bactigs, \$num_fragments: $num_fragments\n";

    %frag_terminus = '';

    # generate a bunch of frags
    for ($icnt = 0; $icnt < $num_fragments; $icnt++)
    {
	my $frag_length = 500 + int (10.0 * gasdev());   # get the paper by Waterman and Lander to get std. dev.

	if ($frag_length < 150)
	{
	    print "\$frag_length: $frag_length";
	    exit;
	}

	$frag_begin = int rand ($bac_length - $frag_length - 1);

	# a slight hack, force the frags to end at different places
	# while (($frag_begin + $frag_length > $bac_length) || ($frag_terminus{$frag_begin} != ''))
	while ($frag_terminus{$frag_begin} != '')
	{
	    $frag_begin = int rand ($bac_length - $frag_length - 1);
	}
	push @{$frag_terminus{$frag_begin}}, $frag_begin + $frag_length;
    }

    $icnt = 0;
    foreach $frag_begin (sort numerically keys %frag_terminus)
    {
	while ($temp = pop @{$frag_terminus{$frag_begin}})
	{
	    $frag_start{$icnt} = $frag_begin;
	    $frag_end{$icnt} = $temp;
	    $icnt++;
	    #print "[$frag_start{$icnt}, $frag_end{$icnt}]";
	}
    }

    $bactig_count = 1;
    $bactig_start{$bactig_count} = $frag_start{0};
    $bactig_end{$bactig_count} = $frag_end{0};
    for ($icnt = 1; $icnt < $num_fragments; $icnt++)
    {
	if ($frag_end{$icnt - 1} - $frag_start{$icnt} > 40)  #if they overlap by > 40 then they are in the same bactig
	{
	    $bactig_end{$bactig_count} = $frag_end{$icnt};
	    #print "1 frag $icnt [$frag_start{$icnt}, $frag_end{$icnt}] is in bactig $bactig_count";
	}
	else
	{
	    #print "3 bactig $bactig_count: [$bactig_start{$bactig_count}, $bactig_end{$bactig_count}]";
	    $gap_length{$bactig_count} = $frag_start{$icnt} - $frag_end{$icnt - 1};
	    $bactig_count++;
	    $bactig_start{$bactig_count} = $frag_start{$icnt};
	    $bactig_end{$bactig_count} = $frag_end{$icnt};
	    #print "2 frag $icnt [$frag_start{$icnt}, $frag_end{$icnt}] is in bactig $bactig_count";
	}
    }
    #print "4 bactig $bactig_count: [$bactig_start{$bactig_count}, $bactig_end{$bactig_count}]";

    # note: @bactig_length is 1-based
    $bac_coverage = 0;
    for ($icnt = 1; $icnt <= $bactig_count; $icnt++)
    {
	# print "bactig_start, bactig_end, gap ", $bactig_start{$icnt}, " ", $bactig_end{$icnt}, " ", $gap_length{$icnt};
	$bactig_length{$icnt} = $bactig_end{$icnt} - $bactig_start{$icnt};
	$bac_coverage += $bactig_length{$icnt};
    }

    #print "in build_bactigs, \$bac_coverage: $bac_coverage";
    #print "in build_bactigs, coverage: " . 100.0 * ($bac_coverage/$bac_length) . "%";

    return $bactig_count;
}

sub numerically { $a <=> $b }

# the following subroutine is modelled on gasdev() as given in "Numerical Recipes in C", 2nd Ed.
sub gasdev
{
    $iset = 0;

    if ($iset == 0)
    {
	do
	{
	    $v1 = 2.0 * rand(1.0) - 1.0;
	    $v2 = 2.0 * rand(1.0) - 1.0;
	    $rsq = $v1 * $v1 + $v2 * $v2;
	} while ($rsq >= 1.0 || $rsq == 0.0);
	$fac = sqrt( - 2.0 * log( $rsq) / $rsq);
	$gset = $v1 * $fac;
	$iset = 1;
	return $v2 * $fac;
    }
    else
    {
	$iset = 0;
	return $gset;
    }
}
