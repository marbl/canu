#!/usr/local/bin/perl

# $Id: caqc.pl,v 1.27 2007-12-17 19:30:52 brianwalenz Exp $
#
# This program reads a Celera .asm file and produces aggregate information
# about the assembly
#
# Written by Mihai Pop, Martin Shumway

#
#

#use warnings;
use strict;
use Getopt::Long;
use IO::File;
use File::Basename;
use Statistics::Descriptive;
use File::Copy;
use Math::BigFloat;
use FindBin qw($Bin);

my $MY_VERSION = "caqc Version 2.13 (Build " . (qw/$Revision: 1.27 $/)[1] . ")";

# Constants
my $MINQUAL    = 20;
my $MINCONTIG  = 10000;
my @help_lines = undef;
my $help_file  = $Bin . '/caqc_help.ini';

my $MY_HELPTEXT = qq~
Generate quality statistics from the specified Celera assembly .asm file.  

  caqc  <prefix>  [options]

    <prefix>   caqc reads the <prefix>.asm, which is the output file of the Celera Assembler. 
               caqc requires that the file be in the current directory.
    
    options:
      -minqual   <n>   Minimum quality value threshhold to report as bad 
                       quality (default $MINQUAL)
      -mincontig <n>   Minimum contig size to report as a big contig
                       (default $MINCONTIG) 
      -d <c>           Specify the record delimiter tag<c>value
      -s <c>           Specify the list delimiter el[<c>el ...]
      -h               Display help information and exit
      -h <term>        Display help on output qc term. Use 'all' for all qc terms.
      -html            Output qc help as html
      -silent          Do not emit output to the stdout
      -euid            Print EUIDs
      -g <n>           Genome size used in the calculation of N50 numbers
                       (default: TotalBasesInContigs for contigs and TotalBasesInScaffolds for scaffolds)
      -debug <level>   Set the debug <level>.
      -l <logfile>     Specify the logfile.  Only created if debug level set.
      -t <count>       Specify the number of top scaffolds and contigs to use.
                       (default: 5)
      -frg             Computes frag related statistics: BasesCount and ClearRangeLengthFRG. 
                       This option requires that a <prefix>.frg file exist in the current directory.
      -metrics         Option to output a <prefix>.qc.metrics file which is inputted to ametrics.  This 
                       option requires that a <prefix>.frg file exist in the current directory. 
                       [Note: -frg is a subset of -metrics]

		       
    output files:
	<prefix>.qc		Described below.
	<prefix>.qc.metrics 	(optional) Input file for ametrics


caqc produces output on stdout as well as a <prefix>.qc file in the 
current directory.  The file output is in a tag-value format organized in a
hierarchical manner compatible with the INI style for config files.  This
file can be read using TIGR::ConfigFile (run perldoc TIGR::ConfigFile). 
Delimiters can be controlled by using the tag-value separator specified in
the -d option, and the minor delimiter separating list values by using the 
-s option.  By default both are set to whitespace for human readability
in the stdout output, and to = and comma respectively in the .qc file.  The 
optional <prefix>.qc.metrics file is a copy of the <prefix>.qc file with additional 
metrics information for ametrics, which uploads it into the ASDB database.

NOTE: caqc has backwards compatibility to support usage of:
    caqc <prefix>.asm [options]

See also:
  http://intranet.tigr.org/software_docs/CeleraAssembler.shtml
  TIGR::ConfigFile 
~;

# Global delimiters
my $DEFAULT_TAG_DELIM   = '=';
my $DEFAULT_FIELD_DELIM = ',';
my $DEFAULT_TOP_COUNT   = 5;
my $d                   = $DEFAULT_TAG_DELIM;      # tag-value delimiter
my $s                   = $DEFAULT_FIELD_DELIM;    # list separator
my $silent              = 0;
my $euid                = 0;

##############################################

=over

=item B<my $rec = getCARecord(\*STDIN);>

Reads from stdin the text between "extreme" { and } .

 for example:

  {A
    {B
    }
   }
 
Returns the whole: {A{B}}

=cut

sub getCARecord {
    my $file = shift;

    my $level = 0;
    my $block = "";

    while (<$file>) {
        if (/^\s*\{/) {
            $level++;
        }
        if (/^\s*\}/) {
            $level--;
        }
        $block .= $_;
        if ( $level == 0 ) {
            last;
        }
    }

    if ( $level != 0 ) {
        die("end of file reached before end of block\n");
    }

    if ( $block ne "" ) {
        return $block;
    }
    else {
        return undef;
    }
}    # getCARecord

######################################################3

=item B<my($id, $fields, $recs) = parseCARecord($rec);>

Parses a record and returns a triplet consisting of
   - record type
   - hash of fields and values
   - array of sub-records

=cut

sub parseCARecord {
    my $record = shift;

    my @lines = split( '\n', $record );

    my $type;
    my %fields;
    my @recs;

    # get record type
    $lines[0] =~ /\{(\w+)/;
    if ( !defined $1 ) {
        die("Wierd start of record: $record\n");
    }
    $type = $1;

    if ( $lines[$#lines] !~ /^\s*\}/ ) {
        die("Wierd end of record: $record\n");
    }

    my $level = 0;
    my $fieldname;
    for ( my $i = 1 ; $i < $#lines ; $i++ ) {
        if ( $lines[$i] =~ /^(\w+):(\S+)$/ ) {    # simple field
            $fields{$1} = $2;
        }    # simple field
        if ( $lines[$i] =~ /^(\w+):$/ ) {    # complex field
            $fieldname = $1;
            $fields{$fieldname} = "";
            $i++;
            while ( $i < $#lines && ( $lines[$i] !~ /^\.$/ ) ) {
                $fields{$fieldname} .= "$lines[$i]\n";
                $i++;
            }
        }    # complex field
        if ( $lines[$i] =~ /^\s*\{/ ) {    # subrecord
            my $level = 1;

            my $thisrec = ++$#recs;

            $recs[$thisrec] .= "$lines[$i]\n";
            $i++;
            while ( $level > 0 && $i < $#lines ) {
                if ( $lines[$i] =~ /^\s*\{/ ) {
                    $level++;
                }
                if ( $lines[$i] =~ /^\s*\}/ ) {
                    $level--;
                }
                $recs[$thisrec] .= "$lines[$i]\n";
                if ( $level == 0 ) {
                    last;
                }
                else {
                    $i++;
                }
            }
            if ( $level != 0 ) {
                die("Error parsing sub_record in: $record\n");
            }
        }    # subrecord
    }    # for $i...

    return ( $type, \%fields, \@recs );
}    # parseCARecord

################################################

=item B<my($id) = getBothCAIds($CAid);>

Obtains both IDs from a "paired" id, that is, converts (10, 1000) into 10,1000.
If the Id is not a pair in parantheses, it returns the input.
Thus, getBothCAIds('(10, 1000)') returns 10,1000 while getBothCAIds("abba") returns "abba".

=cut

sub getBothCAIds {
    my $string = shift;

    if ( $string =~ /\((\d+),(\d+)\)/ ) {
        return ( $1, $2 );
    }
    else {
        return $string;    # just in case we have a real ID
    }
}    # getBothCAIds

################################################

=item B<my($id) = getCAId($CAid);>

Obtains the ID from a "paired" id, that is, converts (10, 1000) into 10.
If the Id is not a pair in parantheses, it returns the input.
Thus, getCAId('(10, 1000)') returns 10 while getCAId("abba") returns "abba".

=cut

sub getCAId {
    my $string = shift;

    if ( $string =~ /\((\d+),(\d+)\)/ ) {
        return $1;
    }
    else {
        return $string;    # just in case we have a real ID
    }
}    # getCAId

################################################

# calculate the gc content of contigs in scaffolds by getting the
# the percentage of g and c bases in the contig consensus for all the
# scaffold contigs
sub calculateGC($$$) {
    my $scaffcontigs_ref = shift;
    my $gcs_ref          = shift;
    my $sizes_ref        = shift;
    my %scaffcontigs     = %$scaffcontigs_ref;
    my %gcs              = %$gcs_ref;
    my %sizes            = %$sizes_ref;
    my @contigids        = keys(%scaffcontigs);
    my $gc_tot           = 0;
    my $cons_tot         = 0;
    foreach my $contig_id (@contigids) {
        my $gc_num      = $gcs{$contig_id};
        my $contig_size = $sizes{$contig_id};
        $gc_tot   += $gc_num;
        $cons_tot += $contig_size;
    }
    my $gc_perc = 0;
    if ( $cons_tot != 0 ) {
        $gc_perc = ( $gc_tot / $cons_tot );
    }
    return $gc_perc;
}

# Emit the results to the terminal, and if prefix specified, to the
# prefix.qc file, as a config-readable output.
#
sub printl($$$) {
    my $tag      = shift;
    my $rh_value = shift;
    my $rfh      = shift;

    my $s1 = sprintf( '%-32s', $tag );
    print( STDERR sprintf( "%s%-10s", $s1, $$rh_value{$tag} ) . "\n" )
      if ( !$silent );
    $rfh->print( "$tag=" . $$rh_value{$tag} . "\n" );
}

# Emit the results to the terminal, and if prefix specified, to the
# prefix.qc file, as a config-readable output.  Integer form.
#
sub printld($$$) {
    my $tag      = shift;
    my $rh_value = shift;
    my $rfh      = shift;

    my $s1 = sprintf( "%-32s", $tag );
    print( STDERR sprintf( "%s%-d",
            $s1, Math::BigFloat->new( $$rh_value{$tag} )->ffround(0) )
          . "\n"
      )
      if ( !$silent );
    $rfh->print( "$tag="
          . sprintf( "%d", Math::BigFloat->new( $$rh_value{$tag} )->ffround(0) )
          . "\n" );
}

# Emit the results to the terminal, and if prefix specified, to the
# prefix.qc file, as a config-readable output.  Float form.
#
sub printlf($$$) {
    my $tag      = shift;
    my $rh_value = shift;
    my $rfh      = shift;

    my $s1 = sprintf( "%-32s", $tag );
    print( STDERR sprintf( "%s%-0.2f", $s1, $$rh_value{$tag} ) . "\n" )
      if ( !$silent );
    $rfh->print( "$tag=" . sprintf( "%0.2f", $$rh_value{$tag} ) . "\n" );
}

# output percentage of total reads

sub printp($$$$) {
    my ( $tag, $rh_value, $rfh, $tot ) = @_;

    my $val = sprintf "%s(%.2f%%)", $$rh_value{$tag},
      $$rh_value{$tag} / $tot * 100;
    printf( STDERR "%-32s%-10s\n", $tag, $val ) if ( !$silent );
    $rfh->print("$tag=$val\n");
}

sub max($$) {
    my ( $a, $b ) = @_;
    return ( $a > $b ) ? $a : $b;
}

sub min($$) {
    my ( $a, $b ) = @_;
    return $a if ( $a != 0 && $b == 0 );
    return $b if ( $a == 0 && $b != 0 );
    return ( $a < $b ) ? $a : $b;
}

sub calcNStats($$$$$) {
    my ( $type, $gsz, $lens, $sortContig, $Results ) = @_;
    my $prevSz   = 0;
    my $prevIncr = 0;
    my $incrSize = $gsz / 10;
    my @N        = ( .25, .5, .75 );
    my $sum      = 0;
    if ( $gsz < 10000000 ) {    # 10M
        $incrSize = 1000000;
    }
    elsif ( $gsz < 100000000 ) {    # 100M
        $incrSize = 10000000;
    }
    elsif ( $gsz < 1000000000 ) {    # 1G
        $incrSize = 100000000;
    }
    elsif ( $gsz < 10000000000 ) {    # 10G
        $incrSize = 1000000000;
    }
    foreach my $cc (@$sortContig) {
        my $len = $lens->{$cc};
        $sum += $len;

  # use while loops instead of if's since big contig/scaffold can span divisions
        while ( @N && $sum > $gsz * $N[0] ) {
            push @{ $Results->{"N${type}Bases"} }, [ shift(@N), $len ];
        }
        while ( int( $sum / $incrSize ) > $prevIncr ) {
            $prevIncr++;
            if ( $len != $prevSz ) {
                push @{ $Results->{"Incr${type}Bases"} },
                  [ $prevIncr * $incrSize, $len ];
                $prevSz = $len;
            }
        }
    }
}

my %bad_mates = ();
my @MDI_Rec   = ();

sub trackBadMates($$$$$) {
    my ( $type, $uid, $size, $mid, $sfs ) = @_;
    my $orient = '+';
    my ( $l, $r ) = split ',', $$sfs{'pos'};
    if ( $l > $r ) {
        ( $l, $r ) = ( $r, $l );
        $orient = '-';
    }
    my $status = $bad_mates{$mid}[0];
    my $libIID = $bad_mates{$mid}[1];
    my $mateExtent;
    if ( @{ $bad_mates{$mid} } > 2 ) {
        $status     = $bad_mates{$mid}[-2];
        $mateExtent = $bad_mates{$mid}[-1];
    }
    else {
        my $stddev = $MDI_Rec[$libIID]{std};
        my $mean   = $MDI_Rec[$libIID]{mea};
        if ( !$stddev ) {
            print STDERR "undef stddev for lib $libIID, frag $mid\n";
        }

        # calculate max range withing 5 stddevs of mean
        if ( $orient eq '+' ) {
            $mateExtent = int( $l + $mean + 5 * $stddev );
        }
        else {
            $mateExtent = int( $r - $mean - 5 * $stddev );
        }
    }
    $bad_mates{$mid} =
      [ "$type:$uid", $size, $l, $r, $orient, $status, $mateExtent ];
}
my @LibIIDRange = ();

sub getLibIID($) {
    my ($frgIID) = @_;

    # ranges should be ordered, so first end > fragIID should give correct lib
    for ( my $i = 1 ; $i < @LibIIDRange ; $i++ ) {
        return $i if $frgIID <= $LibIIDRange[$i];
    }
    return -1;
}

sub outputHelp($) {
    my $helpflag = shift;
    if ( lc($helpflag) eq 'all' ) {
        foreach (@help_lines) {
            print;
        }
    }
    elsif ( $helpflag eq '' ) {
        print STDERR $MY_HELPTEXT;
    }
    else {
        my %help_text_headers    = ();
        my %help_text_headers_lc = ();
        my %help_text_map        = ();
        my %help_text_map_lc     = ();

        my $current_header = undef;
        my @current_array  = undef;
        foreach my $line (@help_lines) {
            if ( $line =~ /^\[(.*)\]/ ) {    #Section
                if ( defined $current_header ) {
                    $help_text_headers{$current_header} = [@current_array];
                    my @fields = @{ $help_text_headers{$current_header} };
                    @current_array = ();
                }
                $current_header = $1;
            }
            elsif ( $line =~ /^(\S+)\s*-\s*(.*)$/ ) {
                $help_text_map{$1} = $2;
                if ( defined $current_header ) {
                    push @current_array, $1;
                }
            }
        }
        foreach my $key ( keys %help_text_map ) {
            $help_text_map_lc{ lc $key } = $key;
        }
        foreach my $key ( keys %help_text_headers ) {
            $help_text_headers_lc{ lc $key } = $key;
        }

        if ( exists $help_text_headers_lc{ lc $helpflag } ) {
            $helpflag = $help_text_headers_lc{ lc $helpflag };
            print STDOUT "[$helpflag]\n";
            my @fields = @{ $help_text_headers{$helpflag} };
            foreach my $field (@fields) {
                my $value = undef;
                $value = $help_text_map{$field} if ( defined $field );
                print STDOUT "$field - $value\n" if ( defined $value );
            }
        }
        elsif ( exists $help_text_map_lc{ lc $helpflag } ) {
            $helpflag = $help_text_map_lc{ lc $helpflag };
            my $help_text = $help_text_map{$helpflag};
            print STDOUT "$helpflag - $help_text\n";
        }
        else {
            print STDOUT
"Term '$helpflag' is not a valid help search.  To see all terms type use '-help all'.\n";
        }
    }
}

sub outputHtml() {
    my $HTML_HEAD = qq~
    <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN" "http://www.w3.org/TR/REC-html40/loose.dtd">
    <html>
    <head>
      <meta http-equiv="Content-Type"
     content="text/html; CHARSET=iso-iso-8859-1">
      <style type="text/css"><!--
    BODY { font-family: serif }
    H1 { font-family: sans-serif }
    H2 { font-family: sans-serif }
    --></style>
    </head>
    <body>
    <center>
    OUTPUT GUIDE
    <h1>CAQC</h1>
    A description of the outputs of qc output of the <a href="../../projects/CeleraAssembler/">Celera Assembler</a>.
    </center>
    <hr>
    <ul>  
    ~;

    my $HTML_TAIL = qq~
    </ul>
    </body>
    </html>
    ~;

    print STDOUT $HTML_HEAD;
    foreach my $line (@help_lines) {
        if ( $line =~ /^\[.*\]/ ) {    #Section
            print STDOUT "<p><b>$line</b></p><br>";
        }
        elsif ( $line =~ /^(\S+)\s*-\s*(.*)$/ ) {
            print STDOUT "<li><b>$1</b> - $2</li><br>";
        }
    }
    print STDOUT $HTML_TAIL;
}

sub readHelpFile() {
    my $help_fh = new IO::File("< $help_file")
      or die "Unable to read $help_file. Error: '$!'";
    @help_lines = <$help_fh>;
}

MAIN:
{
    my $infile;
    my $delim = undef;
    my $sep   = undef;
    my $genomesize;
    my $topCount = undef;
    my $frg      = undef;
    my $metrics  = undef;
    my $helpflag = undef;
    my $version  = undef;
    my $html     = undef;

    Getopt::Long::Configure('no_ignore_case');
    my $err = GetOptions(
        'i=s'         => \$infile,
        'minqual=i'   => \$MINQUAL,
        'mincontig=i' => \$MINCONTIG,
        'd=s'         => \$delim,
        's=s'         => \$sep,
        'silent!'     => \$silent,
        'euid!'       => \$euid,
        'g=i'         => \$genomesize,
        't=i'         => \$topCount,
        'metrics'     => \$metrics,
        'frg'         => \$frg,
        'help:s'      => \$helpflag,
        'V|v'         => \$version,
        'html'        => \$html,

    );

    if ( $err == 0 ) {
        die "Try -h for help.\n";
    }

    if ($version) {
        print STDERR $MY_VERSION . "\n";
        exit(1);
    }

    if ( defined $helpflag ) {
        readHelpFile();
        outputHelp($helpflag);
        exit(1);
    }

    if ($html) {
        readHelpFile();
        outputHtml();
        exit(1);
    }

    $d = ( defined $delim ) ? $delim : $d;
    $s = ( defined $sep )   ? $sep   : $s;
    $topCount = $DEFAULT_TOP_COUNT if ( !defined $topCount );

    if ( !defined $infile ) {
        if ( $#ARGV < 0 ) {
            print STDERR $MY_HELPTEXT;
            die("\nMust specify an input file name.\n");
        }
        else {
            $infile = $ARGV[0];
        }
    }
    my $mateLinkRangeFile = 'mateLinkIIDRanges.txt';
    if ( -s $mateLinkRangeFile ) {
        open( MLRF, "<$mateLinkRangeFile" )
          || die "Read of $mateLinkRangeFile failed.";

        # Lines look like: Distance 1 has IID range 1 to 5675
        # as created from mateLinkIIDRanges.rb
        while (<MLRF>) {
            chop;
            my ( undef, $libNum, undef, undef, undef, $b, undef, $e ) = split;
            $LibIIDRange[$libNum] = $e;
        }
        close MLRF;
    }

   #If input is just the prefix i.e. "gbaf" it will append the ".asm",
   #otherwise, if input already has the ".asm" extension, it will not be touched
    if ( $infile !~ /\.asm$/ ) {
        $infile .= '.asm';
    }
    open( IN, $infile ) or die("Cannot open $infile: ($!)");

    my ( $prefix, $path, $suffix ) = fileparse( $infile, '.asm' );

    my $frgfile    = "$path/$prefix.frg";
    my $catmapfile = "$path/$prefix.catmap";
    my $frgOpen    = undef;
    my $catOpen    = undef;
    if ( defined $frg or defined $metrics ) {
        $frgOpen = open( FRG, $frgfile )    or die("Cannot open $frgfile");
        $catOpen = open( CAT, $catmapfile ) or warn("Cannot open $catmapfile");
    }

    my $record;

    my %lens;               # contig lengths wihtout gaps
    my %sizes;              # contig sizes with gaps
    my %gcs;                # g,c content of the contig consensus
    my %seqs;               # contig number of reads
    my %scaflen;            # sum of contig sizes in scaffold
    my %scafcontig;         # number of contigs in scaffold with key scaffold_id
    my %gc_scaffcontigs;    # contigs in scaffolds with key contigid
    my %inscaff;            # contigs in scaffolds
    my %adjscaflen;         # length of scaffold including inner-scaffold gaps
    my %readlens;           # sum of lengths of reads in contig
    my %utgFrags;           # frags in surrogate unitigs
    my %libs_initial;       # the libraries for a project before assembly
    my %libs_final;         # the libraries for a project after assembly

    my $totalCLRReadLengthASM =
      0; # Aggregate length of CLR of input reads used in assembly from ASM File
    my $totalSeqs = 0;    # Total number of reads used in assembly

    my $ncontigs                = 0;
    my $nunitigs                = 0;
    my $uunitigs                = 0;     #Daniela Nov 16th 2005
    my $numSingletons           = 0;     #VI
    my $numVarRecords           = 0;     #VI
    my $lenPlacedSurrogates     = 0;     #VI
    my $numPlacedSurrogates     = 0;     #VI
    my $numPlacedSurrogateFrags = 0;     #VI
    my $lenPlacedSurrogateFrags = 0;     #VI
    my $utgSeqs                 = 0;
    my $uutgSeqs                = 0;     #Daniela Nov 16th 2005
    my $surrReadLen             = 0;     # reads in surrogates
    my $contigReadLen           = 0;     # reads in contigs
    my $degenReadLen            = 0;     # reads in degenerates
    my $singReadLen             = 0;     # reads in singletons
    my %utglens                 = ();
    my $utglensTotal            = 0;
    my %uutglens                = ();    #Daniela Nov 16th 2005
    my $utgs    = Statistics::Descriptive::Sparse->new();
    my $uutgs   = Statistics::Descriptive::Sparse->new(); #Daniela Nov 16th 2005
    my %Results = ();
    $Results{ReadsWithChaffMate}      = 0;
    $Results{ReadsWithBadShortMate}   = 0;
    $Results{ReadsWithBothDegenMate}  = 0;
    $Results{ReadsWithDegenMate}      = 0;
    $Results{ReadsWithDiffScafMate}   = 0;
    $Results{ReadsWithGoodMate}       = 0;
    $Results{ReadsWithBothChaffMate}  = 0;
    $Results{ReadsWithBadLongMate}    = 0;
    $Results{ReadsWithNoMate}         = 0;
    $Results{ReadsWithOuttieMate}     = 0;
    $Results{ReadsWithSurrogateMate}  = 0;
    $Results{ReadsWithSameOrientMate} = 0;
    $Results{ReadsWithBothSurrMate}   = 0;
    $Results{ReadsWithUnassignedMate} = 0;

    my %readLen;
    my $scaffLinkWeights  = 0;    # sum of links connecting scaffolds
    my $Kb_scaff_count    = 0;
    my $Kb_scaff_tot_span = 0;

    my %cats = ();

    if ($catOpen) {
        while ( defined( my $cat_line = <CAT> ) ) {
            chomp($cat_line);
            my @maps = split( / /, $cat_line );
            $cats{ $maps[0] } = $maps[1];
        }
    }

    my $totalScaffolds          = 0;
    my $totalContigsInScaffolds = 0;
    my $intraScaffoldGaps       = 0;
    my $minContigsPerScaffold   = 0;
    my $maxContigsPerScaffold   = 0;
    my $totalSpanOfScaffolds    = 0;
    my $minScaffoldSpan         = 0;
    my $maxScaffoldSpan         = 0;
    my $totalBasesInScaffolds   = 0;
    my $minBasesInScaffolds     = 0;
    my $maxBasesInScaffolds     = 0;

    while ( $record = getCARecord( \*IN ) ) {

        my $type;
        if ( $record =~ m/^\{?(\w+)/ ) {
            $type = $1;
        }

        if ( $type eq 'CCO' ) {
            my ( $type, $fields, $recs ) = parseCARecord($record);
            $ncontigs++;
            my $contiglen = $$fields{'len'};
            my $nreads    = $$fields{'npc'};
            my $qual      = $$fields{'qlt'};
            my ( $contig_euid, $contigid ) = getBothCAIds( $$fields{'acc'} );

            $sizes{$contig_euid} = $contiglen;
            $lens{$contig_euid}  = $contiglen;
            my $ngaps    = $$fields{'cns'} =~ tr/-/-/;
            my $n_gbases = $$fields{'cns'} =~ tr/G/G/;
            my $n_cbases = $$fields{'cns'} =~ tr/C/C/;

            $lens{$contig_euid} -= $ngaps;

            $seqs{$contig_euid} = $nreads;
            $gcs{$contig_euid}  = $n_gbases + $n_cbases;

            foreach my $rec (@$recs) {
                my ( $sid, $sfs, $srecs ) = parseCARecord($rec);
                if ( $sid eq 'MPS' ) {

                    my $mid  = $$sfs{'mid'};
                    my $rlen = $readLen{$mid};    #abs($r - $l);
                    $readlens{$contig_euid} += $rlen;
                    if ( exists $utgFrags{$mid} ) {
                        $numPlacedSurrogateFrags++;
                        $lenPlacedSurrogateFrags += $rlen;
                    }
                    if ( exists $bad_mates{$mid} ) {
                        trackBadMates( 'CCO', $contig_euid, $contiglen, $mid,
                            $sfs );
                    }
                }
                elsif ( $sid eq 'UPS' ) {
                    my $typ = $$sfs{'typ'};
                    if ( $typ eq 'S' ) {
                        my ( $ub, $ue ) = split ',', $$sfs{'pos'};
                        $lenPlacedSurrogates += abs( $ue - $ub );
                        $numPlacedSurrogates++;
                    }
                }
                elsif ( $sid eq 'VAR' ) {
                    $numVarRecords++;
                }
            }
        }    # if $type = CCO

        elsif ( $type eq 'UTG' ) {
            my ( $type, $fields, $recs ) = parseCARecord($record);

            #Daniela Nov 16th 2005
            if ( $$fields{'sta'} eq 'U' )    # Daniela Jan 24th 2006
            {
                $uunitigs++;
                $uutgs->add_data( $$fields{'len'} );
                $uutgSeqs += $$fields{'nfr'};
                $uutglens{ getCAId( $$fields{'acc'} ) } = $$fields{'len'};
            }
            elsif ( $$fields{'sta'} ne 'N' ) {

                $nunitigs++;
                my $utgLen = $$fields{'len'};
                $utgs->add_data($utgLen);
                $utgSeqs += $$fields{'nfr'};
                $utglens{ getCAId( $$fields{'acc'} ) } = $utgLen;
                $utglensTotal += $utgLen;
                for ( my $rd = 0 ; $rd <= $#$recs ; $rd++ ) {
                    my ( $sid, $sfs, $srecs ) = parseCARecord( $$recs[$rd] );
                    if ( $sid eq 'MPS' ) {

                        my $mid = $$sfs{'mid'};
                        $surrReadLen += $readLen{$mid};
                        $utgFrags{$mid} = 1;
                        if ( exists $bad_mates{$mid} ) {
                            my ( $uuid, $uiid ) =
                              getBothCAIds( $$fields{'acc'} );
                            trackBadMates( 'UTG', $uuid, $utgLen, $mid, $sfs );
                        }
                    }
                }
            }

        }

        elsif ( $type eq 'SCF' ) {
            my ( $type, $fields, $recs ) = parseCARecord($record);
            my ( $scf_euid, $scfid ) = getBothCAIds( $$fields{'acc'} );
            my %scaffcontigs = ();
            my $adj          = 0;
            for ( my $i = 0 ; $i <= $#$recs ; $i++ ) {
                my ( $lrec, $lfield, $lrecs ) = parseCARecord( $$recs[$i] );

                $scaffcontigs{ $$lfield{'ct1'} }++;
                $scaffcontigs{ $$lfield{'ct2'} }++;
                $gc_scaffcontigs{ $$lfield{'ct1'} }++;
                $gc_scaffcontigs{ $$lfield{'ct2'} }++;
                $inscaff{ $$lfield{'ct1'} } = 1;
                $inscaff{ $$lfield{'ct2'} } = 1;
                $adj += $$lfield{'mea'};
            }

            my $contigsInScaffCount = 0;
            my $scaflenTmp          = 0;
            while ( my ( $ct, $freq ) = each %scaffcontigs ) {
                $contigsInScaffCount++;
                $scaflenTmp += $lens{$ct};
            }
            $scafcontig{$scf_euid} = $contigsInScaffCount;
            $scaflen{$scf_euid}    = $scaflenTmp;

            $adjscaflen{$scf_euid} = $scaflenTmp + int($adj);

            if ( ( $adjscaflen{$scf_euid} ) >= 2000 ) {
                $Kb_scaff_count++;
                $Kb_scaff_tot_span =
                  $Kb_scaff_tot_span + $adjscaflen{$scf_euid};
            }

            $totalScaffolds++;

            $totalContigsInScaffolds += $contigsInScaffCount;
            $intraScaffoldGaps       += $contigsInScaffCount - 1;

            $minContigsPerScaffold =
              min( $minContigsPerScaffold, $contigsInScaffCount );
            $maxContigsPerScaffold =
              max( $maxContigsPerScaffold, $contigsInScaffCount );
            $totalSpanOfScaffolds += $adjscaflen{$scf_euid};

            $minScaffoldSpan = min( $minScaffoldSpan, $adjscaflen{$scf_euid} );
            $maxScaffoldSpan = max( $maxScaffoldSpan, $adjscaflen{$scf_euid} );

            $totalBasesInScaffolds += $scaflenTmp;

            $minBasesInScaffolds = min( $minBasesInScaffolds, $scaflenTmp );
            $maxBasesInScaffolds = max( $maxBasesInScaffolds, $scaflenTmp );

        }    # if $type = SCF

        elsif ( $type eq 'SLK' ) {    # scaffold link
            my ( $type, $fields, $recs ) = parseCARecord($record);
            $Results{TotalScaffoldLinks}++;
            $scaffLinkWeights += $$fields{num};
        }

        elsif ( $type eq 'AFG' ) {
            my ( $type, $fields, $recs ) = parseCARecord($record);
            my ( $frag_id, $fragIID ) = getBothCAIds( $$fields{'acc'} );
            my $clrs = $$fields{'clr'};
            my ( $clrl, $clrr ) = split /,/, $clrs;
            my $length = $clrr - $clrl;
            $totalCLRReadLengthASM += $length;
            $readLen{$frag_id} = $length;
            $totalSeqs++;

            my $mstField = $$fields{'mst'};
            $bad_mates{$frag_id} = [ $mstField, getLibIID($fragIID) ];
            if ( $mstField eq 'A' ) {
                $Results{ReadsWithChaffMate}++;
            }
            elsif ( $mstField eq 'C' ) {
                $Results{ReadsWithBadShortMate}++;
            }
            elsif ( $mstField eq 'D' ) {
                $Results{ReadsWithBothDegenMate}++;
            }
            elsif ( $mstField eq 'E' ) {
                $Results{ReadsWithDegenMate}++;
            }
            elsif ( $mstField eq 'F' ) {
                $Results{ReadsWithDiffScafMate}++;
            }
            elsif ( $mstField eq 'G' ) {
                $Results{ReadsWithGoodMate}++;
                delete $bad_mates{$frag_id};
            }
            elsif ( $mstField eq 'H' ) {
                $Results{ReadsWithBothChaffMate}++;
            }
            elsif ( $mstField eq 'L' ) {
                $Results{ReadsWithBadLongMate}++;
            }
            elsif ( $mstField eq 'N' ) {
                $Results{ReadsWithNoMate}++;
                delete $bad_mates{$frag_id};
            }
            elsif ( $mstField eq 'O' ) {
                $Results{ReadsWithOuttieMate}++;
            }
            elsif ( $mstField eq 'R' ) {
                $Results{ReadsWithSurrogateMate}++;
            }
            elsif ( $mstField eq 'S' ) {
                $Results{ReadsWithSameOrientMate}++;
            }
            elsif ( $mstField eq 'U' ) {
                $Results{ReadsWithBothSurrMate}++;
            }
            elsif ( $mstField eq 'Z' ) {
                $Results{ReadsWithUnassignedMate}++;
            }
            else {
                $Results{ReadsWithUnassignedMate}++;
            }

            $numSingletons++ if $$fields{'cha'} == 1;
        }

        elsif ( $type eq 'MDI' ) {
            my ( $type, $fields, $recs ) = parseCARecord($record);
            my ( $id, $iid ) = getBothCAIds( $$fields{'ref'} );
            my $mea  = $$fields{'mea'};
            my $std  = $$fields{'std'};
            my $min  = $$fields{'min'};
            my $max  = $$fields{'max'};
            my $buc  = $$fields{'buc'};
            my $hist = $$fields{'his'};
            $hist =~ s/\s*$//g;
            $hist =~ s/\n/,/g;

            $MDI_Rec[$iid] = { 'mea' => $mea, 'std' => $std };

            if ($metrics) {
                my $lib_id = undef;

                if ( defined $cats{$id} ) {
                    $lib_id = $cats{$id};
                }
                elsif ( !-e ($catmapfile) ) {
                    $lib_id = $id;
                }

                if ( defined $lib_id ) {
                    $libs_final{$lib_id}->{'mea'}  = $mea;
                    $libs_final{$lib_id}->{'std'}  = $std;
                    $libs_final{$lib_id}->{'min'}  = $min;
                    $libs_final{$lib_id}->{'max'}  = $max;
                    $libs_final{$lib_id}->{'buc'}  = $buc;
                    $libs_final{$lib_id}->{'hist'} = $hist;
                }
            }
        }
    }    # while $record

    my $totalBases            = 'NA';
    my $totalCLRReadLengthFRG = 'NA'
      ;  # Aggregate length of CLR of input reads used in assembly from FRG File
    my $totalReadsInFRG =
      'NA';    # Aggregate number of reads used in assembly from FRG File

    if ( $frgOpen or $metrics ) {

        $totalBases            = 0;
        $totalCLRReadLengthFRG = 0;
        $totalReadsInFRG       = 0;

        while ( $record = getCARecord( \*FRG ) ) {
            my $type;
            if ( $record =~ m/^\{?(\w+)/ ) {
                $type = $1;
            }

            if ( $type eq 'DST' and $metrics ) {
                my ( $type, $fields, $recs ) = parseCARecord($record);
                my $id     = $$fields{'acc'};
                my $lib_id = undef;

                if ( defined $cats{$id} ) {
                    $lib_id = $cats{$id};
                }
                elsif ( !-e ($catmapfile) ) {
                    $lib_id = $id;
                }
                if ( defined $lib_id ) {
                    my $mea = $$fields{'mea'};
                    my $std = $$fields{'std'};
                    $libs_initial{$lib_id}->{'mea'} = $mea;
                    $libs_initial{$lib_id}->{'std'} = $std;
                }
            }
            elsif ( $type eq 'FRG' ) {
                my ( $type, $fields, $recs ) = parseCARecord($record);
                my $seq = $$fields{'seq'};
                if (!defined $seq) {
                    print STDERR "Seq not defined for record:$record\n";
                    next;
                }
                my $seqLen = $seq =~ tr/a-zA-Z//;
                $totalBases += $seqLen;
                my $clrStr = $$fields{'clr'};
                if (!defined $clrStr) {
                    print STDERR "clr not defined for record:$record\n";
                    next;
                }
                if ( $clrStr =~ /^(\d+),(\d+)$/ ) {
                    my $clearRangeLen = $2 - $1 + 1;
                    $totalCLRReadLengthFRG += $clearRangeLen;
                }
                $totalReadsInFRG++;
            }
        }
    }

    # Initialize the Results hash so that we never get blanks for values
    $Results{BasesCount}               = $totalBases;
    $Results{ClearRangeLengthFRG}      = $totalCLRReadLengthFRG;
    $Results{ClearRangeLengthASM}      = $totalCLRReadLengthASM;
    $Results{SurrogateBaseLength}      = $surrReadLen;
    $Results{ContigBaseLength}         = 0;
    $Results{DegenBaseLength}          = 0;
    $Results{SingletonBaseLength}      = 0;
    $Results{Contig_SurrBaseLength}    = 0;
    $Results{TotalScaffolds}           = 0;
    $Results{TotalContigsInScaffolds}  = 0;
    $Results{MeanContigsPerScaffold}   = 0.0;
    $Results{MinContigsPerScaffold}    = 0;
    $Results{MaxContigsPerScaffold}    = 0;
    $Results{TotalBasesInScaffolds}    = 0;
    $Results{MeanBasesInScaffolds}     = 0.0;
    $Results{MinBasesInScaffolds}      = 0;
    $Results{MaxBasesInScaffolds}      = 0;
    $Results{TotalSpanOfScaffolds}     = 0;
    $Results{MeanSpanOfScaffolds}      = 0.0;
    $Results{MinScaffoldSpan}          = 0;
    $Results{MaxScaffoldSpan}          = 0;
    $Results{IntraScaffoldGaps}        = 0;
    $Results{MeanSequenceGapLength}    = 0.0;
    $Results{'2KbScaffolds'}           = 0;
    $Results{'2KbScaffoldsSpan'}       = 0;
    $Results{MeanContigLength}         = 0.0;
    $Results{MinContigLength}          = 0;
    $Results{MaxContigLength}          = 0;
    $Results{TotalBigContigs}          = 0;
    $Results{BigContigLength}          = 0;
    $Results{MeanBigContigLength}      = 0.0;
    $Results{MinBigContig}             = 0;
    $Results{MaxBigContig}             = 0;
    $Results{BigContigsPercentBases}   = 0.0;
    $Results{TotalSmallContigs}        = 0;
    $Results{SmallContigLength}        = 0;
    $Results{MeanSmallContigLength}    = 0.0;
    $Results{MinSmallContig}           = 0;
    $Results{MaxSmallContig}           = 0;
    $Results{SmallContigsPercentBases} = 0.0;
    $Results{TotalDegenContigs}        = 0;
    $Results{DegenContigLength}        = 0;
    $Results{MeanDegenContigLength}    = 0.0;
    $Results{MinDegenContig}           = 0;
    $Results{MaxDegenContig}           = 0;
    $Results{DegenPercentBases}        = 0.0;
    $Results{TotalUsableReads}         = 0;
    $Results{TotalReadsInput}          = $totalReadsInFRG;
    $Results{AvgClearRange}            =
      ( $totalSeqs > 0 ) ? $totalCLRReadLengthASM / $totalSeqs : 0;
    $Results{ContigReads}               = 0;
    $Results{BigContigReads}            = 0;
    $Results{SmallContigReads}          = 0;
    $Results{DegenContigReads}          = 0;
    $Results{SurrogateReads}            = 0;
    $Results{SingletonReads}            = 0;
    $Results{SingletonPercentBases}     = 0.0;
    $Results{ChaffReads}                = $numSingletons;
    $Results{TotalVarRecords}           = $numVarRecords;
    $Results{AllReads}                  = 0.0;
    $Results{ContigsOnly}               = 0.0;
    $Results{Contigs_Surrogates}        = 0.0;
    $Results{Contigs_Degens_Surrogates} = 0.0;
    $Results{TotalSurrogates}           = $nunitigs;
    $Results{SurrogatesPercentBases}    = 0.0;
    $Results{SurrogateInstances}        = $numPlacedSurrogates;
    $Results{SurrogateLength}           = $utgs->sum();
    $Results{SurrogateInstanceLength}   = $lenPlacedSurrogates;
    $Results{UnPlacedSurrReadLen}  = $surrReadLen - $lenPlacedSurrogateFrags;
    $Results{PlacedSurrReadLen}    = $lenPlacedSurrogateFrags;
    $Results{PlacedSurrogateReads} = $numPlacedSurrogateFrags;
    $Results{MeanSurrogateLength}  = $utgs->mean();
    $Results{MinSurrogateLength}   = $utgs->min();
    $Results{MaxSurrogateLength}   = $utgs->max();
    $Results{SDSurrogateLength}    = $utgs->standard_deviation();
    $Results{TotalUUnitigs}        = $uunitigs;
    $Results{MeanUUnitigLength}    = $uutgs->mean();
    $Results{MinUUnitigLength}     = $uutgs->min();
    $Results{MaxUUnitigLength}     = $uutgs->max();
    $Results{SDUUnitigLength}      = $uutgs->standard_deviation();

    if ( !exists $Results{ReadsWithNoMate} ) {
        $Results{ReadsWithNoMate} = 0;
    }
    if ( !exists $Results{ReadsWithBadMate} ) {
        $Results{ReadsWithBadMate} = 0;
    }
    if ( !exists $Results{ReadsWithGoodMate} ) {
        $Results{ReadsWithGoodMate} = 0;
    }
    if ( !exists $Results{TotalScaffoldLinks} ) {
        $Results{TotalScaffoldLinks}     = 0;
        $Results{MeanScaffoldLinkWeight} = 0;
    }
    else {
        $Results{MeanScaffoldLinkWeight} =
          $scaffLinkWeights / $Results{TotalScaffoldLinks};
    }

#Begin Marwan 10-17-06:
#Statistics::Descriptive::Sparse->new() functions return undef if no value exists
    if ( !defined $Results{MeanSurrogateLength} ) {
        $Results{MeanSurrogateLength} = 0;
    }
    if ( !defined $Results{MinSurrogateLength} ) {
        $Results{MinSurrogateLength} = 0;
    }
    if ( !defined $Results{MaxSurrogateLength} ) {
        $Results{MaxSurrogateLength} = 0;
    }
    if ( !defined $Results{SDSurrogateLength} ) {
        $Results{SDSurrogateLength} = 0;
    }

    if ( !defined $Results{MeanUUnitigLength} ) {
        $Results{MeanUUnitigLength} = 0;
    }
    if ( !defined $Results{MinUUnitigLength} ) {
        $Results{MinSurrogateLength} = 0;
    }
    if ( !defined $Results{MaxUUnitigLength} ) {
        $Results{MaxSurrogateLength} = 0;
    }
    if ( !defined $Results{SDUUnitigLength} ) {
        $Results{SDSurrogateLength} = 0;
    }

    #End Marwan 10-17-06

    my $totlen  = 0;
    my $totseqs = 0;

    my $biglen  = 0;
    my $bigseqs = 0;
    my $nbigs   = 0;

    my $maxcontig = 0;

    for ( my $i = 0 ; $i < $ncontigs ; $i++ ) {
        my ( $contig_id, $len ) = each(%lens);

        my $number;
        my $min;
        my $max;
        my $totlen;
        my $nseq;

        if ( exists $inscaff{$contig_id} ) {    # the good guys
            if ( $len >= $MINCONTIG ) {         # the good big guys
                $number = 'TotalBigContigs';
                $min    = 'MinBigContig';
                $max    = 'MaxBigContig';
                $totlen = 'BigContigLength';
                $nseq   = 'BigContigReads';
            }
            else {                              # the small good guys
                $number = 'TotalSmallContigs';
                $min    = 'MinSmallContig';
                $max    = 'MaxSmallContig';
                $totlen = 'SmallContigLength';
                $nseq   = 'SmallContigReads';
            }
            $contigReadLen += $readlens{$contig_id};
        }
        else {                                  # the chaff
            $number = 'TotalDegenContigs';
            $min    = 'MinDegenContig';
            $max    = 'MaxDegenContig';
            $totlen = 'DegenContigLength';
            $nseq   = 'DegenContigReads';
            $degenReadLen += $readlens{$contig_id};
        }

        $Results{$max} = max( $len, $Results{$max} );
        $Results{$min} = min( $len, $Results{$min} );
        $Results{$number}++;
        $Results{$totlen} += $len;
        $Results{$nseq}   += $seqs{$contig_id};
    }    # for $i < $ncontigs

    #This is the same thing as TotalContigsInScaffolds
    $Results{'TotalContigs'} =
      0 + $Results{'TotalBigContigs'} + $Results{'TotalSmallContigs'};

    #This is the same thing as TotalBasesInScaffolds
    $Results{'TotalBasesInContigs'} =
      0 + $Results{'BigContigLength'} + $Results{'SmallContigLength'};
    $Results{'MeanContigLength'} =
      ( $Results{'TotalContigs'} > 0 )
      ? $Results{'TotalBasesInContigs'} * 1.0 / $Results{'TotalContigs'}
      : 0;

    $Results{'ContigReads'} =
      0 + $Results{'BigContigReads'} + $Results{'SmallContigReads'};

    #How could this ever be MinBigContig?
    $Results{'MinContigLength'} =
      min( $Results{'MinSmallContig'}, $Results{'MinBigContig'} );

    #How could this ever be MaxBigContig?
    $Results{'MaxContigLength'} =
      max( $Results{'MaxBigContig'}, $Results{'MaxSmallContig'} );

    $Results{'MeanBigContigLength'} =
      ( $Results{'TotalBigContigs'} > 0 )
      ? $Results{'BigContigLength'} * 1.0 / $Results{'TotalBigContigs'}
      : 0;

    $Results{'MeanSmallContigLength'} =
      ( $Results{'TotalSmallContigs'} )
      ? $Results{'SmallContigLength'} * 1.0 / $Results{'TotalSmallContigs'}
      : 0;

    $Results{'BigContigsPercentBases'} =
      ( $Results{'TotalBasesInContigs'} > 0 )
      ? $Results{'BigContigLength'} * 100.0 / $Results{'TotalBasesInContigs'}
      : 0;

    $Results{'SmallContigsPercentBases'} =
      ( $Results{'TotalBasesInContigs'} > 0 )
      ? $Results{'SmallContigLength'} * 100.0 / $Results{'TotalBasesInContigs'}
      : 0;

    $Results{'MeanDegenContigLength'} =
      ( $Results{'TotalDegenContigs'} )
      ? $Results{'DegenContigLength'} * 1.0 / $Results{'TotalDegenContigs'}
      : 0;

    $Results{'DegenPercentBases'} =
      ( $Results{'TotalBasesInContigs'} > 0 )
      ? $Results{'DegenContigLength'} * 100.0 / $Results{'TotalBasesInContigs'}
      : 0;

    # compute N50 values and top $topCount guys (contigs)
    my @sortContig = sort { $lens{$b} <=> $lens{$a} } ( keys %inscaff );
    my $top5contig;
    my $sum       = 0;
    my $reads_tot = 0;
    my $bases_tot = 0;
    for ( my $cc = 0 ; $cc <= $#sortContig && $cc < $topCount ; $cc++ ) {
        $top5contig .= "$cc$d$seqs{$sortContig[$cc]}$s$lens{$sortContig[$cc]}";
        $top5contig .= "$s$sortContig[$cc]" if ($euid);
        $top5contig .= "\n";
        $reads_tot += $seqs{ $sortContig[$cc] };
        $bases_tot += $lens{ $sortContig[$cc] };
    }
    $top5contig .= "total$d$reads_tot$s$bases_tot\n";

    my $gsz;
    if ( !defined $genomesize ) {
        $gsz = $Results{'TotalBasesInContigs'};
    }
    else {
        $gsz = $genomesize;
    }
    calcNStats( 'Contig', $gsz, \%lens, \@sortContig, \%Results );

    $Results{'TotalScaffolds'}          = $totalScaffolds;
    $Results{'TotalContigsInScaffolds'} = $totalContigsInScaffolds;
    $Results{'IntraScaffoldGaps'}       = $intraScaffoldGaps;
    $Results{'MinContigsPerScaffold'}   = $minContigsPerScaffold;
    $Results{'MaxContigsPerScaffold'}   = $maxContigsPerScaffold;
    $Results{'TotalSpanOfScaffolds'}    = $totalSpanOfScaffolds;
    $Results{'MinScaffoldSpan'}         = $minScaffoldSpan;
    $Results{'MaxScaffoldSpan'}         = $maxScaffoldSpan;
    $Results{'TotalBasesInScaffolds'}   = $totalBasesInScaffolds;
    $Results{'MinBasesInScaffolds'}     = $minBasesInScaffolds;
    $Results{'MaxBasesInScaffolds'}     = $maxBasesInScaffolds;

    $Results{'2KbScaffolds'}    = $Kb_scaff_count;
    $Results{'2KbScaffoldSpan'} = $Kb_scaff_tot_span;

    $Results{'MeanContigsPerScaffold'} =
      ( $Results{'TotalScaffolds'} > 0 )
      ? $Results{'TotalContigsInScaffolds'} * 1.0 / $Results{'TotalScaffolds'}
      : 0.0;
    $Results{'MeanBasesInScaffolds'} =
      ( $Results{'TotalScaffolds'} > 0 )
      ? $Results{'TotalBasesInScaffolds'} * 1.0 / $Results{'TotalScaffolds'}
      : 0.0;
    $Results{'MeanSpanOfScaffolds'} =
      ( $Results{'TotalScaffolds'} > 0 )
      ? $Results{'TotalSpanOfScaffolds'} * 1.0 / $Results{'TotalScaffolds'}
      : 0.0;
    $Results{'MeanSequenceGapLength'} =
        ( $Results{'IntraScaffoldGaps'} > 0 )
      ? ( $Results{'TotalSpanOfScaffolds'} - $Results{'TotalBasesInScaffolds'} )
      * 1.0 / $Results{'IntraScaffoldGaps'}
      : 0.0;

    # compute N50 values and top $topCount guys (scaffolds)
    my @sortScaff = sort { $scaflen{$b} <=> $scaflen{$a} } ( keys %scaflen );
    my $top5scaff;
    $sum = 0;

    my $contigs_tot   = 0;
    my $gaps_tot      = 0;
    my $size_tot      = 0;
    my $span_tot      = 0;
    my $avgContig_tot = 0;
    my $avgGap_tot    = 0;

    for ( my $ss = 0 ; $ss <= $#sortScaff && $ss < $topCount ; $ss++ ) {
        my $scf = $sortScaff[$ss];
        $top5scaff .=
          "$ss$d$scafcontig{$scf}$s$scaflen{$scf}$s$adjscaflen{$scf}$s"
          . sprintf( '%d',
            Math::BigFloat->new( $scaflen{$scf} * 1.0 / $scafcontig{$scf} )
              ->ffround(0) )
          . $s
          . sprintf(
            '%d',
            ( $scafcontig{$scf} - 1 > 0 )
            ? (
                Math::BigFloat->new(
                    ( $adjscaflen{$scf} - $scaflen{$scf} ) * 1.0 /
                      ( $scafcontig{$scf} - 1 )
                  )->ffround(0)
              )
            : 0
          );

        $top5scaff .= "$s$scf" if ($euid);
        $top5scaff .= "\n";

        $contigs_tot += $scafcontig{$scf};
        $gaps_tot    += $scafcontig{$scf} - 1;
        $size_tot    += $scaflen{$scf};
        $span_tot    += $adjscaflen{$scf};
        my $avgContig_size =
          sprintf( '%.2f', $scaflen{$scf} * 1.0 / $scafcontig{$scf} );
        my $avgGap_size = sprintf(
            '%.2f',
            ( $scafcontig{$scf} - 1 > 0 )
            ? ( ( $adjscaflen{$scf} - $scaflen{$scf} ) * 1.0 /
                  ( $scafcontig{$scf} - 1 ) )
            : 0.0
        );

        $avgContig_tot += ( $avgContig_size * $scafcontig{$scf} );
        $avgGap_tot += ($avgGap_size) * ( $scafcontig{$scf} - 1 );
    }
    if ( $contigs_tot != 0 ) {
        $avgContig_tot = $avgContig_tot / $contigs_tot;
    }
    else {
        $avgContig_tot = 0;
    }

    if ( $gaps_tot != 0 ) {
        $avgGap_tot = $avgGap_tot / $gaps_tot;
    }
    else {
        $avgGap_tot = 0;
    }
    my $top5scaff_place_holder =
        "total$d$contigs_tot$s$size_tot$s$span_tot$s"
      . sprintf( '%d', Math::BigFloat->new("$avgContig_tot")->ffround(0) )
      . $s
      . sprintf( '%d', Math::BigFloat->new("$avgGap_tot")->ffround(0) ) . "\n";
    $top5scaff .= $top5scaff_place_holder;
    if ( !defined $genomesize ) {
        $gsz = $Results{'TotalBasesInScaffolds'};
    }
    else {
        $gsz = $genomesize;
    }
    calcNStats( 'Scaffold', $gsz, \%scaflen, \@sortScaff, \%Results );

    $Results{'AllReads'} =
      ( $gsz > 0 )
      ? sprintf( '%0.2f', $totalCLRReadLengthASM / $gsz )
      : 0.0;

    #Coverage:
    $Results{'ContigsOnly'} =
        ( $Results{'TotalBasesInScaffolds'} > 0 )
      ? ( ($contigReadLen) / $Results{'TotalBasesInScaffolds'} )
      : 0.0;

#$base->logLocal( "surrReadLen = $surrReadLen\ncontigReadLen = $contigReadLen\ndegenReadLen = $degenReadLen\ntotalReadLen = $totalCLRReadLengthASM", 2 );

    # Bases
    $Results{ContigBaseLength}      = $contigReadLen;
    $Results{DegenBaseLength}       = $degenReadLen;
    $Results{Contig_SurrBaseLength} =
      $Results{'UnPlacedSurrReadLen'} + $contigReadLen;
    $Results{SingletonBaseLength} =
      $totalCLRReadLengthASM - $degenReadLen - $contigReadLen -
      $Results{'UnPlacedSurrReadLen'};

    $Results{'Contigs_Degens_Surrogates'} =
      ( $Results{'TotalBasesInScaffolds'} + $Results{'DegenContigLength'} > 0 )
      ? ( ( $Results{'UnPlacedSurrReadLen'} + $contigReadLen + $degenReadLen ) /
          ( $Results{'TotalBasesInScaffolds'} + $Results{'DegenContigLength'} )
      )
      : 0.0;

    $Results{'Contigs_Surrogates'} =
      ( $Results{'TotalBasesInScaffolds'} > 0 )
      ? (
        $Results{'Contig_SurrBaseLength'} / $Results{'TotalBasesInScaffolds'} )
      : 0.0;

    # Surrogates
    $Results{'SurrogateReads'} = $utgSeqs;

    # Coverage
    $Results{'TotalUsableReads'} = $totalSeqs;
    $Results{'SingletonReads'}   =
      $Results{'TotalUsableReads'} - $Results{'BigContigReads'} -
      $Results{'SmallContigReads'} - $Results{'DegenContigReads'} -
      $Results{'SurrogateReads'} +
      $Results{'PlacedSurrogateReads'};    # Should be in contigs

    # Parameter settings
    $Results{'MinBigContigSizeParm'} = $MINCONTIG;
    $Results{'MinQualityParm'}       = $MINQUAL;

    # gc Content of contigs in scaffolds(placed contigs)
    # Daniela's change Nov 15th 2005
    $Results{'Content'} =
      calculateGC( \%gc_scaffcontigs, \%gcs, \%sizes ) * 100;

    # Emit the results
    my $fh = new IO::File("> $prefix.qc")
      or die("Could not open $prefix.qc ($!)");

    print STDERR "[Scaffolds]\n" if ( !$silent );
    $fh->print("[Scaffolds]\n");
    printl( 'TotalScaffolds',          \%Results, $fh );
    printl( 'TotalContigsInScaffolds', \%Results, $fh );
    printlf( 'MeanContigsPerScaffold', \%Results, $fh );
    printl( 'MinContigsPerScaffold', \%Results, $fh );
    printl( 'MaxContigsPerScaffold', \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");
    printl( 'TotalBasesInScaffolds', \%Results, $fh );
    printld( 'MeanBasesInScaffolds', \%Results, $fh );
    printl( 'MinBasesInScaffolds', \%Results, $fh );
    printl( 'MaxBasesInScaffolds', \%Results, $fh );

    for my $val ( @{ $Results{'NScaffoldBases'} } ) {
        printf $fh "N%2d%s=%d\n", $val->[0] * 100, 'ScaffoldBases', $val->[1];
        printf STDERR "N%2d%-29s%d\n", $val->[0] * 100, 'ScaffoldBases',
          $val->[1];
    }
    for my $val ( @{ $Results{'IncrScaffoldBases'} } ) {
        print $fh "ScaffoldAt$val->[0]=$val->[1]\n";
        printf STDERR "ScaffoldAt%-22s%d\n", $val->[0], $val->[1];
    }
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");
    printl( 'TotalSpanOfScaffolds', \%Results, $fh );
    printld( 'MeanSpanOfScaffolds', \%Results, $fh );
    printl( 'MinScaffoldSpan',   \%Results, $fh );
    printl( 'MaxScaffoldSpan',   \%Results, $fh );
    printl( 'IntraScaffoldGaps', \%Results, $fh );
    printl( '2KbScaffolds',      \%Results, $fh );
    printl( '2KbScaffoldSpan',   \%Results, $fh );
    printld( 'MeanSequenceGapLength', \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    my $top5scaff_header =
"[Top${topCount}Scaffolds${d}contigs${s}size${s}span${s}avgContig${s}avgGap";

    if ( !$silent ) {
        if ($euid) {
            $top5scaff_header .= "${s}EUID]\n";
        }
        else {
            $top5scaff_header .= "]\n";
        }
        print STDERR $top5scaff_header;
    }
    $fh->print($top5scaff_header);

    if ( !$silent ) {
        print STDERR $top5scaff;
    }

    $fh->print($top5scaff);
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[Contigs]\n" if ( !$silent );
    $fh->print("[Contigs]\n");
    printl( 'TotalContigsInScaffolds', \%Results, $fh );
    printl( 'TotalBasesInScaffolds',   \%Results, $fh );
    printl( 'TotalVarRecords',         \%Results, $fh );
    printld( 'MeanContigLength', \%Results, $fh );
    printl( 'MinContigLength', \%Results, $fh );
    printl( 'MaxContigLength', \%Results, $fh );
    for my $val ( @{ $Results{'NContigBases'} } ) {
        printf $fh "N%2d%s=%d\n", $val->[0] * 100, 'ContigBases', $val->[1];
        printf STDERR "N%2d%-29s%d\n", $val->[0] * 100, 'ContigBases',
          $val->[1];
    }
    for my $val ( @{ $Results{'IncrContigBases'} } ) {
        print $fh "ContigAt$val->[0]=$val->[1]\n";
        printf STDERR "ContigAt%-24s%d\n", $val->[0], $val->[1];
    }
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[BigContigs_greater_$MINCONTIG]\n" if ( !$silent );
    $fh->print("[BigContigs_greater_$MINCONTIG]\n");
    printl( 'TotalBigContigs', \%Results, $fh );
    printl( 'BigContigLength', \%Results, $fh );
    printld( 'MeanBigContigLength', \%Results, $fh );
    printl( 'MinBigContig', \%Results, $fh );
    printl( 'MaxBigContig', \%Results, $fh );
    printlf( 'BigContigsPercentBases', \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[SmallContigs]\n" if ( !$silent );
    $fh->print("[SmallContigs]\n");
    printl( 'TotalSmallContigs', \%Results, $fh );
    printl( 'SmallContigLength', \%Results, $fh );
    printld( 'MeanSmallContigLength', \%Results, $fh );
    printl( 'MinSmallContig', \%Results, $fh );
    printl( 'MaxSmallContig', \%Results, $fh );
    printlf( 'SmallContigsPercentBases', \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[DegenContigs]\n" if ( !$silent );
    $fh->print("[DegenContigs]\n");
    printl( 'TotalDegenContigs', \%Results, $fh );
    printl( 'DegenContigLength', \%Results, $fh );
    printld( 'MeanDegenContigLength', \%Results, $fh );
    printl( 'MinDegenContig', \%Results, $fh );
    printl( 'MaxDegenContig', \%Results, $fh );
    printlf( 'DegenPercentBases', \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    my $top5contig_header = "[Top${topCount}Contigs${d}reads${s}bases";
    if ( !$silent ) {
        if ($euid) {
            $top5contig_header .= "${s}EUID]\n";
        }
        else {
            $top5contig_header .= "]\n";
        }
        print STDERR $top5contig_header;
    }

    $fh->print($top5contig_header);
    if ( !$silent ) {
        print STDERR $top5contig;
    }
    $fh->print($top5contig);
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[UniqueUnitigs]\n" if ( !$silent );
    $fh->print("[UniqueUnitigs]\n");
    printl( 'TotalUUnitigs',    \%Results, $fh );
    printl( 'MinUUnitigLength', \%Results, $fh );
    printl( 'MaxUUnitigLength', \%Results, $fh );
    printld( 'MeanUUnitigLength', \%Results, $fh );
    printld( 'SDUUnitigLength',   \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[Surrogates]\n" if ( !$silent );
    $fh->print("[Surrogates]\n");
    printl( 'TotalSurrogates',         \%Results, $fh );
    printl( 'SurrogateInstances',      \%Results, $fh );
    printl( 'SurrogateLength',         \%Results, $fh );
    printl( 'SurrogateInstanceLength', \%Results, $fh );
    printl( 'UnPlacedSurrReadLen',     \%Results, $fh );
    printl( 'PlacedSurrReadLen',       \%Results, $fh );
    printl( 'MinSurrogateLength',      \%Results, $fh );
    printl( 'MaxSurrogateLength',      \%Results, $fh );
    printld( 'MeanSurrogateLength', \%Results, $fh );
    printld( 'SDSurrogateLength',   \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[Mates]\n" if ( !$silent );
    $fh->print("[Mates]\n");
    printp( 'ReadsWithNoMate',         \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithGoodMate',       \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithBadShortMate',   \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithBadLongMate',    \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithSameOrientMate', \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithOuttieMate',     \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithBothChaffMate',  \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithChaffMate',      \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithBothDegenMate',  \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithDegenMate',      \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithBothSurrMate',   \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithSurrogateMate',  \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithDiffScafMate',   \%Results, $fh, $totalSeqs );
    printp( 'ReadsWithUnassignedMate', \%Results, $fh, $totalSeqs );
    printl( 'TotalScaffoldLinks', \%Results, $fh );
    printlf( 'MeanScaffoldLinkWeight', \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[Reads]\n" if ( !$silent );
    $fh->print("[Reads]\n");
    printl( 'TotalReadsInput',  \%Results, $fh );
    printl( 'TotalUsableReads', \%Results, $fh );
    printld( 'AvgClearRange', \%Results, $fh );
    printp( 'ContigReads',          \%Results, $fh, $totalSeqs );
    printp( 'BigContigReads',       \%Results, $fh, $totalSeqs );
    printp( 'SmallContigReads',     \%Results, $fh, $totalSeqs );
    printp( 'DegenContigReads',     \%Results, $fh, $totalSeqs );
    printp( 'SurrogateReads',       \%Results, $fh, $totalSeqs );
    printp( 'PlacedSurrogateReads', \%Results, $fh, $totalSeqs );
    printp( 'SingletonReads',       \%Results, $fh, $totalSeqs );
    printp( 'ChaffReads',           \%Results, $fh, $totalSeqs );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[Coverage]\n" if ( !$silent );
    $fh->print("[Coverage]\n");
    printlf( 'ContigsOnly',               \%Results, $fh );
    printlf( 'Contigs_Surrogates',        \%Results, $fh );
    printlf( 'Contigs_Degens_Surrogates', \%Results, $fh );
    printlf( 'AllReads',                  \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[TotalBaseCounts]\n" if ( !$silent );
    $fh->print("[TotalBaseCounts]\n");
    printl( 'BasesCount',            \%Results, $fh );
    printl( 'ClearRangeLengthFRG',   \%Results, $fh );
    printl( 'ClearRangeLengthASM',   \%Results, $fh );
    printl( 'SurrogateBaseLength',   \%Results, $fh );
    printl( 'ContigBaseLength',      \%Results, $fh );
    printl( 'DegenBaseLength',       \%Results, $fh );
    printl( 'SingletonBaseLength',   \%Results, $fh );
    printl( 'Contig_SurrBaseLength', \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    print STDERR "[gcContent]\n" if ( !$silent );
    $fh->print("[gcContent]\n");
    printlf( 'Content', \%Results, $fh );
    print STDERR "\n" if ( !$silent );
    $fh->print("\n");

    $fh->close() or die("Could not close $prefix.qc ($!)");

    if ($metrics) {
        die("Could not copy $prefix.qc to $prefix.qc.metrics")
          if ( !copy( "$prefix.qc", "$prefix.qc.metrics" ) );

        my $fm = new IO::File(">>$prefix.qc.metrics")
          or die("Could not open $prefix.qc.metrics ($!)");

        my @lib_ids = keys %libs_initial;

        foreach my $lib_id (@lib_ids) {
            $fm->print("[library_$lib_id\_initial]\n");
            my $mea = $libs_initial{$lib_id}->{'mea'};
            $fm->print("mean=$mea\n");
            my $std = $libs_initial{$lib_id}->{'std'};
            $fm->print("sd=$std\n\n");

            $fm->print("[library_$lib_id\_final]\n");
            $mea = $libs_final{$lib_id}->{'mea'};
            $fm->print("mean=$mea\n");
            $std = $libs_final{$lib_id}->{'std'};
            $fm->print("sd=$std\n");
            my $min = $libs_final{$lib_id}->{'min'};
            $fm->print("min=$min\n");
            my $max = $libs_final{$lib_id}->{'max'};
            $fm->print("max=$max\n");
            my $buc = $libs_final{$lib_id}->{'buc'};
            $fm->print("buc=$buc\n");
            my $his = $libs_final{$lib_id}->{'hist'};
            $fm->print("hist=$his\n\n");
        }

        $fm->print("[ContigHistogram]\n");
        my @len_arr        = values(%sizes);
        my @sorted_len_arr = sort { $a <=> $b } @len_arr;
        my $len_list       = join( ',', @sorted_len_arr );
        $fm->print("hist=$len_list\n");
        $fm->close() or die("Could not close $prefix.qc.metrics ($!)");

    }
    my $badMateFile = 'badMateFragmentIDs.txt';

    # sort on UID of UTG, or CCO, then begin position, then end
    # UTGs should only be unplaced surrogates
    open( BADMATES, "| sort -k2.5,2 -k4,4n -k5,5n > $badMateFile" )
      || die "Couldn't write to $badMateFile.";
    local $\ = "\n";
    local $, = "\t";
    while ( my ( $badId, $coord ) = each %bad_mates ) {
        if ($coord) {
            print BADMATES $badId, @$coord;
        }
        else {
            print BADMATES $badId, 'Chaff';
        }
    }
    close BADMATES;

    exit(0);
}

