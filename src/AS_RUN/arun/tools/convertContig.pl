#!/usr/local/bin/perl -w

#
#  Make Perl behave better
#
use strict;

#
#  Use local modules first
#
use FindBin qw($Bin);
use lib "$Bin";

#
#  Use the XContig Suite
#
use TIGR::XContig;
use TIGR::XContig::Util;
use TIGR::XContig::Logger;

#
#  Other Imported Modules
#

#
#  Set up the Helptext
#
my $HELP = '
convertContig - Convert tiling data files to other formats.

USAGE: convertContig [options] <input> <ouput> [<output>...]

Any of the following formats are considered valid input, and will
be auto-selected based on the given filename patterns:

Contig files (*.contig)
XContig files (*.xcontig, *.xml)
Coninfo files (*.coninfo)
XIF files (*.xif)

Any of the following formats are considered valid output, and will
be auto-selected based on the given filename patterns:

Contig files (*.contig)
XContig files (*.xcontig, *.xml)
Tasm files (*.tasm)
Fasta files (*.fasta) [Writes tiling consensuses only]
Sequence files (*.seq)
Quality files (*.qual)
XIF Files (*.xif)

Any file of a valid input format can be converted to any valid output
format.  A single input file can be written to multiple output files.

EXAMPLES
----------------------------

Convert a .contig file to an .xcontig file

  convertContig in.contig out.xcontig

Convert an .xcontig file to a .contig file and a .xif file

  convertContig a.xcontig b.contig x.xif

Convert a .coninfo file to a full working set of files

  convertContig import.coninfo a.xcontig s.seq q.qual

Normalize a .contig file

  convertContig in.contig clean.contig

GENERAL OPTIONS
----------------------------
--help (-h)                     :  Display this message
--quiet (-q)                    :  Do not write status messages to the console
--input-type [type]             :  Interpret the input file as the given type

SUPPORTING INPUT
----------------------------
--seq [sequence file]           :  Specify a sequence file to load sequences
                                   from (more than one can be specified)
--gi [id]:[file]                :  Load the consensus sequence from the given
                                   file and associate it with the given GI.
--acc [id]:[file]               :  Load the consensus sequence from the given
                                   file and associate it with the given ACCESSION
--force-circ [id]               :  Force the given contig id to be circular

FORMAT-SPECIFIC OPTIONS
----------------------------
--comment [text]                :  Sets a comment string for each of the contigs
                                   in the output file, if the format supports it
--com-name [text]               :  Sets a com_name string for each of the contigs
                                   in the output file, if the format supports it
';

#
#  External Executables
#

#
#  Set up the Program Version
#
my $REVISION = (qw/$Revision: 1.1 $/)[1];
my $VERSION_FULL = " Version $XCONTIG_VERSION (Build " . $REVISION . ")";

#
#  Note Program Dependancies
#
my @DEPEND =
    (
     "TIGR::Foundation",
     "TIGR::XContig",
    );

#
#  Create a TIGR Foundation Object
#
my $tf = new TIGR::Foundation;

#
#  Set the Program Info in the Foundation Object
#
$tf->addDependInfo(@DEPEND);
$tf->setVersionInfo($VERSION_FULL);
$tf->setHelpInfo($HELP);

#
#  Use the XContig Suite for logging
#
TIGR::XContig::Logger::setTIGRFoundation($tf);
TIGR::XContig::Logger::setLogFile();

# ############################################################################
#
#   Global Variables
#
# ############################################################################

# Standard Options
my $quiet = 0;

my $input_type = undef;
my $comment = undef;
my $com_name = undef;
my $force_circ_spec = undef;
my @seqs = ();
my @gis = ();
my @accs = ();

# ############################################################################
#
#   Default State Setup
#
# ############################################################################

#
#  Set the messaging behavior
#
TIGR::XContig::Logger::setQuiet($quiet);

# ############################################################################
#
#   Command Line Processing
#
# ############################################################################

#
#  Get the Command Line Options
#
my $opt_result = $tf->TIGR_GetOptions(
                                      # Standard Options
                                      'q|quiet' => \$quiet,
                                      'input-type=s' => \$input_type,
                                      'seq=s@' => \@seqs,
                                      'gi=s@' => \@gis,
                                      'acc=s@' => \@accs,
                                      'comment=s' => \$comment,
                                      'com-name=s' => \$com_name,
                                      'force-circ=s' => \$force_circ_spec,
                                      );
#
#  Die unless the command line parsing went okay.
#
$tf->bail("Failure parsing command line options.") unless ($opt_result);

#
#  Fix the messaging behavior
#
TIGR::XContig::Logger::setQuiet($quiet);

# ############################################################################
#
#   Main Program
#
# ############################################################################

#
#  Load GI and ACCESSION Data
#
foreach my $gi (@gis)
{
    my ($gi_id, $file) = $gi =~ /^(.*?)\:(.*)/;
    if (defined $gi_id and defined $file)
    {
        TIGR::XContig::FastaData::loadFastaFile('GI', $gi_id, $file);
    }
}
foreach my $acc (@accs)
{
    my ($acc_id, $file) = $acc =~ /^(.*?)\:(.*)/;
    if (defined $acc_id and defined $file)
    {
        TIGR::XContig::FastaData::loadFastaFile('ACCESSION', $acc_id, $file);
    }
}

#
#  The input will be first argument on the command line.
#
my $input = shift;

#
#  The rest of the command line should be a list of outputs
#
my @outputs = @ARGV;

if (scalar(@outputs) < 1)
{
    logError(1, "No outputs specified.");
}

my $asm = readContigData($input, $input_type);

unless (defined $asm)
{
  if (defined $input)
  {
      logError(1, "Failed to load input file '$input'.");
  }
  else
  {
      logError(1, "No input file given.");
  }
}

foreach my $seqfile (@seqs)
{
    TIGR::XContig::SeqData::loadSequenceFile($seqfile);
}

#
#  Attempt to read associated files
#
if (TIGR::XContig::SeqData::getSeqCount() < 1)
{
    # Guess the sequence file
    my $seqfile = TIGR::XContig::Util::replaceExtension($input, "seq");

    if (-r $seqfile)
    {
        logMessage(1, "Adopted sequence file: $seqfile");
        logConsole("Adopted sequence file: $seqfile\n");
        TIGR::XContig::SeqData::loadSequenceFile($seqfile);
    }
}

#
#  Force circularity on listed tilings, if supplied
#
if (defined $force_circ_spec)
{
    foreach my $tiling (@{$asm->getTilingList($force_circ_spec)})
    {
        $tiling->setConformation("CIRCULAR");
    }
}

#
#  Set the comment, if supplied
#

if (defined $comment)
{
    foreach my $tiling (@{$asm->{'tiling'}})
    {
        $tiling->{'comment'} = $comment;
    }
}

#
#  Set the com_name, if supplied
#

if (defined $com_name)
{
    foreach my $tiling (@{$asm->{'tiling'}})
    {
        $tiling->{'com_name'} = $com_name;
    }
}

#
#  Write out the requested files
#

OUTPUT: foreach my $output (@outputs)
{
    chomp($output);
    if ($output =~ /^(([A-Za-z]+):\/\/)?(.+?)\s*$/)
    {
        my $format = $2;
        my $file = $3;
        undef $format if (defined $format and length($format) < 1);

        writeContigData($asm, $file, $format);
    }
    else
    {
        logError(0, "Could not parse output specification: $output");
        logConsole("Unidentifiable output: $output\n");
    }
}

#
#  Clean up and finish up
#
END
{

}

