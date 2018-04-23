
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2017-FEB-15
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Grid_Cloud;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(fileExists           fileExistsShellCode
             fetchFile            fetchFileShellCode
             stashFile            stashFileShellCode
             fetchSeqStore        fetchSeqStoreShellCode   fetchSeqStorePartitionShellCode
             fetchOvlStore        fetchOvlStoreShellCode
             stashSeqStore
             stashSeqStorePartitions
             stashOvlStore);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path qw(make_path);
use File::Basename;

use canu::Defaults;
use canu::Execution;

#use canu::Grid "formatAllowedResources";



#  This file contains most of the magic needed to access an object store.  Two flavors of each
#  function are needed: one that runs in the canu.pl process (rooted in the base assembly directory,
#  where the 'correction', 'trimming' and 'unitigging' directories exist) and one that is
#  used in shell scripts (rooted where the shell script is run from).


#  Convert a/path/to/file to ../../../..
sub pathToDots ($) {
    return(join("/", map("..", (1..scalar(split '/', $_[0])))));
}

#  True if we're using an object store.
sub isOS () {
    return(getGlobal("objectStore"));
}


#
#  fileExists() returns true if the file exists on disk or in the object store.  It does not fetch
#  the file.  It returns undef if the file doesn't exist.  The second argument to
#  fileExistsShellCode() is an optional indent level (a whitespace string).
#
#  The shellCode version should emit the if test for file existence, but nothing else (not even the
#  endif).
#

sub fileExists ($) {
    my $file   = shift @_;
    my $exists = "";

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return(1)   if (-e $file);           #  If file exists, it exists.

    if    (isOS() eq "TEST") {
        $exists = `$client describe --name $ns/$file`;
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $exists = "";
    }

    $exists =~ s/^\s+//;
    $exists =~ s/\s+$//;

    return(($exists ne "") ? 1 : undef);
}



sub fileExistsShellCode ($@) {
    my $file   = shift @_;
    my $indent = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $code   = "";

    if    (isOS() eq "TEST") {
        $code .= "${indent}if [ ! -e $file ] ; then\n";
        $code .= "${indent}  exists=`$client describe --name $ns/$file`\n";
        $code .= "${indent}fi\n";
        $code .= "${indent}if [ -e $file -o x\$exists != x ] ; then\n";
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $code .= "${indent}if [ -e $file ]; then\n";
    }

    return($code);
}



#
#  fetchFile() and stashFile() both expect to be called from the assembly root directory, and have
#  the path to the file passed in, e.g., "correction/0-mercounts/whatever.histogram".
#
#  The shellCode versions expect the same, but need the path from the assembly root to the location
#  the shell script is running split.  A meryl script would give "correction/0-mercounts" for the
#  first arg, and could give "some/directory/file" for the file.
#

sub fetchFile ($) {
    my $file   = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return   if (-e $file);   #  If it exists, we don't need to fetch it.

    if    (isOS() eq "TEST") {
        print STDERR "fetchFile()-- '$file'\n";
        make_path(dirname($file));
        runCommandSilently(".", "$client download --output $file $ns/$file", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        #  Nothing we can be obnoxious about here, I suppose we could log...
    }
}



sub fetchFileShellCode ($$$) {
    my $path   = shift @_;
    my $dots   = pathToDots($path);
    my $file   = shift @_;
    my $indent = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $code   = "";

    #  We definitely need to be able to fetch files from places that are
    #  parallel to us, e.g., from 0-mercounts when we're in 1-overlapper.
    #
    #  To get a file, we first go up to the assembly root, then check if the
    #  file exists, and fetch it if not.
    #
    #  The call needs to be something like:
    #    stashFileShellCode("correction/0-mercounts", "whatever", "");

    if    (isOS() eq "TEST") {
        $code .= "${indent}if [ ! -e $dots/$path/$file ] ; then\n";
        $code .= "${indent}  mkdir -p $dots/$path\n";
        $code .= "${indent}  cd       $dots/$path\n";
        $code .= "${indent}  $client download --output $file $ns/$path/$file\n";
        $code .= "${indent}  cd -\n";
        $code .= "${indent}fi\n";
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $code .= "#  File must exist: $file\n";
    }

    return($code);
}



sub stashFile ($) {
    my $file   = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return   if (! -e $file);

    if    (isOS() eq "TEST") {
        print STDERR "stashFile()-- '$file'\n";
        runCommandSilently(".", "$client upload --path $ns/$file $file", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        #  Nothing we can be obnoxious about here, I suppose we could log...
    }

}



sub stashFileShellCode ($$$) {
    my $path   = shift @_;
    my $dots   = pathToDots($path);
    my $file   = shift @_;
    my $indent = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $code   = "";

    #  Just like for fetching, we allow stashing files from parallel
    #  directories (even though that should never happen).

    if    (isOS() eq "TEST") {
        $code .= "${indent}if [ -e $dots/$path/$file ] ; then\n";
        $code .= "${indent}  cd $dots/$path\n";
        $code .= "${indent}  $client upload --path $ns/$path/$file $file\n";
        $code .= "${indent}  cd -\n";
        $code .= "${indent}fi\n";
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $code .= "#  File is important: $file\n";
    }

    return($code);
}



### SEQUENCE

#
#  Given $base/$asm.seqStore, fetch or stash it.
#
#  The non-shell versions are assumed to be running in the assembly directory, that is, where
#  $base/$asm.seqStore would exist naturally.  This is consistent with canu.pl - it runs in the
#  assembly directory, and then chdir to subdirectories to run binaries.
#
#  The shell versions usually run within a subdirectory (e.g., in correction/0-mercounts).  They
#  need to know this location, so they can go up to the assembly directory to fetch and unpack the
#  store.  After fetching, they chdir back to the subdirectory.
#

sub fetchSeqStore ($) {
    my $asm    = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return   if (-e "./$asm.seqStore/info");
    return   if (! fileExists("$asm.seqStore.tar"));

    if    (isOS() eq "TEST") {
        print STDERR "fetchStore()-- Retrieving store '$asm.seqStore'\n";
        runCommandSilently(".", "$client download --output - $ns/$asm.seqStore.tar | tar -xf -", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
    }
}


sub fetchSeqStoreShellCode ($$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $indent = shift @_;           #

    my $base   = dirname($path);     #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $code   = "";

    if    (isOS() eq "TEST") {
        $code .= "${indent}if [ ! -e $root/$asm.seqStore/info ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching $ns/$asm.seqStore.tar, unzipping in '$root'\n";
        $code .= "${indent}  $client download --output - $ns/$asm.seqStore.tar | tar -C $root -xf -\n";
        $code .= "${indent}fi\n";
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $code .= "#  Store must exist: $root/$asm.seqStore\n";
    }

    return($code);
}



sub stashSeqStore ($) {
    my $asm    = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return   if (! -e "./$asm.seqStore/info");

    if    (isOS() eq "TEST") {
        print STDERR "stashSeqStore()-- Saving sequence store '$asm.seqStore'\n";
        runCommandSilently(".", "tar -cf - ./$asm.seqStore* ./*/$asm.seqStore* | $client upload --path $ns/$asm.seqStore.tar -", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
    }
}



sub fetchSeqStorePartitionShellCode ($$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $indent = shift @_;           #

    my $base   = dirname($path);     #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $code   = "";

    my $storePath = "$base/$asm.\${tag}Store";
    my $storeName = "partitionedReads.seqStore";

    if    (isOS() eq "TEST") {
        $code .= "${indent}if [ ! -e $root/$storePath/$storeName/info ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching $ns/$storePath.$storeName.tar, unzipping in '$root/$storePath'\n";
        $code .= "${indent}  $client download --output - $ns/$storePath.$storeName.tar | tar -C $root/$storePath -xf -\n";
        $code .= "${indent}fi\n";
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $code .= "#  Store must exist: $root/$storePath/$storeName\n";
    }

    return($code);
}



sub stashSeqStorePartitions ($$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $base   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $tag    = shift @_;           #  Which tigs are the partitions for

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    my $storePath = "$base/$asm.\${tag}Store";
    my $storeName = "partitionedReads.seqStore";

    return   if (! -e "$storePath/$storeName/info");

    if    (isOS() eq "TEST") {
        print STDERR "stashPartitionedSeqStore()-- Saving partitioned sequence store '$storePath/$storeName'\n";
        runCommandSilently($storePath, "tar -cf ./$storeName* | $client upload --path $ns/$storePath.$storeName.tar -", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
    }
}



### OVERLAPS



sub fetchOvlStore ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return   if (-e "./$base/$asm.ovlStore/index");
    return   if (! fileExists("$base/$asm.ovlStore.tar"));

    if    (isOS() eq "TEST") {
        print STDERR "fetchStore()-- Retrieving store '$base/$asm.ovlStore'\n";
        runCommandSilently($base, "$client download --output - $ns/$base/$asm.ovlStore.tar | tar -xf -", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
    }
}



sub stashOvlStore ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return   if (! -e "./$base/$asm.ovlStore/index");

    if    (isOS() eq "TEST") {
        print STDERR "stashOvlStore()-- Saving overlap store '$base/$asm.ovlStore'\n";
        runCommandSilently($base, "tar -cf - ./$asm.ovlStore* | $client upload --path $ns/$base/$asm.ovlStore.tar -", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
    }
}



sub fetchOvlStoreShellCode ($$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $indent = shift @_;           #

    my $base   = dirname($path);    #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $code   = "";

    if    (isOS() eq "TEST") {
        $code .= "${indent}if [ ! -e $root/$base/$asm.ovlStore/index ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching $ns/$base/$asm.ovlStore.tar, unzipping in '$root/$base'\n";
        $code .= "${indent}  $client download --output - $ns/$base/$asm.ovlStore.tar | tar -C $root/$base -xf -\n";
        $code .= "${indent}fi\n";
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $code .= "#  Store must exist: $root/$base/$asm.ovlStore\n";
    }

    return($code);
}
