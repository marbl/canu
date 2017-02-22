
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
 #    Brian P. Walenz beginning on 2017-JAN-17
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Grid_Cloud;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(fileExists
             fileExistsShellCode
             fetchFile
             fetchFileShellCode
             stashFile
             stashFileShellCode
             fetchStore
             fetchStoreShellCode
             stashStore
             stashStoreShellCode);

use strict;

use File::Path qw(make_path);
use File::Basename;

use Cwd qw(getcwd);

use canu::Defaults;
use canu::Grid;
use canu::Execution qw(runCommand runCommandSilently);


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
    my $code   = "";
    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

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
    my $code   = "";
    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

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
    my $code   = "";
    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

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



#
#  Given $base/$asm.gkpStore, fetch or stash it.
#
#  The non-shell versions are assumed to be running in the assembly directory, that is, where
#  $base/$asm.gkpStore would exist naturally.  This is consistent with canu.pl - it runs in the
#  assembly directory, and then chdir to subdirectories to run binaries.
#
#  The shell versions usually run within a subdirectory (e.g., in correction/0-mercounts).  They
#  need to know this location, so they can go up to the assembly directory to fetch and unpack the
#  store.  After fetching, they chdir back to the subdirectory.
#

sub fetchStore ($) {
    my $store  = shift @_;                           #  correction/asm.gkpStore
    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return   if (-e "$store/info");                  #  Store exists on disk
    return   if (! fileExists("$store.tar"));        #  Store doesn't exist in object store

    if    (isOS() eq "TEST") {
        runCommandSilently(".", "$client download --output - $ns/$store.tar | tar -xf -", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
    }
}



sub stashStore ($) {
    my $store  = shift @_;                         #  correction/asm.gkpStore
    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    return   if (! -e "$store/info");              #  Store doesn't exist on disk

    if    (isOS() eq "TEST") {
        runCommandSilently(".", "tar -cf - $store | $client upload --path $ns/$store.tar -", 1);
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
    }
}



sub fetchStoreShellCode ($$@) {
    my $store  = shift @_;           #  correction/asm.gkpStore - store we're trying to get
    my $root   = shift @_;           #  correction/1-overlapper - place the script is running in
    my $indent = shift @_;           #
    my $base   = dirname($store);    #  correction
    my $basep  = pathToDots($root);  #  ../..
    my $name   = basename($store);   #             asm.gkpStore
    my $code;
    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    if    (isOS() eq "TEST") {
        $code .= "${indent}if [ ! -e $basep/$store/info ] ; then\n";
        $code .= "${indent}  echo Fetching $ns/$store\n";
        $code .= "${indent}  $client download --output - $ns/$store.tar | tar -C $basep -xf -\n";
        $code .= "${indent}fi\n";
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $code .= "#  Store must exist: $store\n";
    }

    return($code);
}



sub stashStoreShellCode ($$@) {
    my $store  = shift @_;           #  correction/asm.gkpStore - store we're trying to get
    my $root   = shift @_;           #  correction/1-overlapper - place the script is running in
    my $indent = shift @_;           #
    my $base   = dirname($store);    #  correction
    my $basep  = pathToDots($root);  #  ../..
    my $name   = basename($store);   #             asm.gkpStore
    my $code;
    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");

    if    (isOS() eq "TEST") {
        $code .= "${indent}if [ -e $basep/$store/info ] ; then\n";
        $code .= "${indent}  echo Stashing $ns/$store\n";
        $code .= "${indent}  tar -C $basep -cf - $store | $client upload --path $ns/$store.tar -\n";
        $code .= "${indent}fi\n";
    }
    elsif (isOS() eq "DNANEXUS") {
    }
    else {
        $code .= "#  Store is important: $store\n";
    }

    return($code);
}

