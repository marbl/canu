
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
 #    Sergey Koren beginning on 2018-MAY-09
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Grid_Cloud;

#  This file contains most of the magic needed to access an object store.  Two flavors of each
#  function are needed: one that runs in the canu.pl process (rooted in the base assembly directory,
#  where the 'correction', 'trimming' and 'unitigging' directories exist) and one that is
#  used in shell scripts (rooted where the shell script is run from).

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(configureCloud
             fileExists           fileExistsShellCode
             objectStoreFileExists
             renameStashedFile
             removeStashedFile
             fetchFile            fetchFileShellCode
             fetchObjectStoreFile fetchObjectStoreFileShellCode
             stashFile            stashFileShellCode
             fetchSeqStore        fetchSeqStoreShellCode   fetchSeqStorePartitionShellCode
             fetchOvlStore        fetchOvlStoreShellCode
             stashSeqStore
             stashSeqStorePartitions
             stashOvlStore        stashOvlStoreShellCode
             stashMeryl           stashMerylShellCode
             fetchMeryl           fetchMerylShellCode);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path qw(make_path);
use File::Basename;

use canu::Defaults;
use canu::Execution;

#use canu::Grid "formatAllowedResources";

my $showWork = 0;


#  Convert a/path/to/file to ../../../..
sub pathToDots ($) {
    return(join("/", map("..", (1..scalar(split '/', $_[0])))));
}

#  True if we're using an object store.
sub isOS () {
    my $os = getGlobal("objectStore");

    #  TEST mode is just like DNANEXUS, except where jobs run; see Execution.pm.
    $os = "DNANEXUS"   if ($os eq "TEST");

    return($os);
}



sub configureCloud ($) {
    my $asm = shift @_;

    setGlobalIfUndef("objectStoreNameSpace", $asm);

    #  The seqStore and ovlStore use these to pull data files directly
    #  from cloud storage.  Right now, it's hard coded to use a 'dx' like
    #  command.  Anything else will need a third variable to pass in the
    #  type of object store in use.

    $ENV{"CANU_OBJECT_STORE_CLIENT"}    = getGlobal("objectStoreClient");
    $ENV{"CANU_OBJECT_STORE_NAMESPACE"} = getGlobal("objectStoreNameSpace");
    $ENV{"CANU_OBJECT_STORE_PROJECT"}   = getGlobal("objectStoreProject");
}



#
#  fileExists() returns true if the file exists on disk or in the object store.
#  It does not fetch the file.  It returns undef if the file doesn't exist.
#
#  If a second parameter is supplied, this only tests if the file exists
#  in the object store.  It is not intended to be used outside this module.
#

sub fileExists ($@) {
    my $file     = shift @_;
    my $nonLocal = shift @_;

    my $exists = "";

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    return(1)   if ((-e $file) && (!defined($nonLocal)));   #  If file exists, it exists.

    if    (isOS() eq "DNANEXUS") {
        $exists = `$client describe --name \"$pr:$ns/$file\"`;
    }

    $exists =~ s/^\s+//;
    $exists =~ s/\s+$//;

    return(($exists ne "") ? 1 : undef);
}



#  Test if the file exists in object storage.  This is just fileExists(),
#  removing the canu namespace and path.
sub objectStoreFileExists ($) {
    my $path     = shift @_;

    my $exists   = "";

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    if    (isOS() eq "DNANEXUS") {
        print STDERR "$client describe --name \"$path\"\n";
        $exists = `$client describe --name \"$path\"`;
    }

    $exists =~ s/^\s+//;
    $exists =~ s/\s+$//;

    return(($exists ne "") ? 1 : undef);
}



sub fileExistsShellCode ($$$@) {
    my $var    = shift @_;
    my $path   = shift @_;
    my $file   = shift @_;
    my $indent = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    #  NOTE that this does NOT emit the closing fi.

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $file ] ; then\n";
        $code .= "${indent}  $var=`$client describe --name \"$pr:$ns/$path/$file\"`\n";
        $code .= "${indent}fi\n";
        $code .= "${indent}if [ -e $file -o x\$$var != x ] ; then\n";
        $code .= "${indent}  $var=true\n";
        $code .= "${indent}else\n";
        $code .= "${indent}  $var=false\n";
        $code .= "${indent}fi\n";
    }
    else {
        $code .= "\n";
        $code .= "${indent}if [ -e $file ]; then\n";
        $code .= "${indent}  $var=true\n";
        $code .= "${indent}else\n";
        $code .= "${indent}  $var=false\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



sub renameStashedFile ($$) {
    my $oldname = shift @_;
    my $newname = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    print STDERR "renameStashedFile()-- \"$ns/$oldname\" -> \"$ns/$newname\"\n"   if ($showWork);

    if    (isOS() eq "DNANEXUS") {
        runCommandSilently(".", "$client mv \"$pr:$ns/$oldname\" \"$pr:$ns/$newname\"", 1);
    }
}



sub removeStashedFile ($) {
    my $name = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    return   if (! fileExists("$name", 1));

    print STDERR "removeStashedFile()-- $ns/$name'\n"   if ($showWork);

    if    (isOS() eq "DNANEXUS") {
        runCommandSilently(".", "$client rm --recursive \"$pr:$ns/$name\"", 1);
    }
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
    my $pr     = getGlobal("objectStoreProject");

    return   if (-e $file);                 #  If it exists, we don't need to fetch it.
    return   if (! fileExists($file, 1));   #  If it doesn't exist in the store, we don't fetch it either.  Because it doesn't exist.

    print STDERR "fetchFile()-- '$file' from '$ns/$file'\n"   if ($showWork);

    make_path(dirname($file));

    if    (isOS() eq "DNANEXUS") {
        runCommandSilently(".", "dirname \"$file\" | xargs -n 1 mkdir -p", 1);
        runCommandSilently(".", "$client download --output \"$file\" \"$pr:$ns/$file\"", 1);
    }
}



sub fetchFileShellCode ($$$) {
    my $path   = shift @_;
    my $dots   = pathToDots($path);
    my $file   = shift @_;
    my $indent = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    my $code   = "";

    #  We definitely need to be able to fetch files from places that are
    #  parallel to us, e.g., from 0-mercounts when we're in 1-overlapper.
    #
    #  To get a file, we first go up to the assembly root, then check if the
    #  file exists, and fetch it if not.
    #
    #  The call needs to be something like:
    #    stashFileShellCode("correction/0-mercounts", "whatever", "");

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $dots/$path/$file ] ; then\n";
        $code .= "${indent}  mkdir -p $dots/$path\n";
        $code .= "${indent}  cd       $dots/$path\n";
        $code .= "${indent}  dirname \"$file\" | xargs -n 1 mkdir -p \n";
        $code .= "${indent}  $client download --output \"$file\" \"$pr:$ns/$path/$file\"\n";
        $code .= "${indent}  cd -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



#  Fetch an object store file using a full path name.
sub fetchObjectStoreFile ($$) {
    my $path   = shift @_;    #  Path of file to fetch.
    my $file   = shift @_;    #  Name of file to write in current directory.

    my $client = getGlobal("objectStoreClient");

    my $exists = objectStoreFileExists($path);

    return   if (-e $file);       #  If it exists, we don't need to fetch it.
    return   if (! $exists);      #  If it doesn't exist in the store, we don't fetch it either.  Because it doesn't exist.

    print STDERR "fetchObjectStoreFile()-- '$file' from '$path'\n"   if ($showWork);

    if    (isOS() eq "DNANEXUS") {
        runCommandSilently(".", "$client download --output \"$file\" \"$path\"", 1);
    }
}



sub fetchObjectStoreFileShellCode ($$$) {
    my $path   = shift @_;    #  Path of file to fetch.
    my $file   = shift @_;    #  Name of file to write in current directory.
    my $indent = shift @_;

    my $client = getGlobal("objectStoreClient");

    my $code   = "";

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $$file ] ; then\n";
        $code .= "${indent}  $client download --output \"$file\" \"$path\"\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



sub stashFile ($) {
    my $file   = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    return   if (! -e $file);

    removeStashedFile($file);

    print STDERR "stashFile()-- '$file' to '$ns/$file'\n"   if ($showWork);

    if    (isOS() eq "DNANEXUS") {
        runCommandSilently(".", "$client upload --wait --parents --path \"$pr:$ns/$file\" \"$file\"", 1);
    }
}



sub stashFileShellCode ($$$) {
    my $path   = shift @_;
    my $dots   = pathToDots($path);
    my $file   = shift @_;
    my $indent = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    #  Just like for fetching, we allow stashing files from parallel
    #  directories (even though that should never happen).

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ -e \"$dots/$path/$file\" ] ; then\n";
        $code .= "${indent}  cd $dots/$path\n";
        $code .= fileExistsShellCode("exists", "$path", "$file", "$indent  ");
        $code .= "${indent}  if [ \$exists = true ] ; then\n";
        $code .= "${indent}    $client rm --recursive \"$pr:$ns/$path/$file\"\n";
        $code .= "${indent}  fi\n";
        $code .= "${indent}  $client upload --wait --parents --path \"$pr:$ns/$path/$file\" \"$file\"\n";
        $code .= "${indent}  cd -\n";
        $code .= "${indent}else\n";
        $code .= "${indent}  # Could not find file that we are meant to stash.  We should exit with an error.\n";
        $code .= "${indent}  exit 1\n";
        $code .= "${indent}fi\n";
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
    my $pr     = getGlobal("objectStoreProject");

    return   if (-e "./$asm.seqStore/info");
    return   if (! fileExists("$asm.seqStore.tar", 1));

    print STDERR "fetchStore()-- Retrieving store '$asm.seqStore'\n"   if ($showWork);

    if    (isOS() eq "DNANEXUS") {
        runCommandSilently(".", "$client download --output - $pr:$ns/$asm.seqStore.tar | tar -xf -", 1);
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
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $root/$asm.seqStore/info ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching $ns/$asm.seqStore.tar, unzipping in '$root'\n";
        $code .= "${indent}  $client download --output - $pr:$ns/$asm.seqStore.tar | tar -C $root -xf -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



sub stashSeqStore ($) {
    my $asm    = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    return   if (! -e "$asm.seqStore/info");

    if    (isOS() eq "DNANEXUS") {
        my $cmd;

        #  Stash the store metadata.

        $cmd  = "tar -cf - ";
        $cmd .= "./$asm.seqStore.err";
        $cmd .= " ./$asm.seqStore.ssi";
        $cmd .= " ./$asm.seqStore/errorLog";
        $cmd .= " ./$asm.seqStore/info";
        $cmd .= " ./$asm.seqStore/info.txt";
        $cmd .= " ./$asm.seqStore/libraries";
        $cmd .= " ./$asm.seqStore/libraries.txt";
        $cmd .= " ./$asm.seqStore/load.dat";
        $cmd .= " ./$asm.seqStore/readNames.txt";
        $cmd .= " ./$asm.seqStore/readlengths*";
        $cmd .= " ./$asm.seqStore/reads";
        $cmd .= " ./$asm.seqStore/version*" if (-e "./$asm.seqStore/version.001");
        $cmd .= " | $client upload --wait --parents --path $pr:$ns/$asm.seqStore.tar -";

        removeStashedFile("$asm.seqStore.tar");

        print STDERR "stashSeqStore()-- Saving sequence store '$asm.seqStore'\n"   if ($showWork);

        runCommandSilently(".", $cmd, 1);

        #  Stash the store data files.

        for (my $bIdx="0000"; (-e "./$asm.seqStore/blobs.$bIdx"); $bIdx++) {
            if (! fileExists("$asm.seqStore/blobs.$bIdx", 1)) {
                runCommandSilently(".", "$client upload --wait --parents --path $pr:$ns/$asm.seqStore/blobs.$bIdx $asm.seqStore/blobs.$bIdx", 1);
            }
        }
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
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    my $storePath = "$base/$asm.\${tag}Store";
    my $storeName = "partitionedReads.seqStore";

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $root/$storePath/$storeName/partitions/blobs.\$jobid ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching $ns/$storePath.$storeName.\$jobid.tar, unzipping in '$root/$storePath'\n";
        $code .= "${indent}  $client download --output - $pr:$ns/$storePath.$storeName.\$jobid.tar | tar -C $root/$storePath -xf -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



sub stashSeqStorePartitions ($$$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $base   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $tag    = shift @_;           #  Which tigs are the partitions for
    my $nJobs  = shift @_;           #  Number of partitions

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    my $storePath = "$base/$asm.${tag}Store";
    my $storeName = "partitionedReads.seqStore";

    return   if (! -e "$storePath/$storeName/info");

    if    (isOS() eq "DNANEXUS") {
        my $jName = "0001";

        for (my $job=1; $job <= $nJobs; $job++) {
            my $cmd;

            $cmd  = "tar -cf - ";
            $cmd .= " ./$storeName/info";
            $cmd .= " ./$storeName/info.txt";
            $cmd .= " ./$storeName/libraries";
            $cmd .= " ./$storeName/partitions/map";
            $cmd .= " ./$storeName/partitions/blobs.$jName";
            $cmd .= " ./$storeName/partitions/reads.$jName";
            $cmd .= " | $client upload --wait --parents --path $pr:$ns/$storePath.$storeName.$jName.tar -";

            removeStashedFile("$storePath.$storeName.$jName.tar");

            print STDERR "stashPartitionedSeqStore()-- Saving partitioned sequence store '$storePath/$storeName.$jName'\n"   if ($showWork);

            runCommandSilently($storePath, $cmd, 1);

            $jName++;
        }
    }
}



### OVERLAPS



sub fetchOvlStore ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    return   if (-e "./$base/$asm.ovlStore/index");
    return   if (! fileExists("$base/$asm.ovlStore.tar"));

    print STDERR "fetchStore()-- Retrieving store '$base/$asm.ovlStore'\n"   if ($showWork);

    if    (isOS() eq "DNANEXUS") {
        runCommandSilently($base, "$client download --output - $pr:$ns/$base/$asm.ovlStore.tar | tar -xf -", 1);
    }
}



sub stashOvlStore ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");

    return   if (! -e "./$base/$asm.ovlStore/index");

    if    (isOS() eq "DNANEXUS") {
        my $cmd;

        #  Stash the store metadata.

        $cmd  = "tar -cf - ";
        $cmd .= " ./$asm.ovlStore/info";
        $cmd .= " ./$asm.ovlStore/index";
        $cmd .= " ./$asm.ovlStore/statistics";
        $cmd .= " ./$asm.ovlStore.config";
        $cmd .= " ./$asm.ovlStore.config.txt";
        $cmd .= " | $client upload --wait --parents --path $pr:$ns/$base/$asm.ovlStore.tar -";

        removeStashedFile("$base/$asm.ovlStore.tar");

        print STDERR "stashOvlStore()-- Saving overlap store '$base/$asm.ovlStore'\n"   if ($showWork);

        runCommandSilently($base, $cmd, 1);

        #  Stash the store data files.

        for (my $bIdx="0001"; (-e "./$base/$asm.ovlStore/$bIdx<001>");   $bIdx++) {
        for (my $sIdx="001";  (-e "./$base/$asm.ovlStore/$bIdx<$sIdx>"); $sIdx++) {
            if (! fileExists("$base/$asm.ovlStore/$bIdx<$sIdx>", 1)) {
                runCommandSilently(".", "$client upload --wait --parents --path \"$pr:$ns/$base/$asm.ovlStore/$bIdx<$sIdx>\" \"./$base/$asm.ovlStore/$bIdx<$sIdx>\"", 1);
            }
        }
        }
    }
}



sub stashOvlStoreShellCode ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";


    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "#  If the store doesn't exist, ovStoreBuild failed, and we should just quit.\n";
        $code .= "\n";
        $code .= "if [ ! -e ./$asm.ovlStore ] ; then\n";
        $code .= "  exit\n";
        $code .= "fi\n";
        $code .= "\n";
        $code .= "#\n";
        $code .= "#  Upload the metadata files.  These shouldn't exist, so we don't bother trying to remove before uploading.\n";
        $code .= "#\n";
        $code .= "\n";
        $code .= "tar -cf - \\\n";
        $code .= " ./$asm.ovlStore/info \\\n";
        $code .= " ./$asm.ovlStore/index \\\n";
        $code .= " ./$asm.ovlStore/statistics \\\n";
        $code .= " ./$asm.ovlStore.config \\\n";
        $code .= " ./$asm.ovlStore.config.txt \\\n";
        $code .= "| \\\n";
        $code .= "$client upload --wait --parents --path $pr:$ns/$base/$asm.ovlStore.tar -\n";
        $code .= "\n";
        $code .= "#\n";
        $code .= "#  Upload data files.\n";
        $code .= "#\n";
        $code .= "\n";
        $code .= "for ff in `ls $asm.ovlStore/????\\<???\\>` ; do\n";
        $code .= "  $client upload --wait --parents --path \"$pr:$ns/$base/\$ff\" \"./\$ff\"\n";
        $code .= "done\n";
        $code .= "\n";
    }

    return($code);
}



sub fetchOvlStoreShellCode ($$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $indent = shift @_;           #

    my $base   = dirname($path);     #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $root/$base/$asm.ovlStore/index ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching $ns/$base/$asm.ovlStore.tar, unzipping in '$root/$base'\n";
        $code .= "${indent}  $client download --output - $pr:$ns/$base/$asm.ovlStore.tar | tar -C $root/$base -xf -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



sub stashMeryl ($$$) {
}



sub fetchMeryl ($$$) {
}



sub stashMerylShellCode ($$$) {
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $name   = shift @_;           #  The name of the meryl directory
    my $indent = shift @_;           #

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ -e ./$name/merylIndex ] ; then\n";
        $code .= "${indent}  echo In `pwd`, Storing './$name' to '$pr:$ns/$path/$name.tar'\n";
        $code .= "${indent}  tar -cf - ./$name \\\n";
        $code .= "${indent}  | \\\n";
        $code .= "${indent}  $client upload --wait --parents --path $pr:$ns/$path/$name.tar -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



sub fetchMerylShellCode ($$$) {
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $name   = shift @_;           #  The name of the meryl directory
    my $indent = shift @_;           #

    my $client = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    if    (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e ./$name/merylIndex ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching '$pr:$ns/$path/$name.tar', unzipping to './$name'\n";
        $code .= "${indent}  $client download --output - $pr:$ns/$path/$name.tar \\\n";
        $code .= "${indent}  | \\\n";
        $code .= "${indent}  tar -xf -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



1;
