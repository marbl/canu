
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
    fetchTigStore        fetchTigStoreShellCode
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

my $showWork = 1;


#  Convert a/path/to/file to ../../../..
sub pathToDots ($) {
    my $c =  1 + ($_[0] =~ tr!/!/!);        #  Count the number of components in the path.

    return(join("/", map("..", (1..$c))));  #  Return a string with dots for each component.
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
    $ENV{"CANU_OBJECT_STORE_CLIENT_UA"} = getGlobal("objectStoreClientUA");
    $ENV{"CANU_OBJECT_STORE_CLIENT_DA"} = getGlobal("objectStoreClientDA");
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
    my $dx       = getGlobal("objectStoreClient");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");

    return(1)   if ((-e $file) && (!defined($nonLocal)));   #  If file exists, it exists.

    my $exists = "";

    if (isOS() eq "DNANEXUS") {
        $exists = `$dx describe --name \"$pr:$ns/$file\"`;
    }

    $exists =~ s/^\s+//;
    $exists =~ s/\s+$//;

    return(($exists ne "") ? 1 : undef);
}



#  Test if the file exists in object storage.  This is just fileExists(),
#  removing the canu namespace and path.
sub objectStoreFileExists ($) {
    my $path     = shift @_;
    my $dx       = getGlobal("objectStoreClient");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");

    my $exists   = "";

    if (isOS() eq "DNANEXUS") {
        $exists = `$dx describe --name \"$path\"`;
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
    my $dx     = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    if (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $file ] ; then\n";
        $code .= "${indent}  $var=`$dx describe --name \"$pr:$ns/$path/$file\"`\n";
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
    my $dx      = getGlobal("objectStoreClient");
    my $ns      = getGlobal("objectStoreNameSpace");
    my $pr      = getGlobal("objectStoreProject");

    if (isOS() eq "DNANEXUS") {
        print STDERR "renameStashedFile()-- \"$ns/$oldname\" -> \"$ns/$newname\"\n"   if ($showWork);

        if (runCommandSilently(".", "$dx mv \"$pr:$ns/$oldname\" \"$pr:$ns/$newname\"", 1)) {
            caExit("failed to rename object store file", undef);
        }
    }
}



sub removeStashedFile ($) {
    my $name    = shift @_;
    my $dx      = getGlobal("objectStoreClient");
    my $ns      = getGlobal("objectStoreNameSpace");
    my $pr      = getGlobal("objectStoreProject");

    return   if (! fileExists("$name", 1));

    if (isOS() eq "DNANEXUS") {
        print STDERR "removeStashedFile()-- $ns/$name'\n"   if ($showWork);

        if (runCommandSilently(".", "$dx rm --recursive \"$pr:$ns/$name\"", 1)) {
            caExit("failed to remove object store file", undef);
        }
    }
}



#  Runs from the assembly root directory.
sub stashFile ($) {
    my $pathname = shift @_;
    my $path     = dirname($pathname);
    my $name     = basename($pathname);
    my $dx       = getGlobal("objectStoreClient");
    my $ua       = getGlobal("objectStoreClientUA");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");

    my $retries  = 5;    #  Try a few times to upload the file.   (also set in stashFileShellCode())
    my $delay    = 10;   #  Wait a little bit before retrying.

    return   if (! -e $pathname);

    if (isOS() eq "DNANEXUS") {
        print STDERR "stashFile()-- '$pathname' to project '$pr' namespace '$ns' path '$path' name '$name'.\n"   if ($showWork);

        if (runCommandSilently(".", "$dx rm --recursive \"$pr:$ns/$path/$name\"", 1)) {
            caExit("failed to remove object store file", undef);
        }

        #  Try a couple of times to upload the file.  If the UA fails, delay a bit and retry.
        while (($retries > 0) &&
               (runCommandSilently(".", "$ua --do-not-compress --wait-on-close --project \"$pr\" --folder \"$ns/$path/\" --name \"$name\" \"$pathname\"", 0))) { 
            $retries--;
            print STDERR "stashFile()-- Failed to stash file '$pathname', wait $delay seconds and try again ($retries times left).\n";
            sleep($delay);
        }

        if ($retries == 0) {
            caExit("failed to upload file '$pathname' to object store '$pr:$ns/$path/$name'", undef);
        }
    }
}


#  Runs from the assembly root directory.
sub fetchFile ($) {
    my $pathname = shift @_;
    my $path     = dirname($pathname);
    my $name     = basename($pathname);
    my $dx       = getGlobal("objectStoreClient");
    my $da       = getGlobal("objectStoreClientDA");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");

    return   if (-e $pathname);                 #  If it exists, we don't need to fetch it.
    return   if (! fileExists($pathname, 1));   #  If it doesn't exist in the store, we don't fetch it either.  Because it doesn't exist.

    if (isOS() eq "DNANEXUS") {
        print STDERR "fetchFile()-- from project '$pr' path '$path' name '$name' -> '$pathname'\n"   if ($showWork);
        make_path($path);

        if (runCommandSilently(".", "$da download --output \"$pathname\" \"$pr:$ns/$path/$name\"", 1)) {
            caExit("failed to download file from object store", undef);
        }
    }
}


#  Runs from the directory where $file exists.  $file shouldn't have any directory components,
#  but if it does, they're moved to $path.
#
#  $path is the path to the file in object storage.
#
sub stashFileShellCode ($$$) {
    my $path   = shift @_;     #  Path, relative to assembly root, of where we are called from.
    my $file   = shift @_;     #  Name of the file to stash, in current directory.
    my $indent = shift @_;
    my $dx     = getGlobal("objectStoreClient");
    my $ua     = getGlobal("objectStoreClientUA");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    my $retries  = 5;    #  Try a few times to upload the file.   (also set in stashFile())
    my $delay    = 10;   #  Wait a little bit before retrying.

    if (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ -e \"./$file\" ] ; then\n";
        $code .= "${indent}  fold=`dirname \"$path/$file\"`\n";
        $code .= "${indent}  name=`basename \"$file\"`\n";
        $code .= "${indent}\n";
        $code .= "${indent}  $dx rm --recursive \"$pr:$ns/$path/$file\"\n";
        $code .= "${indent}\n";
        $code .= "${indent}  retries=$retries\n";
        $code .= "${indent}  while [ \$retries -gt 0 ] && \\\n";
        $code .= "${indent}        ! $ua --do-not-compress --wait-on-close --project \"$pr\" --folder \"$ns/\$fold/\" --name \"\$name\" \"./$file\" ; do\n";
        $code .= "${indent}    retries=`expr \$retries - 1`\n";
        $code .= "${indent}    echo \"Failed to stash file '$file', wait $delay seconds and try again ($retries times left).\"\n";
        $code .= "${indent}    sleep $delay\n";
        $code .= "${indent}  done\n";
        $code .= "${indent}  if [ \$retries -eq 0 ] ; then\n";
        $code .= "${indent}    echo Failed to stash file '$file', removing incomplete copy.\n";
        $code .= "${indent}    $dx rm --recursive \"$pr:$ns/$path/$file\"\n";
        $code .= "${indent}    exit 1\n";
        $code .= "${indent}  fi\n";
        $code .= "${indent}\n";
        $code .= "${indent}else\n";
        $code .= "${indent}  # Could not find file that we are meant to stash.  We should exit with an error.\n";
        $code .= "${indent}  echo \"Can't upload: failed to find '$file' in '$path'.\"\n";
        $code .= "${indent}  exit 1\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}


#  Runs from the directory where 'file' exists; $file CANNOT have path components.
sub fetchFileShellCode ($$$) {
    my $path   = shift @_;     #  Path, relative to assembly root, of where we are called from.
    my $name   = shift @_;     #  Name of the file to stash, in current directory.
    my $indent = shift @_;
    my $dx     = getGlobal("objectStoreClient");
    my $da     = getGlobal("objectStoreClientDA");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $code   = "";

    if (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e \"./$name\" ] ; then\n";
        $code .= "${indent}  $da download --output \"./$name\" \"$pr:$ns/$path/$name\"\n";
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

sub stashSeqStore ($) {
    my $asm    = shift @_;
    my $files;

    return   if (! -e "$asm.seqStore/info");

    if (defined(isOS())) {
        #  Tar up the store metadata, then upload.

        $files  = "./$asm.seqStore.err ";
        $files .= "./$asm.seqStore.ssi ";
        $files .= "./$asm.seqStore/errorLog ";
        $files .= "./$asm.seqStore/info ";
        $files .= "./$asm.seqStore/info.txt ";
        $files .= "./$asm.seqStore/libraries ";
        $files .= "./$asm.seqStore/libraries.txt ";
        $files .= "./$asm.seqStore/load.dat ";
        $files .= "./$asm.seqStore/readNames.txt ";
        $files .= "./$asm.seqStore/readlengths* ";
        $files .= "./$asm.seqStore/reads ";
        $files .= "./$asm.seqStore/version*"      if (-e "./$asm.seqStore/version.001");

        if (runCommandSilently(".", "tar -cf - $files | gzip -1c > ./$asm.seqStore.tar.gz", 1)) {
            caExit("failed to tar and compress sequence store files", undef);
        }

        stashFile("$asm.seqStore.tar.gz");

        unlink "./$asm.seqStore.tar.gz";

        #  Now stash the individual data files.  If they exist in the store, don't restash.

        for (my $bIdx="0000"; (-e "./$asm.seqStore/blobs.$bIdx"); $bIdx++) {
            if (! fileExists("$asm.seqStore/blobs.$bIdx", 1)) {
                stashFile("$asm.seqStore/blobs.$bIdx");
            }
        }
    }
}


sub fetchSeqStore ($) {
    my $asm    = shift @_;

    return   if (-e "./$asm.seqStore/info");
    return   if (! fileExists("$asm.seqStore.tar.gz", 1));

    if (defined(isOS())) {
        print STDERR "\n"                                                  if ($showWork);
        print STDERR "fetchStore()-- Retrieving store '$asm.seqStore'\n"   if ($showWork);
        print STDERR "\n"                                                  if ($showWork);

        fetchFile("$asm.seqStore.tar.gz");

        if (runCommandSilently(".", "gzip -dc ./$asm.seqStore.tar.gz | tar -xf -", 1)) {
            caExit("failed to uncompress and untar sequence store files", undef);
        }

        unlink "./$asm.seqStore.tar.gz";
    }
}


sub fetchSeqStoreShellCode ($$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $indent = shift @_;           #
    my $base   = dirname($path);     #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'
    my $code   = "";

    if (defined(isOS())) {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $root/$asm.seqStore/info ] ; then\n";
        $code .= "${indent}  echo \"\"\n";
        $code .= "${indent}  echo In '`pwd`', fetching '$asm.seqStore' into '$root'.\n";
        $code .= "${indent}  echo \"\"\n";
        $code .= "${indent}  cd $root\n";
        $code .= fetchFileShellCode(".", "$asm.seqStore.tar.gz", "${indent}  ");
        $code .= "${indent}  tar -xf ./$asm.seqStore.tar.gz\n";
        $code .= "${indent}  rm -f ./$asm.seqStore.tar.gz\n";
        $code .= "${indent}  cd -\n";
        $code .= "${indent}fi\n";
    }
    #  The data files are fetched on demand.

    return($code);
}




sub stashSeqStorePartitions ($$$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $base   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $tag    = shift @_;           #  Which tigs are the partitions for
    my $nJobs  = shift @_;           #  Number of partitions

    my $storePath = "$base/$asm.${tag}Store";
    my $storeName = "partitionedReads.seqStore";

    return   if (! -e "$storePath/$storeName/info");

    if (defined(isOS())) {
        #  Tar up store metadata and data for each partition, and upload.

        my $jName = "0001";

        for (my $job=1; $job <= $nJobs; $job++) {
            my $files;

            $files .= "./$storeName/info ";
            $files .= "./$storeName/info.txt ";
            $files .= "./$storeName/libraries ";
            $files .= "./$storeName/partitions/map ";
            $files .= "./$storeName/partitions/blobs.$jName ";
            $files .= "./$storeName/partitions/reads.$jName";

            if (runCommandSilently($storePath, "tar -cf - $files | gzip -1c > ./$storeName.$jName.tar.gz", 1)) {
                caExit("failed to tar and compress partitioned seqStore files", undef);
            }

            stashFile("$storePath/$storeName.$jName.tar.gz");

            unlink "./$storeName.$jName.tar.gz";

            $jName++;
        }
    }
}


sub fetchSeqStorePartitionShellCode ($$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $indent = shift @_;           #
    my $base   = dirname($path);     #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'
    my $code   = "";

    my $storePath = "$base/$asm.\${tag}Store";
    my $storeName = "partitionedReads.seqStore";

    if (defined(isOS())) {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $root/$storePath/$storeName/partitions/blobs.\$jobid ] ; then\n";
        $code .= "${indent}  cd $root/$storePath\n";
        $code .= fetchFileShellCode($storePath, "$storeName.\$jobid.tar.gz", "${indent}  ");
        $code .= "${indent}  gzip -dc ./$storeName.\$jobid.tar.gz | tar -xf -\n";
        $code .= "${indent}  rm -f ./$storeName.\$jobid.tar.gz\n";
        $code .= "${indent}  cd -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}






sub stashOvlStore ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;

    return   if (! -e "./$base/$asm.ovlStore/index");

    if (defined(isOS())) {
        print STDERR "stashOvlStore()-- Saving overlap store '$base/$asm.ovlStore'\n"   if ($showWork);

        #  Stash the store metadata.

        my $files;
        $files .= " ./$asm.ovlStore/info";
        $files .= " ./$asm.ovlStore/index";
        $files .= " ./$asm.ovlStore/statistics";
        $files .= " ./$asm.ovlStore.config";
        $files .= " ./$asm.ovlStore.config.txt";

        if (runCommandSilently($base, "tar -cf - $files | gzip -1c > ./$asm.ovlStore.tar.gz", 1)) {
            caExit("failed to tar and compress overlap store files", undef);
        }

        stashFile("$base/$asm.ovlStore.tar.gz");
        unlink "$base/$asm.ovlStore.tar.gz";

        #  Stash data files.

        for (my $bIdx="0001"; (-e "./$base/$asm.ovlStore/$bIdx<001>");   $bIdx++) {
        for (my $sIdx="001";  (-e "./$base/$asm.ovlStore/$bIdx<$sIdx>"); $sIdx++) {
            if (! fileExists("$base/$asm.ovlStore/$bIdx<$sIdx>", 1)) {
                stashFile("$base/$asm.ovlStore/$bIdx<$sIdx>");
            }
        }
        }
    }
}


sub stashOvlStoreShellCode ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;
    my $code   = "";

    if (defined(isOS())) {
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
        $code .= "gzip -1c > $base/$asm.ovlStore.tar.gz\n";
        $code .= "\n";
        $code .= stashFileShellCode("$base", "$asm.ovlStore.tar.gz", "");
        $code .= "#\n";
        $code .= "#  Upload data files.\n";
        $code .= "#\n";
        $code .= "\n";
        $code .= "for ff in `ls $asm.ovlStore/????\\<???\\>` ; do\n";
        $code .= stashFileShellCode("$base", "\$ff", "  ");
        $code .= "done\n";
        $code .= "\n";
    }

    return($code);
}



sub fetchOvlStore ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;

    return   if (-e "./$base/$asm.ovlStore/index");
    return   if (! fileExists("$base/$asm.ovlStore.tar.gz"));

    if (defined(isOS())) {
        print STDERR "fetchStore()-- Retrieving store '$base/$asm.ovlStore'\n"   if ($showWork);

        fetchFile("$base/$asm.ovlStore.tar.gz");

        if (runCommandSilently($base, "gzip -dc ./$asm.ovlStore.tar.gz | tar -xf -", 1)) {
            caExit("failed to uncompress and untar overlap store files", undef);
        }

        unlink "$base/$asm.ovlStore.tar.gz";
    }
}


sub fetchOvlStoreShellCode ($$$) {
    my $asm    = shift @_;           #  The name of the assembly.
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $indent = shift @_;           #
    my $base   = dirname($path);     #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'
    my $code   = "";

    if (defined(isOS())) {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $root/$base/$asm.ovlStore/index ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching '$base' / '$asm.ovlStore.tar.gz', unzipping in '$root/$base'\n";
        $code .= fetchFileShellCode($base, "$asm.ovlStore.tar.gz", "${indent}  ");
        $code .= "${indent}  gzip -dc $asm.ovlStore.tar.gz | tar -C $root/$base -xf -\n";
        $code .= "${indent}  rm -f $asm.ovlStore.tar.gz\n";
        $code .= fileExistsShellCode("exists", "unitigging", "$asm.ovlStore/evalues", "${indent}  ");
        $code .= "${indent}  if [ \$exists = false ] ; then\n";
        $code .= "${indent}    cd $root/$base/$asm.ovlStore\n";
        $code .= fetchFileShellCode("$base/$asm.ovlStore", "evalues", "${indent}    ");
        $code .= "${indent}    cd -\n";
        $code .= "${indent}  fi\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



#  So we can create a directory and fetch the file into that directory,
#  functions exist to fetch a tig store.
#
#  There's no corresponding stashTigStore because it's only two files.
#
sub fetchTigStore ($$$$) {
    my $base   = shift @_;
    my $asm    = shift @_;           #  The name of the assembly.
    my $type   = shift @_;           #  Type of store: 'corStore' or 'ctgStore' or 'utgStore'.
    my $vers   = shift @_;           #  Version to fetch, '001' or '002'.

    return   if (-e "./$base/$asm.$type/seqDB.v$vers.tig");
    return   if (! fileExists("$base/$asm.$type/seqDB.v$vers.tig"));

    make_path("$base/$asm.$type");

    fetchFile("$base/$asm.$type/seqDB.v$vers.dat");
    fetchFile("$base/$asm.$type/seqDB.v$vers.tig");
}


sub fetchTigStoreShellCode ($$$$$) {
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging', etc.
    my $asm    = shift @_;           #  The name of the assembly.
    my $type   = shift @_;           #  Type of store: 'corStore' or 'ctgStore' or 'utgStore'.
    my $vers   = shift @_;           #  Version to fetch, '001' or '002'.
    my $indent = shift @_;           #
    my $base   = dirname($path);     #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'
    my $code   = "";

    if (defined(isOS())) {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $root/$base/$asm.$type/seqDB.v$vers.dat ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching '$base' / '$asm.$type/seqDBv$vers'\n";
        $code .= "${indent}  mkdir -p $root/$base/$asm.$type\n";
        $code .= "${indent}  cd $root/$base/$asm.$type\n";
        $code .= fetchFileShellCode("$base/$asm.$type", "seqDB.v$vers.dat", "${indent}  ");
        $code .= fetchFileShellCode("$base/$asm.$type", "seqDB.v$vers.tig", "${indent}  ");
        $code .= "${indent}  cd -\n";
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
    my $code   = "";

    if (defined(isOS())) {
        $code .= "\n";
        $code .= "${indent}if [ -e ./$name/merylIndex ] ; then\n";
        $code .= "${indent}  echo In `pwd`, Storing '$path' '$name'\n";
        $code .= "\n";
        $code .= "${indent}  tar -cf - ./$name/*Index | gzip -1c > ./$name.tar.gz\n";
        $code .= stashFileShellCode($path, "$name.tar.gz", "${indent}  ");
        $code .= "${indent}  rm -f ./$name.tar.gz\n";
        $code .= "\n";
        $code .= "${indent}  cd $name\n";
        $code .= "\n";
        $code .= "${indent}  for ff in 0x000000 0x000001 0x000010 0x000011 0x000100 0x000101 0x000110 0x000111 \\\n";
        $code .= "${indent}            0x001000 0x001001 0x001010 0x001011 0x001100 0x001101 0x001110 0x001111 \\\n";
        $code .= "${indent}            0x010000 0x010001 0x010010 0x010011 0x010100 0x010101 0x010110 0x010111 \\\n";
        $code .= "${indent}            0x011000 0x011001 0x011010 0x011011 0x011100 0x011101 0x011110 0x011111 \\\n";
        $code .= "${indent}            0x100000 0x100001 0x100010 0x100011 0x100100 0x100101 0x100110 0x100111 \\\n";
        $code .= "${indent}            0x101000 0x101001 0x101010 0x101011 0x101100 0x101101 0x101110 0x101111 \\\n";
        $code .= "${indent}            0x110000 0x110001 0x110010 0x110011 0x110100 0x110101 0x110110 0x110111 \\\n";
        $code .= "${indent}            0x111000 0x111001 0x111010 0x111011 0x111100 0x111101 0x111110 0x111111 ; do\n";
        $code .= stashFileShellCode("$path/$name", "\$ff.merylData", "${indent}    ");
        $code .= "${indent}  done\n";
        $code .= "${indent}  cd -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



sub fetchMerylShellCode ($$$) {
    my $path   = shift @_;           #  The subdir we're running in; 'unitigging/4-unitigger', etc.
    my $name   = shift @_;           #  The name of the meryl directory
    my $indent = shift @_;           #
    my $code   = "";

    if (defined(isOS())) {
        $code .= "\n";
        $code .= "${indent}if [ ! -e ./$name/merylIndex ] ; then\n";
        $code .= "${indent}  echo In `pwd`, fetching '$path' '$name.tar.gz'\n";
        $code .= "\n";
        $code .= fetchFileShellCode($path, "$name.tar.gz", "${indent}  ");
        $code .= "${indent}  gzip -dc ./$name.tar.gz | tar -xf -\n";
        $code .= "${indent}  rm -f ./$name.tar.gz\n";
        $code .= "\n";
        $code .= "${indent}  cd $name\n";
        $code .= "\n";
        $code .= "${indent}  for ff in 0x000000 0x000001 0x000010 0x000011 0x000100 0x000101 0x000110 0x000111 \\\n";
        $code .= "${indent}            0x001000 0x001001 0x001010 0x001011 0x001100 0x001101 0x001110 0x001111 \\\n";
        $code .= "${indent}            0x010000 0x010001 0x010010 0x010011 0x010100 0x010101 0x010110 0x010111 \\\n";
        $code .= "${indent}            0x011000 0x011001 0x011010 0x011011 0x011100 0x011101 0x011110 0x011111 \\\n";
        $code .= "${indent}            0x100000 0x100001 0x100010 0x100011 0x100100 0x100101 0x100110 0x100111 \\\n";
        $code .= "${indent}            0x101000 0x101001 0x101010 0x101011 0x101100 0x101101 0x101110 0x101111 \\\n";
        $code .= "${indent}            0x110000 0x110001 0x110010 0x110011 0x110100 0x110101 0x110110 0x110111 \\\n";
        $code .= "${indent}            0x111000 0x111001 0x111010 0x111011 0x111100 0x111101 0x111110 0x111111 ; do\n";
        $code .= fetchFileShellCode("$path/$name", "\$ff.merylData", "${indent}    ");
        $code .= "${indent}  done\n";
        $code .= "${indent}  cd -\n";
        $code .= "${indent}fi\n";
    }

    return($code);
}



1;
