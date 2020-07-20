
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
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
             fetchFileFromLink    fetchFileFromLinkShellCode
             stashFile            stashFileShellCode
                                  stashFilesShellCode
             fetchSeqStore        fetchSeqStoreShellCode
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

my $retryCount = 5;    #  Try a few times to upload the file.
my $retryDelay = 10;   #  Wait a little bit before retrying.


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

sub sanitizeName ($) {
    my $name = shift @_;

    $name =~ s!/\./!/!g;
    $name =~ s!//!/!g;

    return($name);
}



sub configureCloud ($$) {
    my $asm = shift @_;
    my $dir = shift @_;

    #  If no object store defined, nothing to do here.

    if (!defined(getGlobal("objectStore"))) {
        return;
    }

    #  The 'namespace' is just our output directory.

    setGlobalIfUndef("objectStoreNameSpace", $dir);

    #  Pull any object store environment out of the environment.  These
    #  explicitly overwrite any current setting as its what the last
    #  execution of canu said they were.
    #
    #  These are set in canu-scripts/canu.##.sh.

    if (defined($ENV{"CANU_OBJECT_STORE_CLIENT"}))      { setGlobal("objectStoreClient",    $ENV{"CANU_OBJECT_STORE_CLIENT"});    }
    if (defined($ENV{"CANU_OBJECT_STORE_CLIENT_UA"}))   { setGlobal("objectStoreClientUA",  $ENV{"CANU_OBJECT_STORE_CLIENT_UA"}); }
    if (defined($ENV{"CANU_OBJECT_STORE_CLIENT_DA"}))   { setGlobal("objectStoreClientDA",  $ENV{"CANU_OBJECT_STORE_CLIENT_DA"}); }
    if (defined($ENV{"CANU_OBJECT_STORE_PROJECT"}))     { setGlobal("objectStoreProject",   $ENV{"CANU_OBJECT_STORE_PROJECT"});   }
    if (defined($ENV{"CANU_OBJECT_STORE_NAMESPACE"}))   { setGlobal("objectStoreNameSpace", $ENV{"CANU_OBJECT_STORE_NAMESPACE"}); }

    #  Fix up parameters, check for missing and/or invalid settings.

    if ((getGlobal("objectStore") ne "TEST") &&
        (getGlobal("objectStore") ne "DNANEXUS")) {
        addCommandLineError("ERROR:  Invalid 'objectStore' specified (" . getGlobal("objectStore") . "); must be unset or 'DNANEXUS'\n");
    }

    if (!defined(getGlobal("objectStoreClient"))) {
        addCommandLineError("ERROR:  objectStoreClient must be specified if objectStore is specified\n");
    }

    #  If no upload or download agents, default to the client.

    if (!defined(getGlobal("objectStoreClientUA"))) {
        setGlobal("objectStoreClientUA", getGlobal("objectStoreClient"));
    }
    if (!defined(getGlobal("objectStoreClientDA"))) {
        setGlobal("objectStoreClientDA", getGlobal("objectStoreClient"));
    }

    #  DNAnexus needs a project too.

    if (( getGlobal("objectStore") eq "DNANEXUS") &&
        (!defined(getGlobal("objectStoreProject")))) {
        addCommandLineError("ERROR:  objectStoreProject must be specified if objectStore=DNANEXUS is specified\n");
    }

    #  The seqStore and ovlStore use these to pull data files directly
    #  from cloud storage.  Right now, it's hard coded to use a 'dx' like
    #  command.  Anything else will need a third variable to pass in the
    #  type of object store in use.

    $ENV{"CANU_OBJECT_STORE_CLIENT"}    = getGlobal("objectStoreClient");
    $ENV{"CANU_OBJECT_STORE_CLIENT_UA"} = getGlobal("objectStoreClientUA");
    $ENV{"CANU_OBJECT_STORE_CLIENT_DA"} = getGlobal("objectStoreClientDA");
    $ENV{"CANU_OBJECT_STORE_PROJECT"}   = getGlobal("objectStoreProject");
    $ENV{"CANU_OBJECT_STORE_NAMESPACE"} = getGlobal("objectStoreNameSpace");

    print STDERR "--\n";
    print STDERR "-- Object storage enabled!  Using:\n";
    print STDERR "--   CANU_OBJECT_STORE_CLIENT     = '", getGlobal("objectStoreClient"),    "'\n";
    print STDERR "--   CANU_OBJECT_STORE_CLIENT_UA  = '", getGlobal("objectStoreClientUA"),  "'\n";
    print STDERR "--   CANU_OBJECT_STORE_CLIENT_DA  = '", getGlobal("objectStoreClientDA"),  "'\n";
    print STDERR "--   CANU_OBJECT_STORE_PROJECT    = '", getGlobal("objectStoreProject"),   "'\n";
    print STDERR "--   CANU_OBJECT_STORE_NAMESPACE  = '", getGlobal("objectStoreNameSpace"), "'\n";
}



#
#  fileExists() returns true if the file exists on disk or in the object store.
#  It does not fetch the file.  It returns undef if the file doesn't exist.
#
#  If a second parameter is supplied, this only tests if the file exists
#  in the object store.  It is not intended to be used outside this module.
#

sub fileExists ($@) {
    my $dx       = getGlobal("objectStoreClient");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");
    my $file     = shift @_;
    my $link     = sanitizeName("$pr:$ns/$file");
    my $nonLocal = shift @_;

    return(1)   if ((-e $file) && (!defined($nonLocal)));   #  If file exists, it exists.

    my $exists = "";

    if (isOS() eq "DNANEXUS") {
        $exists = `$dx describe --name \"$link\" 2> /dev/null`;
    }

    $exists =~ s/^\s+//;
    $exists =~ s/\s+$//;

    return(($exists ne "") ? 1 : undef);
}



#  Test if the file exists in object storage.  This is just fileExists(),
#  removing the canu namespace and path.
sub objectStoreFileExists ($) {
    my $dx       = getGlobal("objectStoreClient");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");
    my $path     = santizeName(shift @_);

    my $exists   = "";

    if (isOS() eq "DNANEXUS") {
        $exists = `$dx describe --name \"$path\" 2> /dev/null`;
    }

    $exists =~ s/^\s+//;
    $exists =~ s/\s+$//;

    return(($exists ne "") ? 1 : undef);
}



sub fileExistsShellCode ($$$@) {
    my $dx     = getGlobal("objectStoreClient");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $var    = shift @_;
    my $path   = shift @_;
    my $file   = shift @_;
    my $indent = shift @_;
    my $code   = "";

    if (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $file ] ; then\n";
        $code .= "${indent}  $var=`$dx describe --name \"$pr:$ns/$path/$file\" 2> /dev/null`\n";
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
    my $dx      = getGlobal("objectStoreClient");
    my $ns      = getGlobal("objectStoreNameSpace");
    my $pr      = getGlobal("objectStoreProject");
    my $oldname = shift @_;
    my $newname = shift @_;
    my $oldlink = sanitizeName("$pr:$ns/$oldname");
    my $newlink = sanitizeName("$pr:$ns/$newname");

    if (isOS() eq "DNANEXUS") {
        print STDERR "-- Rename storage object '$oldlink' -> '$newlink'\n"   if ($showWork);

        if (runCommandSilently(".", "$dx mv \"$oldlink\" \"$newlink\"", 1)) {
            caExit("failed to rename object store file", undef);
        }
    }
}



sub removeStashedFile ($) {
    my $dx      = getGlobal("objectStoreClient");
    my $ns      = getGlobal("objectStoreNameSpace");
    my $pr      = getGlobal("objectStoreProject");
    my $name    = shift @_;
    my $link    = sanitizeName("$pr:$ns/$name");

    return   if (! fileExists("$name", 1));

    if (isOS() eq "DNANEXUS") {
        print STDERR "-- Remove storage object '$link'\n"   if ($showWork);

        if (runCommandSilently(".", "$dx rm --recursive \"$link\"", 1)) {
            caExit("failed to remove object store file", undef);
        }
    }
}



#  Runs from the assembly root directory.
sub stashFile ($) {
    my $dx       = getGlobal("objectStoreClient");
    my $ua       = getGlobal("objectStoreClientUA");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");
    my $pathname = shift @_;
    my $path     = dirname($pathname);
    my $name     = basename($pathname);
    my $link;

    if ($path eq ".") {
        $link = sanitizeName("$pr:$ns/$name");
        $path = sanitizeName("$ns");
    } else {
        $link = sanitizeName("$pr:$ns/$path/$name");
        $path = sanitizeName("$ns/$path");
    }

    return   if (! -e $pathname);

    my $retries  = $retryCount;    #  Try a few times to upload the file.
    my $delay    = $retryDelay;    #  Wait a little bit before retrying.

    if (isOS() eq "DNANEXUS") {
        print STDERR "-- Store '$pathname' to storage object '$link'.\n"   if ($showWork);

        #  If the filename exists, remove it before uploading.
        if (fileExists($pathname, 1)) {
            if (runCommandSilently(".", "$dx rm --recursive \"$link\"", 1)) {
                caExit("failed to remove object store file", undef);
            }
        }

        #  Try a couple of times to upload the file.  If the UA fails, delay a bit and retry.
        while (($retries > 0) &&
               (runCommandSilently(".", "$ua --do-not-compress --wait-on-close --project \"$pr\" --folder \"$path/\" --name \"$name\" \"$pathname\"", 0))) {
            $retries--;
            print STDERR "-- WARNING:\n";
            print STDERR "-- WARNING:  Failed to store '$pathname'.  Wait $delay seconds and try again ($retries times left).\n";
            print STDERR "-- WARNING:\n";
            sleep($delay);
        }

        if ($retries == 0) {
            caExit("failed to upload file '$pathname' to object store '$link'", undef);
        }
    }
}


#  Fetches a file from object storage with '$pathname' relative to the
#  assembly root directory.  Expects to be called from that same directory -
#  no big deal, since most of not all of Canu is run from that directory.
#
#  This function only fails if the download fails.
#
#  It succeeds if the file does not exist in object storage.  Client code can
#  therefore attempt to fetch a file from the store, and create/compute it if
#  it then still doesn't exist locally.
#
sub fetchFile ($) {
    my $dx       = getGlobal("objectStoreClient");
    my $da       = getGlobal("objectStoreClientDA");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");
    my $pathname = shift @_;
    my $path     = dirname($pathname);
    my $name     = basename($pathname);
    my $link;

    if ($path eq ".") {
        $link = sanitizeName("$pr:$ns/$name");
    } else {
        $link = sanitizeName("$pr:$ns/$path/$name");
    }

    return   if (-e $pathname);                 #  If it exists, we don't need to fetch it.
    return   if (! fileExists($pathname, 1));   #  If it doesn't exist in the store, we don't fetch it either.  Because it doesn't exist.

    if (isOS() eq "DNANEXUS") {
        print STDERR "-- Fetch '$pathname' from storage object '$link'.\n"   if ($showWork);

        if (! -d "$path") {
            print STDERR "--   Make path '$path'\n";
            make_path($path);
        }

        if (runCommandSilently(".", "$da download --output \"$pathname\" \"$link\"", 1)) {
            caExit("failed to download file from object store", undef);
        }
    }
}

#  Runs from the assembly root directory.
sub fetchFileFromLink ($$) {
    my $dx       = getGlobal("objectStoreClient");
    my $da       = getGlobal("objectStoreClientDA");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");
    my $link     = shift @_;
    my $file     = shift @_;
    my $path     = dirname($file);

    return   if (-e $file);                 #  If it exists, we don't need to fetch it.
    return   if (! fileExists($link, 1));   #  If it doesn't exist in the store, we don't fetch it either.  Because it doesn't exist.

    #  'dx describe --name file-Ff12PPQ0JbZVfKvK98zYKX17' will return the name too.

    if (isOS() eq "DNANEXUS") {
        print STDERR "-- Fetch '$file' from storage object '$link'.\n"   if ($showWork);
        make_path($path);

        if (runCommandSilently(".", "$da download --output \"$file\" \"$link\"", 1)) {
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
    my $dx     = getGlobal("objectStoreClient");
    my $ua     = getGlobal("objectStoreClientUA");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $path   = shift @_;     #  Path, relative to assembly root, of where we are called from.
    my $file   = shift @_;     #  Name of the file to stash, in current directory.
    my $indent = shift @_;
    my $code   = "";

    #  Not currently used by grid jobs, untested.
    #
    #if ($path eq ".") {
    #    $link = sanitizeName("$pr:$ns/$name");
    #    $path = sanitizeName("$ns");
    #} else {
    #    $link = sanitizeName("$pr:$ns/$path/$name");
    #    $path = sanitizeName("$ns/$path");
    #}

    if (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ -e \"./$file\" ] ; then\n";
        $code .= "${indent}  fold=`dirname \"$path/$file\"`\n";
        $code .= "${indent}  name=`basename \"$file\"`\n";
        $code .= "${indent}\n";
        $code .= "${indent}  if $dx describe --name \"$pr:$ns/$path/$file\" > /dev/null 2>&1 ; then\n";
        $code .= "${indent}    $dx rm --recursive \"$pr:$ns/$path/$file\"\n";
        $code .= "${indent}  fi\n";
        $code .= "${indent}\n";
        $code .= "${indent}  retries=$retryCount\n";
        $code .= "${indent}  while [ \$retries -gt 0 ] && \\\n";
        $code .= "${indent}        ! $ua --do-not-compress --wait-on-close --project \"$pr\" --folder \"$ns/\$fold/\" --name \"\$name\" \"./$file\" ; do\n";
        $code .= "${indent}    retries=`expr \$retries - 1`\n";
        $code .= "${indent}    echo \"Failed to stash file '$file', wait $retryDelay seconds and try again (\$retries times left).\"\n";
        $code .= "${indent}    sleep $retryDelay\n";
        $code .= "${indent}  done\n";
        $code .= "${indent}  if [ \$retries -eq 0 ] ; then\n";
        $code .= "${indent}    echo \"Failed to stash file '$file', removing incomplete copy in '$pr:$ns/$path/$file'.\"\n";
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


sub stashFilesShellCode ($$$$) {
    my $dx     = getGlobal("objectStoreClient");
    my $ua     = getGlobal("objectStoreClientUA");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $folder = shift @_;     #  Path to store files in.
    my $files  = shift @_;     #  Glob of files to store, in current directory.
    my $indent = shift @_;
    my $purge  = shift @_;
    my $code   = "";

    #  Upload a set of files to a directory.
    #
    #  If the fourth argument is "purge", The directory, and all files in it,
    #  is REMOVED from the object store before new files are uploaded.

    #  DO NOT quote "$files" below.  It could be a glob, or a list of
    #  multiple files.  If a glob, the quotes will prevent it from expanding.
    #  If multiple files, the quotes will make it one string.

    if (isOS() eq "DNANEXUS") {
        if ($purge eq "purge") {
            $code .= "\n";
            $code .= "${indent}  $dx rm --recursive \"$pr:$ns/$folder/\"\n";
            $code .= "${indent}  $dx mkdir          \"$pr:$ns/$folder/\"\n";
        }
        $code .= "\n";
        $code .= "${indent}retries=$retryCount\n";
        $code .= "${indent}while [ \$retries -gt 0 ] && \\\n";
        $code .= "${indent}      ! $ua --do-not-compress --wait-on-close --project \"$pr\" --folder \"$ns/$folder/\" $files ; do\n";
        $code .= "${indent}  retries=`expr \$retries - 1`\n";
        $code .= "${indent}  echo \"Failed to stash files '$files', wait $retryDelay seconds and try again (\$retries times left).\"\n";
        $code .= "${indent}  sleep $retryDelay\n";
        $code .= "${indent}done\n";
        $code .= "${indent}if [ \$retries -eq 0 ] ; then\n";
        $code .= "${indent}  echo \"Failed to stash files '$files', removing incomplete copy in '$pr:$ns/$folder/'.\"\n";
        $code .= "${indent}  $dx rm --recursive \"$pr:$ns/$folder/\"\n";
        $code .= "${indent}  exit 1\n";
        $code .= "${indent}fi\n";
        $code .= "\n";
    }

    return($code);
}


#  Fetch an objet store file '$path/$name' and place it at './$name'.
#
#  It is up to the caller to set the working directory BEFORE calling this
#  function.  Usually, the current directory is correctly set as is.
#
#  Occasionally, for example in overlapper, the native directory is
#  'trimming/1-overlapper/` but we need to fetch file
#  `trimming/0-mercounts/assembly.ms22.ignore.gz`:
#    cd ../0-mercounts
#    fetchFileShellCode("trimming/0-mercounts", "assembly.ms22.ignore.gz")
#    cd -
#
#  Or we can fetch files in subdirectories.  The 'blocks/' directory will be
#  created if needed.
#    fetchFileShellCode("correction/1-overlapper", "blocks/0001.dat")
#
#  Not that this DOES NOT fail if the file doesn't exist in object storage.
#  It is therefore possible to blindly fetch a file, then do something
#  depending on if the file now exists or not.
#
sub fetchFileShellCode ($$$) {
    my $dx     = getGlobal("objectStoreClient");
    my $da     = getGlobal("objectStoreClientDA");
    my $ns     = getGlobal("objectStoreNameSpace");
    my $pr     = getGlobal("objectStoreProject");
    my $path   = shift @_;     #  Path, relative to assembly root, of where we are called from.
    my $name   = shift @_;     #  Name of the file to fetch, in current directory.
    my $link;
    my $indent = shift @_;
    #my $code   = "";

    if ($path eq ".") {
        $link = sanitizeName("$pr:$ns/$name");
    } else {
        $link = sanitizeName("$pr:$ns/$path/$name");
    }

    return(fetchFileFromLinkShellCode($link, $name, $indent));

    #if (isOS() eq "DNANEXUS") {
    #    $code .= "\n";
    #    $code .= "${indent}if [ ! -e \"./$name\" ] ; then\n";
    #    $code .= "${indent}  mkdir -p `dirname \"./$name\"`\n";
    #    $code .= "${indent}  $da download --output \"./$name\" \"$link\"\n";
    #    $code .= "${indent}fi\n";
    #}

    #return($code);
}



sub fetchFileFromLinkShellCode ($$$) {
    my $dx       = getGlobal("objectStoreClient");
    my $da       = getGlobal("objectStoreClientDA");
    my $ns       = getGlobal("objectStoreNameSpace");
    my $pr       = getGlobal("objectStoreProject");
    my $link     = shift @_;
    my $name     = shift @_;
    my $indent   = shift @_;
    my $code     = "";

    if (isOS() eq "DNANEXUS") {
        $code .= "\n";
        $code .= "${indent}if [ ! -e \"./$name\" ] ; then\n";
        $code .= "${indent}  mkdir -p `dirname \"./$name\"`\n";
        $code .= "${indent}  $da download --output \"./$name\" \"$link\"\n";
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

        #  Stash the individual data files.  If they exist in the store, don't restash.

        for (my $bIdx="0000"; (-e "./$asm.seqStore/blobs.$bIdx"); $bIdx++) {
            if (! fileExists("$asm.seqStore/blobs.$bIdx", 1)) {
                stashFile("$asm.seqStore/blobs.$bIdx");
            }
        }

        #  Tar up the store metadata, then upload.  If we failed to upload any blob file,
        #  we won't upload the metadata, and the store will need to be regenerated.

        $files  = "./$asm.seqStore.err ";
        $files .= "./$asm.seqStore.sh ";
        $files .= "./$asm.seqStore/errorLog ";
        $files .= "./$asm.seqStore/info ";
        $files .= "./$asm.seqStore/info.txt ";
        $files .= "./$asm.seqStore/libraries ";
        $files .= "./$asm.seqStore/libraries.txt ";
        $files .= "./$asm.seqStore/readNames.txt ";
        $files .= "./$asm.seqStore/readlengths* ";
        $files .= "./$asm.seqStore/reads ";
        $files .= "./$asm.seqStore/reads-corc ";
        $files .= "./$asm.seqStore/reads-coru ";
        $files .= "./$asm.seqStore/reads-rawc ";
        $files .= "./$asm.seqStore/reads-rawu ";
        $files .= "./$asm.seqStore/version*"      if (-e "./$asm.seqStore/version.001");

        if (runCommandSilently(".", "tar -cf - $files | gzip -1c > ./$asm.seqStore.tar.gz", 1)) {
            caExit("failed to tar and compress sequence store files", undef);
        }

        stashFile("$asm.seqStore.tar.gz");

        unlink "./$asm.seqStore.tar.gz";
    }
}


sub fetchSeqStore ($) {
    my $asm    = shift @_;

    return   if (-e "./$asm.seqStore/info");
    return   if (! fileExists("$asm.seqStore.tar.gz", 1));

    if (defined(isOS())) {
        print STDERR "--\n"                                                   if ($showWork);
        print STDERR "-- Retrieving '$asm.seqStore' from object storage.\n"   if ($showWork);
        print STDERR "--\n"                                                   if ($showWork);

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
        $code .= "${indent}  echo \"In '`pwd`', fetching '$asm.seqStore' into '$root'.\"\n";
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



sub stashOvlStore ($$) {
    my $asm    = shift @_;
    my $base   = shift @_;

    return   if (! -e "./$base/$asm.ovlStore/index");

    if (defined(isOS())) {
        print STDERR "--\n"                                                   if ($showWork);
        print STDERR "-- Saving '$base/$asm.ovlStore' to object storage.\n"   if ($showWork);
        print STDERR "--\n"                                                   if ($showWork);

        #  Stash data files.

        for (my $bIdx="0001"; (-e "./$base/$asm.ovlStore/$bIdx<001>");   $bIdx++) {
        for (my $sIdx="001";  (-e "./$base/$asm.ovlStore/$bIdx<$sIdx>"); $sIdx++) {
            if (! fileExists("$base/$asm.ovlStore/$bIdx<$sIdx>", 1)) {
                stashFile("$base/$asm.ovlStore/$bIdx<$sIdx>");
            }
        }
        }

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
        $code .= "#  Upload data files.\n";
        $code .= "#\n";
        $code .= "\n";
        $code .= "cd $asm.ovlStore\n";
        $code .= "\n";
        $code .= stashFilesShellCode("$base/$asm.ovlStore", "????\\<???\\>", "", "purge");
        $code .= "\n";
        $code .= "cd ..\n";
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
        print STDERR "--\n"                                                         if ($showWork);
        print STDERR "-- Retrieving '$base/$asm.ovlStore' from object storage.\n"   if ($showWork);

        fetchFile("$base/$asm.ovlStore.tar.gz");

        if (runCommandSilently($base, "gzip -dc ./$asm.ovlStore.tar.gz | tar -xf -", 1)) {
            caExit("failed to uncompress and untar overlap store files", undef);
        }

        unlink "$base/$asm.ovlStore.tar.gz";

        print STDERR "--\n"                                                         if ($showWork);
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
        $code .= "${indent}  echo \"In `pwd`, fetching '$base' / '$asm.ovlStore.tar.gz', unzipping in '$root/$base'.\"\n";
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
    my $type   = shift @_;           #  Type of store: 'corStore' or 'ctgStore'.
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
    my $type   = shift @_;           #  Type of store: 'corStore' or 'ctgStore'.
    my $vers   = shift @_;           #  Version to fetch, '001' or '002'.
    my $indent = shift @_;           #
    my $base   = dirname($path);     #  'unitigging'
    my $root   = pathToDots($path);  #  '../..'
    my $code   = "";

    if (defined(isOS())) {
        $code .= "\n";
        $code .= "${indent}if [ ! -e $root/$base/$asm.$type/seqDB.v$vers.dat ] ; then\n";
        $code .= "${indent}  echo \"In `pwd`, fetching '$base' / '$asm.$type/seqDBv$vers'/\"\n";
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
        $code .= "${indent}  echo \"In `pwd`, storing '$path' '$name'.\"\n";
        $code .= "\n";
        $code .= "${indent}  cd $name\n";
        $code .= "\n";
        $code .= stashFilesShellCode("$path/$name", "0x??????.merylData", "${indent}  ", "purge");
        $code .= "\n";
        $code .= "${indent}  cd -\n";
        $code .= "\n";
        $code .= "${indent}  tar -cf - ./$name/*Index | gzip -1c > ./$name.tar.gz\n";
        $code .= stashFileShellCode($path, "$name.tar.gz", "${indent}  ");
        $code .= "${indent}  rm -f ./$name.tar.gz\n";
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
        $code .= "${indent}  echo \"In `pwd`, fetching '$path' '$name.tar.gz'.\"\n";
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
