package ca3g::OverlapStore;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(createOverlapStore);

use strict;

use ca3g::Defaults;
use ca3g::Execution;




sub createOverlapStoreSequential ($$$@) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $files        = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= " -o $wrk/$asm.ovlStore.BUILDING \\\n";
    $cmd .= " -g $wrk/$asm.gkpStore \\\n";
    $cmd .= " -M " . getGlobal("ovlStoreMemory") . " \\\n";
    $cmd .= " -L $files \\\n";
    $cmd .= " > $wrk/$asm.ovlStore.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("failed to create the overlap store", "$wrk/$asm.ovlStore.err");
    }

    rename "$wrk/$asm.ovlStore.BUILDING", "$wrk/$asm.ovlStore";
}




sub createOverlapStoreParallel ($$$@) {
    die "Not yet.";
}



sub createOverlapStore ($$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $seq          = shift @_;
    my $path         = "$wrk/1-overlapper";

    goto alldone if (-d "$wrk/$asm.ovlStore");
    goto alldone if (-d "$wrk/$asm.tigStore");

    #  Did we _really_ complete?

    caFailure("overlapper claims to be finished, but no job list found in '$path/ovljob.files'", undef)  if (! -e "$path/ovljob.files");
    caFailure("overlapper claims to be finished, but no job list found in '$path/ovljob.files'", undef)  if (  -z "$path/ovljob.files");

    #  Then just build the store!  Simple!

    createOverlapStoreSequential($wrk, $asm, "$path/ovljob.files")  if ($seq eq "sequential");
    createOverlapStoreParallel  ($wrk, $asm, "$path/ovljob.files")  if ($seq eq "parallel");

    goto alldone  if (getGlobal("saveOverlaps"));

    #  Delete the inputs and directories.

    my %directories;

    open(F, "< $path/ovljob.files");
    while (<F>) {
        chomp;
        unlink "$path/$_";

        my @components = split '/', $_;
        pop @components;
        my $dir = join '/', @components;
        
        $directories{$dir}++;
    }
    close(F);
    
    foreach my $dir (keys %directories) {
        rmdir "$path/$dir";
    }
    
    unlink "$path/ovljob.files";

    #  Now all done!
  alldone:
    stopAfter("overlapper");
}

1;
