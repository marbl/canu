
use strict;

#  Global arg list.  New command line options MUST be added here, but
#  there is no mechanism to enforce this rule.
#
#  The value of the hash is the number of extra args needed.
#
my %validArgs = ("-help"            => "0",
                 "-restart"         => "1",
                 "-configure"       => "1",
                 "-searchest"       => "1",
                 "-searchmrna"      => "1",
                 "-filterest"       => "1",
                 "-filtermrna"      => "1",
                 "-polish"          => "1",
                 "-assembleoutput"  => "1",
                 "-mapest"          => "3",
                 "-mapest-nofilter" => "3",
                 "-mapsnp"          => "3",
                 "-mapesttoest"     => "3",
                 "-mapmrna"         => "3",

                 "-snpdelimiter"    => "1",
                 "-snpoffset"       => "1",
                 "-snppostag"       => "1",
                 "-snpsizetag"      => "1",
                 "-abort"           => "1",
                 "-aligns"          => "0",
                 "-alwaysprint"     => "1",
                 "-batchsize"       => "1",
                 "-cdna"            => "1",
                 "-cleanup"         => "0",
                 "-farmpolishes"    => "2",
                 "-farmsearches"    => "2",
                 "-genomic"         => "1",
                 "-hitsortmemory"   => "1",
                 "-interspecies"    => "0",
                 "-localpolishes"   => "1",
                 "-localsearches"   => "1",
                 "-longintron"      => "1",
                 "-maskmers"        => "1",
                 "-maxintron"       => "1",
                 "-memory"          => "1",
                 "-mersize"         => "1",
                 "-mincoverage"     => "1",
                 "-minidentity"     => "1",
                 "-minlength"       => "1",
                 "-minsim4coverage" => "1",
                 "-minsim4identity" => "1",
                 "-minsim4length"   => "1",
                 "-noaligns"        => "0",
                 "-nocleanup"       => "0",
                 "-nomaskmers"      => "0",
                 "-nostats"         => "0",
                 "-numbatches"      => "1",
                 "-prebuild"        => "2",
                 "-relink"          => "1",
                 "-runlater"        => "0",
                 "-searchopts"      => "1",
                 "-searchthreads"   => "1",
                 "-segments"        => "1",
                 "-skiphitfilter"   => "0",
                 "-species"         => "1",
                 "-stats"           => "0",
                 "-verbose"         => "0" );



sub checkArgs {
    my $fail = 0;

    while (scalar(@_)) {
        my $a = shift @_;
        my $n = $validArgs{$a};
        if (defined($n)) {
            while ($n > 0) {
                my $b = shift @_;
                if (!defined($b) || defined($validArgs{$b})) {
                    print "ERROR: Too few parameters given to option '$a'\n";
                    $fail = 1;
                    last;
                }
                $n--;
            }
        } else {
            print "ERROR: Unknown option '$a' supplied.\n";
            $fail = 1;
        }
    }

    exit(1) if ($fail);
}



1;
