use strict;

#  Global arg list.  New command line options MUST be added here, but
#  there is no mechanism to enforce this rule.
#
#  The value of the hash is the number of extra args needed.
#
my %validArgs = ("-help"            => "0",  # 
                 "-restart"         => "1",  # 
                 "-time"            => "1",  # 
                 "-mapest"          => "3",  # 
                 "-mapmrna"         => "3",  # 
                 "-mapsnp"          => "3",  # 

                 "-snpdelimiter"    => "1",  # 
                 "-snpoffset"       => "1",  # 
                 "-snppostag"       => "1",  # 
                 "-snpsizetag"      => "1",  # 
                 "-snpoutformat"    => "1",  #

                 "-abort"           => "1",  # 4
                 "-aligns"          => "0",  # 4
                 "-alwaysprint"     => "1",  # 4
                 "-batchsize"       => "1",  # 4
                 "-cleanup"         => "1",  # 5
                 "-extrahitfilter"  => "1",  # 3, undocumented
                 "-filternone"      => "0",  # 3
                 "-hitsortmemory"   => "1",  # 23
                 "-interspecies"    => "0",  # 4
                 "-localpolishes"   => "1",  # 4
                 "-localsearches"   => "1",  # 2
                 "-lsfjobname"      => "1",  # 234
                 "-lsfproject"      => "1",  # 234
                 "-lsfsearchqueue"  => "1",  # 2
                 "-lsffilterqueue"  => "1",  # 23
                 "-lsfpolishqueue"  => "1",  # 4
                 "-lsffinishqueue"  => "1",  # 45
                 "-maskmers"        => "1",  # 2, undocumented
                 "-mersize"         => "1",  # 2
                 "-mincoverage"     => "1",  # 45
                 "-minidentity"     => "1",  # 45
                 "-minlength"       => "1",  # 45
                 "-minsim4coverage" => "1",  # 4
                 "-minsim4identity" => "1",  # 4
                 "-minsim4length"   => "1",  # 4
                 "-noaligns"        => "0",  # 4
                 "-nocleanup"       => "0",  # 5
                 "-nostats"         => "0",  # 24
                 "-numbatches"      => "1",  # 4
                 "-relink"          => "1",  # 4
                 "-runlater"        => "0",  # 24
                 "-savetemporary"   => "0",  # 5
                 "-searchopts"      => "1",  # 2
                 "-searchthreads"   => "1",  # 2
                 "-sge"             => "1",  # 2, takes name
                 "-sgeaccount"      => "1",  # 2, takes account string
                 "-sgepriority"     => "1",  # 2, takes integer
                 "-species"         => "1",  # 2
                 "-stats"           => "0",  # 24
                 "-verbose"         => "0"); # 123



sub checkArgs {
    my $fail = 0;

    my %argFound;

    my $locl = 0;
    my $farm = 0;

    while (scalar(@_)) {
        my $a = shift @_;
        my $n = $validArgs{$a};

        $locl++ if ($a =~ m/^-local/);
        $farm++ if ($a =~ m/^-lsf/);

        if (defined($n)) {
            $argFound{$a} = 1;

            while ($n > 0) {
                my $b = shift @_;
                if (!defined($b) || defined($validArgs{$b})) {
                    print STDERR "ERROR: Too few parameters given to option '$a'\n";
                    $fail = 1;
                    last;
                }
                $n--;
            }
        } else {
            print STDERR "ERROR: Unknown option '$a' supplied.\n";
            $fail = 1;
        }
    }

    exit(1) if ($fail);

    if (($locl > 0) &&
        ($farm > 0)) {
        print STDERR "ERROR: Both local and farm modes requested.\n";
        exit(1);
    }

    if (($farm > 0) && ($farm != 6)) {
        print STDERR "ERROR: Incomplete LSF specification.\n";
        print STDERR "         missing -lsfjobname\n"     if (!defined($argFound{"-lsfjobname"}));
        print STDERR "         missing -lsfproject\n"     if (!defined($argFound{"-lsfproject"}));
        print STDERR "         missing -lsfsearchqueue\n" if (!defined($argFound{"-lsfsearchqueue"}));
        print STDERR "         missing -lsffilterqueue\n" if (!defined($argFound{"-lsffilterqueue"}));
        print STDERR "         missing -lsfpolishqueue\n" if (!defined($argFound{"-lsfpolishqueue"}));
        print STDERR "         missing -lsffinishqueue\n" if (!defined($argFound{"-lsffinishqueue"}));
        exit(1);
    }
}



1;
