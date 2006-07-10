#!/usr/bin/perl

my $genomeDir = "/bioinfo/assembly/walenz/GENOMES";
my $id1;
my $id2;

while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-g/) {
        $genomeDir = shift @ARGV;
    } elsif ($arg =~ m/^-1/) {
        $id1 = shift @ARGV;
    } elsif ($arg =~ m/^-2/) {
        $id2 = shift @ARGV;
    }

}

die "usage: $0 [-genomedir path] -1 id1 -2 id2\n" if (!defined($id1) || !defined($id2));

die "No bin dir?\n" if (! -e "/bioinfo/assembly/walenz/src/genomics/linux64/bin");
die "No bin?\n"     if (! -x "/bioinfo/assembly/walenz/src/genomics/linux64/bin/snapper2");

$genomeDir = "$ENV{'PWD'}/$genomeDir" if ($genomeDir !~ m!^/!);

my $name = "${id1}vs${id2}";
my $cmd;

$cmd .= "perl /bioinfo/assembly/walenz/src/genomics/atac-driver/briatac.pl ";
$cmd .= "  -dir $name ";
$cmd .= "  -id1 $id1 -id2 $id2 ";
$cmd .= "  -genomedir $genomeDir ";
$cmd .= "  -meryldir  $genomeDir ";
$cmd .= "  -bindir    /bioinfo/assembly/walenz/src/genomics/linux64/bin ";
$cmd .= "  -merylthreads 4 ";
$cmd .= "  -numsegments 2 ";
$cmd .= "  -numthreads  4 ";
$cmd .= "  -samespecies";
print "$cmd\n";
system($cmd) and die "Failed to briatac.pl!\n";


$cmd  = "cd $name && ln -s $name.k20.u1.f20.g0.matches.sorted.extended $name.atac && ";
$cmd .= "time sh /bioinfo/assembly/walenz/src/genomics/atac/atacdriver.sh $name.atac && ";
$cmd .= "grep ^M $name.atac.ckpLast | cut -d' ' -f 1-12 | sort -k5,5 -k6n > $name.atac.ckpLast.sorted && ";
$cmd .= "/bioinfo/assembly/walenz/src/genomics/atac-driver/clumpMaker/clumpMaker -c 5000 -2 -S -f $name.atac.ckpLast.sorted > $name.atac.ckpLast.clumps";
system($cmd) and die "Failed to atacdriver.sh!\n";


