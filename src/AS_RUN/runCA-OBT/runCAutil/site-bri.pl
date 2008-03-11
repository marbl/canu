

sub localDefaults() {

    #  Hmm.  We need to set the pathMap before runCA calls
    #  setParameters(), which is done before localSetup().
    #
    #  localDefaults() is called before any command line processing
    #  occurs, and we don't know where the assembly directory is, so
    #  we can't make the pathMap there.

    my $siteName = "/home/work/w/wgs/siteName";

    if (! -e $siteName) {
        open(F, "> $siteName") or die;
        print F "twobyfour.home  /home/work/w/wgs/FreeBSD-6.2-RELEASE-amd64/bin\n";
        print F "node5.home      /home/work/w/wgs/FreeBSD-7.0-RELEASE-amd64/bin\n";
        print F "node6.home      /home/work/w/wgs/FreeBSD-7.0-RELEASE-amd64/bin\n";
        print F "node7.home      /home/work/w/wgs/FreeBSD-7.0-RELEASE-amd64/bin\n";
        print F "node8.home      /home/work/w/wgs/FreeBSD-7.0-RELEASE-amd64/bin\n";
        close(F);
    }

    setGlobal("pathMap", $siteName);
}

sub localOption ($@) {
    my $arg  = shift @_;
    my @ARGV = @_;

    return($arg, @ARGV);
}

sub localSetup($) {
    my $numSteps = shift @_;
}

sub localStart ($) {
    my $commandToRun = shift @_;
}

sub localFinish ($) {
    my $commandThatRan = shift @_;
}

sub localFailure ($) {
    my $errorMessage = shift @_;
}

sub localPostTerminator($) {
    my $terminatorDirectory = shift @_;
}

sub localFinalize() {
}

1;
