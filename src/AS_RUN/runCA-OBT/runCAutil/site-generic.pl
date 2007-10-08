
#  This file allows customizations based on your local environment.
#  You can set site-wide default options, like SGE accounts,
#  priorities, queues.  You can do something before and after each
#  major component runs -- JCVI likes to update a database showing
#  progress, users can use the web to display this info.  Finally, you
#  can post-process the assembly.


sub localDefaults() {
}

sub localOption ($@) {
    my $arg  = shift @_;
    my @ARGV = @_;
    my $found = 1;

    #  This function is given the chance to process command line
    #  options.  If it finds an option it likes, it can do whatever it
    #  wants to process the option, including modifying the @ARGV
    #  array (for example, by shifting off the next value).
    #
    if ($arg =~ m/^-thisIsNotAnOption/ ) {
        #  $local_value needs to be declared ("my $local_value")
        #  somewhere in this file.  It can then be used in the other
        #  local*() functions.
        #
        #$local_value = shift @ARGV;
    } else {
        $found = 0;
    }

    #  The else block that sets $found to zero is VERY IMPORTANT!  If
    #  we processed an option, we are responsible for returning the
    #  next value.  It's a bit of a burden, and unfortunately, the
    #  best method I could think of.
    #
    if ($found) {
        $arg = shift @ARGV;
    }

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

1;
