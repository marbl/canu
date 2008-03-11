use ExtUtils::MakeMaker;

my %AS;
my (@del, $i);
$i = 0;
for my $s (@ARGV) {
    if($s =~ /AS_(.*?)=(.*)/) {
        $AS{$1} = $2;
        push(@del, $i);
    }
    $i++;
}
for $i (@del) {
    splice(@ARGV, $i, 1);
}

unless($AS{"BASE"} || ($AS{"LIB"} && $AS{"SRC"})) {
    die("The path to celera assembler must be specified. E.g. AS_BASE=/...");
}
unless($AS{"OS"}) {
    my $ostype = `uname`;
    my $machinetype = `uname -m`;
    chomp $ostype;
    chomp $machinetype;
    $machinetype = "amd64" if($machinetype eq "x86_64");
    $machinetype = "ppc" if($machinetype eq "Power Macintosh");
    $AS{"OS"} = "${ostype}-${machinetype}";
}
$AS{"LIB"} = $AS{"BASE"} . "/" . $AS{"OS"} . "/lib" unless($AS{"LIB"});
$AS{"SRC"} = $AS{"BASE"} . "/" . "/src" unless($AS{"SRC"});

my %PM;
my $OBJECT = "MSG.o";
my $CLEAN  = "swigperlrun.h";


# Generate module from swig
my @SWIGS = qw(MSG);
sub MY::postamble {
    return ".PRECIOUS: %_wrap.c %.pm
%_wrap.d: %.swig
\tswig -perl5 \$(INC) -MM \$< | sed -e 's/: / \$\@:/g' > \$\@

MMCCMD = \$(CC) -MM -MG \$(PASTHRU_INC) \$(INC) \\
        \$(CCFLAGS) \$(OPTIMIZE) \\
        \$(PERLTYPE) \$(MPOLLUTE) \$(DEFINE_VERSION) \\
        \$(XS_DEFINE_VERSION) \$(CCCDLFLAGS) \"-I\$(PERL_INC)\" \\
        \$(PASTHRU_DEFINE) \$(DEFINE)

%.d: %.c
\t\$(MMCCMD) \$< > \$\@

swigperlrun.h: \$(addsuffix _wrap.c, @{SWIGS})
\tswig -perl5 -external-runtime \$\@

%_wrap.c %.pm: %.swig
\tswig -perl5 -shadow \$(INC) \$<


include \$(addsuffix _wrap.d, @{SWIGS})
include \$(addsuffix .d, @{SWIGS})
"

}

foreach my $pm (@SWIGS) {
    $PM{$pm . '.pm'} = '$(INST_LIB)/AS/' . $pm . '.pm';
    $OBJECT .= ' ' . $pm . '_wrap$(OBJ_EXT)';
    $CLEAN .= ' ' . $pm . '_wrap.c' . ' ' . $pm . '_wrap.d' . ' ' . $pm . '.pm' . ' ' . $pm . '.d';
}
$OBJECT .= " $AS{'LIB'}/libCA.a";
$PM{"Parser.pm"} = '$(INST_LIB)/AS/MSG/Parser.pm';

WriteMakefile('NAME'            => 'AS::MSG',
              'VERSION'         => '0.0.1',
              'LIBS'            => ['-lm'],
              'PM'              => \%PM,
              'INC'             => "-I$AS{'SRC'} -I$AS{'SRC'}/AS_UTL -I/usr/include",
              'OBJECT'          => $OBJECT,
              'CC'              => 'gcc',
              'DEFINE'          => '-DX86_GCC_LINUX',
              'clean'           => { FILES => $CLEAN },
    );
