#!/usr/local/bin/perl


sub splitFile {
    my ($name) = @_;
    my $len = 0;
    my $cnt = 0;
    my $idx = "00";

    open(F, "< $name");
    while (!eof(F)) {

        print "bzip2 -9v $name.$idx &\n";
        open(Z, "> $name.$idx");

        $_ = <F>;
        while ((!eof(F)) && ($len < 200000000)) {
            $len += 1 + length $_;
            $cnt++;
            print Z $_;
            $_ = <F>;
        }
        print Z $_;

        #
        #  Read until the next sim4end
        #

        if (!eof(F)) {
            $_ = <F>;
            while ((!eof(F)) && (! m/sim4end/)) {
                print Z $_;
                $_ = <F>;
            }
            print Z $_;
        }

        $len = 0;

        close(Z);

        $idx++;
    }
    close(F);
}

splitFile("/work/assembly/walenzbp/dbEST-20020331/polishes-good");

