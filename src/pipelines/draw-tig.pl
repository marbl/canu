#!/usr/bin/env perl

use SVG;       #  https://metacpan.org/pod/SVG, textproc/p5-SVG
use strict;

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }

#  -tig
#  -read
#  -pad
#  -region
#  -layout

my  $input;
my  $desiredTig;
my  %desiredRead;
my  %desiredEnd;

#  Coordinates of what to show.
#  displayPad will extend the region, if that region was not otherwise set.
my ($displayBgn, $displayEnd, $displayFixed, $displayPad) = (undef, undef, 0, 1000);

#  Draw a highlight around specific regions.
#  More than one highlight is allowed.
#  This is a list of '$bgn\0$end'
my  @highlights;

my  $layout;
my  $bogart;
my  %bogartClass;
my  $noContained = 0;

while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    #  Limit the display to a specific tig (if you give it all layouts for
    #  the whole assembly), or limit to a specific set of reads.
    #  Reads can be specified as "ID-3" or "ID-5" to have that read end highlighted,
    #  or as "ID-h" or "ID-l" to have the high/low coordinate highlighted.
    if    ($arg eq "-tig") {
        $desiredTig = shift @ARGV;
    }

    elsif ($arg eq "-read") {
        while ($ARGV[0] =~ m/^(\d+)-*(.*)$/) {
            my $r = $1;
            my $e = $2;

            $desiredRead{$r}  = 1;

            $desiredEnd{$r}  .= "\0$e"   if ( defined($desiredEnd{$r}) && (defined($e)));
            $desiredEnd{$r}   = "$e"     if (!defined($desiredEnd{$r}) && (defined($e)));

            shift @ARGV;
        }
    }

    #  A region to highlight, will pull in any read that intersects it and
    #  display (at least) this much.
    elsif ($arg eq "-highlight") {
        while ($ARGV[0] =~ m/^(\d+)-(\d+)$/) {
            my $b = $1;
            my $e = $2;

            push @highlights, "$b\0$e";

            shift @ARGV;
        }
    }

    #  A specific point to highlight.
    elsif ($arg eq "-flag") {
    }

    #  A region to display.  If not set, will pick an appropriate region
    #  based on -highlight, or -read or -tig.
    elsif ($arg eq "-display") {
        if ($ARGV[0] =~ m/^(\d+)-(\d+)$/) {
            $displayBgn   = $1;
            $displayEnd   = $2;
            $displayFixed = 1;

            shift @ARGV;
        }
    }

    elsif ($arg eq "-pad") {
        $displayPad = shift @ARGV;
    }

    #  Some bogart related options.  Load bogart annotations or turn off contained reads.
    elsif ($arg eq "-bogart") {
        $bogart = shift @ARGV;
    }

    elsif ($arg eq "-no-contained") {
        $noContained = 1;
    }

    #  Your input file!
    elsif ($arg eq "-layout") {
        $layout = shift @ARGV;
    }

    else {
        print STDERR "Unknown option '$arg'.\n";
    }
}

if (!defined($layout)) {
    print "usage: $0 ....\n";
    print "  -tig  T               limit drawing to just tig T\n";
    print "  -read R[-e]           limit drawing to just read R, marking end e\n";
    print "                          (any number of -read options can be supplied)\n";
    print "                          ('e' can be any of '3', '5' - to flag a specific end\n";
    print "                                          or 'l', 'h' - to flag a specific coordinate)\n";
    print "  -highlight bgn-end    draw a highlight box around coordinates bgn..end\n";
    print "                          (any number of highlight boxes can be drawn)\n";
    print "  -display bgn-end      explicitly specify the region to show.  if not\n";
    print "                        set, will be determined from -read and -highlight.\n";
    print "  -pad P                extend the display area by P bases on each end\n";
    print "  -bogart B             annotate reads with bogart status\n";
    print "  -no-contawined        do not draw contained reads (needs -bogart)\n";
    print "  -layout L             read the layout from tgStoreDump file L\n";
    print "\n";
    print "Output is SVG to stdout.\n";
    print "\n";
    exit(0);
}

#  Adjust or set displayBgn/displayEnd based on the highlight regions.
#  But only if the display region was not specified already.

if (($displayFixed == 0) &&
    (scalar(%desiredRead) == 0)) {
    foreach my $h (@highlights) {
        my ($hb, $he) = split '\0', $h;

        $displayBgn = $hb   if (!defined($displayBgn));
        $displayEnd = $he   if (!defined($displayEnd));

        $displayBgn = ($hb < $displayBgn) ? $hb : $displayBgn;
        $displayEnd = ($displayEnd < $he) ? $he : $displayEnd;
    }
}

#  Read the layout, and save any reads that either intersect a highlight region
#  or are a read flagged for interest.

my $tig;

open(LAY, "< $layout") or die "Failed to open '$layout' for reading: $!\n";
while (<LAY>) {
    my @v = split '\s+', $_;

    $tig = $v[1]   if ($v[0] eq "tig");

    next if (defined($desiredTig) && ($tig ne $desiredTig));   #  Skip lines for tigs we don't care about
    next if ($v[0] ne "read");                                 #  Skip non-read lines

    my $readID =  $v[1];
    my $bgn    = ($v[8] < $v[9]) ? $v[8] : $v[9];   #  Pick the lower coordinate
    my $end    = ($v[8] < $v[9]) ? $v[9] : $v[8];   #  Pick the higher coordinate

    #  If the read is listed as desired, check that all reads are from the
    #  desired tig (setting it if it isn't set) and then update the region to
    #  display to show that read.

    foreach my $r (keys %desiredRead) {
        next   if ($r ne $readID);

        if (!defined($desiredTig)) {
            $desiredTig = $tig;
        }

        if ($displayFixed == 0) {
            $displayBgn = $bgn   if (!defined($displayBgn));
            $displayEnd = $end   if (!defined($displayEnd));

            $displayBgn = ($bgn < $displayBgn) ? $bgn : $displayBgn;
            $displayEnd = ($displayEnd < $end) ? $end : $displayEnd;
        }

        print STDERR "  read $readID $bgn $end -> region $displayBgn $displayEnd\n";
    }

    #  If there are no desired reads listed, all we can do is check/set
    #  desired tig and update regions.

    if (scalar(%desiredRead) == 0) {
        if (!defined($desiredTig)) {
            $desiredTig = $tig;
        }

        if ($displayFixed == 0) {
            $displayBgn = $bgn   if (!defined($displayBgn));
            $displayEnd = $end   if (!defined($displayEnd));

            $displayBgn = ($bgn < $displayBgn) ? $bgn : $displayBgn;
            $displayEnd = ($displayEnd < $end) ? $end : $displayEnd;
        }

        print STDERR "  read $readID $bgn $end -> region $displayBgn $displayEnd\n";
    }
}

#  Extend the display region by any highlight region.

#foreach my $h (@highlights) {
#    my ($hb, $he) = split '\0', $h;
#
#    $displayBgn = ($hb < $displayBgn) ? $hb : $displayBgn;
#    $displayEnd = ($displayEnd < $he) ? $he : $displayEnd;
#}

#
#  Report what we're going to be showing.
#

print STDERR "Showing tig $desiredTig from ", $displayBgn - $displayPad, " to ", $displayEnd + $displayPad, "\n";

#
#  Load the bogart classifications
#

if (defined($bogart)) {
    print STDERR "Loading bogart classification from '$bogart'\n";
    open(B, "< $bogart") or die "Failed to open '$bogart': $!\n";
    $_ = <B>;    #  Header
    $_ = <B>;    #  Header
    while (<B>) {
        my @v = split '\s+', $_;

        $bogartClass{$v[0]} = $v[2];
    }
    close(B);
}

#
#  Read the layout again, loading all the tig positions.
#

my $minBgn;   #  Region that is actually drawn, includes reads that intersect
my $maxEnd;   #  the displayBgn/displayEnd as picked above.
my @reads;

open(LAY, "< $layout") or die "Failed to open '$layout' for reading: $!\n";
while (<LAY>) {
    my @v = split '\s+', $_;

    $tig = $v[1]   if ($v[0] eq "tig");

    next if ($tig ne $desiredTig);
    next if ($v[0] ne "read");

    my $readID    = $v[1];
    my $anchorID  = $v[3];
    my $ahang     = $v[5];
    my $bhang     = $v[6];
    my $bgn       = ($v[8] < $v[9]) ? $v[8] : $v[9];   #  Pick the lower coordinate
    my $end       = ($v[8] < $v[9]) ? $v[9] : $v[8];   #  Pick the higher coordinate
    my $orient    = ($v[8] < $v[9]) ? 'f'   : 'r';
    my $class     = (exists($bogartClass{$readID})) ? $bogartClass{$readID} : undef;

    next   if ($end < $displayBgn - $displayPad);   #  Ignore reads before the region we care about.
    next   if ($displayEnd + $displayPad < $bgn);   #  Ignore reads after  the region we care about.
    next   if (($noContained) && ($class eq "C----"));

    $minBgn = $bgn   if (!defined($minBgn));
    $maxEnd = $end   if (!defined($maxEnd));

    $minBgn = min($minBgn, $bgn);
    $maxEnd = max($maxEnd, $end);

    push @reads, "$readID\0$anchorID\0$ahang\0$bhang\0$bgn\0$end\0$orient\0$class";
}
close(LAY);


print STDERR "Found ", scalar(@reads), " reads intersecting $displayBgn-$displayEnd.\n";
print STDERR "Will draw $minBgn [ $displayBgn-$displayEnd ] $maxEnd\n";

my $nReads = scalar(@reads);

my $xBorder  = 50;
my $yBorder  = 50;

my $ySpacing = 20;

my $xPage    = $xBorder + 1000 + $xBorder;
my $yPage    = $yBorder + ($nReads - 1) * $ySpacing + $yBorder;

my $xScale   = ($xPage - 2 * $xBorder) / ($maxEnd - $minBgn);

my $yPos     = $yBorder;


#  Create an SVG object and make some drawing groups.

my $svg = SVG->new(width  => "${xPage}",
                   height => "${yPage}");

#  Draw any boxes around regions of interest

foreach my $h (@highlights) {
    my ($hb, $he) = split '\0', $h;

    my $xb = ($hb - $minBgn) * $xScale + $xBorder;
    my $xe = ($he - $minBgn) * $xScale + $xBorder;

    my $y = $svg->group(id    => "highlight$hb$he",
                        style => { 'stroke'       => 'rgb(0,0,0)',
                                   'fill'         => 'rgb(237,97,43)',
                                   'fill-opacity' => '255' });

    #print STDERR "RECTANGLE $xb to $xe from 0 to $yPage height\n";

    $y->rectangle(id => "highlightbox$hb$he",
                  x => $xb,         width  => $xe - $xb,
                  y => $yBorder-10, height => 10 + ($nReads-1) * $ySpacing + 10);
}

foreach my $read (@reads) {
    my ($readID, $anchorID, $ahang, $bhang, $bgn, $end, $orient, $class) = split '\0', $read;

    my $x1 = ($bgn - $minBgn) * $xScale + $xBorder;
    my $x2 = ($end - $minBgn) * $xScale + $xBorder;
    my $yp = $yPos;

    #print STDERR "read $readID $x1 $x2 at y $yp\n";

    my $y = $svg->group(id    => "read$readID",
                        style => { 'stroke' => 'rgb(0,0,0)',
                                   'fill'   => 'rgb(0,0,0)' });

    $y->line(id => "line$readID",
             x1 => $x1, y1 => $yp,
             x2 => $x2, y2 => $yp);

    #  Draw arrows.
    #
    #  fill stroke stroke-width stroke-opacity fill-opacity
    #  stoke-linecap: butt square round
    #  stroke-linejoin: miter round bevel
    #  stroke-dasharray="5,5"
    #  fill="none"
        
    if ($orient eq "f") {
        my $xv = [ $x2 - $ySpacing/3, $x2, $x2 - $ySpacing/3 ];
        my $yv = [ $yp - $ySpacing/3, $yp, $yp + $ySpacing/3 ];

        my $points = $y->get_path(x => $xv, y => $yv, -type => "polyline", -closed => 'false');
        
        $y->polyline(%$points, id => "arrow$readID",
                     style => {'fill'   => 'none',
                               'stroke' => 'rgb(255,0,0)'});

    } else {
        my $xv = [ $x1 + $ySpacing/3, $x1, $x1 + $ySpacing/3 ];
        my $yv = [ $yp - $ySpacing/3, $yp, $yp + $ySpacing/3 ];

        my $points = $y->get_path(x => $xv, y => $yv, -type => "polyline", -closed => 'false');

        $y->polyline(%$points, id => "arrow$readID",
                     style => {'fill'   => 'none',
                               'stroke' => 'rgb(255,0,0)'});
    }

    #  Label

    {
        my $t;

        if (exists($desiredRead{$readID})) {
            $t = $y->text(id => "label$readID",
                          x => $x2 + 1,
                          y => $yp,
                          style => { 'font'      => 'Courier',
                                     'stroke'    => 'rgb(255,0,0)',
                                     'font-size' => 12 });
        } else {
            $t = $y->text(id => "label$readID",
                          x => $x2 + 1,
                          y => $yp,
                          style => { 'font'      => 'Courier',
                                     'stroke'    => 'rgb(0,0,0)',
                                     'font-size' => 12 });
        }

        if (defined($class)) {
            $t->cdata("$readID $class");
        } else {
            $t->cdata("$readID");
        }
    }

    #  End points

    {
        if (exists($desiredEnd{$readID})) {
            my @tags = split '\0', $desiredEnd{$readID};

            foreach my $tag (@tags) {
                my $xc;

                $xc = $x1   if  ($tag eq "l");
                $xc = $x2   if  ($tag eq "h");
                $xc = $x1   if (($tag eq "5") && ($orient eq "f"));
                $xc = $x1   if (($tag eq "3") && ($orient eq "r"));
                $xc = $x2   if (($tag eq "3") && ($orient eq "f"));
                $xc = $x2   if (($tag eq "5") && ($orient eq "r"));

                $y->circle(cx => $xc, cy => $yp, r => $ySpacing/4, id => "mark$readID$tag");

                $y->line(x1 => $xc, y1 => $yBorder - $ySpacing / 2,
                         x2 => $xc, y2 => $yBorder - $ySpacing / 2 + scalar(@reads) * $ySpacing,
                         id => "endline$readID$tag");
            }
        }
    }

    $yPos += $ySpacing;
}

#  Render the SVG object, implicitly use svg namespace

#open(F, "> test.svg");
print $svg->xmlify;
#close(F);

exit(0);
