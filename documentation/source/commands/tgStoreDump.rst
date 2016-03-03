tgStoreDump
~~~~~~

::

  usage: tgStoreDump -G <gkpStore> -T <tigStore> <v> [opts]
  
  STORE SELECTION (mandatory)
  
    -G <gkpStore>           path to the gatekeeper store
    -T <tigStore> <v>       path to the tigStore, version, to use
  
  TIG SELECTION - if nothing specified, all tigs are reported
                - all ranges are inclusive.
  
    -tig A[-B]              only dump tigs between ids A and B
    -unassembled            only dump tigs that are 'unassembled'
    -bubbles                only dump tigs that are 'bubbles'
    -contigs                only dump tigs that are 'contigs'
    -nreads min max         only dump tigs with between min and max reads
    -length min max         only dump tigs with length between 'min' and 'max' bases
    -coverage c C g G       only dump tigs with between fraction g and G at coverage between c and C
                              example:  -coverage 10 inf 0.5 1.0 would report tigs where half of the
                                        bases are at 10+ times coverage.
  
  DUMP TYPE - all dumps, except status, report on tigs selected as above
  
    -status                 the number of tigs in the store
  
    -tigs                   a list of tigs, and some information about them
  
    -consensus [opts]       the consensus sequence, with options:
                              -gapped           report the gapped (multialignment) consensus sequence
                              -fasta            report sequences in FASTA format (the default)
                              -fastq            report sequences in FASTQ format
  
    -layout [opts]          the layout of reads in each tig
                            if '-o' is supplied, three files are created, otherwise just the layout is printed to stdout
                              -gapped           report the gapped (multialignment) positions
                              -o outputPrefix   write plots to 'outputPrefix.*' in the current directory
  
    -multialign [opts]      the full multialignment, output is to stdout
                              -w width          width of the page
                              -s spacing        spacing between reads on the same line
  
    -sizes [opts]           size statistics
                              -s genomesize     denominator to use for n50 computation
  
    -coverage [opts]        read coverage plots, one plot per tig
                              -o outputPrefix   write plots to 'outputPrefix.*' in the current directory
  
    -depth [opts]           a histogram of depths
                              -single           one histogram per tig
  
    -overlap                read overlaps
                              -thin overlap     report regions where the (thickest) read overlap is less than 'overlap' bases
  
    -overlaphistogram       a histogram of the thickest overlaps used
                              -o outputPrefix   write plots to 'outputPrefix.*' in the current directory
  
  
