tgStoreDump
===========

::

  usage: tgStoreDump -G <gkpStore> -T <tigStore> <v> [opts]
  
    -G <gkpStore>         Path to the gatekeeper store
    -T <tigStore> <v>     Path to the tigStore, version, to use
  
  
    -D <operation>        Dump something about the store
       list               ...a list of the unitigs in the store (NOT IMPLEMENTED)
       properties         ...a list of properties for ALL multialigns in the store (for -E)
  
    -u id[-id]            Unitig to dump (for -d option); if A-B, dump tigs from id A to id B, inclusive
    -c id[-id]            Contig to dump (for -d option); if A-B, dump tigs from id A to id B, inclusive
  
    -U                    Dump ALL unitigs (for -d option)
    -C                    Dump ALL contigs (for -d option)
  
    -nreads min max       Dump tigs with between min and max reads (inclusive)
  
    -d <operation>        Dump something about a multialign (-c or -u) in the store
       properties         ...properties
       frags              ...a list of fragments
       unitigs            ...a list of unitigs
       consensus [C]      ...the consensus sequence
                               if C supplied, only consensus with coverage >= C is output
       consensusgapped    ...the consensus sequence, with gaps as indicated in the multialignment
       layout             ...the layout
       multialign         ...the full multialignment
       sizes              ...an analysis of sizes of the tigs
       coverage           ...an analysis of read coverage of the tigs
       overlap            ...an analysis of read overlaps in the tigs
       fmap               ...a map from fragment IID to unitig IID
  
  
    -E <editFile>         Change properties of multialigns
  
  
    -compress             Move tigs from earlier versions into the specified version.  This removes
                          historical versions of unitigs/contigs, and can save tremendous storage space,
                          but makes it impossible to back up the assembly past the specified versions
  
    For '-d multialign':
    -w width              Width of the page.
    -s spacing            Spacing between reads on the same line.
  
    For '-d coverage':
    -o prefix             Output files will be written to 'prefix.*' in the current directory.
                          (defaults to 'tigStore' (the -t option) if not set.)
  
    For '-d sizes':
    -s genomesize         Denominator to use for n50 computation
