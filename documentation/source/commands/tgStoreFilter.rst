tgStoreFilter
~~~~~~~~~~~~~

::

  usage: tgStoreFilter -g gkpStore -t tigStore version
  
    -G <G>       Mandatory, path G to a gkpStore directory.
    -T <T> <v>   Mandatory, path T to a tigStore, and version V.
  
    -j J         Unitig is not unique if astat is below J (cgbUniqueCutoff)
    -k K         (unused) (cgbDefinitelyUniqueCutoff)
  
    -span F      Unitig is not unique if a single read spans more than fraction F (default 1.0) of unitig
    -lowcov D F  Unitig is not unique if fraction F (default 1.0) of unitig is below read depth D (default 2)
    -reads R     Unitig is not unique if unitig has fewer than R (default 2) reads
                 If R is fractional, the least populous unitigs containing fraction R of reads are marked as repeat
                 Example: unitigs with 9 or fewer reads contain 10% of the reads.  -reads 0.10 would mark these are repeat.
    -long L      Unitig is unique if unitig is at least L (default unlimited) bases long
    -short S     Unitig is not unique if unitig is shorter than S (default 1000) bases long
  
    -o <name>    Prefix for output files.
    -n           Do not update the tigStore.
  
  Algorithm:  The first rule to trigger will mark the unitig.
  
    1)  A unitig with a single read is NOT unique.
    2)  A unitig with fewer than R (-reads) reads is NOT unique.
    3)  A unitig with a single read spanning fraction F (-span) of the unitig is NOT unique.
    4)  A unitig longer than L (-length) bases IS unique.
    5)  A unitig with astat less than J (-j) is NOT unique.
    6)  A unitig with fraction F below coverage D (-lowcov) is NOT unique.
    7)  A unitig shorter than S (-short) bases long is NOT unique.
    8)  Otherwise, the unitig IS unique.
  No gatekeeper store (-G option) supplied.
  No input tigStore (-T option) supplied.
  No output prefix (-o option) supplied.
