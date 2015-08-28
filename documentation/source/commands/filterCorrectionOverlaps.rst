filterCorrectionOverlaps
========================

::

  usage: filterCorrectionOverlaps [options]
  
  Rewrites an ovlStore, filtering overlaps that shouldn't be used for correcting reads.
  
    -G gkpStore     input reads
    -O ovlStore     input overlaps
    -S scoreFile    output scores for each read, binary file
  
    -c coverage     retain at most this many overlaps per read
  
    -l length       filter overlaps shorter than this length
    -e (min-)max    filter overlaps outside this range of fraction error
                      example:  -e 0.20          filter overlaps above 20% error
                      example:  -e 0.05-0.20     filter overlaps below 5% error
                                                              or above 20% error
  
    -logfile L      write detailed per-read logging to file L
  
  The following are not implemented:
  
    -nocontain      filter overlaps that are contained in the target read
    -nodovetail     filter dovetail overlaps
  
    -maxhang h      filter overlaps with more than 'h' bases unaligned on both the
                    target and evidence read, on either end
  ERROR: no gatekeeper store (-G) supplied.
  ERROR: no overlap store (-O) supplied.
  ERROR: no log file name (-f) supplied.
