filterCorrectionOverlaps
~~~~~~

::

  usage: filterCorrectionOverlaps [options]
  
  Rewrites an ovlStore, filtering overlaps that shouldn't be used for correcting reads.
  
    -G gkpStore     input reads
    -O ovlStore     input overlaps
    -S scoreFile    output scores for each read, binary file, to 'scoreFile'
                    per-read logging to 'scoreFile.log' (see -nolog)
                    summary statistics to 'scoreFile.stats' (see -nostats)
  
    -c coverage     retain at most this many overlaps per read
  
    -l length       filter overlaps shorter than this length
    -e (min-)max    filter overlaps outside this range of fraction error
                      example:  -e 0.20          filter overlaps above 20% error
                      example:  -e 0.05-0.20     filter overlaps below 5% error
                                                              or above 20% error
  
    -nolog          don't create 'scoreFile.log'
    -nostats        don't create 'scoreFile.stats'
  ERROR: no gatekeeper store (-G) supplied.
  ERROR: no overlap store (-O) supplied.
  ERROR: no output scoreFile (-S) supplied.
