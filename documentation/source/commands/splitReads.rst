splitReads
~~~~~~

::

  usage: splitReads -G gkpStore -O ovlStore -Ci input.clearFile -Co output.clearFile -o outputPrefix]
  
    -G gkpStore    path to read store
    -O ovlStore    path to overlap store
  
    -o name        output prefix, for logging
  
    -t bgn-end     limit processing to only reads from bgn to end (inclusive)
  
    -Ci clearFile  path to input clear ranges (NOT SUPPORTED)
    -Co clearFile  path to ouput clear ranges
  
    -e erate       ignore overlaps with more than 'erate' percent error
  
    -minlength l   reads trimmed below this many bases are deleted
  
