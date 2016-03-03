trimReads
~~~~~~

::

  usage: trimReads -G gkpStore -O ovlStore -Co output.clearFile -o outputPrefix
  
    -G gkpStore    path to read store
    -O ovlStore    path to overlap store
  
    -o name        output prefix, for logging
  
    -t bgn-end     limit processing to only reads from bgn to end (inclusive)
  
    -Ci clearFile  path to input clear ranges (NOT SUPPORTED)
    -Co clearFile  path to ouput clear ranges
  
    -e erate       ignore overlaps with more than 'erate' percent error
  
    -ol l          the minimum evidence overlap length
    -oc c          the minimum evidence overlap coverage
                     evidence overlaps must overlap by 'l' bases to be joined, and
                     must be at least 'c' deep to be retained
  
    -minlength l   reads trimmed below this many bases are deleted
  
