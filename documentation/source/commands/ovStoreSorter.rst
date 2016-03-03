ovStoreSorter
~~~~~~

::

  usage: ovStoreSorter ...
    -O x.ovlStore    path to overlap store to build the final index for
    -G asm.gkpStore  path to gkpStore for this assembly
  
    -F s             number of slices used in bucketizing/sorting
    -job j m         index of this overlap input file, and max number of files
  
    -M m             maximum memory to use, in gigabytes
  
    -deleteearly     remove intermediates as soon as possible (unsafe)
    -deletelate      remove intermediates when outputs exist (safe)
  
    -force           force a recompute, even if the output exists
  
      DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER
      DANGER                                                DANGER
      DANGER   This command is difficult to run by hand.    DANGER
      DANGER          Use ovStoreCreate instead.            DANGER
      DANGER                                                DANGER
      DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER
  
  ERROR: No overlap store (-O) supplied.
  ERROR: no slice number (-F) supplied.
  ERROR: no max job ID (-job) supplied.
