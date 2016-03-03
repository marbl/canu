ovStoreIndexer
~~~~~~

::

  usage: ovStoreIndexer ...
    -O x.ovlStore    path to overlap store to build the final index for
    -F s             number of slices used in bucketizing/sorting
  
    -t x.ovlStore    explicitly test a previously constructed index
    -f               when testing, also create a new 'idx.fixed' which might
                     resolve rare problems
  
    -nodelete        do not remove intermediate files when the index is
                     successfully created
  
      DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER
      DANGER                                                DANGER
      DANGER   This command is difficult to run by hand.    DANGER
      DANGER          Use ovStoreCreate instead.            DANGER
      DANGER                                                DANGER
      DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER
  
  ERROR: No overlap store (-O) supplied.
  ERROR: One of -F (number of slices) or -t (test a store) must be supplied.
