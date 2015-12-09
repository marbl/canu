ovStoreIndexer
~~~~~~~~~~~~~~

::

  usage: ovStoreIndexer ...
    -o x.ovlStore    path to overlap store to build the final index for
    -F s             number of slices used in bucketizing/sorting
  
    -t x.ovlStore    explicitly test a previously constructed index
    -f               when testing, also create a new 'idx.fixed' which might
                     resolve rare problems
  
    -nodelete        do not remove intermediate files when the index is
                     successfully created
