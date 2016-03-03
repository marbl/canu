ovStoreBuild
~~~~~~

::

  usage: ovStoreBuild -O asm.ovlStore -G asm.gkpStore [opts] [-L fileList | *.ovb.gz]
    -O asm.ovlStore       path to store to create
    -G asm.gkpStore       path to gkpStore for this assembly
  
    -L fileList           read input filenames from 'flieList'
  
    -F f                  use up to 'f' files for store creation
    -M g                  use up to 'g' gigabytes memory for sorting overlaps
                            default 4; g-0.125 gb is available for sorting overlaps
  
    -e e                  filter overlaps above e fraction error
    -l l                  filter overlaps below l bases overlap length (needs gkpStore to get read lengths!)
  
  Non-building options:
    -evalues              input files are evalue updates from overlap error adjustment
    -config out.dat       don't build a store, just dump a binary partitioning file for ovStoreBucketizer
  
  ERROR: No overlap store (-o) supplied.
  ERROR: No gatekeeper store (-g) supplied.
  ERROR: No input overlap files (-L or last on the command line) supplied.
