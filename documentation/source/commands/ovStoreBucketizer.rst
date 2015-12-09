ovStoreBucketizer
~~~~~~~~~~~~~~~~~

::

  usage: ovStoreBucketizer -c asm.ovlStore -g asm.gkpStore -i file.ovb.gz -job j [opts]
    -c asm.ovlStore       path to store to create
    -g asm.gkpStore       path to gkpStore for this assembly
  
    -i file.ovb.gz        input overlaps
    -job j                index of this overlap input file
  
    -F f                  use up to 'f' files for store creation
  
    -obt                  filter overlaps for OBT
    -dup                  filter overlaps for OBT/dedupe
  
    -e e                  filter overlaps above e fraction error
  
    -raw                  write uncompressed buckets
  ERROR: No overlap store (-o) supplied.
  ERROR: No gatekeeper store (-g) supplied.
  ERROR: No input (-i) supplied.
  ERROR: No job index (-job) supplied.
