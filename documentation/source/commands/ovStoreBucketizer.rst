ovStoreBucketizer
~~~~~~

::

  usage: ovStoreBucketizer -O asm.ovlStore -G asm.gkpStore -i file.ovb.gz -job j [opts]
    -O asm.ovlStore       path to store to create
    -G asm.gkpStore       path to gkpStore for this assembly
  
    -C config             path to previously created ovStoreBuild config data file
  
    -i file.ovb.gz        input overlaps
    -job j                index of this overlap input file
  
    -F f                  use up to 'f' files for store creation
  
    -obt                  filter overlaps for OBT
    -dup                  filter overlaps for OBT/dedupe
  
    -e e                  filter overlaps above e fraction error
  
    -raw                  write uncompressed buckets
  
      DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER
      DANGER                                                DANGER
      DANGER   This command is difficult to run by hand.    DANGER
      DANGER          Use ovStoreCreate instead.            DANGER
      DANGER                                                DANGER
      DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER
  
  ERROR: No overlap store (-O) supplied.
  ERROR: No gatekeeper store (-G) supplied.
  ERROR: No input (-i) supplied.
  ERROR: No job index (-job) supplied.
