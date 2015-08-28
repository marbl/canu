ovStoreBuild
============

::

  usage: ovStoreBuild -o asm.ovlStore -g asm.gkpStore [opts] [-L fileList | *.ovb.gz]
    -o asm.ovlStore       path to store to create
    -g asm.gkpStore       path to gkpStore for this assembly
  
    -F f                  use up to 'f' files for store creation
    -M m                  use up to 'm' MB memory for store creation
  
    -e e                  filter overlaps above e fraction error
    -l l                  filter overlaps below l bases overlap length (needs gkpStore to get read lengths!)
  
    -L fileList           read input filenames from 'flieList'
  
    -big iid              handle a large number of overlaps in the last library
                          iid is the first read iid in the last library, from
                          'gatekeeper -dumpinfo *gkpStore'
  
    -evalues              Input files are evalue updates from overlap error adjustment
  
  ERROR: No overlap store (-o) supplied.
  ERROR: No gatekeeper store (-g) supplied.
  ERROR: No input overlap files (-L or last on the command line) supplied.
