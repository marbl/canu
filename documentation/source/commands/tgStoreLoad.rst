tgStoreLoad
~~~~~~

::

  usage: tgStoreLoad -G <gkpStore> -T <tigStore> <v> [input.cns]
  
    -G <gkpStore>         Path to the gatekeeper store
    -T <tigStore> <v>     Path to the tigStore and version to add tigs to
  
    -L <file-of-files>    Load the tig(s) from files listed in 'file-of-files'
  
    -n                    Don't replace, just report what would have happened
  
    The primary operation is to replace tigs in the store with ones in a set of input files.
    The input files can be either supplied directly on the command line or listed in
    a text file (-L).
  
    A new store is created if one doesn't exist, otherwise, whatever tigs are there are
    replaced with those in the -R file.  If version 'v' doesn't exist, it is created.
  
    Even if -n is supplied, a new store is created if one doesn't exist.
  
    To add a new tig, give it a tig id of -1.  New tigs must be added to the latest version.
    To delete a tig, remove all children, and set the number of them to zero.
  
  ERROR:  no gatekeeper store (-G) supplied.
  ERROR:  no tig store (-T) supplied.
  ERROR:  no input tigs (-R) supplied.
