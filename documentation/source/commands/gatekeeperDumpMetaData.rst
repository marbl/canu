gatekeeperDumpMetaData
~~~~~~

::

  usage: gatekeeperDumpMetaData -G gkpStore [p] [...]
  
    -G gkpStore [p]  dump reads from 'gkpStore', restricted to
                     partition 'p', if supplied.
  
    -libs            dump information about libraries
    -reads [-full]   dump information about reads
                       (-full also dumps some storage metadata)
  
    -stats           dump summary statistics on reads
  
    -b id            output starting at read/library 'id'
    -e id            output stopping after read/library 'id'
  
    -r id            output only the single read 'id'
  
  ERROR: no gkpStore (-G) supplied.
