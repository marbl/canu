tgStoreCoverageStat
~~~~~~~~~~~~~~~~~~~

::

  usage: tgStoreCoverageStat -g gkpStore -t tigStore version
  
    -g <G>     Mandatory, path G to a gkpStore directory.
    -t <T> <v> Mandatory, path T to a tigStore, and version V.
    -s <S>     Optional, assume genome size S.
    -o <name>  Recommended, prefix for output files.
    -n         Do not update the tigStore (default = do update).
    -u         Do not estimate based on N50 (default = use N50).
  
    -L         Be leniant; don't require reads start at position zero.
  No gatekeeper store (-g option) supplied.
  No input tigStore (-t option) supplied.
  No output prefix (-o option) supplied.
