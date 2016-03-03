tgStoreCoverageStat
~~~~~~

::

  usage: tgStoreCoverageStat -G gkpStore -T tigStore version -o output-prefix [-s genomeSize] ...
  
    -G <G>     Mandatory, path G to a gkpStore directory.
    -T <T> <v> Mandatory, path T to a tigStore, and version V.
    -o <name>  Mandatory, prefix for output files.
    -s <S>     Optional, assume genome size S.
  
    -n         Do not update the tigStore (default = do update).
    -u         Do not estimate based on N50 (default = use N50).
  
    -L         Be leniant; don't require reads start at position zero.
  
  No gatekeeper store (-G option) supplied.
  No input tigStore (-T option) supplied.
  No output prefix (-o option) supplied.
