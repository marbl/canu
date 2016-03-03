readConsensus
~~~~~~

::

  usage: readConsensus ...
    -G gkpStore     Mandatory, path to gkpStore
  
  Inputs can come from either an overlap or a tig store.
    -O ovlStore     
    -T tigStore tigVers      
  
  If from an ovlStore, the range of reads processed can be restricted.
    -b bgnID        
    -e endID        
  
  Outputs will be written as the full multialignment and the final consensus sequence
    -c output.cns   
    -f output.fastq 
  
    -erate e        Overlaps are computed at 'e' fraction error; must be larger than the original erate
    -memory m       Use up to 'm' GB of memory
  
    -t n            Use up to 'n' cores
  
  ERROR: no gatekeeper (-G) supplied.
  ERROR: no inputs (-O or -T) supplied.
  ERROR: no outputs (-c or -f) supplied.
