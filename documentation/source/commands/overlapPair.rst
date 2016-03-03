overlapPair
~~~~~~

::

  usage: overlapPair ...
    -G gkpStore     Mandatory, path to gkpStore
  
  Inputs can come from either a store or a file.
    -O ovlStore     
    -O ovlFile      
  
  If from an ovlStore, the range of reads processed can be restricted.
    -b bgnID        
    -e endID        
  
  Outputs will be written to a store or file, depending on the input type
    -o ovlStore     
    -o ovlFile      
  
    -erate e        Overlaps are computed at 'e' fraction error; must be larger than the original erate
    -partial        Overlaps are 'overlapInCore -G' partial overlaps
    -memory m       Use up to 'm' GB of memory
  
    -t n            Use up to 'n' cores
  
  Advanced options:
  
    -invert         Invert the overlap A <-> B before aligning (they are not re-inverted before output)
  
