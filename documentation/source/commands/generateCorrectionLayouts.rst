generateCorrectionLayouts
~~~~~~

::

  usage: generateCorrectionLayouts -G gkpStore -O ovlStore [ -T tigStore | -F ] ...
    -G gkpStore   mandatory path to gkpStore
    -O ovlStore   mandatory path to ovlStore
  
    -S file       global score (binary) input file
  
    -T corStore   output layouts to tigStore corStore
    -F            output falconsense-style input directly to stdout
  
    -p  name      output prefix name, for logging and summary
  
    -b  bgnID     
    -e  endID     
  
    -rl file      
  
    -L  length    minimum length of evidence overlaps
    -E  erate     maxerror rate of evidence overlaps
  
    -C  coverage  maximum coverage of evidence reads to emit
    -M  length    minimum length of a corrected read
  
  ERROR: no gkpStore input (-G) supplied.
  ERROR: no ovlStore input (-O) supplied.
