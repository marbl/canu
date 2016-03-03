overlapImport
~~~~~~

::

  usage: overlapImport [options] ascii-ovl-file-input.[.gz]
  
  Required:
    -G name.gkpStore   path to valid gatekeeper store
  
  Output options:
    -o file.ovb        output file name
    -O name.ovlStore   output overlap store
  Format options:
    -legacy            'CA8 overlapStore -d' format
    -coords            'overlapConvert -coords' format (not implemented)
    -hangs             'overlapConvert -hangs' format (not implemented)
    -raw               'overlapConvert -raw' format
  
  Input file can be stdin ('-') or a gz/bz2/xz compressed file.
  
  ERROR: need to supply a gkpStore (-G).
  ERROR: need to supply a format type (-legacy, -coords, -hangs, -raw).
