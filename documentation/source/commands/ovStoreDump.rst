ovStoreDump
~~~~~~~~~~~

::

  usage: ovStoreDump -G gkpStore -O ovlStore [-b bgnID] [-e endID] ...
  
  There are three modes of operation:
    -d         dump a store (range selected with -b and -e)
    -q a b     report the a,b overlap, if it exists.
    -p a       dump a picture of overlaps to fragment 'a'.
  
    FORMAT (for -d)
  
    -coords    dump overlap showing coordinates in the reads (default)
    -hangs     dump overlap showing dovetail hangs unaligned
    -raw       dump overlap showing its raw native format (four hangs)
    -binary    dump overlap as raw binary data
    -counts    dump the number of overlaps per read
  
    MODIFIERS (for -d and -p)
  
    -E erate          Dump only overlaps <= erate fraction error.
    -L length         Dump only overlaps that are larger than L bases (only for -p picture mode).
    -d5               Dump only overlaps off the 5' end of the A frag.
    -d3               Dump only overlaps off the 3' end of the A frag.
    -dC               Dump only overlaps that are contained in the A frag (B contained in A).
    -dc               Dump only overlaps that are containing the A frag (A contained in B).
    -v                Report statistics (to stderr) on some dumps (-d).
  
