ovStoreDump
~~~~~~

::

  usage: ovStoreDump -G gkpStore -O ovlStore ...
  
  There are three modes of operation:
    -d [a[-b]] dump overlaps for reads a to b, inclusive
    -q a b     report the a,b overlap, if it exists.
    -p a       dump a picture of overlaps to fragment 'a'.
  
    FORMAT (for -d)
  
    -coords    dump overlap showing coordinates in the reads (default)
    -hangs     dump overlap showing dovetail hangs unaligned
    -raw       dump overlap showing its raw native format (four hangs)
    -paf       dump overlaps in miniasm/minimap format
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
    -unique           Report only overlaps where A id is < B id, do not report both A to B and B to A overlap
  
    -best prefix      Annotate picture with status from bogart outputs prefix.edges, prefix.singletons, prefix.edges.suspicious
    -noc              With -best data, don't show overlaps to contained reads.
    -nos              With -best data, don't show overlaps to suspicious reads.
  
  ERROR: no operation (-d, -q or -p) supplied.
  ERROR: no input gkpStore (-G) supplied.
  ERROR: no input ovlStore (-O) supplied.
