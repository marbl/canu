utgcns
~~~~~~

::

  usage: utgcns [opts]
  
    INPUT
      -G g            Load reads from gkStore 'g'
      -T t v p        Load unitigs from tgStore 't', version 'v', partition 'p'.
                        Expects reads will be in gkStore partition 'p' as well
                        Use p='.' to specify no partition
      -t file         Test the computation of the unitig layout in 'file'
                        'file' can be from:
                          'tgStoreDump -d layout' (human readable layout format)
                          'utgcns -L'             (human readable layout format)
                          'utgcns -O'             (binary multialignment format)
  
    OUTPUT
      -O results      Write computed tigs to binary output file 'results'
      -L layouts      Write computed tigs to layout output file 'layouts'
      -F fastq        Write computed tigs to fastq  output file 'fastq'
  
    TIG SELECTION (if -T input is used)
      -u b            Compute only unitig ID 'b' (must be in the correct partition!)
      -u b-e          Compute only unitigs from ID 'b' to ID 'e'
      -f              Recompute unitigs that already have a multialignment
      -maxlength l    Do not compute consensus for unitigs longer than l bases.
  
    PARAMETERS
      -e e            Expect alignments at up to fraction e error
      -em m           Don't ever allow alignments more than fraction m error
      -l l            Expect alignments of at least l bases
      -maxcoverage c  Use non-contained reads and the longest contained reads, up to
                      C coverage, for consensus generation.  The default is 0, and will
                      use all reads.
  
    LOGGING
      -v              Show multialigns.
      -V              Enable debugging option 'verbosemultialign'.
  
  ERROR:  No gkpStore (-G) supplied.
  ERROR:  No tigStore (-T) OR no test unitig (-t) supplied.
