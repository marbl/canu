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
  
      -p package      Load unitig and read from 'package' created with -P.  This
                      is usually used by developers.
  
  
    ALGORITHM
      -quick          No alignments, just paste read sequence into the unitig positions.
                      This is very fast, but the consensus sequence is formed from a mosaic
                      of read sequences, and there can be large indel.  This is useful for
                      checking intermediate assembly structure by mapping to reference, or
                      possibly for use as input to a polishing step.
      -pbdagcon       Use pbdagcon (https://github.com/PacificBiosciences/pbdagcon).
                      This is fast and robust.  It is the default algorithm.  It does not
                      generate a final multialignment output (the -v option will not show
                      anything useful).
      -utgcns         Use utgcns (the original Celera Assembler consensus algorithm)
                      This isn't as fast, isn't as robust, but does generate a final multialign
                      output.
  
  
    OUTPUT
      -O results      Write computed tigs to binary output file 'results'
      -L layouts      Write computed tigs to layout output file 'layouts'
      -A fasta        Write computed tigs to fasta  output file 'fasta'
      -Q fastq        Write computed tigs to fastq  output file 'fastq'
  
      -P package      Create a copy of the inputs needed to compute the unitigs.  This
                      file can then be sent to the developers for debugging.  The unitig(s)
                      are not processed and no other outputs are created.  Ideally,
                      only one unitig is selected (-u, below).
  
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
  
  ERROR:  No gkpStore (-G) and no package (-p) supplied.
  ERROR:  No tigStore (-T) OR no test unitig (-t) OR no package (-p)  supplied.
