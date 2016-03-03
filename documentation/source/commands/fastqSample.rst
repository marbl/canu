fastqSample
~~~~~~

::

  
  usage: fastqSample [opts]
    Input Specification
      -I NAME  input name (prefix) of the reads
      -T T     total number of mate pairs in the input (if not supplied, will be counted)
      -L L     length of a single read (if not supplied, will be determined)
      -U       reads are unmated, expected in *.u.fastq
  
    Output Specification
      -O NAME  output name (prefix) of the reads (default is same as -I)
      -A       automatically include coverage or number of reads in the output name
      -m L     ignore reads shorter than L bases
      -max     don't sample randomly, pick the longest reads
  
    Method 1: specify desired output coverage:
      -g G     genome size
      -c C     desired coverage in the output reads
  
    Method 2: specify desired number of output pairs
      -p N     for mated reads, output 2N reads, or N pairs of reads
               for unmated reads, output N reads
  
    Method 3: specify a desired fraction of the input:
      -f F     output F * T pairs of reads (T as above in -t option)
               0.0 < F <= 1.0
  
    Method 4: specify a desired total length
      -b B     output reads/pairs until B bases is exceeded
  
  
  Samples reads from paired Illumina reads NAME.1.fastq and NAME.2.fastq and outputs:
      NAME.Cx.1.fastq and N.Cx.2.fastq (for coverage based sampling)
      NAME.n=N.1.fastq and N.n=N.2.fastq (for coverage based sampling)
  
  If -T is not supplied, the number of reads will be counted for you.
  
  ERROR: no name supplied with -I.
  ERROR: no method supplied with -c, -p, -f or -b
  
