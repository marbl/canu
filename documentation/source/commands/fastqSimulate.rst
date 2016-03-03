fastqSimulate
~~~~~~

::

  usage: fastqSimulate -f reference.fasta -o output-prefix -l read-length ....
    -f ref.fasta    Use sequences in ref.fasta as the genome.
    -o name         Create outputs name.1.fastq and name.2.fastq (and maybe others).
    -l len          Create reads of length 'len' bases.
    -n n            Create 'n' reads (for -se) or 'n' pairs of reads (for -pe and -mp).
    -x read-cov     Set 'np' to create reads that sample the genome to 'read-cov' read coverage.
    -X clone-cov    Set 'np' to create reads that sample the genome to 'clone-cov' clone coverage.
  
    -em err         Reads will contain fraction mismatch  error 'e' (0.01 == 1% error).
    -ei err         Reads will contain fraction insertion error 'e' (0.01 == 1% error).
    -ed err         Reads will contain fraction deletion  error 'e' (0.01 == 1% error).
  
    -seed s         Seed randomness with 32-bit integer s.
  
    -allowgaps      Allow pairs to span N regions in the reference.  By default, pairs
                    are not allowed to span a gap.  Reads are never allowed to cover N's.
  
    -allowns        Allow reads to contain N regions.  Implies -allowgaps
  
    -nojunction     For -mp, do not create chimeric junction reads.  Create only fully PE or
                    fully MP reads.
  
    -normal p       Output a normal-oriented (both forward or both reverse) pair with
                    probability p.  Only for -pe and -mp.
  
    -se
                    Create single-end reads.
  
    -cc junkSize junkStdDev false
                    Create chimeric single-end reads.  The chimer is formed from two uniformly
                    distributed positions in the reference.  Some amount of random junk is inserted
                    at the junction.  With probability 'false' the read is not chimeric, but still
                    the junk bases inserted in the middle.
  
    -pe shearSize shearStdDev
                    Create paired-end reads, from fragments of size 'shearSize +- shearStdDev'.
  
    -mp insertSize insertStdDev shearSize shearStdDev enrichment
                    Create mate-pair reads.  The pairs will be 'insertSize +- insertStdDev'
                    apart.  The circularized insert is then sheared into fragments of size
                    'shearSize +- shearStdDev'.  With probability 'enrichment' the fragment
                    containing the junction is used to form the pair of reads.  The junction
                    location is uniformly distributed through this fragment.
                    Reads are labeled as:
                      tMP - a MP pair
                      fMP - a PE pair
                      aMP - a MP pair with junction in the first read
                      bMP - a MP pair with junction in the second read
                      cMP - a MP pair with junction in both reads (the reads overlap)
  
  Output QV's are the Sanger spec.
  
  ERROR:  No fasta file (-f) supplied.
  ERROR:  No output prefix (-o) supplied.
  ERROR:  No type (-se or -pe or -mp) selected.
