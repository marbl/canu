canu
~~~~~~

::

  -- Detected Java(TM) Runtime Environment '1.8.0_60' (from 'java').
  
  usage: canu [-correct | -trim | -assemble] \
              [-s <assembly-specifications-file>] \
               -p <assembly-prefix> \
               -d <assembly-directory> \
               genomeSize=<number>[g|m|k] \
               errorRate=0.X \
              [other-options] \
              [-pacbio-raw | -pacbio-corrected | -nanopore-raw | -nanopore-corrected] *fastq
  
    By default, all three stages (correct, trim, assemble) are computed.
    To compute only a single stage, use:
      -correct  - generate corrected reads
      -trim     - generate trimmed reads
      -assemble - generate an assembly
  
    The assembly is computed in the (created) -d <assembly-directory>, with most
    files named using the -p <assembly-prefix>.
  
    The genome size is your best guess of the genome size of what is being assembled.
    It is used mostly to compute coverage in reads.  Fractional values are allowed: '4.7m'
    is the same as '4700k' and '4700000'
  
    The errorRate is not used correctly (we're working on it).  Don't set it
    If you want to change the defaults, use the various utg*ErrorRate options.
  
    A full list of options can be printed with '-options'.  All options
    can be supplied in an optional sepc file.
  
    Reads can be either FASTA or FASTQ format, uncompressed, or compressed
    with gz, bz2 or xz.  Reads are specified by the technology they were
    generated with:
      -pacbio-raw         <files>
      -pacbio-corrected   <files>
      -nanopore-raw       <files>
      -nanopore-corrected <files>
  
  Complete documentation at http://canu.readthedocs.org/en/latest/
  
  ERROR:  Assembly name prefix not supplied with -p.
  ERROR:  Directory not supplied with -d.
  ERROR:  Required parameter 'genomeSize' is not set
  
