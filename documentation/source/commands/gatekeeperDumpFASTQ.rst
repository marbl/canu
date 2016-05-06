gatekeeperDumpFASTQ
~~~~~~

::

  usage: gatekeeperDumpFASTQ [...] -o fastq-prefix -g gkpStore
    -G gkpStore
    -o fastq-prefix     write files fastq-prefix.(libname).fastq, ...
                        if fastq-prefix is '-', all sequences output to stdout
                        if fastq-prefix ends in .gz, .bz2 or .xz, output is compressed
  
    -l libToDump        output only read in library number libToDump (NOT IMPLEMENTED)
    -r id[-id]          output only the single read 'id', or the specified range of ids
  
    -c clearFile        clear range file from OBT modules
    -allreads           if a clear range file, lower case mask the deleted reads
    -allbases           if a clear range file, lower case mask the non-clear bases
    -onlydeleted        if a clear range file, only output deleted reads (the entire read)
  
    -fastq              output is FASTQ format (with extension .fastq, default)
    -fasta              output is FASTA format (with extension .fasta)
  
    -nolibname          don't include the library name in the output file name
  
  ERROR: no gkpStore (-G) supplied.
  ERROR: no output prefix (-o) supplied.
