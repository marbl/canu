gatekeeperDumpFASTQ
~~~~~~~~~~~~~~~~~~~

::

  usage: gatekeeperDumpFASTQ [...] -o fastq-prefix -g gkpStore
    -G gkpStore
    -o fastq-prefix     write files fastq-prefix.(libname).fastq, ...
  
    -l libToDump        output only read in library number libToDump (NOT IMPLEMENTED)
    -b id               output starting at read 'id'
    -e id               output stopping after read 'id'
  
    -c clearFile        clear range file from OBT modules
    -allreads           if a clear range file, lower case mask the deleted reads
    -allbases           if a clear range file, lower case mask the non-clear bases
    -onlydeleted        if a clear range file, only output deleted reads (the entire read)
  
    -r id               output only the single read 'id'
  
    -fastq              output is FASTQ format (with extension .fastq, default)
    -fasta              output is FASTA format (with extension .fasta)
  
    -nolibname          don't include the library name in the output file name
  
  ERROR: no gkpStore (-G) supplied.
  ERROR: no output prefix (-o) supplied.
