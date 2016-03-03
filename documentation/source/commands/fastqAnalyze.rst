fastqAnalyze
~~~~~~

::

  usage: fastqAnalyze [-stats] [-o output.fastq] input.fastq
    If no options are given, input.fastq is analyzed and a best guess for the
    QV encoding is output.  Otherwise, the QV encoding is converted to Sanger-style
    using this guess.
  
    In some cases, the encoding cannot be determined.  When this occurs, no guess is
    output.  For conversion, you can force the input QV type with:
  
    -solexa     input QV is solexa
    -illumina   input QV is illumina
    -sanger     input QV is sanger
  
    -o          sanger-style-output.fastq
  
    If -stats is supplied, no QV analysis or conversion is performed, but some simple
    statistics are computed and output to stdout.
  
