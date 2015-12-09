gatekeeperCreate
~~~~~~~~~~~~~~~~

::

  usage: gatekeeperCreate [...] -o gkpStore
    -o gkpStore         create this gkpStore
    
    -minlength L        discard reads shorter than L
    
    
  ERROR: no gkpStore (-g) supplied.
  ERROR: no input files supplied.


.. _gkp-files:

Gatekeeper Input File
~~~~~~~~~~~~~~~~~~~~~

The gatekeeper input file is a text file detaling the library the reads are to be placed in, and
the input read files.  Multiple read files are allowed.

A simple file is::

 !  Comments use the usual
 #  delimiters.
 #
 name    pacbio-library-1
 trim    bestedge
 /data/pacbio/reads/lib.1.fastq
 /data/pacbio/reads/lib.2.fastq
 /data/pacbio/reads/lib.3.fastq

This will create the library names 'pacbio-library-1', with option 'trim' set to 'bestedge'.  The
library will be loaded with reads from three files.

Parameters are defined in ``src/stores/gatekeeperCreate.C``:

preset
  Set all parameters for a given technology.  Defined in ``src/stores/gkStore.C``.
qv
  For FASTA format reads, the QV for every base.
isNonRandom
  This library is from a non-random construction, and should not be used for coverage computations.
trustHomopolymerRuns
  This library has good homopolymer runs.  A leftover from 454 reads.  Currently not used.
removeDuplicateReads
  Not used.
finalTrim
  Select the algorithm used for trimming reads.
removeSpurReads
  If set, remove spurs from the ends of reads.
removeChimericReads
  If set, clean up chimeric reads.
checkForSubReads
  If set, clean up PacBio RS II hairpin adapters.
