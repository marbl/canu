canu
~~~~

::

  -- Detected 12 CPUs and 64 gigabytes of memory.
  -- Detected Sun Grid Engine in '/usr/local/sge/default'.
  
  usage: canu [-correct | -trim | -assemble] \
              [-s <assembly-specifications-file>] \
               -p <assembly-prefix> \
               -d <assembly-directory> \
	       genomeSize~Ng \
               [other-options] \
               [read-type *fastq]
  
By default, all three stages (correct, trim, assemble) are computed. To compute only a single stage, use

-correct
	 generate corrected reads

-trim
	generate trimmed reads

-assemble
	generate an assembly
  
A full list of options can be printed with '-options'.  All options can be supplied in an optional sepc file.
  
Reads can be either FASTA or FASTQ format, uncompressed, or compressed with gz, bz2 or xz::

	-pacbio-raw         
	-pacbio-corrected   
	-nanopore-raw       
	-nanopore-corrected 
  
Complete documentation at http://canu.github.io/
