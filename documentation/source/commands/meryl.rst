meryl
~~~~~~

::

  usage: meryl [personality] [global options] [options]
  
  where personality is:
          -P -- compute parameters
          -B -- build table
          -S -- scan table
          -M -- "math" operations
          -D -- dump table
  
  -P:  Given a sequence file (-s) or an upper limit on the
       number of mers in the file (-n), compute the table size
       (-t in build) to minimize the memory usage.
          -m #          (size of a mer; required)
          -c #          (homopolymer compression; optional)
          -p            (enable positions)
          -s seq.fasta  (seq.fasta is scanned to determine the number of mers)
          -n #          (compute params assuming file with this many mers in it)
  
       Only one of -s, -n need to be specified.  If both are given
       -s takes priority.
  
  
  -B:  Given a sequence file (-s) and lots of parameters, compute
       the mer-count tables.  By default, both strands are processed.
          -f            (only build for the forward strand)
          -r            (only build for the reverse strand)
          -C            (use canonical mers, assumes both strands)
          -L #          (DON'T save mers that occur less than # times)
          -U #          (DON'T save mers that occur more than # times)
          -m #          (size of a mer; required)
          -c #          (homopolymer compression; optional)
          -p            (enable positions)
          -s seq.fasta  (sequence to build the table for)
          -o tblprefix  (output table prefix)
          -v            (entertain the user)
  
       By default, the computation is done as one large sequential process.
       Multi-threaded operation is possible, at additional memory expense, as
       is segmented operation, at additional I/O expense.
  
       Threaded operation: Split the counting in to n almost-equally sized
       pieces.  This uses an extra h MB (from -P) per thread.
          -threads n    (use n threads to build)
  
       Segmented, sequential operation: Split the counting into pieces that
       will fit into no more than m MB of memory, or into n equal sized pieces.
       Each piece is computed sequentially, and the results are merged at the end.
       Only one of -memory and -segments is needed.
          -memory mMB   (use at most m MB of memory per segment)
          -segments n   (use n segments)
  
       Segmented, batched operation: Same as sequential, except this allows
       each segment to be manually executed in parallel.
       Only one of -memory and -segments is needed.
          -memory mMB     (use at most m MB of memory per segment)
          -segments n     (use n segments)
          -configbatch    (create the batches)
          -countbatch n   (run batch number n)
          -mergebatch     (merge the batches)
       Initialize the compute with -configbatch, which needs all the build options.
       Execute all -countbatch jobs, then -mergebatch to complete.
         meryl -configbatch -B [options] -o file
         meryl -countbatch 0 -o file
         meryl -countbatch 1 -o file
         ...
         meryl -countbatch N -o file
         meryl -mergebatch N -o file
       Batched mode can run on the grid.
          -sge        jobname      unique job name for this execution.  Meryl will submit
                                   jobs with name mpjobname, ncjobname, nmjobname, for
                                   phases prepare, count and merge.
          -sgebuild "options"    any additional options to sge, e.g.,
          -sgemerge "options"    "-p -153 -pe thread 2 -A merylaccount"
                                   N.B. - -N will be ignored
                                   N.B. - be sure to quote the options
  
  -M:  Given a list of tables, perform a math, logical or threshold operation.
       Unless specified, all operations take any number of databases.
  
       Math operations are:
          min       count is the minimum count for all databases.  If the mer
                    does NOT exist in all databases, the mer has a zero count, and
                    is NOT in the output.
          minexist  count is the minimum count for all databases that contain the mer
          max       count is the maximum count for all databases
          add       count is sum of the counts for all databases
          sub       count is the first minus the second (binary only)
          abs       count is the absolute value of the first minus the second (binary only)
  
       Logical operations are:
          and       outputs mer iff it exists in all databases
          nand      outputs mer iff it exists in at least one, but not all, databases
          or        outputs mer iff it exists in at least one database
          xor       outputs mer iff it exists in an odd number of databases
  
       Threshold operations are:
          lessthan x            outputs mer iff it has count <  x
          lessthanorequal x     outputs mer iff it has count <= x
          greaterthan x         outputs mer iff it has count >  x
          greaterthanorequal x  outputs mer iff it has count >= x
          equal x               outputs mer iff it has count == x
       Threshold operations work on exactly one database.
  
          -s tblprefix  (use tblprefix as a database)
          -o tblprefix  (create this output)
          -v            (entertain the user)
  
       NOTE:  Multiple tables are specified with multiple -s switches; e.g.:
                meryl -M add -s 1 -s 2 -s 3 -s 4 -o all
       NOTE:  It is NOT possible to specify more than one operation:
                meryl -M add -s 1 -s 2 -sub -s 3
              will NOT work.
  
  
  -D:  Dump the table (not all of these work).
  
       -Dd        Dump a histogram of the distance between the same mers.
       -Dt        Dump mers >= a threshold.  Use -n to specify the threshold.
       -Dc        Count the number of mers, distinct mers and unique mers.
       -Dh        Dump (to stdout) a histogram of mer counts.
       -s         Read the count table from here (leave off the .mcdat or .mcidx).
  
  
