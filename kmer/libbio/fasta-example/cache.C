#include <stdio.h>
#include <stdlib.h>

#include "fasta.H"
#include "fasta-cache.H"

//  Quick-n-dirty example of using the quick-n-dirty fasta cache
//
//  cxx -I../libbri -I../fasta -o cache cache.C ../fasta/libfasta.a ../libbri/libbri.a


int
main(int argc, char **argv) {

  //  Build a fasta-cache of whatever file is on the command line
  //    The 32 is a useless parameter here...
  //    'true' says to really load all the sequences at the start
  //
  FastACache  cache(argv[1], 32, true);

  //  Grab the second sequence from the cache
  //
  FastASequenceInCore  *seq = cache.getSequence(1);

  //  And do something with the sequence
  //
  fprintf(stdout, "%s\n%s\n", seq->header(), seq->sequence());

  //  Don't delete seq, it's owned by the cache!
}
