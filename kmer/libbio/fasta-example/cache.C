#include <stdio.h>
#include <stdlib.h>

#include "fasta.H"
#include "fasta-cache.H"

//  Quick-n-dirty example of using the quick-n-dirty fasta cache
//
//  cxx -I../libbri -I../fasta -o cache cache.C ../fasta/libfasta.a ../libbri/libbri.a


int
main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "usage: %s idx-to-get fasta-file\n", argv[0]);
    exit(1);
  }

  //  Build a fasta-cache of whatever file is on the command line
  //    The 32 is a useless parameter here...
  //    'true' says to really load all the sequences at the start
  //
  FastACache  cache(argv[2], 32, true);

  //  Grab the second sequence from the cache
  //
  FastASequenceInCore  *seq = cache.getSequence(atoi(argv[1]));

  //  And do something with the sequence
  //
  fprintf(stdout, "%s\n%s\n", seq->header(), seq->sequence());

  //  Don't delete seq, it's owned by the cache!
}
