#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "bio++.H"

//  Build a merStreamFile using small mers, read it back using bigger mers.
//
//  construct a fasta sequence:
//    small sequence
//    big sequence, multiple of mersize
//    small sequence
//    big sequence
//    etc
//    small sequence
//
//  Then we reconstruct the sequences using mers.  All three merstream
//  sources are tested (the character string source is not tested).
//  The merStreamFile is tested both forwards (nextMer(), via the
//  merstream interface) and backwards (setIterationStart()).

#define BUILD_SIZE        88
#define TEST_SIZE        403
#define MERS_PER_SEQ      37
#define TEST_ITERATIONS  300

#define MSF_FILENAME    "junk.bigmer"
#define FASTA_FILENAME  "junk.bigmer.fasta"

//  construct a multi-fasta sequence, alternating short and long
//  sequences.  Short sequences are less than TEST_SIZE long (most,
//  not all, longer than BUILD_SIZE).  Long sequences are exactly
//  TEST_SIZE * MERS_PER_SEQ long -- this lets us dump the mers to
//  reconstruct the sequence.
//
void
buildFastA(void) {
  mt_s  *mtctx  = mtInit(time(0L));
  char  *seq    = new char [TEST_SIZE * MERS_PER_SEQ + 1];
  char   dna[4] = { 'A', 'C', 'G', 'T' };

  FILE *F = fopen(FASTA_FILENAME, "w");

  for (u32bit i=0; i<TEST_ITERATIONS; i++) {
    u32bit len;

    fprintf(F, ">"u32bitFMT"short\n", i);
    len = mtRandom32(mtctx) % (TEST_SIZE-1) + 1;
    for (u32bit s=0; s<len; s++)
      seq[s] = dna[ mtRandom32(mtctx) % 4 ];
    seq[len] = 0;
    fprintf(F, "%s\n", seq);

    fprintf(F, ">"u32bitFMT"long\n", i);
    len = TEST_SIZE * MERS_PER_SEQ;
    for (u32bit s=0; s<len; s++)
      seq[s] = dna[ mtRandom32(mtctx) % 4 ];
    seq[len] = 0;
    fprintf(F, "%s\n", seq);
  }

  fclose(F);
}


//  Uses the merStreamFile directly to read mers, chains mers
//  into a sequence, compares against the correct sequence.
//
void
test1(u32bit style) {
  seqCache   *fasta = new seqCache(FASTA_FILENAME);
  seqInCore  *sseq  = fasta->getSequenceInCore();
  seqInCore  *lseq  = fasta->getSequenceInCore();

  char                  mseq[TEST_SIZE * MERS_PER_SEQ + 1];

  //  Construct a reader, and load the first mer.

  merStream           *MS = 0L;
  merStreamFileReader *RD = 0L;
  chainedSequence     *CS = 0L;

  switch (style) {
    case 0:
      fprintf(stderr, "test1(0)-- Testing merStreamFileReader -> merStream\n");
      RD = new merStreamFileReader(MSF_FILENAME, TEST_SIZE);
      MS = new merStream(RD);
      break;
    case 1:
      fprintf(stderr, "test1(2)-- Testing chainedSequence -> merStream\n");
      CS = new chainedSequence();
      CS->setSource(FASTA_FILENAME);
      CS->finish();
      MS = new merStream(TEST_SIZE, CS);
      break;
    case 2:
      fprintf(stderr, "test1(3)-- Testing merStreamFileReader (backwards)\n");
      RD = new merStreamFileReader(MSF_FILENAME, TEST_SIZE);
      break;
    default:
      break;
  }


  for (u32bit s=0; fasta->eof() == false; s++) {
    for (u32bit i=0; i<TEST_SIZE * MERS_PER_SEQ + 1; i++)
      mseq[i] = 0;

    switch (style) {
      case 0:
      case 1:
        //  Fill the sequence using non-overlapping mers, skipping
        //  intermediate mers (there aren't intermediate mers if we're the
        //  last mer in the sequence!)
        //
        MS->nextMer();
        for (u32bit i=0; i<MERS_PER_SEQ; i++) {
          MS->theFMer().merToString(mseq + i * TEST_SIZE);
          if (i != MERS_PER_SEQ-1)
            MS->nextMer(TEST_SIZE - 1);
        }
        break;
      case 2:
        //  Same thing, but read the mers backwards -- we could read
        //  the sequences backwards, too, but that doesn't gain us
        //  anything (we still seek to every location).
        //
        for (u32bit i=MERS_PER_SEQ; i--; ) {
          char  copy[TEST_SIZE + 1];
          RD->setIterationStart(s * (MERS_PER_SEQ * TEST_SIZE - TEST_SIZE + 1) + i * (TEST_SIZE));
          RD->nextMer();
          RD->theFMer().merToString(copy);
          strncpy(mseq + i * TEST_SIZE, copy, TEST_SIZE);

          //  Aww, what the hell!  Test reverse complement stuff too!
          //
          kMer  f = RD->theFMer();
          kMer  r = RD->theRMer();
          f.reverseComplement();

          if (f != r) {
            char str[1025];
            fprintf(stderr, "Reverse Complement mismatch:\n");
            fprintf(stderr, "  reversed fwd = '%s'\n", f.merToString(str));
            fprintf(stderr, "           rev = '%s'\n", r.merToString(str));
            exit(1);
          }

          f = RD->theFMer();
          r = RD->theRMer();
          r.reverseComplement();

          if (f != r) {
            char str[1025];
            fprintf(stderr, "Reverse Complement mismatch:\n");
            fprintf(stderr, "           fwd = '%s'\n", f.merToString(str));
            fprintf(stderr, "  reversed rev = '%s'\n", r.merToString(str));
            exit(1);
          }


        }
        mseq[MERS_PER_SEQ * TEST_SIZE] = 0;
        break;
      default:
        break;
    }

    //  Compare our mer-constructed sequence to the long sequence in
    //  the file
    //
    if (strcmp(mseq, lseq->sequence()) != 0) {
      fprintf(stderr, "FAIL:  seq="u32bitFMT"\nmseq=%s\nlseq=%s\n", s, mseq, lseq->sequence());
      exit(1);
    }

    delete sseq;
    delete lseq;

    sseq = fasta->getSequenceInCore();
    lseq = fasta->getSequenceInCore();
  }

  delete sseq;
  delete lseq;

  delete CS;
  delete RD;
  delete MS;

  fprintf(stderr, "  OK!\n");
}




int
main(int argc, char **argv) {

  //  Minimum KMER_WORDS is 13 -- mersizes up to 416 bases
  if (KMER_WORDS < 13) {
    fprintf(stderr, "I need at least KMER_WORDS == 13; test not run.\n");
    exit(0);
  }

  buildFastA();

  merStreamFileBuilder   *B = new merStreamFileBuilder(BUILD_SIZE, FASTA_FILENAME, MSF_FILENAME);
  B->build(true);
  delete B;

  test1(0);
  test1(1);
  test1(2);

  unlink(FASTA_FILENAME);
  unlink(FASTA_FILENAME "idx");
  unlink(MSF_FILENAME ".merStream");

  exit(0);
}
