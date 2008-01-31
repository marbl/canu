#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

//  Driver for the positionDB creation.  Reads a sequence.fasta, builds
//  a positionDB for the mers in the file, and then writes the internal
//  structures to disk.
//
//  The positionDB constructor is smart enough to read either a pre-built
//  image or a regular multi-fasta file.


#define MERSIZE 20

char const *usage =
"usage: %s [args]\n"
"       -mersize k         The size of the mers, default=20.\n"
"       -merskip k         The skip between mers, default=0\n"
"       -use a-b,c         Specify which sequences to use, default=all\n"
"       -merbegin b        Build on a subset of the mers, starting at mer #b, default=all mers\n"
"       -merend e          Build on a subset of the mers, ending at mer #e, default=all mers\n"
"       -sequence s.fasta  Input sequences.\n"
"       -output p.posDB    Output filename.\n"
"\n"
"       To dump information about an image:\n"
"         -dump datafile\n"
"\n"
"       To run sanity tests:\n"
"         -buildonly [build opts] sequence.fasta\n"
"           --  just builds a table and exits\n"
"         -existence [build opts] sequence.fasta\n"
"           --  builds (or reads) a table reports if any mers\n"
"               in sequence.fasta cannot be found\n"
"         -extra [build opts] sequence.fasta\n"
"           --  builds (or reads) a table reports if any mers\n"
"               NOT in sequence.fasta are be found\n"
"         -test1 sequence.fasta\n"
"           --  Tests if each and every mer is found in the\n"
"               positionDB.  Reports if it doesn't find a mer\n"
"               at the correct position.  Doesn't report if table\n"
"               has too much stuff.\n"
"         -test2 db.fasta sequence.fasta\n"
"           --  Builds a positionDB from db.fasta, then searches\n"
"               the table for each mer in sequence.fasta.  Reports\n"
"               all mers it finds.\n"
"            -- This is a silly test and you shouldn't do it.\n";



int
test1(char *filename) {
  seqStream         *C       = new seqStream(filename, true);
  kMerBuilder        KB(MERSIZE);
  merStream         *T       = new merStream(&KB, C);
  positionDB        *M       = new positionDB(T, MERSIZE, 0, 0L, 0L, 0, true);
  u64bit            *posn    = new u64bit [1024];
  u64bit             posnMax = 1024;
  u64bit             posnLen = u64bitZERO;
  u32bit             missing = u32bitZERO;
  u32bit             failed  = u32bitZERO;
  char               str[33];

  if (T->rewind() == false) {
    fprintf(stderr, "test 1 failed to rewind the merStream.\n");
    return(true);
  }

  while (T->nextMer()) {
    if (M->get(T->theFMer(),
                posn,
                posnMax,
                posnLen)) {

      missing = u32bitZERO;
      for (u32bit i=0; i<posnLen; i++)
        if (posn[i] == T->thePositionInStream())
          missing++;

      if (missing != 1) {
        failed++;

        fprintf(stdout, "%s @ "u64bitFMT"/"u64bitFMT": Found "u64bitFMT" table entries, and "u32bitFMT" matching positions (",
                T->theFMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen, missing);

        for (u32bit i=0; i<posnLen; i++) {
          fprintf(stdout, u64bitFMT, posn[i]);
          if (i < posnLen - 1)
            fprintf(stdout, " ");
          else
            fprintf(stdout, ")\n");
        }
      }
    } else {
      failed++;

      fprintf(stdout, "Found no matches for mer=%s at pos="u64bitFMT"\n",
              T->theFMer().merToString(str), T->thePositionInStream());
    }
  }

  delete M;
  delete T;
  delete C;

  return(failed != 0);
}



int
test2(char *filename, char *query) {
  seqStream         *C       = new seqStream(filename, true);
  kMerBuilder        KB(MERSIZE);
  merStream         *T       = new merStream(&KB, C);
  positionDB        *M       = new positionDB(T, MERSIZE, 0, 0L, 0L, 0, true);
  u64bit            *posn    = new u64bit [1024];
  u64bit             posnMax = 1024;
  u64bit             posnLen = u64bitZERO;
  char               str[33];

  delete T;
  delete C;

  C = new seqStream(query, true);
  T = new merStream(&KB, C);

  while (T->nextMer()) {
    if (M->get(T->theFMer(),
                posn,
                posnMax,
                posnLen)) {
      fprintf(stdout, "Got a F match for mer=%s at "u64bitFMT"/"u64bitFMT" (in mers), numMatches="u64bitFMT"\n",
              T->theFMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen);
    }

    if (M->get(T->theRMer(),
                posn,
                posnMax,
                posnLen)) {
      fprintf(stdout, "Got a R match for mer=%s at "u64bitFMT"/"u64bitFMT" (in mers), numMatches="u64bitFMT"\n",
              T->theRMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen);
    }
  }

  delete C;
  delete M;
  delete T;

  return(0);
}



//  Builds a positionDB possibly using a subset of the file.
//
//  Subset on entire sequences:
//    -use x-y,a,b  
//
//  Subset on a range of mers, in this case, use only the 1000th
//  through 1999th (inclusive) mer:
//    -merbegin 1000 -merend 2000
//
//  Or do both, use the first 1000 mers from the 3rd sequence:
//    -use 3 -merbegin 0 -merend 1000



int
main(int argc, char **argv) {
  u32bit           mersize = 20;
  u32bit           merskip = 0;

  char            *maskF = 0L;
  char            *onlyF = 0L;

  seqStream        SS;
  u64bit           merBegin = ~u64bitZERO;
  u64bit           merEnd   = ~u64bitZERO;

  merStream       *MS = 0L;
  seqStore        *seqsto = 0L;

  char            *outputFile = 0L;

  if (argc < 3) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-mersize") == 0) {
      mersize = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-merskip") == 0) {
      merskip = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-mask") == 0) {
      maskF = argv[++arg];
    } else if (strcmp(argv[arg], "-only") == 0) {
      onlyF = argv[++arg];

    } else if (strcmp(argv[arg], "-use") == 0) {
      SS.parse(argv[++arg]);
    } else if (strcmp(argv[arg], "-merbegin") == 0) {
      merBegin = strtou64bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-merend") == 0) {
      merEnd = strtou64bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-sequence") == 0) {
      SS.setFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-output") == 0) {
      outputFile = argv[++arg];

    } else if (strcmp(argv[arg], "-dump") == 0) {
      positionDB *e = new positionDB(argv[argc-1], false);
      e->printState(stdout);
      delete e;
      exit(0);
    } else if (strcmp(argv[arg], "-test1") == 0) {
      exit(test1(argv[arg+1]));
    } else if (strcmp(argv[arg], "-test2") == 0) {
      exit(test2(argv[arg+1], argv[arg+2]));
    } else {
      fprintf(stderr, "ERROR: unknown arg '%s'\n", argv[arg]);
      exit(1);
    }
    arg++;
  }

  //  Exit quickly if the output file exists.
  //
  if (fileExists(outputFile)) {
    fprintf(stderr, "Output file '%s' exists already!\n", outputFile);
    exit(0);
  }

  SS.finish();

  //  Approximate the number of mers in the sequences.
  //
  u64bit     numMers = SS.lengthOfSequences();

  //  Reset the limits.
  //
  //  XXX: If the user somehow knows how many mers are in the input
  //  file, and specifies an end between there and the amount of
  //  sequence, we'll pointlessly still make a merStreamFile, even
  //  though we shouldn't.
  //
  if (merBegin == ~u64bitZERO)   merBegin = 0;
  if (merEnd   == ~u64bitZERO)   merEnd   = SS.lengthOfSequences();

  if (merBegin >= merEnd) {
    fprintf(stderr, "ERROR: merbegin="u64bitFMT" and merend="u64bitFMT" are incompatible.\n",
            merBegin, merEnd);
    exit(1);
  }

  kMerBuilder        KB(mersize);

  //  We need a merStreamFile to implement -merbegin (using
  //  seekToMer()) and -merend (using setIterationLimit()).
  //
  if ((merBegin > 0) || (merEnd < numMers)) {
    seqsto = new seqStore(outputFile, &SS);
    MS     = new merStream(&KB, seqsto);

    if (MS->setRange(merBegin, merEnd) == false) {
      fprintf(stderr, "ERROR:  Couldn't set the range to "u64bitFMT","u64bitFMT".\n", merBegin, merEnd);
      exit(1);
    }
  } else {
    MS = new merStream(&KB, &SS);
  }

  existDB *maskDB = 0L;
  if (maskF) {
    fprintf(stderr, "Building maskDB from '%s'\n", maskF);
    maskDB = new existDB(maskF, mersize, existDBnoFlags, 0, ~u32bitZERO);
  }

  existDB *onlyDB = 0L;
  if (onlyF) {
    fprintf(stderr, "Building onlyDB from '%s'\n", onlyF);
    onlyDB = new existDB(onlyF, mersize, existDBnoFlags, 0, ~u32bitZERO);
  }

  fprintf(stderr, "Building table with merSize "u32bitFMT", merSkip "u32bitFMT"\n", mersize, merskip);

  positionDB *positions = new positionDB(MS, mersize, merskip, maskDB, onlyDB, 0, true);

  fprintf(stderr, "Dumping positions table to '%s'\n", outputFile);

  positions->saveState(outputFile);

  delete MS;
  delete seqsto;
  delete positions;

  exit(0);
}
