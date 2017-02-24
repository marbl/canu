
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-AUG-14 to 2003-SEP-18
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-30 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAY-19 to 2014-APR-11
 *      are Copyright 2005,2007-2008,2011,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

//  Driver for the positionDB creation.  Reads a sequence.fasta, builds
//  a positionDB for the mers in the file, and then writes the internal
//  structures to disk.
//
//  The positionDB constructor is smart enough to read either a pre-built
//  image or a regular multi-fasta file.


#define MERSIZE 20


int
test1(char *filename) {
  merStream         *T       = new merStream(new kMerBuilder(MERSIZE), new seqStream(filename), true, true);
  positionDB        *M       = new positionDB(T, MERSIZE, 0, 0L, 0L, 0L, 0, 0, 0, 0, true);
  uint64            *posn    = new uint64 [1024];
  uint64             posnMax = 1024;
  uint64             posnLen = uint64ZERO;
  uint64             count   = uint64ZERO;
  uint32             missing = uint32ZERO;
  uint32             failed  = uint32ZERO;
  char               str[33];

  T->rewind();

  while (T->nextMer()) {
    if (M->getExact(T->theFMer(),
                    posn,
                    posnMax,
                    posnLen,
                    count)) {

      missing = uint32ZERO;
      for (uint32 i=0; i<posnLen; i++)
        if (posn[i] == T->thePositionInStream())
          missing++;

      if (missing != 1) {
        failed++;

        fprintf(stdout, "%s @ "uint64FMT"/"uint64FMT": Found "uint64FMT" table entries, and "uint32FMT" matching positions (",
                T->theFMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen, missing);

        for (uint32 i=0; i<posnLen; i++) {
          fprintf(stdout, uint64FMT, posn[i]);
          if (i < posnLen - 1)
            fprintf(stdout, " ");
          else
            fprintf(stdout, ")\n");
        }
      }
    } else {
      failed++;

      fprintf(stdout, "Found no matches for mer=%s at pos="uint64FMT"\n",
              T->theFMer().merToString(str), T->thePositionInStream());
    }
  }

  delete M;
  delete T;

  return(failed != 0);
}



int
test2(char *filename, char *query) {
  merStream         *T       = new merStream(new kMerBuilder(MERSIZE), new seqStream(filename), true, true);
  positionDB        *M       = new positionDB(T, MERSIZE, 0, 0L, 0L, 0L, 0, 0, 0, 0, true);
  uint64            *posn    = new uint64 [1024];
  uint64             posnMax = 1024;
  uint64             posnLen = uint64ZERO;
  uint64             count   = uint64ZERO;
  char               str[33];

  delete T;

  T = new merStream(new kMerBuilder(MERSIZE), new seqStream(query), true, true);

  while (T->nextMer()) {
    if (M->getExact(T->theFMer(),
                    posn,
                    posnMax,
                    posnLen,
                    count)) {
      fprintf(stdout, "Got a F match for mer=%s at "uint64FMT"/"uint64FMT" (in mers), numMatches="uint64FMT"\n",
              T->theFMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen);
    }

    if (M->getExact(T->theRMer(),
                    posn,
                    posnMax,
                    posnLen,
                    count)) {
      fprintf(stdout, "Got a R match for mer=%s at "uint64FMT"/"uint64FMT" (in mers), numMatches="uint64FMT"\n",
              T->theRMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen);
    }
  }

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
  uint32           mersize = 20;
  uint32           merskip = 0;

  char            *maskF = 0L;
  char            *onlyF = 0L;

  uint64           merBegin = ~uint64ZERO;
  uint64           merEnd   = ~uint64ZERO;

  char            *sequenceFile = 0L;
  char            *outputFile   = 0L;

  if (argc < 3) {
    fprintf(stderr, "usage: %s [args]\n", argv[0]);
    fprintf(stderr, "       -mersize k         The size of the mers, default=20.\n");
    fprintf(stderr, "       -merskip k         The skip between mers, default=0\n");
    fprintf(stderr, "       -use a-b,c         Specify which sequences to use, default=all\n");
    fprintf(stderr, "       -merbegin b        Build on a subset of the mers, starting at mer #b, default=all mers\n");
    fprintf(stderr, "       -merend e          Build on a subset of the mers, ending at mer #e, default=all mers\n");
    fprintf(stderr, "       -sequence s.fasta  Input sequences.\n");
    fprintf(stderr, "       -output p.posDB    Output filename.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       To dump information about an image:\n");
    fprintf(stderr, "         -dump datafile\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       To run sanity tests:\n");
    fprintf(stderr, "         -buildonly [build opts] sequence.fasta\n");
    fprintf(stderr, "           --  just builds a table and exits\n");
    fprintf(stderr, "         -existence [build opts] sequence.fasta\n");
    fprintf(stderr, "           --  builds (or reads) a table reports if any mers\n");
    fprintf(stderr, "               in sequence.fasta cannot be found\n");
    fprintf(stderr, "         -extra [build opts] sequence.fasta\n");
    fprintf(stderr, "           --  builds (or reads) a table reports if any mers\n");
    fprintf(stderr, "               NOT in sequence.fasta are be found\n");
    fprintf(stderr, "         -test1 sequence.fasta\n");
    fprintf(stderr, "           --  Tests if each and every mer is found in the\n");
    fprintf(stderr, "               positionDB.  Reports if it doesn't find a mer\n");
    fprintf(stderr, "               at the correct position.  Doesn't report if table\n");
    fprintf(stderr, "               has too much stuff.\n");
    fprintf(stderr, "         -test2 db.fasta sequence.fasta\n");
    fprintf(stderr, "           --  Builds a positionDB from db.fasta, then searches\n");
    fprintf(stderr, "               the table for each mer in sequence.fasta.  Reports\n");
    fprintf(stderr, "               all mers it finds.\n");
    fprintf(stderr, "            -- This is a silly test and you shouldn't do it.\n");
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-mersize") == 0) {
      mersize = strtouint32(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-merskip") == 0) {
      merskip = strtouint32(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-mask") == 0) {
      maskF = argv[++arg];
    } else if (strcmp(argv[arg], "-only") == 0) {
      onlyF = argv[++arg];

    } else if (strcmp(argv[arg], "-merbegin") == 0) {
      merBegin = strtouint64(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-merend") == 0) {
      merEnd = strtouint64(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-sequence") == 0) {
      sequenceFile = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      outputFile = argv[++arg];

    } else if (strcmp(argv[arg], "-dump") == 0) {
      positionDB *e = new positionDB(argv[++arg], 0, 0, 0, false);
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


  merStream *MS = new merStream(new kMerBuilder(MERSIZE),
                                new seqStream(sequenceFile),
                                true, true);

  //  Approximate the number of mers in the sequences.
  //
  uint64     numMers = MS->approximateNumberOfMers();

  //  Reset the limits.
  //
  //  XXX: If the user somehow knows how many mers are in the input
  //  file, and specifies an end between there and the amount of
  //  sequence, we'll pointlessly still make a merStreamFile, even
  //  though we shouldn't.
  //
  if (merBegin == ~uint64ZERO)   merBegin = 0;
  if (merEnd   == ~uint64ZERO)   merEnd   = numMers;

  if (merBegin >= merEnd) {
    fprintf(stderr, "ERROR: merbegin="uint64FMT" and merend="uint64FMT" are incompatible.\n",
            merBegin, merEnd);
    exit(1);
  }

  if ((merBegin > 0) || (merEnd < numMers))
    MS->setBaseRange(merBegin, merEnd);

  existDB *maskDB = 0L;
  if (maskF) {
    fprintf(stderr, "Building maskDB from '%s'\n", maskF);
    maskDB = new existDB(maskF, mersize, existDBnoFlags, 0, ~uint32ZERO);
  }

  existDB *onlyDB = 0L;
  if (onlyF) {
    fprintf(stderr, "Building onlyDB from '%s'\n", onlyF);
    onlyDB = new existDB(onlyF, mersize, existDBnoFlags, 0, ~uint32ZERO);
  }

  fprintf(stderr, "Building table with merSize "uint32FMT", merSkip "uint32FMT"\n", mersize, merskip);

  positionDB *positions = new positionDB(MS, mersize, merskip, maskDB, onlyDB, 0L, 0, 0, 0, 0, true);

  fprintf(stderr, "Dumping positions table to '%s'\n", outputFile);

  positions->saveState(outputFile);

  delete MS;
  delete positions;

  exit(0);
}
