#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"
//#include "fasta.H"
#include "sim4.H"

using namespace std;
#include <vector>

//
// Current ESTmapper generated summary is:
//
// GOOD: >= 95% identity, >= 80% composite, >= 0 bp
// cDNA-genomic matches  28715039 matches (24921387 different cDNA and 81 genomic)
// Matches per cDNA      1.1522 matches/cDNA
// Matches per genomic   354506.6543 matches/genomic
//
// cDNA COUNTS:
// cDNA:            27440540
// cDNA-good:       24921387 ( 90.8196%)
// cDNA-missing:       26071 (  0.0950%)
// cDNA-zero:        2493082 (  9.0854%)
//
//
// New format / summary should be
//
//
// X% identity     coverage:  50  55  60  65  70  75  80  85  90  95 100
// sequence-genomic matches  %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u
// Unique sequences          %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u
// Matches per sequence      %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u
// Unique genomic            %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u
// Matches per genomic       %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u
//
// usage:
//
// Report raw numbers at 90, 95, 99 percent identity, 50, 60, 70, 80, 90,
// 100 percent coverage:
//   summarizePolishes -i 90 95 99 -c 50 60 70 80 90 100 -p polishes-file
//
// Report percentages at same
//   summarizePolishes -i 90 95 99 -c 50 60 70 80 90 100 -nf cdna-file -p polishes-file
//   summarizePolishes -i 90 95 99 -c 50 60 70 80 90 100 -n  num-seqs  -p polishes-file
//
// Read from stdin, default to 95 percent identity, 50 percent coverage:
//   summarizePolishes -p -
//
// Buckets?  Cumulative?  Both?  If we do buckets with size 1, we'll
// use lots of space, but be fast.  Compute correct bucket sizes on
// output.  101*101 entries, 6 million sequences -> 190GB.
//
// So keep sorted list of values, find first bucket that is <= the
// match we have.  792MB for the example below (3 %i, 11 %c, 6 million
// seqs).
//
// Read in all %i,%c.  Compute each identity x coverage pair
// separately.  48MB for scores + 24MB for a pair.  Memory efficient,
// maybe not compute efficient.
//

struct match {
  u32bit   _estid;
  u32bit   _genid;
  u32bit   _identity;
  u32bit   _coverage;
};



void
readMatches(char           *filename,
            vector<match>  &matches) {
  FILE     *F = stdin;

  if ((filename != 0L) && (strcmp(filename, "-") != 0)) {
    errno = 0;
    F = fopen(filename, "r");
    if (errno) {
      fprintf(stderr, "ERROR: Can't open '%s'.\n%s\n", filename, strerror(errno));
      exit(1);
    }
    fprintf(stderr, "Reading matches from '%s'\n", filename);
  } else {
    fprintf(stderr, "Reading matches from 'stdin'\n");
  }

  matches.clear();

  sim4polish *p;

  while ((p = s4p_readPolish(F)) != 0L) {
    match  m;

    m._estid    = p->estID;
    m._genid    = p->genID;
    m._identity = p->percentIdentity;
    m._coverage = p->querySeqIdentity;

    matches.push_back(m);

    s4p_destroyPolish(p);
  }

  if ((filename != 0L) && (strcmp(filename, "-") == 0)) {
    fclose(F);
  }

  fprintf(stderr, "read %d matches.\n", matches.size());
}



int
main(int argc, char **argv) {
  char      *polishesFile = 0L;
  u32bit     numSeqs      = 0;
  char      *sequenceFile = 0L;;
  u32bit     idLen        = 0;
  u32bit     id[101]      = { 0 };
  u32bit     cvLen        = 0;
  u32bit     cv[101]      = { 0 };
  bool       formatExcel  = false;

  if (argc == 1) {
    fprintf(stderr, "usage: %s [-excel] [-p polishes-file] [-n num-seqs | -nf seq-file] [-i val ...] [-c val ...]\n", argv[0]);
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-polishes", 2) == 0) {
      polishesFile = argv[++arg];
    } else if (strncmp(argv[arg], "-n", 3) == 0) {
      numSeqs = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-nf", 3) == 0) {
      sequenceFile = argv[++arg];
    } else if (strncmp(argv[arg], "-identity", 2) == 0) {
      arg++;
      while ((argv[arg]) && (argv[arg][0] != '-'))
        id[idLen++] = atoi(argv[arg++]);
      arg--;
    } else if (strncmp(argv[arg], "-coverage", 2) == 0) {
      arg++;
      while ((argv[arg]) && (argv[arg][0] != '-'))
        cv[cvLen++] = atoi(argv[arg++]);
      arg--;
    } else if (strncmp(argv[arg], "-excel", 2) == 0) {
      formatExcel=true;
    }
    arg++;
  }

  if (polishesFile == 0L) {
    fprintf(stderr, "ERROR:  No polishes file specified!\n");
    exit(1);
  }

  if (idLen == 0) {
    fprintf(stderr, "WARNING:  Defaulting to 95%% identity.\n");
    id[idLen++] = 95;
  }

  if (cvLen == 0) {
    fprintf(stderr, "WARNING:  Defaulting to 50%% coverage.\n");
    cv[cvLen++] = 50;
  }


  fprintf(stderr, "Polishes:      %s\n", polishesFile);
  fprintf(stderr, "numSeqs:       "u32bitFMT"\n", numSeqs);
  fprintf(stderr, "sequenceFile:  %s\n", sequenceFile);
  fprintf(stderr, "ids:           "u32bitFMT" -- ", idLen);
  for (u32bit i=0; i<idLen; i++)
    fprintf(stderr, " "u32bitFMT"", id[i]);
  fprintf(stderr, "\n");
  fprintf(stderr, "cvs:           "u32bitFMT" -- ", cvLen);
  for (u32bit i=0; i<cvLen; i++)
    fprintf(stderr, " "u32bitFMT"", cv[i]);
  fprintf(stderr, "\n");

  vector<match>  matches;

  readMatches(polishesFile, matches);

  //  Find the largest cDNA and genomic idx
  //
  u32bit  estmax = 0;
  u32bit  genmax = 0;
  for (u32bit i=0; i<matches.size(); i++) {
    if (estmax < matches[i]._estid)
      estmax = matches[i]._estid;
    if (genmax < matches[i]._genid)
      genmax = matches[i]._genid;
  }

  estmax++;
  genmax++;

  //  Allocate space for statistics
  //
  u32bit  *estcounts = new u32bit [estmax];
  u32bit  *gencounts = new u32bit [genmax];
  u32bit   mapped;
  u32bit   notmapped;
  u32bit   uniqest;
  u32bit   uniqgen;

  if (formatExcel) {
    fprintf(stdout, "identity\tcoverage\tmapped\tnotmapped\tuniqest\tuniqgen\n");
    fflush(stdout);
  }

  //  Foreach identity and each coverage, find how many things
  //  are above that level.
  //
  for (u32bit i=0; i<idLen; i++) {
    for (u32bit c=0; c<cvLen; c++) {
      mapped    = 0;
      notmapped = 0;
      for (u32bit z=0; z<estmax; z++)
        estcounts[z] = 0;
      for (u32bit z=0; z<genmax; z++)
        gencounts[z] = 0;

      for (u32bit z=0; z<matches.size(); z++) {
        if ((id[i] <= matches[z]._identity) &&
            (cv[c] <= matches[z]._coverage)) {
          mapped++;
          estcounts[ matches[z]._estid ]++;
          gencounts[ matches[z]._genid ]++;
        } else {
          notmapped++;
        }
      }

      uniqest = 0;
      uniqgen = 0;

      for (u32bit z=0; z<estmax; z++)
        if (estcounts[z])
          uniqest++;
      for (u32bit z=0; z<genmax; z++)
        if (gencounts[z])
          uniqgen++;

      if (formatExcel) {
        fprintf(stdout, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n", id[i], cv[c], mapped, notmapped, uniqest, uniqgen);
        fflush(stdout);
      } else {
        fprintf(stdout, u32bitFMTW(3)" "u32bitFMTW(3)": mapped="u32bitFMTW(8)" notmapped="u32bitFMTW(8)"  est="u32bitFMTW(8)" gen="u32bitFMTW(8)"\n", id[i], cv[c], mapped, notmapped, uniqest, uniqgen);
        fflush(stdout);
      }
    }
  }
}

