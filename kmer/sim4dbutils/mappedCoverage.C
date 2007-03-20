#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "bio.h"
#include "sim4.H"

//  Reports the amount of sequence covered by ALL matches for that
//  sequence.  Example: if sequence iid 4 has two matches, one
//  covering the first 30% and the second covering the last 30%, this
//  will report that sequence iid 4 is covered 60%.
//
//  Takes no options, reads from stdin, writes to stdout.

void
usage(char *name) {
  fprintf(stderr, "usage: %s [-mask in.fasta] [-cov dat] [-raw | -blast] < sim4db-results\n", name);
  fprintf(stderr, "       -mask    Read sequences from in.fasta, lower-case mask\n");
  fprintf(stderr, "                any base with an alignment, write to out.fasta\n");
  fprintf(stderr, "       -cov     Write coverage statistics to 'dat' instead of stdout\n");
  fprintf(stderr, "       -raw     If present, assume the 'sim4db-results' are\n");
  fprintf(stderr, "                a space-separated list of 'iid begin end', one per line\n");
  fprintf(stderr, "       -blast   Same idea as raw, expects 'UID.IID' for query id,\n");
  fprintf(stderr, "                blast format (-m) 9.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Output on stdout is the masked sequence if -mask is specified,\n");
  fprintf(stderr, "otherwise, it is the coverage statistics.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "-mask is almost a required option - we need it to get the length.\n");
  fprintf(stderr, "of sequences with no mapping (100%% uncovered) and to get the\n");
  fprintf(stderr, "number of sequences.\n");
  fprintf(stderr, "\n");

  if (isatty(fileno(stdin)))
    fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");
}


int
main(int argc, char **argv) {
  u32bit                covMax    = 0;
  intervalList        **cov       = 0L;
  u32bit               *len       = 0L;

  u32bit                lastIID   = 0;

  bool                  isRaw     = false;
  bool                  isBlast   = false;

  char                 *fastaname = 0L;
  char                 *covname   = 0L;

  FastABase            *F = 0L;
  FastASequenceInCore  *S = 0L;

  FILE                 *C = stdout;

  if (isatty(fileno(stdin)))
    usage(argv[0]), exit(1);

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-mask") == 0) {
      fastaname = argv[++arg];
    } else if (strcmp(argv[arg], "-cov") == 0) {
      covname   = argv[++arg];
    } else if (strcmp(argv[arg], "-raw") == 0) {
      isRaw = true;
    } else if (strcmp(argv[arg], "-blast") == 0) {
      isBlast = true;
    } else {
      fprintf(stderr, "unknown arg: '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (fastaname) {
    C = 0L;
    F = new FastAFile(fastaname);
    F->openIndex();
  }

  if (covname) {
    errno = 0;
    C     = fopen(covname, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for write: %s\n", covname, strerror(errno)), exit(1);
  }

  covMax   = 1024 * 1024;
  if (F)
    covMax = F->getNumberOfSequences();
  cov      = new intervalList * [covMax];
  len      = new u32bit [covMax];

  fprintf(stderr, "Found "u32bitFMT" sequences in the input file.\n", covMax);

  for (u32bit i=0; i<covMax; i++) {
    cov[i] = 0L;
    len[i] = 0;
  }

  if (isRaw || isBlast) {
    char          inLine[1024];
    splitToWords  S;

    while (!feof(stdin)) {
      fgets(inLine, 1024, stdin);
      S.split(inLine);

      u32bit  iid=0, beg=0, end=0;

      if (isRaw) {
        iid = strtou32bit(S[0], 0L);
        beg = strtou32bit(S[1], 0L) - 1;   //  Convert to space-based
        end = strtou32bit(S[2], 0L);
      }
      if (isBlast) {
        char *iii = S[0];
        while ((*iii != '.') && (*iii))
          iii++;
        iii++;
        if (*iii == 0)
          fprintf(stderr, "UID.IID error: '%s'\n", S[0]);

        iid = strtou32bit(iii, 0L);
        beg = strtou32bit(S[6], 0L) - 1;   //  Convert to space-based
        end = strtou32bit(S[7], 0L);
      }

      if (iid >= covMax) {
        fprintf(stderr, "ERROR:  Found iid "u32bitFMT", but only allocated "u32bitFMT" places!\n",
                iid, covMax);
        exit(1);
      }
      if (cov[iid] == 0L) {
        cov[iid] = new intervalList;
        len[iid] = 0;
      }
      if (iid >= lastIID) {
        lastIID = iid + 1;
      }
      cov[iid]->add(beg, end-beg);
    }
  } else {
    sim4polish  *p = 0L;

    while ((p = s4p_readPolish(stdin)) != 0L) {
      if (p->estID > covMax)
        fprintf(stderr, "DIE!  You have more sequences in your polishes than in your source!\n"), exit(1);

      if (p->estID >= covMax) {
        fprintf(stderr, "ERROR:  Found iid "u32bitFMT", but only allocated "u32bitFMT" places!\n",
                p->estID, covMax);
        exit(1);
      }
      if (cov[p->estID] == 0L) {
        cov[p->estID] = new intervalList;
        len[p->estID] = p->estLen;
      }
      if (p->estID >= lastIID) {
        lastIID = p->estID + 1;
      }

      for (u32bit e=0; e<p->numExons; e++) {
        p->exons[e].estFrom--;        //  Convert to space-based

        if (p->matchOrientation == SIM4_MATCH_FORWARD)
          cov[p->estID]->add(p->exons[e].estFrom,
                             p->exons[e].estTo - p->exons[e].estFrom);
        else
          cov[p->estID]->add(p->estLen - p->exons[e].estTo,
                             p->exons[e].estTo - p->exons[e].estFrom);
      }

      s4p_destroyPolish(p);
    }
  }


  //  Scan the list of intervalLists, compute the amount covered, print.
  //
  for (u32bit iid=0; iid<lastIID; iid++) {

    //  Argh!  If there are no intervals, we need to report the whole
    //  sequence is uncovered!

    u32bit  numRegions  = 0;
    u32bit  sumLengths  = 0;
    u32bit  l, h;

    //  Save the number of regions and the sum of their lengths,
    //  then merge regions
    //
    if (cov[iid]) {
      numRegions = cov[iid]->numberOfIntervals();
      sumLengths = cov[iid]->sumOfLengths();
      cov[iid]->merge();
    }

    if (F) {
      F->find(iid);
      S = F->getSequence();

      if (len[iid] == 0)
        len[iid] = S->sequenceLength();

      assert(len[iid] == S->sequenceLength());

      char   *seq = new char [len[iid] + 1];
      strcpy(seq, S->sequence());

      for (u32bit p=0; p<len[iid]; p++)
        seq[p] = toUpper[seq[p]];

      if (cov[iid]) {
        for (u32bit c=0; c<cov[iid]->numberOfIntervals(); c++) {
          l = cov[iid]->lo(c);
          h = cov[iid]->hi(c);

          if (h > len[iid]) {
            fprintf(stderr, "ERROR:  range "u32bitFMT"-"u32bitFMT" out of bounds (seqLen = "u32bitFMT")\n",
                    l, h, len[iid]);
            assert(h <= len[iid]);
          }

          for (u32bit p=l; p<h; p++)
            //seq[p] = toLower[seq[p]];
            seq[p] = 'N';
        }
      }

      fprintf(stdout, "%s\n%s\n", S->header(), seq);

      delete [] seq;
      delete    S;
    }

    if (C) {
      double  percentCovered = 0.00;

      if (cov[iid])
        percentCovered = cov[iid]->sumOfLengths() / (double)len[iid];

      fprintf(C, u32bitFMT"\t"u32bitFMT"\t%5.3f\t"u32bitFMT"\t"u32bitFMT"\n",
              iid,
              len[iid],
              percentCovered,
              numRegions,
              sumLengths);
    }
  }
}

