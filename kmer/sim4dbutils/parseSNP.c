#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "sim4reader.h"

//
//  Writes things with mappings that don't contain the snp itself to a
//  failure file.  Otherwise, if the mapping is above the threshold, a
//  line describing the snp is output.
//


//  Define this if the input SNPs are directly from dbSNP -- the
//  defline is different, and the SNP position is base-based.
//
//#define SNPS_ARE_RS


char *usage =
"usage: %s <options> < polishes > unprocessed-polishes\n"
"\n"
"  -i I      Require that a mapping have at least I % identity overall (not just the exon with the SNP)\n"
"  -c C      Require that a mapping have at least C % coverage\n"
"  -F file   Print failed mappings to 'file'\n"
"            A failed mapping is one that doesn't contain the SNP itself\n"
"  -O file   Save parsed SNP mappings to 'file'\n"
"\n"
"  Mappings that are below the quality thresholds are printed to stdout.\n"
"\n";



FILE *multiMultiFile   = 0L;  //  multiple hits, at least one is multiple exon
FILE *multiSingleFile  = 0L;  //  multiple hits, all are single exon
FILE *singleMultiFile  = 0L;  //  single hit, and it has more than one exon
FILE *singleSingleFile = 0L;  //  single hit, single exon

int   sm = 0;
int   ss = 0;
int   mm = 0;
int   ms = 0;

int   failedsnps    = 0;
int   failedmatches = 0;

FILE *validSNPMap   = 0L;
FILE *failedSNPMap  = 0L;
FILE *lowQualityMap = 0L;

char *
findSNPid(char *defline) {
  char *ret = 0L;
  int   sta = 0;
  int   len = 0;
  int   i = 0;

#if 0
  for (len=1; defline[len] && !isspace(defline[len]); len++)
    ;
#endif

  for (len=1; defline[len] && defline[len] != '_'; len++)
    ;
  for (sta=len-1; sta > 0 && defline[sta] != '|'; sta--)
    ;

  errno = 0;
  ret = (char *)malloc(sizeof(char) * (len+1));
  if (errno) {
    fprintf(stderr, "malloc() problem: %s\n", strerror(errno));
    exit(1);
  }

  for (i=sta; i<len-1; i++)
    ret[i-sta] = defline[i+1];

  ret[len-sta-1] = 0;

  return(ret);
}

char *
findGENid(char *defline) {
  char *ret = 0L;
  int   len = 0;
  int   i = 0;

  for (len=1; defline[len] && !isspace(defline[len]); len++)
    ;

  errno = 0;
  ret = (char *)malloc(sizeof(char) * (len+1));
  if (errno) {
    fprintf(stderr, "malloc() problem: %s\n", strerror(errno));
    exit(1);
  }

  for (i=0; i<len-1; i++)
    ret[i] = defline[i+1];

  ret[len-1] = 0;

  return(ret);
}



int
findPosition(char *defline) {
  int   i=0;

#ifdef SNPS_ARE_RS
  while ((defline[i+4] != 0) && ((defline[i]   != 'E') ||
                                 (defline[i+1] != 'P') ||
                                 (defline[i+2] != 'O') ||
                                 (defline[i+3] != 'S') ||
                                 (defline[i+4] != '=')))
    i++;
#else
  while ((defline[i+4] != 0) && ((defline[i]   != '/') ||
                                 (defline[i+1] != 'p') ||
                                 (defline[i+2] != 'o') ||
                                 (defline[i+3] != 's') ||
                                 (defline[i+4] != '=')))
    i++;
#endif

  if (defline[i] == 0) {
    fprintf(stderr, "pos not found in defline '%s'!\n", defline);
    exit(1);
  }

#ifdef SNPS_ARE_RS
  return(atoi(defline+i+5) - 1);
#else
  return(atoi(defline+i+5));
#endif
}


int
findSize(char *defline) {
  int   i=0;

#ifdef SNPS_ARE_RS
  return(1);
#else
  while ((defline[i+5] != 0) && ((defline[i]   != '/') ||
                                 (defline[i+1] != 's') ||
                                 (defline[i+2] != 'i') ||
                                 (defline[i+3] != 'z') ||
                                 (defline[i+4] != 'e') ||
                                 (defline[i+5] != '=')))
    i++;

  if (defline[i] == 0) {
    fprintf(stderr, "size not found in defline '%s'!\n", defline);
    exit(1);
  }

  return(atoi(defline+i+6));
#endif
}



int
printSNP(FILE *F, sim4polish *p) {

  char *SNPid = findSNPid(p->estDefLine);
  char *GENid = findGENid(p->genDefLine);
  int   pos   = findPosition(p->estDefLine);
  int   siz   = findSize(p->estDefLine);
  int   exonWithSNP   = -1;
  int   i             = 0;
  int   seqOffset     = 0;
  int   bpToExamine   = 0;
  int   examinePos    = 0;
  int   genPosition   = 0;

  //  Find the exon with the SNP
  //
  for (i=0; i<p->numExons; i++) {
    if (((p->exons[i].estFrom-1) <= pos) && (pos <= (p->exons[i].estTo-1)))
      exonWithSNP = i;
  }

  if (exonWithSNP == -1)
    return(1);

  //  If the match is complement, then the alignment is printed using
  //  the reverse complemented SNP sequence, and so we need to find
  //  the offset at the end of the sequence (not always the same as
  //  the offset at the start of the sequence).
  //
  seqOffset  = pos;
  if (p->matchOrientation == MATCH_COMPLEMENT)
    seqOffset = p->estLen - pos - siz;

  //  Now, we examine the alignment strings to decide exactly
  //  where the SNP is located in the genomic.
  //
  bpToExamine = seqOffset - p->exons[exonWithSNP].estFrom + 1;
  examinePos  = 0;
  genPosition = p->genLo + p->exons[exonWithSNP].genFrom - 1;

  while (bpToExamine > 0) {

    //  If the SNP alignment eats up a base pair, decrement
    //  the number of bp left to examine.
    //
    if (p->exons[exonWithSNP].estAlignment[examinePos] != '-')
      bpToExamine--;

    //  If the the genomic alignment is not a gap, increment the
    //  position.
    //
    if (p->exons[exonWithSNP].genAlignment[examinePos] != '-')
      genPosition++;

    examinePos++;
  }

  fprintf(F, "%s %s %d %s global[%d %d] exon[%d %d %d %d]\n",
          SNPid, GENid,
          genPosition,
          (p->matchOrientation == MATCH_FORWARD) ? "forward" : "complement",
          p->percentIdentity,
          p->querySeqIdentity,
          p->numExons,
          exonWithSNP,
          p->exons[exonWithSNP].percentIdentity,
          (int)floor(100.0 * (double)p->exons[exonWithSNP].numMatches / (double)p->estLen));

  free(SNPid);
  free(GENid);

  return(0);
}





void
pickBestSlave(sim4polish **p, int pNum) {
  int   numMulti  = 0;
  int   numFailed = 0;
  int   i;

  //  Count the number of matches that have more than one exon
  //
  for (i=0; i<pNum; i++)
    if (p[i]->numExons > 1)
      numMulti++;


  if (pNum == 1) {
    if (numMulti == 0) {
      ss++;

      if (singleSingleFile)
        printPolish(singleSingleFile, p[0]);

      if (validSNPMap)
        if (printSNP(validSNPMap, p[0]))
          if (failedSNPMap) {
            numFailed++;
            printPolish(failedSNPMap, p[0]);
          }
    } else {
      sm++;

      if (singleMultiFile)
        printPolish(singleMultiFile, p[0]);

      if (validSNPMap)
        if (printSNP(validSNPMap, p[0]))
          if (failedSNPMap) {
            numFailed++;
            printPolish(failedSNPMap, p[0]);
          }
    }
  } else {
    if (numMulti == 0) {
      ms++;

      if (multiSingleFile)
        for (i=0; i<pNum; i++)
          printPolish(multiSingleFile, p[i]);

      if (validSNPMap)
        for (i=0; i<pNum; i++)
          if (printSNP(validSNPMap, p[i]))
            if (failedSNPMap) {
              numFailed++;
              printPolish(failedSNPMap, p[i]);
            }
    } else {
      mm++;

      if (multiMultiFile)
        for (i=0; i<pNum; i++)
          printPolish(multiMultiFile, p[i]);

      if (validSNPMap)
        for (i=0; i<pNum; i++)
          if (printSNP(validSNPMap, p[i]))
            if (failedSNPMap) {
              numFailed++;
              printPolish(failedSNPMap, p[i]);
            }
    }
  }


  if (numFailed == pNum)
    failedsnps++;

  failedmatches += numFailed;
}




//  Just a wrapper around the real best picker, so that we can easily
//  destroy polishes when we're done.
//
void
pickBest(sim4polish **p, int pNum) {
  int i;

  pickBestSlave(p, pNum);

  for (i=0; i<pNum; i++)
    destroyPolish(p[i]);
}


int
main(int argc, char **argv) {
  int          arg    = 1;
  int          pNum   = 0;
  int          pAlloc = 8388608;
  sim4polish **p      = 0L;
  sim4polish  *q      = 0L;
  int          estID  = ~0;

  int          percentID = 0;
  int          percentCO = 0;

  validSNPMap   = 0L;
  failedSNPMap  = 0L;
  lowQualityMap = 0L;

  while (arg < argc) {
    if        (strncmp(argv[arg], "-i", 2) == 0) {
      percentID = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-c", 2) == 0) {
      percentCO = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-F", 2) == 0) {
      errno = 0;
      failedSNPMap = fopen(argv[++arg], "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
    } else if (strncmp(argv[arg], "-O", 2) == 0) {
      errno = 0;
      validSNPMap = fopen(argv[++arg], "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
    } else if (strncmp(argv[arg], "-D", 2) == 0) {
      char   name[1025];

      errno = 0;

      arg++;

      sprintf(name, "%s-multi-multi", argv[arg]);
      multiMultiFile = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", name, strerror(errno));
        exit(1);
      }

      sprintf(name, "%s-multi-single", argv[arg]);
      multiSingleFile = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", name, strerror(errno));
        exit(1);
      }

      sprintf(name, "%s-single-multi", argv[arg]);
      singleMultiFile = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", name, strerror(errno));
        exit(1);
      }

      sprintf(name, "%s-single-single", argv[arg]);
      singleSingleFile = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", name, strerror(errno));
        exit(1);
      }
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }


  //  Show help if we don't get an output file
  //
  if (validSNPMap == 0L) {
    fprintf(stderr, "usage: %s [options]\n", argv[0]);
    fprintf(stderr, "             -i min-identity     filter matches on percent identity\n");
    fprintf(stderr, "             -c min-coverage     filter matches on percent coverage\n");
    fprintf(stderr, "             -F failed           save matches that do not contain the\n");
    fprintf(stderr, "                                 to the file 'failed'\n");
    fprintf(stderr, "             -O output           save the parsed SNPs to the file\n");
    fprintf(stderr, "                                 'output'\n");
    fprintf(stderr, "             -D prefix           report debugging stuff into files\n");
    fprintf(stderr, "                                 prefixed with 'prefix'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "             only -O is required.  Input is read from stdin.\n");
    fprintf(stderr, "\n");
    exit(1);
  }



  //  Read polishes, picking the best when we see a change in
  //  the estID.

  p = (sim4polish **)malloc(sizeof(sim4polish *) * pAlloc);

  while ((q = readPolish(stdin)) != 0L) {

    //printSNP(stdout, q);

    if (q->estID < estID) {
      fprintf(stderr, "ERROR:  Polishes not sorted by SNP idx!\n");
      exit(1);
    }

    if ((q->estID != estID) && (pNum > 0)) {
      pickBest(p, pNum);
      pNum  = 0;
    }

    //  Reallocate pointers?
    //
    if (pNum >= pAlloc) {
      p = (sim4polish **)realloc(p, sizeof(sim4polish *) * (pAlloc *= 2));
      if (p == 0L) {
        fprintf(stderr, "Out of memory: Couldn't allocate space for polish pointers.\n");
        exit(1);
      }
    }

    estID     = q->estID;

    if ((q->percentIdentity >= percentID) &&
        (q->querySeqIdentity >= percentCO)) {
      p[pNum++] = q;
    } else {
      if (lowQualityMap)
        printPolish(lowQualityMap, q);
      destroyPolish(q);
    }
  }

  if (pNum > 0)
    pickBest(p, pNum);

  fprintf(stderr, "SNPs with:\n");
  fprintf(stderr, "  single hit, single exon:        %6d\n", ss);
  fprintf(stderr, "  single hit, multiple exons:     %6d\n", sm);
  fprintf(stderr, "  multiple hits, single exon:     %6d\n", ms);
  fprintf(stderr, "  multiple hits, multiple exons:  %6d\n", mm);
  fprintf(stderr, "\n");
  fprintf(stderr, "SNPs that failed:                 %6d\n", failedsnps);
  fprintf(stderr, "matches that failed:              %6d\n", failedmatches);

  fclose(validSNPMap);
  fclose(failedSNPMap);
  //fclose(lowQualityMap);  //  It's stdout!

  fclose(multiMultiFile);
  fclose(multiSingleFile);
  fclose(singleMultiFile);
  fclose(singleSingleFile);

  return(0);
}

