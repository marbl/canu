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


char *usage =
"usage: %s [options]\n"
"             -i min-identity     filter matches on percent identity\n"
"             -c min-coverage     filter matches on percent coverage\n"
"             -F failed           save matches that do not contain the\n"
"                                 to the file 'failed'\n"
"             -O output           save the parsed SNPs to the file\n"
"                                 'output'\n"
"             -D prefix           report debugging stuff into files\n"
"                                 prefixed with 'prefix'\n"
"             -d delimiter        Use the single character delimiter as\n"
"                                 the end of the defline ID field.  The\n"
"                                 default is to split on any whitespace.\n"
"\n"
"             -s sizeTag          Use this tag as the size of the snp.\n"
"             -p posTag           Use this tag as the position of the snp.\n"
"\n"
"                                 TAGS: The number immediately after the first\n"
"                                 occurance of the tag will be used.\n"
"                                 The defaults are \"-s /size= -p /pos=\"\n"
"\n"
"             -o offset           An additive offset to the SNP position.\n"
"                                 The default is 0.\n"
"\n"
"             only -O is required.  Input is read from stdin.\n"
"\n"
"             NOTE!  Sizes and sizeTag is NOT IMPLEMENTED!\n"
"                    All SNPs are of size == 1\n"
"\n"
"             If you parse base-based SNPs, the result is returned base-based.\n"
"             You should use an ofset of 0.\n"
"\n"
"             If you parse space-based SNPs, the result is returned base-based.\n"
"             You should use an offset of 1.\n"
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

char  fieldDelimiter = 255;
char *sizeTag        = "/size=";
char *posTag         = "/pos=";
int   positionOffset = 0;

char *
findSNPid(char *defline) {
  char *ret = 0L;
  int   sta = 0;
  int   len = 0;
  int   i = 0;

  if (fieldDelimiter == 255) {
    for (len=1; defline[len] && !isspace(defline[len]); len++)
      ;
  } else {
    for (len=1; defline[len] && defline[len] != fieldDelimiter; len++)
      ;
  }

#if 0
  //  This was used for a set of SNPs with a non-standard defline
  //  structure.  It returns the field between the first '|' and the
  //  next '_'.
  //
  for (len=1; defline[len] && defline[len] != '_'; len++)
    ;
  for (sta=len-1; sta > 0 && defline[sta] != '|'; sta--)
    ;
#endif

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
  int   i = 0;
  char *p = 0L;

  p = strstr(defline, posTag);

#if 0
  while ((defline[i+4] != 0) && ((defline[i]   != '/') ||
                                 (defline[i+1] != 'p') ||
                                 (defline[i+2] != 'o') ||
                                 (defline[i+3] != 's') ||
                                 (defline[i+4] != '=')))
    i++;
#endif

  if (p == 0L) {
    fprintf(stderr, "posTag '%s' not found in defline '%s'!\n", posTag, defline);
    exit(1);
  }

  while (*p && !isdigit(*p))
    p++;

  if (*p == 0) {
    fprintf(stderr, "Found posTag '%s' in defline '%s', but didn't find any numbers!\n", posTag, defline);
    exit(1);
  }

  return(atoi(p) + positionOffset);
}


int
findSize(char *defline) {
  int   i=0;

  //  XXX: We only handle size 1 SNPs.  Probably should be extended
  //  like findPosition.

  return(1);
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
  //  XXX:  Previous version had this as "- siz"  I don't understand
  //  if this is a "size" or just a "+1" thing.
  //
  seqOffset  = pos;
  if (p->matchOrientation == MATCH_COMPLEMENT)
    seqOffset = p->estLen - pos + siz;

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

  fprintf(F, "%s %s %d %c/%c %s global[%d %d] exon[%d %d %d %d]\n",
          SNPid,
          GENid,
          genPosition,
          p->exons[exonWithSNP].estAlignment[examinePos-1],
          p->exons[exonWithSNP].genAlignment[examinePos-1],
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
    } else if (strncmp(argv[arg], "-d", 2) == 0) {
      arg++;
      fieldDelimiter = argv[arg][0];
    } else if (strncmp(argv[arg], "-p", 2) == 0) {
      arg++;
      posTag = argv[arg];
    } else if (strncmp(argv[arg], "-s", 2) == 0) {
      arg++;
      sizeTag = argv[arg];
    } else if (strncmp(argv[arg], "-o", 2) == 0) {
      arg++;
      positionOffset = atoi(argv[arg]);
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }


  //  Show help if we don't get an output file
  //
  if (validSNPMap == 0L) {
    fprintf(stderr, usage, argv[0]);
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
      destroyPolish(q);
    }
  }

  if (pNum > 0)
    pickBest(p, pNum);

  fprintf(stdout, "SNPs with:\n");
  fprintf(stdout, "  single hit, single exon:        %6d\n", ss);
  fprintf(stdout, "  single hit, multiple exons:     %6d\n", sm);
  fprintf(stdout, "  multiple hits, single exon:     %6d\n", ms);
  fprintf(stdout, "  multiple hits, multiple exons:  %6d\n", mm);
  fprintf(stdout, "SNPs that failed:                 %6d\n", failedsnps);
  fprintf(stdout, "matches that failed:              %6d\n", failedmatches);

  fclose(validSNPMap);
  fclose(failedSNPMap);

  fclose(multiMultiFile);
  fclose(multiSingleFile);
  fclose(singleMultiFile);
  fclose(singleSingleFile);

  return(0);
}

